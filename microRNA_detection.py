import argparse, os, platform, re, regex, shutil, subprocess, sys, concurrent.futures
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalOmegaCommandline

work_dir = os.getcwd() + "/"

parser = argparse.ArgumentParser(description='Detects and validates microRNAs from small RNA reads')
parser.add_argument("-o","--output", default = work_dir + "candiates/", help="Output directory for detected candidates (default: ./candidates)")
parser.add_argument("-s","--smallrnas", help="FASTA file containing small RNA reads")
parser.add_argument("-g","--genome", required=True, help="FASTA file containing the genome assembly")
parser.add_argument("-c", "--conservation", required=True, help="CSV file containing miRNA family information")
parser.add_argument("-t", "--threshold", default=0.775, help="Conservation threshold for rejecting mature miRNA candidates")
parser.add_argument("-p", "--ps-structures", action="store_true", help="Produces PostScript visualisations for best graded pre-miRNA candidates")
parser.add_argument("--conservation-dir", help="Folder containing miRNA family sequences and alignments")
parser.add_argument("-f", "--family-list", metavar='family', nargs='+', help="Restrict list of families to those stated")
parser.add_argument("-n", "--num-threads", default=1, help="Number of threads to run in parallel (default: 1)")
parser.add_argument("-j", "--job-id", help="Cluster job ID to select family for current run")
args = parser.parse_args()

if not args.output:
    candidates_dir = work_dir + "output/"
elif not args.output.startswith("/"):
    candidates_dir = work_dir + args.output + "/"
else:
    candidates_dir = args.output + "/"

if args.smallrnas and not args.smallrnas.startswith("/"):
    smallRNA_reads = work_dir + args.smallrnas

if not args.genome.startswith("/"):
    genome_file = work_dir + args.genome
else:
    genome_file = args.genome

if not args.conservation.startswith("/"):
    conservation_data = work_dir + args.conservation
else:
    conservation_data = args.conservation

if not args.conservation_dir:
    conservation_dir = conservation_data.replace(".csv", "/")

if args.conservation_dir and not args.conservation_dir.startswith("/"):
    conservation_dir = work_dir + args.conservation_dir

if args.smallrnas:
    smallRNA_sequences = list(SeqIO.parse(smallRNA_reads, "fasta"))

genome_sequences = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

mirna_list = []
with open(conservation_data, "r") as f:
    for line in f:
        if len(line) < 2 or line.startswith("family") or (args.family_list and line.split(",")[0] not in args.family_list and line.split(",")[0].lower() not in args.family_list):
            continue
        mirna_list.append(line.strip("\n").split(","))

if not os.path.exists(candidates_dir):
    os.makedirs(candidates_dir)

def reverse(sequence):
    # reverse nuceleotide sequences for reverse complement checks
    reverse_sequence = ""
    for nucleotide in sequence:
        if nucleotide.upper() == "A":
            if "U" in sequence:
                reverse_sequence += "U"
            else:
                reverse_sequence += "T"
        elif nucleotide.upper() == "C":
            reverse_sequence += "G"
        elif nucleotide.upper() == "G":
            reverse_sequence += "C"
        elif nucleotide.upper() == "T" or nucleotide.upper() == "U":
            reverse_sequence += "A"
        else:
            reverse_sequence += nucleotide
    return reverse_sequence[::-1]

def detectHairpin(mirna_position, strand, structure):
    # evaluate quality of RNA folding structure according to a hairpin template
    if strand == "3p":
        structure = structure[::-1]
        open_stretch = ")"
        close_stretch = "("
    else:
        open_stretch = "("
        close_stretch = ")"
    
    if  structure.find(open_stretch) == -1: # circle
        return "g"
    
    # calculate number of arms
    arms = []
    arm_closing = True
    for i in range(len(structure)):
        if structure[i] == open_stretch and arm_closing:
            arm_closing = False
            arms.append([i,-1])
        elif structure[i] == close_stretch:
            arm_closing = True
            arms[-1][1] = i
            
            
    arms[-1][1] = len(structure) - 1
    
    if len(arms) != 1: # bifurcations in structure
        return "g"
    
    before_hairpin = True
    bulge_sizes = []
    loop_size = 0
    
    for i in range(len(structure)):
        
        if len(bulge_sizes) < 1:
            if structure[i] != open_stretch:
                continue
            
            if i >= 10:
                return "g"
            
            bulge_sizes.append(0)
            continue
        
        if structure[i] == ".":
            bulge_sizes[-1] += 1
        
        if before_hairpin and structure[i] == close_stretch:
            before_hairpin = False
            loop_size = bulge_sizes[-1]
            bulge_sizes[-1] = 0
            if (i - 5 - loop_size) < 20: # arm too short to contain mature miRNA
                return "g"
        elif (before_hairpin and structure[i] == open_stretch) or (not before_hairpin and structure[i] == close_stretch):
            bulge_sizes.append(0)
        elif not before_hairpin and structure[i] == open_stretch:
            bulge_sizes[-1] = 0
    
    bulge_sizes[-1] = 0
    bulge_size = max(bulge_sizes)
    
    # deduct 1 per max unpaired nucleotides > 4 in stem region
    # deduct 1 per 2 max unpaired nucleotides > 10 in loop region
    return 1 + (0 if bulge_size < 4 else bulge_size - 3) + (0 if loop_size < 12 else int((loop_size - 10)/2))

def detect_mirna(mirna_tuple, max_candidates=10000):
    
    mirna, strand, seed, min_identity, avg_identity, mature_length, hairpin_length = mirna_tuple
    seed = seed.replace("U", "T")
    min_identity = float(min_identity)
    avg_identity = float(avg_identity)
    mature_length = int(mature_length)
    hairpin_length = int(hairpin_length)
    
    min_identity = float(args.threshold)
    if min_identity > 1.0:
        min_identity /= 100
    avg_identity = max(avg_identity, 0)
    
    log_file = candidates_dir + mirna + "-" + strand + ".log"
    with open(log_file, "w") as f:
        f.write("STARTING microRNA detection of {}-{}\nSeed sequence: {}\n".format(mirna, strand, seed))
    
    reference = None
    member_sequences = list(SeqIO.parse(conservation_dir + mirna + "-" + strand + ".fasta", "fasta"))
    for sequence in member_sequences:
        species = sequence.id[0:3]
        if species in ["cel", "dme", "hsa"]:
            reference = sequence
            break
    
    if not reference:
        reference = member_sequences[0]
    
    ref_sequence = str(reference.seq).replace("U", "T")
    span_ref_seed = regex.search(seed+"{e<1}", ref_sequence)
    index_ref_seed = span_ref_seed.start()
    
    # find candidates amongst small RNAs
    mature_file = log_file.replace(".log", ".mature")
    if os.path.exists(mature_file):
        with open(log_file, "a") as f:
            f.write("SKIPPING searching for candidates: {} already exists!\n".format(mature_file))
    else:
        if not args.smallrnas:
            with open(log_file, "a") as f:
                f.write("ABORT searching for candidates: No small RNAs file given!\n".format(mature_file))
            return "{}-{}: ABORTED!".format(mirna, strand)
        with open(log_file, "a") as f:
            f.write("Reference sequence: {} - {}\n".format(reference.id, ref_sequence))
        
        mature_sequences = list()
        for sequence in smallRNA_sequences:
            span_seed = regex.search(seed+"{e<1}", str(sequence.seq))
            if len(sequence.seq) < 15 or not span_seed or span_seed.start() > 4:
                continue
            index_seed = span_seed.start()
            
            if len(sequence.seq) - index_seed < 15:
                continue # skip candidate sequences that are shorter than 16 nucleotides
            
            identity = 0.0
            for i in range(0, len(ref_sequence)):
                if index_seed - (index_ref_seed - i) < 0: # skipping until start of candidate sequence
                    continue
                if index_seed + (i - index_ref_seed) >= len(sequence.seq): # stopping when end of candidate sequence
                    break
                if sequence.seq[index_seed - index_ref_seed + i] == ref_sequence[i]:
                    identity += 1
            
            identity = identity / (i + 1)
            if identity > 1:
                print(sequence.seq, seed, ref_sequence, identity)
                1/0
            if identity >= min_identity:
                new_id = "{}_{}".format(str(100*identity).split(".")[0], sequence.id)
                new_sequence = SeqRecord(sequence.seq, id=new_id, description = "")
                mature_sequences.append(new_sequence)
        
        if len(mature_sequences) < 1:
            with open(log_file, "a") as f:
                f.write("TERMINATING: No candidates found for {}-{} (min. identity {})\n".format(mirna, strand, min_identity))
            return "{}-{}: No mature miRNA candidates!".format(mirna, strand)
        
        with open(log_file, "a") as f:
            f.write("Found {} candidate sequences in small RNA reads file\n".format(len(mature_sequences)))
        
        
        sorted_mature = sorted(mature_sequences, key=lambda sequence: int(sequence.id.split("_")[0]), reverse=True)
        if len(mature_sequences) > max_candidates:
            max_identity = sorted_mature[0].id.split("_")[1]
            i = 0
            while sorted_mature[i].id.split("_")[0] == max_identity or i < 10:
                if sorted_mature[i].id.split("_")[0] != max_identity:
                    max_identity = sorted_mature[i].id.split("_")[0]
                i += 1
                if i == len(sorted_mature):
                    break
            SeqIO.write(sorted_mature[0:(i-1)], mature_file, "fasta")
            with open(log_file, "a") as f:
                f.write(" > Filtered for best {} candidate sequences\n".format(i))
        else:
            SeqIO.write(sorted_mature, mature_file, "fasta")
    
    mature_sequences = SeqIO.to_dict(SeqIO.parse(mature_file, "fasta"))
    
    # find DNA stretches matching candidates
    pre_mirnas_file = log_file.replace(".log", ".pre_mirnas")
    
    if os.path.exists(pre_mirnas_file):
        with open(log_file, "a") as f:
            f.write("SKIPPING scanning genome for matching sequences: {} already exists!\n".format(pre_mirnas_file))
        pre_mirnas = list(SeqIO.parse(pre_mirnas_file, "fasta"))
    else:
        genome_indices = dict()
        pre_mirnas = list()
        genome_keys = set()
        pre_mirnas_pos = dict()
        for candidate_id in mature_sequences:
            candidate_sequence = mature_sequences[candidate_id]
            for seq_id in genome_sequences:
                if candidate_sequence.seq in genome_sequences[seq_id].seq:
                    starts = [m.start() for m in re.finditer('(?='+str(candidate_sequence.seq)+')', str(genome_sequences[seq_id].seq))]
                    for start in starts:
                        if strand == "5p":
                            begin = max(0, start - 10)
                            end = min(start - 10 + hairpin_length, len(genome_sequences[seq_id].seq))
                        else:
                            begin = max(0, start - hairpin_length + mature_length + 10)
                            if len(candidate_sequence) > mature_length + 10:
                                end = min(start + len(candidate_sequence), len(genome_sequences[seq_id].seq))
                            else:
                                end = min(start + mature_length + 10, len(genome_sequences[seq_id].seq))
                        if seq_id not in pre_mirnas_pos:
                            pre_mirnas_pos[seq_id] = []
                        already_detected = False
                        for position in pre_mirnas_pos[seq_id]:
                            if abs(begin - position) < 5:
                                already_detected = True
                                break
                        if already_detected:
                            continue
                        pre_mirnas_pos[seq_id].append(begin)
                        gene_seq_id = "{}_{}_{}".format(candidate_sequence.description, seq_id, begin)
                        if gene_seq_id not in genome_keys:
                            pre_mirnas.append(SeqRecord(Seq(str(genome_sequences[seq_id].seq)[begin : end]), id=gene_seq_id.replace("|", "-"), description=""))
                            genome_keys.add(gene_seq_id)
                        
                reverse_candidate = reverse(str(candidate_sequence.seq))
                if reverse_candidate in genome_sequences[seq_id].seq:
                    reverse_sequence = reverse(str(genome_sequences[seq_id].seq))
                    starts = [m.start() for m in re.finditer('(?='+str(candidate_sequence.seq)+')', reverse_sequence)]
                    for start in starts:
                        if strand == "5p":
                            begin = max(0, start - 10)
                            end = min(start - 10 + hairpin_length, len(reverse_sequence))
                        else:
                            begin = max(0, start - hairpin_length + mature_length + 10)
                            if len(reverse_candidate) > mature_length + 10:
                                end = min(start + len(reverse_candidate) + 10, len(reverse_sequence))
                            else:
                                end = min(start + mature_length + 10, len(reverse_sequence))
                        if seq_id + "-r" not in pre_mirnas_pos:
                            pre_mirnas_pos[seq_id + "-r"] = []
                        already_detected = False
                        for position in pre_mirnas_pos[seq_id + "-r"]:
                            if abs(begin - position) < 5:
                                already_detected = True
                                break
                        if already_detected:
                            continue
                        pre_mirnas_pos[seq_id + "-r"].append(begin)
                        gene_seq_id = "{}_{}_{}-reversed".format(candidate_sequence.description, seq_id, begin)
                        if gene_seq_id not in genome_keys:
                            pre_mirnas.append(SeqRecord(Seq(reverse_sequence[begin : end]), id=gene_seq_id.replace("|", "_"), description=""))
                            genome_keys.add(gene_seq_id)
        
        if len(pre_mirnas) < 1:
            with open(log_file, "a") as f:
                f.write("TERMINATING: No candidates for {}-{} matched to the genome\n".format(mirna, strand))
            return "{}-{}: No pre-miRNA candidates in genome!".format(mirna, strand)
        
        with open(log_file, "a") as f:
            f.write("Found {} genomic sequences matching the candidates\n".format(len(pre_mirnas)))
        
        if len(pre_mirnas) > max_candidates:
            sorted_pre_mirnas = sorted(pre_mirnas, key=lambda sequence: int(sequence.id.split("_")[0]), reverse=True)
            max_identity = sorted_pre_mirnas[0].id.split("_")[0]
            i = 0
            while sorted_pre_mirnas[i].id.split("_")[0] == max_identity or i < 10:
                if sorted_pre_mirnas[i].id.split("_")[0] != max_identity:
                    max_identity = sorted_pre_mirnas[i].id.split("_")[0]
                i += 1
                if i == len(sorted_pre_mirnas):
                    break
            SeqIO.write(sorted_pre_mirnas[0:(i-1)], pre_mirnas_file, "fasta")
            with open(log_file, "a") as f:
                f.write(" > Filtered for best {} genomic sequences\n".format(i))
        else:
            SeqIO.write(pre_mirnas, pre_mirnas_file, "fasta")
    
    
    # detect hairpin structures
    graded_file = "{}{}-{}.best_graded".format(candidates_dir, mirna, strand)
    if os.path.exists(graded_file):
        with open(log_file, "a") as f:
            f.write("SKIPPING computing secondary structures: {} already exists!\n".format(graded_file))
    else:
        fold_dir = candidates_dir + mirna + "-" + strand + "/"
        if not os.path.exists(fold_dir):
            os.makedirs(fold_dir)
        best_grade = -1
        highest_conservation = 0
        best_graded = []
        for seq in pre_mirnas:
            mature_id = "_".join(seq.id.split("_")[0:2])
            mature_sequence = str(mature_sequences[mature_id].seq).upper()
            
            sequence = str(seq.seq).upper()
            short_sequences = ">{}-{}bp\n{}\n".format(seq.id, len(seq.seq), sequence)
            for i in range(1, len(sequence) - 43):
                if strand == "5p":
                    if mature_sequence not in sequence[:(len(sequence) - i)]:
                           break
                    short_sequences += ">{}-{}bp\n{}\n".format(seq.id, len(sequence)-i, sequence[:(len(sequence) - i)])
                else:
                    if mature_sequence not in sequence[i:]:
                           break
                    short_sequences += ">{}-{}bp\n{}\n".format(seq.id, len(sequence)-i, sequence[i:])
            
            with open(graded_file.replace(".best_graded", ".shortened"), "w") as f:
                f.write(short_sequences)
            
            p = subprocess.Popen(["RNAfold", "--noPS"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            output = p.communicate(input=short_sequences.encode())[0].decode().split("\n")
            
            for i in range(0, len(output), 3):
                if len(output[i]) < 1:
                    continue
                name = output[i].strip("> \n")
                sequence = output[i+1].strip(" \n").upper().replace("U", "T")
                structure = output[i+2].strip(" \n").split()[0]
                energy = output[i+2].strip(" \n").split()[1]
                
                mirna_start = sequence.find(mature_sequence)
                mirna_position = [mirna_start, mirna_start + len(mature_sequence)]
                if -1 in mirna_position:
                    print("ERROR: mature sequence not in pre-miRNA!\n ", mature_id, seq.id, mature_sequence, sequence)
                    continue
                grade = detectHairpin(mirna_position, strand, structure)
                if grade == "g":
                    continue
                if best_grade == -1 or grade < best_grade:
                    best_grade = grade
                    best_graded = []
                    highest_conservation = int(name.split("_")[0])
                if grade == best_grade:
                    best_graded.append(SeqRecord(Seq(sequence), id=str(grade) + "_" + name, description=""))
                    highest_conservation = max(int(name.split("_")[0]), highest_conservation)
                if grade == 1:
                    break
        
        # clean up
        os.remove(graded_file.replace(".best_graded", ".shortened"))
        
        if len(best_graded) < 1:
            with open(log_file, "a") as f:
                f.write("No hairpin loop structures found!\n")
            shutil.rmtree(fold_dir)
            return "{}-{}: No hairpin loop structures found!".format(mirna, strand)
        else:
            SeqIO.write(best_graded, graded_file, "fasta")
            if args.ps_structures: # generate structure files
                os.chdir(fold_dir)
                p = subprocess.Popen(["RNAfold", "-i", graded_file, "-o"])
                p.wait()
                for f in os.listdir(fold_dir):
                    if not f.endswith(".ps"):
                        continue
                    new_f = fold_dir + mirna + "-" + strand + "_" + f
                    p = subprocess.Popen(["mv", fold_dir + f.replace("_ss.ps", ".fold"), new_f.replace("_ss.ps", ".fold")])
                    p.wait()
                    p = subprocess.Popen(["mv", fold_dir + f, new_f])
                    p.wait()
                    if platform.system() == "Linux":
                        p = subprocess.Popen(["ps2pdf", "-dEPSCrop", new_f, new_f.replace(".ps", ".pdf")])
                    else:
                        p = subprocess.Popen(["pstopdf", new_f, "-o", new_f.replace(".ps", ".pdf")])
                    p.wait()
                os.chdir(work_dir)
            else: # delete temporary structures directory
                shutil.rmtree(fold_dir)
            
            with open(log_file, "a") as f:
                f.write("{}-{}: Best graded hairpin structure: {} ({}% conservation)!\n".format(mirna, strand, best_grade, highest_conservation))
            
            return "{}-{}: Best graded hairpin structure: {} ({}% conservation)!".format(mirna, strand, best_grade, highest_conservation)
    
    return "{}-{}: DONE!".format(mirna, strand)

if args.job_id:
    detect_mirna(mirna_list[int(args.job_id) - 1])
else:
    with concurrent.futures.ProcessPoolExecutor(max_workers=int(args.num_threads)) as executor:
        futures = [executor.submit(detect_mirna, mirna) for mirna in mirna_list]

    results = concurrent.futures.wait(futures)
    for rd in results.done:
        print(rd.result())


