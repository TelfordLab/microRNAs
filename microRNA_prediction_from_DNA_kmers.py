import argparse, concurrent.futures, os, platform, re, regex, shutil, subprocess, sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalOmegaCommandline

work_dir = os.getcwd() + "/"

parser = argparse.ArgumentParser(description='Predicts microRNA candidates from genome data')
parser.add_argument("-o","--output", default = work_dir + "candiates/", help="Output directory for detected candidates (default: ./candidates)")
parser.add_argument("-g","--genome", required=True, help="FASTA file containing the genome assembly")
parser.add_argument("-kd","--kmer_dir", required=True, help="directory containing kmer files")
parser.add_argument("-c", "--conservation", required=True, help="CSV file containing miRNA family information")
parser.add_argument("--conservation-dir", help="Folder containing miRNA family sequences and alignments")
parser.add_argument("-t", "--threshold", default=0.775, help="Threshold for keeping predicted mature miRNA candidates")
parser.add_argument("-f", "--family-list", metavar='family', nargs='+', help="Restrict list of families to those stated")
parser.add_argument("-n", "--num-threads", default=1, help="Number of threads to run in parallel (default: 1)")
parser.add_argument("-j", "--job-id", help="Cluster job ID to select family for current run")
args = parser.parse_args()

if args.output and not args.output.startswith("/"):
    candidates_dir = work_dir + args.output + "/"
else:
    candidates_dir = args.output + "/"

log_file = "/" + candidates_dir.strip("/") + ".log"

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

if not args.kmer_dir.startswith("/"):
    kmer_dir = work_dir + args.kmer_dir
else:
    kmer_dir = args.kmer_dir

genome_sequences = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

if not os.path.exists(candidates_dir):
    os.makedirs(candidates_dir)

mirna_list = []
with open(conservation_data, "r") as f:
    for line in f:
        if len(line) < 2 or line.startswith("family") or (args.family_list and line.split(",")[0] not in args.family_list and line.split(",")[0].lower() not in args.family_list):
            continue
        mirna_list.append(line.strip("\n").split(","))

def reverse(sequence):
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

def findMatureCandidateDNA(mirna_tuple):
    mirna, strand, seed, min_identity, avg_identity, mature_length, hairpin_length = mirna_tuple
    seed = seed.replace("U", "T")
    min_identity = float(args.threshold)
    if min_identity > 1.0:
        min_identity /= 100
    avg_identity = float(avg_identity)
    mature_length = int(mature_length)
    hairpin_length = int(hairpin_length)
    
    log_file = candidates_dir + mirna + "-" + strand + ".log"
    
    log_out = "Predicting miRNA candidates for {}\nMinimal conservation threshold: {}\n".format(mirna, min_identity)
    
    mature_file = candidates_dir + mirna + "-" + strand + ".mature"
    if os.path.exists(mature_file):
        with open(log_file, "a") as f:
            f.write("SKIPPING! " + mature_file + " exists!")
        return
    
    if min_identity < 0:
        min_identity = 0.75
    
    if avg_identity < 0:
        avg_identity = 0.75
    
    ref_sequence = None
    sequences = list(SeqIO.parse(conservation_dir + mirna + "-" + strand + ".fasta", "fasta"))
    for sequence in sequences:
        species = sequence.id[0:3]
        if species in ["cel", "dme", "hsa"] or (mirna == "miR-137" and species == "mmu"):
            ref_sequence = sequence
            break
    
    if not ref_sequence:
        ref_sequence = sequences[0]
    
    span_ref_seed = regex.search(seed+"{e<1}", str(ref_sequence.seq).replace("U", "T"))
    if not span_ref_seed:
        return log_out + "ERROR: Seed '{}' not in ref sequence {}: {}!\n\n".format(seed, reference.id, ref_sequence)
    
    index_ref_seed = span_ref_seed.start()
    
    
    mature_candidates = []
    genome_seqs = []
    
    seed_positions = dict()
    with open(kmer_dir + seed + ".pos", "r") as f:
        for line in f:
            seed_positions[line.split(" ", 1)[0]] = [int(s) for s in line.strip(" \n").split(" ")[1:]]
    
    for seq_id in seed_positions:
        starts = seed_positions[seq_id]
        for start in starts:
            if len(genome_sequences[seq_id].seq) - start < 40:
                continue
            mature_candidate = str(genome_sequences[seq_id].seq)[start-1:start-1+mature_length].upper()
            index_seed = 1
            identity = 0.0
            for i in range(0, len(ref_sequence)):
                if index_seed - (index_ref_seed - i) < 0: # skipping until start of candidate sequence
                    continue
                if index_seed + (i - index_ref_seed) >= len(mature_candidate): # stopping when end of candidate sequence
                    break
                if mature_candidate[index_seed - index_ref_seed + i] == ref_sequence.seq[i].replace("U", "T"):
                    identity += 1
            identity = identity / (i + 1)
            if identity > 1:
                print(sequence.seq, seed, str(ref_sequence.seq).replace("U", "T"), identity)
                1/0
            if identity >= min_identity:
                if strand == "5p":
                    begin = max(0, start - 11)
                    end = min(start - 11 + hairpin_length, len(genome_sequences[seq_id].seq))
                else:
                    begin = max(0, start - hairpin_length + mature_length + 10)
                    end = min(start + mature_length + 10, len(genome_sequences[seq_id].seq))
                if end - begin < 10:
                    print(mirna_tuple)
                    print(seed, start, begin, end)
                    print(hairpin_length, mature_length)
                    print(genome_sequences[seq_id].seq)
                    1/0
                mature_candidates.append([mature_candidate, int(100*identity), "{}_{}".format(seq_id, begin + 5)])
                genome_seqs.append([str(genome_sequences[seq_id].seq)[begin:end].upper(), int(100*identity), "{}_{}-{}".format(seq_id, begin, end)])
    
    reverse_seed_positions = dict()
    with open(kmer_dir + reverse(seed) + "_reversed.pos", "r") as f:
        for line in f:
            reverse_seed_positions[line.split(" ", 1)[0]] = [int(s) for s in line.strip(" \n").split(" ")[1:]]
    
    for seq_id in reverse_seed_positions:
        starts = reverse_seed_positions[seq_id]
        
        for start in starts:
            if start < 34:
                continue
            
            mature_candidate = reverse(str(genome_sequences[seq_id].seq)[start+7-mature_length:start+7].upper())
            
            index_seed = 1
            identity = 0.0
            for i in range(0, len(ref_sequence)):
                if index_seed - (index_ref_seed - i) < 0: # skipping until start of candidate sequence
                    continue
                if index_seed + (i - index_ref_seed) >= len(mature_candidate): # stopping when end of candidate sequence
                    break
                if mature_candidate[index_seed - index_ref_seed + i] == ref_sequence.seq[i].replace("U", "T"):
                    identity += 1
            identity = identity / (i + 1)
            
            if identity > 1:
                print(sequence.seq, seed, str(ref_sequence.seq).replace("U", "T"), identity)
                1/0
            
            if identity >= min_identity:
                if strand == "5p":
                    begin = max(0, start - hairpin_length + 17)
                else:
                    begin = max(0, start - mature_length - 3)
                end = min(start + 17, len(genome_sequences[seq_id].seq))
                
                mature_candidates.append([mature_candidate, int(100*identity), "{}_{}-reversed".format(seq_id, begin + 5)])
                genome_seqs.append([reverse(str(genome_sequences[seq_id].seq)[begin:end]), int(100*identity), "{}_{}-{}-reversed".format(seq_id, begin, end)])
    
    if mature_candidates:
        sorted_mature = sorted(mature_candidates, reverse=True, key=lambda mc: mc[1])
        with open(mature_file, "w") as f:
            for i in range(len(sorted_mature)):
                f.write(">{}_DNA{}\n{}\n".format(sorted_mature[i][1], i+1, sorted_mature[i][0].upper()))
        sorted_genome_seqs = sorted(genome_seqs, reverse=True, key=lambda mc: mc[1])
        with open(mature_file.replace(".mature", ".pre_mirnas"), "w") as f:
            for i in range(len(sorted_genome_seqs)):
                if len(sorted_genome_seqs[i][0]) < 5:
                    print(sorted_genome_seqs[i])
                f.write(">{}_DNA{}_{}\n{}\n".format(sorted_genome_seqs[i][1], i+1, sorted_genome_seqs[i][2], sorted_genome_seqs[i][0].upper()))
        with open(log_file, "w") as f:
            f.write(log_out + "{} mature miRNA candidates found for {}-{}!\n\n".format(len(mature_candidates), mirna, strand))
    else:
        with open(mature_file.replace(".mature", ".log"), "w") as f:
            f.write("No mature miRNA candidates found for {}-{} (min. identity {})\n".format(mirna, strand, min_identity))
    
    return "DONE: {}-{}!".format(mirna, strand)

if args.job_id:
    print(findMatureCandidateDNA(mirna_list[int(args.job_id) - 1]))
else:
    with concurrent.futures.ProcessPoolExecutor(max_workers=int(args.num_threads)) as executor:
        futures = [executor.submit(findMatureCandidateDNA, mirna) for mirna in mirna_list]
    
    results = concurrent.futures.wait(futures)
    for rd in results.done:
        print(rd.result())


