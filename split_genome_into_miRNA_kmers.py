import argparse, concurrent.futures, os, subprocess, sys
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Splits genome into 6-mers for efficient microRNA seed scan')
parser.add_argument("-c", "--conservation", required=True, help="CSV file containing miRNA family information")
parser.add_argument("-g","--genome", required=True, help="FASTA file containing the genome assembly")
parser.add_argument("-n", "--num-threads", default=1, help="Number of threads to run in parallel (default: 1)")
args = parser.parse_args()

mirna_families_file = args.conservation
genome_file = args.genome
genome_dir = os.getcwd() + "/"
kmers_dir = genome_dir + genome_file.split("/")[-1].split(".")[0] + "_kmers/"

# split genome into 6-mers
if not os.path.exists(kmers_dir):
    os.makedirs(kmers_dir)
    print("SPLITTING: {}".format(genome_file))    

    genome_sequences = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

    b = set(genome_sequences.keys())
    while len(b) > 0:
        print("Sequences left:", len(b))
        for gs in genome_sequences:
            if os.path.exists(kmers_dir + gs + ".kmers"):
                if gs in b:
                    b.remove(gs)
                    print("Sequences left:", len(b))
                #print("SKIPPING: {}.kmers already exists!".format(gs))
                continue
            
            seq = str(genome_sequences[gs].seq).upper()
            out = ["{} {}\n".format(i, seq[i:i+6]) for i in range(0, len(seq)-6)]
            try:
                with open(kmers_dir + gs + ".kmers", "w") as f:
                    f.write("".join(out))
            except OSError: # Buffer issues when writing 
                print("Caught OSError: Trying to split writing!")
                with open(kmers_dir + gs + ".kmers", "w") as f:
                    f.write("".join(out[:100000000]))
                with open(kmers_dir + gs + ".kmers", "a") as f:
                    f.write("".join(out[100000000:]))
    
    os.makedirs(kmers_dir + "pos_unfiltered")

# create kmer grep bash script
if os.path.exists(kmers_dir + "kmer_grep.sh"):
    print("WARNING: {}kmer_grep.sh already exists!".format(genome_dir))
else:
    with open(kmers_dir + "kmer_grep.sh", "w") as f:
        f.write('if [ "$3" == "r" ]\nthen\nout=$1/pos_unfiltered/$2_reversed.pos\nelse\nout=$1/pos_unfiltered/$2.pos\nfi\nfor i in $1/*.kmers\n\ndo\necho `basename $i .kmers` `grep $2 $i | cut -d " " -f1`\ndone > ${out}')

# grep seeds from split files
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

seeds = []
seeds_reversed = []

with open(mirna_families_file, "r") as f:
    for line in f:
        if line.startswith("family") or len(line) < 2:
            continue
        seeds.append(line.split(",")[2].replace("U", "T"))
        seeds_reversed.append(reverse(line.split(",")[2].replace("U", "T")))

all_seeds = (set(seeds) | set(seeds_reversed))

def grep(seed):
    ret = ""
    if seed in seeds:
        if os.path.exists(kmers_dir + "pos_unfiltered/" + seed + ".pos"):
            ret += seed + ".pos already exists!\n"
        else:
            p = subprocess.Popen(["sh", kmers_dir + "kmer_grep.sh", kmers_dir, seed])
            p.wait()
    if seed in seeds_reversed:
        if os.path.exists(kmers_dir + "pos_unfiltered/" + seed + "_reversed.pos"):
            ret += seed + "_reversed.pos already exists!"
        else:
            p = subprocess.Popen(["sh", kmers_dir + "kmer_grep.sh", kmers_dir, seed, "r"])
            p.wait()
    return ret

with concurrent.futures.ProcessPoolExecutor(max_workers=int(args.num_threads)) as executor:
    futures = [executor.submit(grep, seed) for seed in all_seeds]

results = concurrent.futures.wait(futures)
for rd in results.done:
    print(rd.result())

# create and execute bash script to get all non-empty grep results
with open(kmers_dir + "grep.sh", "w") as f:
    f.write('for s in {}/pos_unfiltered/*.pos; do grep " " $s > {}/`basename $s`; done'.format(kmers_dir, kmers_dir))

p = subprocess.Popen(["sh", kmers_dir + "grep.sh"])
p.wait()

# clean up
os.remove(kmers_dir + "kmer_grep.sh")
os.remove(kmers_dir + "grep.sh")
