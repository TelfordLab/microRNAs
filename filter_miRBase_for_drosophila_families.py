import os, subprocess, sys
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline

if len(sys.argv) != 3:
    print("Use: python3 filter_miRBase_for_drosophila_families.py <miRBase_mature.fa> <miRBase_hairpin.fa>")
    sys.exit(0)

work_dir = os.getcwd() + "/"
conservation_dir = work_dir + "conservation_drosophila/"
miRBase_matures = SeqIO.to_dict(SeqIO.parse(sys.argv[1], "fasta")) # mature.fa from miRBase
miRBase_hairpins = SeqIO.to_dict(SeqIO.parse(sys.argv[2], "fasta")) # hairpin.fa from miRBase

nonbilateria = {"hma", "nve", "aqu", "lco", "sci"}
protostomia = {"isc","rmi","tur","dpu","mja","aae","aga","ame","api","bmo","cqu","dan","der","dgr","dme","dmo","dpe","dps","dse","dsi","dvi","dwi","dya","hme","lmi","mse","ngi","nlo","nvi","pxy","tca","smr","asu","bma","cbn","cbr","cel","crm","hco","ppc","prd","str","cte","gpy","tre","hru","lgi","cla","egr","emu","gsa","sja","sma","sme"}
deuterostomia = {"xbo","bbe","bfl","cin","csa","odi","pma","xla","xtr","gga","tgu","cfa","ocu","aja","eca","efu","mdo","meu","sha","age","lla","sla","mml","mne","pbi","ggo","hsa","ppa","ppy","ptr","ssy","lca","cgr","mmu","rno","bta","chi","oar","tch","ssc","dre","fru","hhi","ipu","ola","pol","ssa","tni","aca","oha","pmi","spu","sko"}
metazoa = (nonbilateria | protostomia | deuterostomia)

flies = {"dan", "der", "dgr", "dme", "dmo", "dpe", "dps", "dse", "dsi", "dvi", "dwi", "dya"}

metazoa_families = dict() # stores all families present in Metazoa
drosophila_families = [] # filtered list for families only present in Drosophila species
miRBase_families = dict() # used to store information about seed and conservation properties

if not os.path.exists(conservation_dir):
    os.mkdir(conservation_dir)

for sequence_name in miRBase_matures:
    species = sequence_name.split("-")[0]
    
    if species not in metazoa:
        continue
    
    if sequence_name.split("-")[1] not in ["miR", "lin", "let", "lsy"]:
        if "miR" in sequence_name.split("-")[1]:
            family_name = sequence_name.split("-")[1].replace("miR", "miR")
        elif sequence_name.split("-")[1] == "bantam":
            family_name = sequence_name.split("-")[1]
        else:
            print(sequence_name)
            1/0
    else:
        family_name = "-".join([sequence_name.split("-")[1], sequence_name.split("-")[2]])
        if family_name != "miR-iab":
            while not family_name[-1].isdigit():
                family_name = family_name[:-1]
    
    if family_name not in metazoa_families:
        metazoa_families[family_name] = set()
    
    metazoa_families[family_name].add(species)

for family in metazoa_families:
    if len(metazoa_families[family]) > 1 and len(metazoa_families[family] - flies) < 1:
        drosophila_families.append(family)

class MicroRNA_family:
    
    def __init__(self, name, work_dir):
        self.name = name
        self.seed = "-"
        self.sequences = []
        self.alignment = None
        self.min_identity = -1.0
        self.avg_identity = -1.0
        self.work_dir = work_dir + "/"
    
    def update_alignment(self):
        fasta_file = self.work_dir + self.name + ".fasta"
        SeqIO.write(self.sequences, fasta_file, "fasta")
        
        if len(self.sequences) < 2:
            return
        
        msa_file = fasta_file.replace(".fasta", ".msa")
        if os.path.exists(msa_file):
            os.remove(msa_file)
        clustalomega_cline = ClustalOmegaCommandline(infile=fasta_file, outfile=msa_file, verbose=True, auto=True)
        stderr, stdout = clustalomega_cline()
        self.alignment = AlignIO.read(msa_file, "fasta")
    
    def compute_seed(self):
        self.update_alignment()
        if len(self.sequences) < 2:
            self.seed = str(self.sequences[0].seq)[1:7]
            return True
        
        self.seed = "-"
        for seed_start in range(1, len(self.alignment[0].seq)):
            counts = dict()
            for sequence in self.alignment:
                if str(sequence.seq)[seed_start:seed_start+6] not in counts:
                    counts[str(sequence.seq)[seed_start:seed_start+6]] = 0.0
                counts[str(sequence.seq)[seed_start:seed_start+6]] += 1
            
            for s in counts:
                if counts[s] >= min(0.9, 1.0 - 1.0/len(self.alignment)) * len(self.alignment):
                    self.seed = s
                    break
            
            if "-" not in self.seed:
                break
        
        return "-" not in self.seed
    
    def compute_conservation(self):
        self.update_alignment()
        if len(self.sequences) < 2:
            return
        msa_file = self.work_dir + self.name + ".msa"
        mview_file = msa_file.replace(".msa", ".mview")
        with open(mview_file, "w") as f:
            p = subprocess.Popen(["mview", "-in", "fasta", msa_file], stdout = f)
            p.wait()
        
        with open(mview_file, "r") as f:
            mview_results = f.readlines()
        
        self.min_identity = 1.0
        self.avg_identity = 0.0
        count = 0
        for mview_result in mview_results:
            if "%" not in mview_result:
                continue
            result = mview_result.strip(" ").split()
            identity = float(result[-2].strip("%"))/100.0
            if identity < self.min_identity:
                self.min_identity = identity
            self.avg_identity += identity
            count += 1
        self.avg_identity /= count

for sequence_name in miRBase_matures:
    
    family_name = "-".join(sequence_name.split("-")[1:3])
    if family_name not in drosophila_families:
        continue
    
    hairpin_name = "-".join(sequence_name.lower().split("-")[0:3])
    
    if sequence_name.endswith("-5p") or sequence_name.endswith("-3p"):
        strand_name = sequence_name.split("-")[-1]
    else:
        if hairpin_name in miRBase_hairpins:
            sequence = str(miRBase_hairpins[hairpin_name].seq)
        elif hairpin_name + "-1" in miRBase_hairpins:
            sequence = str(miRBase_hairpins[hairpin_name + "-1"].seq)
        else:
            print(hairpin_name, hairpin_name + "-1")
            1/0
        
        position = sequence.find(str(miRBase_matures[sequence_name].seq))
        if position < len(sequence) / 2 - 20:
            strand_name = "5p"
        elif position > len(sequence) / 2:
            strand_name = "3p"
        else: # mature miRNA sequence is part of the loop region > manual inspection
            print(sequence_name, hairpin_name)
            print(miRBase_matures[sequence_name].seq)
            print(sequence)
            print(position, len(sequence) / 2)
            1/0
    
    if family_name + "-" + strand_name not in miRBase_families:
        miRBase_families[family_name + "-" + strand_name] = MicroRNA_family(family_name + "-" + strand_name, conservation_dir)
    
    miRBase_families[family_name + "-" + strand_name].sequences.append(miRBase_matures[sequence_name])

# check families for conserved seed region and collect information
csv_out = "family,strand,seed,min_identity,avg_identity,max_mature_size,max_hairpin_size\n"
for family_name in sorted(miRBase_families):
    if not miRBase_families[family_name].compute_seed():
        print("No seed for {}!".format(family_name))
        continue
    
    miRBase_families[family_name].compute_conservation()
    
    if miRBase_families[family_name].min_identity < 0:
        continue
    
    max_mature = 0
    max_hairpin = 0
    for sequence in miRBase_families[family_name].sequences:
        max_mature = max(len(sequence), max_mature)
        
        hairpin_name = "-".join(sequence.id.lower().split("-")[0:3])
        if hairpin_name.endswith("-5p") or hairpin_name.endswith("-3p"):
            hairpin_name = hairpin_name.rsplit("-", 1)[0]
        
        if hairpin_name in miRBase_hairpins:
            max_hairpin = max(len(miRBase_hairpins[hairpin_name].seq), max_hairpin)
        else:
            for sequence_name in miRBase_hairpins:
                if sequence_name.startswith(hairpin_name + "-"):
                    max_hairpin = max(len(miRBase_hairpins[sequence_name].seq), max_hairpin)
    
    csv_out += "{},{},{},{:.3f},{:.3f},{},{}\n".format(miRBase_families[family_name].name.rsplit("-", 1)[0], miRBase_families[family_name].name.rsplit("-", 1)[1], miRBase_families[family_name].seed, miRBase_families[family_name].min_identity, miRBase_families[family_name].avg_identity, max_mature, max_hairpin)

# write miRNA family information to file
conservation_csv = "/" + conservation_dir.strip("/") + ".csv"
with open(conservation_csv, "w") as f:
    f.write(csv_out)

print("Filter {} manually for active strands!".format(conservation_csv))

