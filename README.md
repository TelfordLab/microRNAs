# microRNAs
Pipelines to infer conserved microRNAs found in miRBase, as well as detect and predict miRNA candidates from small RNA and genomic data

## filter_miRBase_for_<phylum>_families.py

Requirements:
- [BioPython](https://biopython.org/)
- [Clustal Omega](http://www.clustal.org/omega/)
- [MView](https://desmid.github.io/mview/)

Use:
`python3 filter_miRBase_for_<phylum>_families.py <miRBase_mature.fa> <miRBase_hairpin.fa>`

Input:
- [miRBase mature miRNA sequences](ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz)
- [miRBase hairpin sequences](ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz)

Output:
- `conservation_<phylum>.csv` - CSV containing information about inferred miRNA families
- `conservation_<phylum>/` - directory of miRNA family sequences and alignments

## microRNA_detection.py

Requirements:
- [BioPython](https://biopython.org/)
- [RNAfold](https://www.tbi.univie.ac.at/RNA/)

Use:
`microRNA_detection.py -s <small_RNA_reads.fa> -g <genome.fa> -c <conservation.csv>`

Input:
- small RNA reads (`-s <small_RNA_reads.fa>`, this is *optional* when using pre-computed steps, e.g. miRNA predictions, see below)
- genome as nucleotide scaffolds (`-g <genome.fa>`)
- conservation information about queried miRNA families (`-c <conservation.csv>`)
- *optional*: apply fixed threshold for all mature miRNA candidates (`-t <threshold>`)
- *optional*: output directory to resume previously computed steps (`-o <output_dir>`)

Output:
- mature miRNA candidates for each miRNA family identified from small RNA reads (`<miRNA-family>.mature`)
- pre-miRNA candidates for each miRNA family identified from genome (`<miRNA-family>.pre_mirnas`)
- best graded pre-miRNA candidates after evaluation of RNA structure (`<miRNA-family>.best_graded`)

## split_genome_into_miRNA_kmers.py

Requirements:
- [BioPython](https://biopython.org/)

Use:
`python3 split_genome_into_miRNA_kmers.py -g <genome.fa> -c <conservation.csv>`

Input:
- genome as nucleotide scaffolds (`-g <genome.fa>`)
- seed information about miRNA families (`-c <conservation.csv>`)

Output:
- directory containing location and position of all miRNA family seeds within the genome (`<seed>.pos`)

## microRNA_prediction_from_DNA_kmers.py

Requirements:
- information about seed locations within genome (see above `split_genome_into_miRNA_kmers.py`)
- [BioPython](https://biopython.org/)

Use:
`python3 microRNA_prediction_from_DNA_kmers.py -g <genome.fa> -c <conservation.csv> -kd <kmers_directory>`

Input:
- genome as nucleotide scaffolds (`-g <genome.fa>`)
- seed information about miRNA families (`-c <conservation.csv>`)
- directory containing miRNA family seed information (`-kd <kmers_directory>`, see above)

Output:
- mature miRNA candidates for each miRNA family predicted from the genome (`<miRNA-family>.mature`)
- pre-miRNA candidates for each miRNA family predicted from genome (`<miRNA-family>.pre_mirnas`)