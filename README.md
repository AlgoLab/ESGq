# ESGq

ESGq is a graph-based approach for the quantification of Alternative Splicing events across two conditions. It works in three steps:
1. simple graphs representing AS events (Event Splicing Graphs) are built from input annotation
2. RNA-Seq reads are aligned to the graphs
3. graph alignments are analyzed to compute 𝜓 and ∆𝜓 values for each event

### Installation
``` sh
git clone https://github.com/AlgoLab/esgq.git
mamba create -c bioconda -n esgq python=3.9 suppa vg biopython gffutils
conda activate esgq
```

### Usage
``` sh
usage: ESGq [-h] [-t THREADS] -1 C1 [C1 ...] -2 C2 [C2 ...] FA GTF WD

required arguments:
  FA              Reference in FASTA format
  GTF             Gene annotation in GTF format
  WD              Output/Working directory
  -1 C1 [C1 ...]  Samples for condition 1 (can be gzipped)
  -2 C2 [C2 ...]  Samples for condition 2 (can be gzipped)

optional arguments:
  -h, --help      show this help message and exit
  -t THREADS      Number of threads to use (default: 1)
```
### Example
``` sh
# Paired-end sample
python3 ESGq.py example/ref.fa example/gene.gtf example/OUT-PE -1 example/A_1.fq,example/A_2.fq example/B_1.fq,example/B_2.fq -2 example/C_1.fq,example/C_2.fq example/D_1.fq,example/D_2.fq
# Single-end sample
python3 ESGq.py example/ref.fa example/gene.gtf example/OUT-SE -1 example/A_1.fq example/B_1.fq -2 example/C_1.fq example/D_1.fq
```
Results are stored in the `events.dpsi` file saved in the output directory. For each event (first column), it reports the 𝜓 value per each replicate and the ∆𝜓 (last column).
