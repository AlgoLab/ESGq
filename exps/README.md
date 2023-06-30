# experiments

### Data
Reference and annotation from [flybase](https://flybase.org/) (r6.51)
``` sh
wget http://ftp.flybase.net/releases/FB2023_02/dmel_r6.51/fasta/dmel-all-chromosome-r6.51.fasta.gz
wget http://ftp.flybase.net/releases/FB2023_02/dmel_r6.51/gtf/dmel-all-r6.51.gtf.gz
gunzip dmel-all-chromosome-r6.51.fasta.gz
samtools faidx dmel-all-chromosome-r6.51.fasta 2L 2R 3L 3R 4 X Y > ref.chroms.fa
for c in 2L 2R 3L 3R 4 X Y ; do zgrep -P "^$c\t" r6.51/dmel-all-r6.51.gtf.gz ; done > genes.chroms.gtf
wget http://ftp.flybase.net/releases/FB2023_02/dmel_r6.51/fasta/dmel-all-transcript-r6.51.fasta.gz
gunzip -c dmel-all-transcript-r6.51.fasta.gz > transcripts.fa
```

RNA-Seq samples from [https://doi.org/10.1016/j.dib.2021.107413](https://doi.org/10.1016/j.dib.2021.107413)
``` sh
mkdir PRJNA718442
cd PRJNA718442
for sraid in SRR14101759 SRR14101760 SRR14101761 SRR14101762 SRR14101763 SRR14101764
do
    fasterq-dump -p -3 $sraid
    mv ${sraid}_1.fastq ${sraid}_l150_1.fq
    mv ${sraid}_2.fastq ${sraid}_l150_2.fq
done
```

##### Trim samples
``` sh
for sraid in SRR14101759 SRR14101760 SRR14101761 SRR14101762 SRR14101763 SRR14101764
do
    seqtk trimfq -b0 -e50 ${sraid}_l150_1.fq > ${sraid}_l100_1.fq
    seqtk trimfq -b0 -e50 ${sraid}_l150_2.fq > ${sraid}_l100_2.fq
    seqtk trimfq -b0 -e100 ${sraid}_l150_1.fq > ${sraid}_l50_1.fq
    seqtk trimfq -b0 -e100 ${sraid}_l150_2.fq > ${sraid}_l50_2.fq
done 
```

##### Simulate single-end dataset
``` sh
for sraid in SRR14101759 SRR14101760 SRR14101761 SRR14101762 SRR14101763 SRR14101764
do
    cat ${sraid}_l150_1.fq ${sraid}_l150_2.fq > ${sraid}_l150.fq
done
```

### Run the experiments
1. Update `config.yaml` and `config-se.yaml`
2. Run Snakemake
``` sh
# Paired-end
snakemake -c32 --use-conda [-p] [-n] # update config and run this for each length
# Single-end
snakemake -c32 --use-conda [-p] [-n] -s Snakefile-se
```
The main results are stored in a `summary.csv` file saved in the working directory `wd`, along with several correlation plots.

To format the results in a simil-latex table, run:
``` sh
python3 exps/scripts/build_corr_table.py /path/to/PE-l50-wd/summary.csv /path/to/PE-l100-wd/summary.csv /path/to/PE-l150-wd/summary.csv /path/to/SE-l150-wd/summary.csv
```