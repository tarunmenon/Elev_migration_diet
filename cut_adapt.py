## Code to rempove primer sequences and heterogeneity spacers from the illuminal reads

#Scripts used:
#cutadapt: python 3.6 (https://cutadapt.readthedocs.io/en/v3.6/)

## you need to do this once for the MiSeq and once for the NovaSeq runs
## make sure you are in the directory with the Illumina reads that you have downloaded from figshare and unzipped completely


# Arthropods
for r1 in *_R1_*.fastq; do
    r2=${r1/_R1_/_R2_}
    cutadapt -g ^AGATATTGGAACWTTATATTTTATTTTTGG...GGAGGATTTGGWAATTGATTAGTW  -G ^WACTAATCAATTWCCAAATCCTCC...CCAAAAATAAAATATAAWGTTCCAATATCT --discard-untrimmed --pair-adapters -o arthropod_trimmed_24/trimmed-${r1} -p arthropod_trimmed_24/trimmed-${r2} ${r1} ${r2} 
done

for r1 in *_R1_*.fastq; do
    r2=${r1/_R1_/_R2_}
    cutadapt -g ^AGATATTGGAACWTTATATTTTATTTTTGG...GGAGGATTTGGWAATTGATTAGTW  -G ^WACTAATCAATTWCCAAATCCTCC...CCAAAAATAAAATATAAWGTTCCAATATCT --discard-untrimmed --pair-adapters -o arthropod_trimmed_23/trimmed-${r1} -p arthropod_trimmed_23/trimmed-${r2} ${r1} ${r2} 
done

# Plants
for r1 in *_R1_*.fastq; do
    r2=${r1/_R1_/_R2_}
    cutadapt -g ^GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGG  -G ^CCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --discard-untrimmed --pair-adapters -o plant_trimmed_24/trimmed-${r1} -p plant_trimmed_24/trimmed-${r2} ${r1} ${r2} 
done

for r1 in *_R1_*.fastq; do
    r2=${r1/_R1_/_R2_}
    cutadapt -g ^GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGG  -G ^CCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --discard-untrimmed --pair-adapters -o plant_trimmed_23/trimmed-${r1} -p plant_trimmed_23/trimmed-${r2} ${r1} ${r2} 
done

## Move to DADA2 in R