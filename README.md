# ACCESS INFORMATION

# DATA & CODE FILE OVERVIEW
This is the repository for the manuscript titled "Resource abundance and dietary specialization predict elevational migration in a hyperdiverse montane bird community". This data repository consists of 2 major datasets and a corresponding set of scripts that allow you to recreate the analysis in the manuscript. 

# DATA

## Metabarcoding dataset
Raw data from the DNA metabarcoding of 261 bird faecal samples in the summer and winter across an elevation gradient in the eastern Himalayas. 
There are 2 zipfiles containing compressed fastq files within. This can be found on zenodo

1) **MiSeq_01122023.zip**: 190 files consisting of the MiSeq 150bp paired end reads of 95 samples

2) **NovaSeq_24012025.zip**: 374 files consisting of the NovaSeq6000 150bp paired end reads of 187 samples


## Stable isotope dataset
Raw data from the Stable Isotope analysis (Carbon and Nitrogen) of 258 bird whole blood samples in the summer and winter and 75 plant samples across an elevation gradient in the eastern Himalayas. This can be found on zenodo

1) **bird_iso.csv**: Bird whole blood samples, 3 columns
Samples: Unique ID for the sample which is combination of the unique number on the ringed bird and the date the sample was collected
D13C: Stable isotope ratio of Carbon (δ13C) 
D15N: Stable isotop ration of Nitrogen (δ15N)

2) **plant_iso.csv**: Leaf samples, 3 columns
Sample name: Unique ID for the sample which is combination of the elevation and nth number of samples collected at that elevation
d13C: Stable isotope ratio of Carbon (δ13C) 
d15N: Stable isotop ration of Nitrogen (δ15N)

## Other Datasets
Datasets required to be used in combination with the above two major datasets to replicate the analysis in the manuscript

1) **arthropod_reference.fasta**: Fasta file containing arthropod reference sequences in ecoPCR format 
2) **sequences_arthropoddb.zip**: ecoPCR taxonomy database for arthropods, unzip and place in the same folder as the reference fasta

3) **plant_reference.fasta**: Fasta file containing plant reference sequences in ecoPCR format 
4) **sequences_plantdb.zip**: ecoPCR taxonomy database for plants, unzip and place in the same folder as the reference fasta

The above four datasets are on zenodo

5) **guilds.csv**: dataset describing the migratory status of the birds in the birds sampled. 2 columns
SPECIES: The common name of the bird species
mig: The migratory status of the bird (LE_res = low elevation resident, ME_res = mid elevation resident, HE_res = High elevation resident, mig = Elevational migrant)

6) **ring_numbers.csv**: Data to map the ring numbers and dates back onto the metabarcoding and isotope data. There are 7 columns
INDID: Unique number on the bird ring
SPECIES: The common name of the bird species
DAY: The day the blood/faecal sample was collected
MONTH: The month the blood/faecal sample was collected
YEAR: The year the blood/faecal sample was collected
ELEV: The elevation at which the blood/faecal sample was collected
DATE: A complete date which is a combination of DAY, MONTH and YEAR

7) **samples_date22.csv**: Samples extracted and amplified for the MiSeq run. The sample names in the MiSeq run did not have dates so this dataset matches the ring number to the date the sample was collected. Its has 2 columns
ID: Unique number on the bird ring
DATE: The date the faecal sample was collected

## Intermediate Datasets
These are datasets that are outputs at various stages of the analysis pipeline in case the user wants to skip any step (Some of the steps in the pipeline require significant computational resources which many users may not have access to)

1) **arthropod2023_dadaoutput.fasta**: Denoised and dereplicated arthropod ASVs from the DADA2 pipeline on the MiSeq run (after using cutadapt to remove primer sequences)

2) **arthropod2024_dadaoutput.fasta**: Denoised and dereplicated arthropod ASVs from the DADA2 pipeline on the NovaSeq run (after using cutadapt to remove primer sequences)

3) **plant2023_dadaoutput.fasta**: Denoised and dereplicated plant ASVs from the DADA2 pipeline on the MiSeq run (after using cutadapt to remove primer sequences)

4) **plant2024_dadaoutput.fasta**: Denoised and dereplicated plant ASVs from the DADA2 pipeline on the NovaSeq run (after using cutadapt to remove primer sequences)

5) **arthropod_ASVs_counts_23.txt**: Counts of each arthropod ASV in each sample, output of the DADA2 pipeline on the MiSeq run (after using cutadapt to remove primer sequences)

6) **arthropod_ASVs_counts_24.txt**: Counts of each arthropod ASV in each sample, output of the DADA2 pipeline on the NovaSeq run (after using cutadapt to remove primer sequences)

7) **plant_ASVs_counts_23.txt**: Counts of each plant ASV in each sample, output of the DADA2 pipeline on the MiSeq run (after using cutadapt to remove primer sequences)

8) **plant_ASVs_counts_24.txt**: Counts of each plant ASV in each sample, output of the DADA2 pipeline on the NovaSeq run (after using cutadapt to remove primer sequences)

9) **arthropod_tag_23.tab**: Taxonomically assigned arthropod sequences from the MiSeq run. Output of the ecotag script from obitools

10) **arthropod_tag_24.tab**: Taxonomically assigned arthropod sequences from the NovaSeq run. Output of the ecotag script from obitools

11) **plant_tag_23.tab**: Taxonomically assigned plant sequences from the MiSeq run. Output of the ecotag script from obitools

12) **plant_tag_24.tab**: Taxonomically assigned plant sequences from the NovaSeq run. Output of the ecotag script from obitools

# Code scripts and workflow
We have a separate workflow and a independent set of code for the metabarcoding analysis and the stable isotope analysis. The workflow uses a combination of (`.R`) and (`.py`) scripts

## Metabarcoding Analysis
Step 1. Use "**cut_adapt.py**" to trim primer sequences and heterogeneity spacers

Step 2. Use "**dada2.R**" to  to denoise and dereplicate the trimmed reads, remove chimaeric sequences and create a list of unique ASVs and table of counts of those ASVs in each sample

Step 3. Use "**tax_assign.py**" to assign taxonomy to the list of ASVs output from dada2

Step 4. Use "**ElevMigrationDiet_metabarcoding.R**" to generate the final set of analyses and visualisation in the manuscript and supplementary files

## Stable Isotope Analysis
Use "**ElevMigrationDiet_StableIsotope.R**" to generate the final set of analyses and visualisation in the manuscript and supplementary files

# SOFTWARE VERSIONS
OBItools 1.12: python 2.7
cutadapt 3.6: python 3.12
R-4.3.3
R packages: dada2 3.20, nicheROVER 1.1.2, tRophicPosition 0.8.0, SIBER 2.1.9, tidyverse 2.0.0, vegan 2.6-4, ggpubr 0.6.0, cowplot 1.2.0, , iNEXT 3.0.1, gridExtra 2.3

