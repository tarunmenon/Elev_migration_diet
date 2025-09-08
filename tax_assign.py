## Using the DADA2 output files we used the ecotag script from the obitools package to assign taxonomy to our sequences 
#OBItools: python 2.7 (https://pythonhosted.org/OBITools/index.html)

# The reference databases were created using a combination existing databases built from BOLD and NCBI and a database we curated using rCRUX (https://github.com/CalCOFI/rCRUX). To work with obitools the reference databases were then brought into the ecopcr format

# Arthropod database: https://github.com/meglecz/mkCOInr
# Plant database: https://ucedna.com/reference-databases-for-metabarcoding

## Make sure you are in the directory that has all the files (dont forget to unzip sequences_arthropoddb)
ecotag -d sequences_arthropoddb  -R arthropod_reference.fasta --sort=count  --uppercase -E 0.5 -M 0.5  -r arthropod_dadaoutput_24.fasta  > arthropod_tag_24.fasta
obitab arthropod_tag_24.fasta > arthropod_tag_24.tab

ecotag -d sequences_arthropoddb  -R arthropod_reference.fasta --sort=count  --uppercase -E 0.5 -M 0.5  -r arthropod_dadaoutput_23.fasta  > arthropod_tag_23.fasta
obitab arthropod_tag_23.fasta > arthropod_tag_23.tab

## Make sure you are in the directory that has all the files (dont forget to unzip sequences_plantdb)
ecotag -d sequences_plantdb  -R plant_reference.fasta --sort=count  --uppercase -E 0.7 -r plant_dadaoutput_24.fasta  > plant_tag_24.fasta
obitab plant_tag_24.fasta > plant_tag_24.tab

ecotag -d sequences_plantdb  -R plant_reference.fasta --sort=count  --uppercase -E 0.7 -r plant_dadaoutput_23.fasta  > plant_tag_23.fasta
obitab plant_tag_23.fasta > plant_tag_23.tab