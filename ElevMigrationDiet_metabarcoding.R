## These scripts are to be run on the output of the dada2 pipeline (per sample ASV counts) and the ecotag script (Taxonomically annotated ASVs) to replicate the faecal DNA metabarcoding related visualisation and analysis in the manuscript

setwd("~/Elev_migration_diet")

library(tidyverse)
library(vegan)
library(ggpubr)
library(gridExtra)
library(iNEXT)
`%ni%` <- Negate(`%in%`)
cbbPalette <- c("#D55E00", "#56B4E9","#000000", "#CC79A7","#009E73", "#E69F00", "#4B0092", "#0072B2","#DC3220" )

## Arthropod reads

## Read in the processed data from the NovaSeq runs  
art_asv_count24<-read.delim("arthropod_ASVs_counts_24.txt")
colnames(art_asv_count24) <- gsub("^[^.]+\\.", "", colnames(art_asv_count24))
#There were 2 samples with 0 reads so setting those to 0 because dada2 couldn't handle zeros so I had to remove these samples at that stage
art_asv_count24$A5191_7124<-0
art_asv_count24$A5200_7124<-0
art_asv_tag24<-read.delim("arthropod_tag_24.tab")

art_asv_count24$id<-row.names(art_asv_count24) # adding rownames

## Read in the processed data from the MiSeq runs
art_asv_count23<-read.delim("arthropod_ASVs_counts_23.txt")
art_asv_tag23<-read.delim("arthropod_tag_23.tab")

art_asv_count23$id<-row.names(art_asv_count23)

# merge the per sample ASV counts with taxonomically annotated ASVs
art_asv_count23<-left_join(art_asv_count23, art_asv_tag23) %>% filter(!is.na(order_name))%>% select(!c(id, match_count.output))
art_asv_count24<-left_join(art_asv_count24, art_asv_tag24) %>% filter(!is.na(order_name)) %>% select(!c(id, match_count.output))

# merge the NovaSeq and MiSeq datasets and add the numbers of sequences in common
art_asv_count <- full_join(art_asv_count23,art_asv_count24) %>% mutate_if(is.numeric,coalesce,0)

# remove ASVs with less than 90% match and where orders hasn't been identified
art_asv_count<-art_asv_count %>% filter(best_identity.output >= 0.90) %>% 
  filter(!is.na(order_name))

# read counts lower or equal to the count of an ASV appearing in the negative controls were replaced with zero
art_asv_count[,c(1:95,117:300)]<-art_asv_count[,c(1:95,117:300)]-art_asv_count$EX.CTRL-art_asv_count$PCR.CTRL-art_asv_count$EX_CTRL1-art_asv_count$EX_CTRL2-art_asv_count$PCR_CTRL
art_asv_count[art_asv_count < 0] <- 0

# removing ASVs with less than 20 total reads and ASVS
art_asv_count$reads<-rowSums(art_asv_count[,c(1:95,117:300)])
art_asv_count<-art_asv_count %>% filter(art_asv_count$reads > 20)
art_asv_count$counts<-rowSums(art_asv_count[,c(1:95,117:300)] != 0) #number of samples with each ASV

# ASVs with less than 10 reads in a sample were considered to be absent
art_asv_count <- art_asv_count %>%
  mutate(across(c(1:95,117:300), ~ ifelse(. < 10, 0, .)))

art_asv_count$ASV<-rownames(art_asv_count) #adding a column with ASV number

# Creating a long table of all possible samples with counts of all possible ASVs 
individual_columns <- names(art_asv_count)[c(1:95,117:300)] #vector of column names
art_non_zero <- art_asv_count %>% pivot_longer(cols = all_of(individual_columns), names_to = "ID", values_to = "Count") %>%
  mutate(Count = replace_na(Count, 0))

# Create the summary table of counts of each arthropod order in each sample
# Fix typos in sample names
result_df_art <- art_non_zero %>%
  group_by(ID, order_name) %>%
  summarise(Count = sum(Count), .groups = 'drop') %>%
  pivot_wider(names_from = order_name, values_from = Count, values_fill = 0) %>% mutate(ID= recode(ID, "Z2662"="A2662", "A2584" = "Z2584", "A1297"="AB1297", "A2591"="Z2591", "A2594.25422" = "A2690", "A2593"="Z2593",  "A2598"="Z2598",  "A2599"="Z2599", "Z2388" = "A2883", "A1739_10124"="AB1739_10124", "A8224_19124" = "A8224_14124", "A8233_16124"="A8233_15124", "A8317_1224"="Z8317_1224","AB1179_16124" = "AB1179_161223",  "Z5093_19124"="Z5092_19124", "Z5019_22124" = "Z2711_23522",  "A3700_22323"="A3700_22523","A3841_22523" =  "Z3841_22523", "A3865_23523" = "Z3865_23523", "A5119_20124" = "Z5119_20124","Z5152_1224" = "Z2513_30422", "AB0525_18523" = "AB0525_19523", "Z3632_6523" = "A3632_6523", "A3606_18523"="Z3806_18523", "A269807_25422" = "A269087_25422", "AB0373_26422" = "AB0373_26423", "Z2604_5622" = "Z2604_5522", "Z2605_6522"= "Z2605_5522", "Z2612_5622"= "Z2612_5522", "Z3428_26422" = "Z3428_26423", "A3110_141222" = "A3110_151222", "A3165_7123" = "A3165_8123", "A5004_121223" = "A5104_121223","A5200_7124" = "A5200_8124"))

# Read in the data that maps the sample names with information on species, elevation and date of sample collection
bird_data<-read.csv("E:/IISC/Global_Change_Lab/Data/ring_numbers.csv")

# Data on migratory status
guilds<-read.csv("guilds.csv")

# The sample names in the MiSeq run did not have dates so this dataset matches the ring number to the date the sample was collected.
samples_date22<-read.csv("samples_date22.csv")

## Merging the all the datasets
bird_data$ID<- paste0(bird_data$INDID, "_",sub("^0", "", sub("^(\\d+)-0?(\\d+)-(\\d{2})(\\d{2})$", "\\1\\2\\4", bird_data$DATE)))
bird_data$DAY<-as.numeric(bird_data$DAY)
bird_data$MONTH<-as.numeric(bird_data$MONTH)
bird_data$season<-ifelse(bird_data$MONTH == 4 | bird_data$MONTH == 5, "Summer", "Winter") #Defining seasons

art_bird<- left_join(result_df_art,samples_date22) %>% mutate(ID = str_remove(ID, "\\..*"))
art_bird$ID <- ifelse(is.na(art_bird$DATE), art_bird$ID, paste0(art_bird$ID, "_",sub("^0", "", sub("^(\\d+)-0?(\\d+)-(\\d{2})(\\d{2})$", "\\1\\2\\4", art_bird$DATE))))
art_bird<-art_bird%>% left_join(bird_data, by = "ID") %>% left_join(guilds)

# Creating a long table of all possible samples with counts of all possible Orders 
z<-art_bird %>% pivot_longer(cols = names(art_bird[2:25]),names_to = "arthropod_order", values_to = "count") %>% drop_na(SPECIES) %>% filter(mig == "HE_res" | mig == "mig")

# Creating a summary table of the frequency of occurrence of each arthropod order in each season at the species level
art_bar_sps <- z %>%
  group_by(season, mig, SPECIES, arthropod_order) %>%
  summarise(reads = sum(count), detections = sum(count > 0), total = n(), .groups = 'drop') %>%
  mutate(proportion = detections/total) %>%
  select(SPECIES, mig, season, proportion, arthropod_order, total, detections, reads)

# Creating a summary table of the frequency of occurrence of each arthropod order in each season at the migratory strategy level
art_bar <- z %>%
  group_by(season, mig, arthropod_order) %>%
  summarise(reads = sum(count), detections = sum(count > 0), total = n(), .groups = 'drop') %>%
  mutate(proportion = detections/total) %>%
  select(season, mig, proportion, arthropod_order, total, detections, reads)

# Plot for Figure 4
art_all<-art_bar %>% filter(arthropod_order %in% c("Araneae", "Coleoptera", "Lepidoptera", "Hemiptera", "Diptera")) %>%
  mutate(mig = recode(mig, "HE_res" = "Resident"))%>%
  mutate(mig = recode(mig, "mig" = "Migrant")) %>% ggplot(aes(x = mig, y = proportion, fill = season)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~arthropod_order ) +
  labs(x = "Migratory Status",
       y = "Frequency of Occurrence",
       fill = "Season") +
  theme_bw(base_size = 15)+
  scale_fill_manual(values = cbbPalette)+
  theme(legend.position = "inside",
        legend.justification = c(0.9,0.2),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

#ggsave(plot = art_all, "Fig4.jpeg", width = 8, height = 5)


## Supplementary plot

art_res<- art_bar_sps %>% filter(arthropod_order %in% c("Araneae", "Coleoptera", "Lepidoptera", "Hemiptera", "Diptera"))%>% 
  filter(mig == "HE_res") %>%
  filter(SPECIES %in% c("Brown-throated Fulvetta", "Rufous-vented Yuhina", "Streak-breasted Scimitar-Babbler", "Stripe-throated Yuhina", 'Rufous-capped Babbler'))%>%
  ggplot(aes(x = arthropod_order, y = proportion, fill = season)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = cbbPalette)+
  facet_wrap(~SPECIES, scales = "free_x", ncol = 1) +
  labs(title = "Residents",
       x = "Orders",
       y = "Frequency of Occurence",
       fill = "Season") +
  #theme_minimal()+
  theme_bw(base_size = 12)+
  theme(plot.title = element_text(size = 13, hjust = 0.5))+
  theme(legend.position="bottom", legend.text = element_text(size=12), legend.title = element_text(size=13))+
  theme(strip.text.x = element_text(size = 12))

art_mig<- art_bar_sps %>% filter(arthropod_order %in% c("Araneae", "Coleoptera", "Lepidoptera", "Hemiptera", "Diptera"))%>%
  filter(mig == "mig") %>%
  filter(SPECIES %in% c("Whistler's Warbler", "Rufous-gorgeted Flycatcher", "Chestnut-headed Tesia", "Rufous-winged Fulvetta", "Rufous-bellied Niltava", "Buff-barred Warbler"))%>%
  ggplot(aes(x = arthropod_order, y = proportion, fill = season)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = cbbPalette)+
  facet_wrap(~SPECIES, scales = "free_x", ncol = 1) +
  labs(title = "Migrants",
       x = "Orders",
       y = "Frequency of Occurence",
       fill = "Season") +
  #theme_minimal()+
  theme_bw(base_size = 12)+
  theme(plot.title = element_text(size = 13, hjust = 0.5))+
  theme(legend.position="bottom", legend.text = element_text(size=12), legend.title = element_text(size=13))+
  theme(strip.text.x = element_text(size = 12))

art_species<-ggarrange(art_res, art_mig, ncol = 2, common.legend = T, legend = "bottom")

#ggsave("FigS1.jpeg", art_species, width = 9, height = 13)



#Proportion tests for each of the 5 orders
art_bar$nondetections<-art_bar$total-art_bar$detections

Lep_res<-prop.test(art_bar %>% filter(arthropod_order == "Lepidoptera" & mig == "HE_res") %>% select(detections,nondetections) %>% as.matrix())
Lep_mig<-prop.test(art_bar %>% filter(arthropod_order == "Lepidoptera" & mig == "mig") %>% select(detections,nondetections) %>% as.matrix())

Hem_res<-prop.test(art_bar %>% filter(arthropod_order == "Hemiptera" & mig == "HE_res") %>% select(detections,nondetections) %>% as.matrix())
Hem_mig<-prop.test(art_bar %>% filter(arthropod_order == "Hemiptera" & mig == "mig") %>% select(detections,nondetections) %>% as.matrix())

Col_res<-prop.test(art_bar %>% filter(arthropod_order == "Coleoptera" & mig == "HE_res") %>% select(detections,nondetections) %>% as.matrix())
Col_mig<-prop.test(art_bar %>% filter(arthropod_order == "Coleoptera" & mig == "mig") %>% select(detections,nondetections) %>% as.matrix())

Ara_res<-prop.test(art_bar %>% filter(arthropod_order == "Araneae" & mig == "HE_res") %>% select(detections,nondetections) %>% as.matrix())
Ara_mig<-prop.test(art_bar %>% filter(arthropod_order == "Araneae" & mig == "mig") %>% select(detections,nondetections) %>% as.matrix())

Dip_res<-prop.test(art_bar %>% filter(arthropod_order == "Diptera" & mig == "HE_res") %>% select(detections,nondetections) %>% as.matrix())
Dip_mig<-prop.test(art_bar %>% filter(arthropod_order == "Diptera" & mig == "mig") %>% select(detections,nondetections) %>% as.matrix())


####### Plant Reads

##### Plant Family

## Read in the processed data from the NovaSeq runs 
pln_asv_count24<-read.delim("plant_ASVs_counts_24.txt")
colnames(pln_asv_count24) <- gsub("^[^.]+\\.", "", colnames(pln_asv_count24))
pln_asv_tag24<-read.delim("plant_tag_24.tab")

pln_asv_count24$id<-row.names(pln_asv_count24)

## Read in the processed data from the MiSeq runs
pln_asv_count23<-read.delim("plant_ASVs_counts_23.txt") %>% rename_with(~ str_remove(., "trimmed\\."))
pln_asv_tag23<-read.delim("plant_tag_23.tab")

pln_asv_count23$id<-row.names(pln_asv_count23)

# merge the per sample ASV counts with taxonomically annotated ASVs
pln_count23_family<-left_join(pln_asv_count23, pln_asv_tag23) %>% filter(!is.na(family_name))%>% select(!c(id, match_count.plant3_unique))
pln_count24_family<-left_join(pln_asv_count24, pln_asv_tag24) %>% filter(!is.na(family_name)) %>% select(!c(id, match_count.plant3_unique))

# merge the NovaSeq and MiSeq datasets and add the numbers of sequences in common
pln_count_family <- full_join(pln_count23_family,pln_count24_family, by = c("sequence", "family_name", "best_identity.plant3_unique")) %>% mutate_if(is.numeric,coalesce,0)

# remove ASVs with less than 90% match and where families hasn't been identified
# One ASV with a very high number of reads had 80% match, we confirmed the sequence manually on NCBI so that it could be included in the analysis
pln_count_family<-pln_count_family %>% filter(best_identity.plant3_unique > 0.90|sequence == "ctcctgctttccaaaagtaagcataaaaaaag") %>% filter(!is.na(family_name))

# read counts lower or equal to the count of an ASV appearing in the negative controls were replaced with zero
pln_count_family[,c(1:95,117:300)]<-pln_count_family[,c(1:95,117:300)]-pln_count_family$EX.CTRL-pln_count_family$PCR.CTRL-pln_count_family$EX_CTRL1-pln_count_family$EX_CTRL2-pln_count_family$PCR_CTRL
pln_count_family[pln_count_family < 0] <- 0

# removing ASVs with less than 20 total reads and ASVS
pln_count_family$reads<-rowSums(pln_count_family[,c(1:95,117:300)])
pln_count_family<-pln_count_family %>% filter(pln_count_family$reads > 20)
pln_count_family$counts<-rowSums(pln_count_family[,c(1:95,117:300)] != 0)

# Samples in which an ASV had <5% of total reads was replaced with 0
# Calculate the threshold for each column (5% of the column sum)
column_sums <- colSums(pln_count_family[,c(1:95,117:300)])
thresholds <- 0.05 * column_sums

# Apply the threshold to each column
# Remove families which are non-fruiting/flowering and not found in the region
pln_count_family <- pln_count_family %>%
  mutate(across(c(1:95,117:300), ~ ifelse(. < thresholds[cur_column()], 0, .))) %>% filter(family_name %ni% c("Dryopteridaceae", "Thelypteridaceae", "Rhamnaceae", "Pteridaceae", "Malphigiaceae", "Hydrangaceae", "Gesneriaceae", "Geraniaceae", "Amaryllidaceae", "Amphorogynaceae", "Atherospermataceae", "Casuarinaceae", "Dioscoreaceae", "Heliconiaceae", "Hymenophyllaceae", "Orobanchaceae", "Canellaceae"))

pln_count_family$ASV<-rownames(pln_count_family) #adding a column with ASV number


# Creating a long table of all possible samples with counts of all possible ASVs
individual_columns <- names(pln_count_family)[c(1:95,117:300)]
plnF_non_zero <- pln_count_family %>% pivot_longer(cols = all_of(individual_columns), names_to = "ID", values_to = "Count") %>%
  mutate(Count = replace_na(Count, 0))

# Create the summary table of counts of each plant family in each sample
# Fix typos in sample names
result_pln_family <- plnF_non_zero %>%
  group_by(ID, family_name) %>%
  summarise(Count = sum(Count), .groups = 'drop') %>%
  pivot_wider(names_from = family_name, values_from = Count, values_fill = 0) %>% mutate(ID= recode(ID, "Z2662"="A2662", "A2584" = "Z2584", "A1297"="AB1297", "A2591"="Z2591", "A2594.25422" = "A2690", "A2593"="Z2593",  "A2598"="Z2598",  "A2599"="Z2599", "Z2388" = "A2883", "A1739_10124"="AB1739_10124", "A8224_19124" = "A8224_14124", "A8233_16124"="A8233_15124", "A8317_1224"="Z8317_1224","AB1179_16124" = "AB1179_161223",  "Z5093_19124"="Z5092_19124", "Z5019_22124" = "Z2711_23522",  "A3700_22323"="A3700_22523","A3841_22523" =  "Z3841_22523", "A3865_23523" = "Z3865_23523", "A5119_20124" = "Z5119_20124","Z5152_1224" = "Z2513_30422", "AB0525_18523" = "AB0525_19523", "Z3632_6523" = "A3632_6523", "A3606_18523"="Z3806_18523", "A269807_25422" = "A269087_25422", "AB0373_26422" = "AB0373_26423", "Z2604_5622" = "Z2604_5522", "Z2605_6522"= "Z2605_5522", "Z2612_5622"= "Z2612_5522", "Z3428_26422" = "Z3428_26423", "A3110_141222" = "A3110_151222", "A3165_7123" = "A3165_8123", "A5004_121223" = "A5104_121223","A5200_7124" = "A5200_8124"))

## Merging the all the datasets
pln_bird_family<- left_join(result_pln_family,samples_date22) %>% mutate(ID = str_remove(ID, "\\..*"))
pln_bird_family$ID <- ifelse(is.na(pln_bird_family$DATE), pln_bird_family$ID, paste0(pln_bird_family$ID, "_",sub("^0", "", sub("^(\\d+)-0?(\\d+)-(\\d{2})(\\d{2})$", "\\1\\2\\4", pln_bird_family$DATE))))
pln_bird_family<-pln_bird_family%>% left_join(bird_data, by = "ID") %>% left_join(guilds)

# Creating a long table of all possible samples with counts of all possible families
y<-pln_bird_family %>% pivot_longer(cols = names(pln_bird_family[2:90]),names_to = "plant_family", values_to = "count") %>% drop_na(SPECIES) %>% filter(mig == "HE_res" | mig == "mig")

# Creating a summary table of the frequency of occurrence of each plant family in each season at the species level
pln_bar_sps <- y %>%
  group_by(season, mig, SPECIES, plant_family) %>%
  summarise(reads= sum(count),detections = sum(count > 0), total = n(), .groups = 'drop') %>%
  mutate(proportion = detections/total) %>%
  select(SPECIES,season, mig, proportion, plant_family, total, detections, reads)

# Creating a summary table of the frequency of occurrence of each plant family in each season at the migratory strategy level
pln_bar <- y %>%
  group_by(season, mig, plant_family) %>%
  summarise(reads= sum(count),detections = sum(count > 0), total = n(), .groups = 'drop') %>%
  mutate(proportion = detections/total) %>%
  select(season, mig, proportion, plant_family, total, detections, reads)

# Subsetting important families (Supplementary table 5)
c<-pln_bar %>% filter(plant_family %in% c("Polygonaceae", "Ericaceae", "Rosaceae", "Caprifoliaceae", "Acanthaceae", "Malvaceae", "Pinaceae", "Poaceae", "Fagaceae" ))

# Proportion test
pln_bar$nondetections<-pln_bar$total-pln_bar$detections

Pol_res<-prop.test(pln_bar %>% filter(plant_family == "Polygonaceae" & mig == "HE_res") %>% select(detections,nondetections) %>% as.matrix())
Pol_mig<-prop.test(pln_bar %>% filter(plant_family == "Polygonaceae" & mig == "mig") %>% select(detections,nondetections) %>% as.matrix())


###### Plant orders
## Analysis of plants also done at the order level for the community level analysis

# merge the per sample ASV counts with taxonomically annotated ASVs
pln_count23_order<-left_join(pln_asv_count23, pln_asv_tag23) %>% filter(!is.na(order_name))%>% select(!c(id, match_count.plant3_unique))
pln_count24_order<-left_join(pln_asv_count24, pln_asv_tag24) %>% filter(!is.na(order_name)) %>% select(!c(id, match_count.plant3_unique))

# merge the NovaSeq and MiSeq datasets and add the numbers of sequences in common
pln_count_order <- full_join(pln_count23_order,pln_count24_order, by = c("sequence", "order_name", "best_identity.plant3_unique")) %>% mutate_if(is.numeric,coalesce,0)

# remove ASVs with less than 90% match and where order hasn't been identified
# One ASV with a very high number of reads had 80% match, we confirmed the sequence manually on NCBI so that it could be included in the analysis
pln_count_order<-pln_count_order %>% filter(best_identity.plant3_unique > 0.90|sequence == "ctcctgctttccaaaagtaagcataaaaaaag") %>% filter(!is.na(order_name))

# read counts lower or equal to the count of an ASV appearing in the negative controls were replaced with zero
pln_count_order[,c(1:95,117:300)]<-pln_count_order[,c(1:95,117:300)]-pln_count_order$EX.CTRL-pln_count_order$PCR.CTRL-pln_count_order$EX_CTRL1-pln_count_order$EX_CTRL2-pln_count_order$PCR_CTRL
pln_count_order[pln_count_order < 0] <- 0

# removing ASVs with less than 20 total reads and ASVS
pln_count_order$reads<-rowSums(pln_count_order[,c(1:95,117:300)])
pln_count_order<-pln_count_order %>% filter(pln_count_order$reads > 20)
pln_count_order$counts<-rowSums(pln_count_order[,c(1:95,117:300)] != 0)

# Samples in which an ASV had <5% of total reads was replaced with 0
# Calculate the threshold for each column (5% of the column sum)
column_sums <- colSums(pln_count_order[,c(1:95,117:300)])
thresholds <- 0.05 * column_sums

# Apply the threshold to each column
# Remove families which are non-fruiting/flowering and not found in the region
pln_count_order <- pln_count_order %>%
  mutate(across(c(1:95,117:300), ~ ifelse(. < thresholds[cur_column()], 0, .))) %>% filter(family_name.x %ni% c("Dryopteridaceae", "Thelypteridaceae", "Rhamnaceae", "Pteridaceae", "Malphigiaceae", "Hydrangaceae", "Gesneriaceae", "Geraniaceae", "Amaryllidaceae", "Amphorogynaceae", "Atherospermataceae", "Casuarinaceae", "Dioscoreaceae", "Heliconiaceae", "Hymenophyllaceae", "Orobanchaceae", "Canellaceae")) %>% filter(family_name.y %ni% c("Dryopteridaceae", "Thelypteridaceae", "Rhamnaceae", "Pteridaceae", "Malphigiaceae", "Hydrangaceae", "Gesneriaceae", "Geraniaceae", "Amaryllidaceae", "Amphorogynaceae", "Atherospermataceae", "Casuarinaceae", "Dioscoreaceae", "Heliconiaceae", "Hymenophyllaceae", "Orobanchaceae", "Canellaceae"))

individual_columns <- names(pln_count_order)[c(1:95,117:300)] #adding a column with ASV number

# Creating a long table of all possible samples with counts of all possible ASVs
plnO_non_zero <- pln_count_order %>% pivot_longer(cols = all_of(individual_columns), names_to = "ID", values_to = "Count") %>%
  mutate(Count = replace_na(Count, 0))

# Create the summary table of counts of each plant order in each sample
# Fix typos in sample names
result_pln_order <- plnO_non_zero %>%
  group_by(ID, order_name) %>%
  summarise(Count = sum(Count), .groups = 'drop') %>%
  pivot_wider(names_from = order_name, values_from = Count, values_fill = 0) %>% mutate(ID= recode(ID, "Z2662"="A2662", "A2584" = "Z2584", "A1297"="AB1297", "A2591"="Z2591", "A2594.25422" = "A2690", "A2593"="Z2593",  "A2598"="Z2598",  "A2599"="Z2599", "Z2388" = "A2883", "A1739_10124"="AB1739_10124", "A8224_19124" = "A8224_14124", "A8233_16124"="A8233_15124", "A8317_1224"="Z8317_1224","AB1179_16124" = "AB1179_161223",  "Z5093_19124"="Z5092_19124", "Z5019_22124" = "Z2711_23522",  "A3700_22323"="A3700_22523","A3841_22523" =  "Z3841_22523", "A3865_23523" = "Z3865_23523", "A5119_20124" = "Z5119_20124","Z5152_1224" = "Z2513_30422", "AB0525_18523" = "AB0525_19523", "Z3632_6523" = "A3632_6523", "A3606_18523"="Z3806_18523", "A269807_25422" = "A269087_25422", "AB0373_26422" = "AB0373_26423", "Z2604_5622" = "Z2604_5522", "Z2605_6522"= "Z2605_5522", "Z2612_5622"= "Z2612_5522", "Z3428_26422" = "Z3428_26423", "A3110_141222" = "A3110_151222", "A3165_7123" = "A3165_8123", "A5004_121223" = "A5104_121223","A5200_7124" = "A5200_8124"))

## Merging the all the datasets
pln_bird_order<- left_join(result_pln_order,samples_date22) %>% mutate(ID = str_remove(ID, "\\..*"))
pln_bird_order$ID <- ifelse(is.na(pln_bird_order$DATE), pln_bird_order$ID, paste0(pln_bird_order$ID, "_",sub("^0", "", sub("^(\\d+)-0?(\\d+)-(\\d{2})(\\d{2})$", "\\1\\2\\4", pln_bird_order$DATE))))
pln_bird_order<-pln_bird_order%>% left_join(bird_data, by = "ID") %>% left_join(guilds)


####### Community analysis, Ordination, Permanova, Permdisp
## Creating a dataframe of all orders (plant and arthropod) in present in a sample
bird_diet<-left_join(art_bird,pln_bird_order, by = c("INDID","ID","SPECIES","mig","season","ELEV"))

## Keep only high elevation residents and elevational migrants
bird_diet<-bird_diet %>% drop_na(SPECIES) %>% filter(mig == "HE_res" | mig == "mig") 

## Create 2 more columns which are a combination of migratory strategy and season and species and season
bird_diet$migsea<-paste0(bird_diet$mig,"_",bird_diet$season)
bird_diet$specsea<-paste0(bird_diet$SPECIES,"_",bird_diet$season)

## Keep only species with sufficient data
bird_diet<-bird_diet%>% filter(SPECIES %in% c("Brown-throated Fulvetta", "Rufous-vented Yuhina", "Streak-breasted Scimitar-Babbler", "Stripe-throated Yuhina", "Rufous-capped Babbler", "Whistler's Warbler", "Rufous-gorgeted Flycatcher", "Chestnut-headed Tesia", "Rufous-winged Fulvetta", "Rufous-bellied Niltava", "Buff-barred Warbler"))

## Create dissimmilarity matrix with jaccard index
x.dist<-vegdist(bird_diet[,-c(1,26:35,76:82)], method = "jaccard", binary = T)
x.dist<-as.matrix(x.dist)

## PCoA plot
pcoa <- cmdscale(x.dist, eig=TRUE)
pcoa_PLOT<- data.frame(pcoa$points, mig = as.factor(bird_diet$mig), migsea = as.factor(bird_diet$migsea), season = as.factor(bird_diet$season)) %>% mutate(migsea = recode(migsea,  "HE_res_Summer" = "Resident_Summer"))%>% mutate(migsea = recode(migsea, "mig_Summer" = "Migrant_Summer"))  %>% mutate(migsea = recode(migsea, "HE_res_Winter" = "Resident_Winter"))%>% mutate(migsea = recode(migsea, "mig_Winter" = "Migrant_Winter"))

## Ellipse parameters
plot(pcoa$points)
ord_se<-ordiellipse(pcoa, pcoa_PLOT$migsea, display = "sites", 
                    kind = "sd", label = T)

## Function for creating ellipses in ggplot
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame()

for(g in levels(pcoa_PLOT$migsea)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(pcoa_PLOT[pcoa_PLOT$migsea==g,],
                                                   veganCovEllipse(ord_se[[g]]$cov,ord_se[[g]]$center,ord_se[[g]]$scale)))
                                ,migsea=g))
}

df_ell$migsea<-as.factor(df_ell$migsea) 

## Plotting the PCoA - Figure 5
PCoA_all<-pcoa_PLOT %>% ggplot(aes(X1, X2)) + geom_point(aes(colour = migsea)) + 
  geom_polygon(data = df_ell, aes(x=Dim1, y=Dim2, colour = migsea), fill = NA , linewidth=1)+
  scale_color_manual(values = cbbPalette)+
  xlab("PCoA 1") + ylab("PCoA2")+
  theme_bw()+
  guides(color = guide_legend(title = "Migratory Behaviour \nand Season"), shape = "none", linetype = "none")

#ggsave(plot = PCoA_all, "Fig5.jpeg", width = 7, height = 4)


## PCoA at the species level
# Migratory species
mig_diet<-bird_diet %>% filter(SPECIES %in% c("Whistler's Warbler", "Rufous-gorgeted Flycatcher", "Chestnut-headed Tesia", "Rufous-winged Fulvetta", "Rufous-bellied Niltava", "Buff-barred Warbler"))

mig.dist<-vegdist(mig_diet[,-c(1,26:35,76:82)], method = "jaccard", binary = T)
mig.dist<-as.matrix(mig.dist)
pcoa_mig <- cmdscale(mig.dist, eig=TRUE)
pcoa_migPLOT<- data.frame(pcoa_mig$points, sp = as.factor(mig_diet$SPECIES), season = as.factor(mig_diet$season), specsea = as.factor(mig_diet$specsea))

plot(pcoa_mig$points)
ord_se<-ordiellipse(pcoa_mig, pcoa_migPLOT$specsea, display = "sites", 
                    kind = "sd", label = T)

df_ell_mig <- data.frame()

for(g in levels(pcoa_migPLOT$specsea)){
  df_ell_mig <- rbind(df_ell_mig, cbind(as.data.frame(with(pcoa_migPLOT[pcoa_migPLOT$specsea==g,],
                                                           veganCovEllipse(ord_se[[g]]$cov,ord_se[[g]]$center,ord_se[[g]]$scale)))
                                        ,specsea=g))
}

df_ell_mig <-
  df_ell_mig %>%
  separate(col = specsea, 
           into = c('sp','season'),
           sep = '_',
           remove = F)

df_ell_mig$sp<-as.factor(df_ell_mig$sp)

## PCoA plot of migratory species
migplot<-pcoa_migPLOT %>% ggplot(aes(X1, X2)) + geom_point(aes(colour = sp, shape = season)) + 
  geom_polygon(data = df_ell_mig, aes(x=Dim1, y=Dim2, colour = sp, linetype = season), fill = NA , linewidth=1)+
  scale_color_manual(values = cbbPalette)+xlim(-0.5,0.5)+
  xlab("PCoA 1") + ylab("PCoA2")+
  guides(color = guide_legend(title = "Migratory Species"), shape = "none", linetype = "none")

## PCoA at the species level
# Resident species
res_diet<-bird_diet %>% filter(SPECIES %in% c("Brown-throated Fulvetta", "Rufous-vented Yuhina", "Streak-breasted Scimitar-Babbler", "Stripe-throated Yuhina", "Rufous-capped Babbler"))

res.dist<-vegdist(res_diet[,-c(1,26:35,76:82)], method = "jaccard", binary = T)
res.dist<-as.matrix(res.dist)
pcoa_res <- cmdscale(res.dist, eig=TRUE)
pcoa_resPLOT<- data.frame(pcoa_res$points, sp = as.factor(res_diet$SPECIES), season = as.factor(res_diet$season), specsea = as.factor(res_diet$specsea))

plot(pcoa_res$points)
ord_se<-ordiellipse(pcoa_res, pcoa_resPLOT$specsea, display = "sites", 
                    kind = "sd", label = T)

df_ell_res <- data.frame()

for(g in levels(pcoa_resPLOT$specsea)){
  df_ell_res <- rbind(df_ell_res, cbind(as.data.frame(with(pcoa_resPLOT[pcoa_resPLOT$specsea==g,],
                                                           veganCovEllipse(ord_se[[g]]$cov,ord_se[[g]]$center,ord_se[[g]]$scale)))
                                        ,specsea=g))
}

df_ell_res <-
  df_ell_res %>%
  separate(col = specsea, 
           into = c('sp','season'),
           sep = '_',
           remove = F)

df_ell_res$sp<-as.factor(df_ell_res$sp)

## PCoA plot of resident species
resplot<-pcoa_resPLOT %>% ggplot(aes(X1, X2)) + geom_point(aes(colour = sp, shape = season)) + 
  geom_polygon(data = df_ell_res, aes(x=Dim1, y=Dim2, colour = sp, linetype = season), fill = NA , linewidth=1)+
  scale_color_manual(values = cbbPalette) +xlim(-0.5,0.5)+
  xlab("PCoA 1") + ylab("PCoA2")+
  guides(color = guide_legend(title = "Resident Species"), shape = "none", linetype = "none")

## Plotting the species level PCoA - Supplementary figure 3
PCoA_sp<-ggarrange(resplot, migplot, nrow = 2)

#ggsave(plot = PCoA_sp, "FigS3.jpeg", width = 9, height = 6)

#########

#PERMANOVA for each species

permres<-adonis2(bird_diet %>% filter(mig == "HE_res") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                   bird_diet %>% filter(mig == "HE_res") %>% .$season , method = "jaccard")

permmig<-adonis2(bird_diet %>% filter(mig == "mig") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                   bird_diet %>% filter(mig == "mig") %>% .$season , method = "jaccard")

permBTFU<-adonis2(bird_diet %>% filter(SPECIES == "Brown-throated Fulvetta") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Brown-throated Fulvetta") %>% .$season , method = "jaccard")


permCHTE<-adonis2(bird_diet %>% filter(SPECIES == "Chestnut-headed Tesia") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Chestnut-headed Tesia") %>% .$season , method = "jaccard")

permRCBA<-adonis2(bird_diet %>% filter(SPECIES == "Rufous-capped Babbler") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Rufous-capped Babbler") %>% .$season , method = "jaccard")

permRBNI<-adonis2(bird_diet %>% filter(SPECIES == "Rufous-bellied Niltava") %>%select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Rufous-bellied Niltava") %>% .$season , method = "jaccard")

permRVYU<-adonis2(bird_diet %>% filter(SPECIES == "Rufous-vented Yuhina") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Rufous-vented Yuhina") %>% .$season , method = "jaccard")

permRGFL<-adonis2(bird_diet %>% filter(SPECIES == "Rufous-gorgeted Flycatcher") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Rufous-gorgeted Flycatcher") %>% .$season , method = "jaccard")

permSBSB<-adonis2(bird_diet %>% filter(SPECIES == "Streak-breasted Scimitar-Babbler") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Streak-breasted Scimitar-Babbler") %>% .$season , method = "jaccard")

permRWFU<-adonis2(bird_diet %>% filter(SPECIES == "Rufous-winged Fulvetta") %>%select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Rufous-winged Fulvetta") %>% .$season , method = "jaccard")

permSTYU<-adonis2(bird_diet %>% filter(SPECIES == "Stripe-throated Yuhina") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Stripe-throated Yuhina") %>% .$season , method = "jaccard")

permWHWA<-adonis2(bird_diet %>% filter(SPECIES == "Whistler's Warbler") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Whistler's Warbler") %>% .$season , method = "jaccard")

permBBWA<-adonis2(bird_diet %>% filter(SPECIES == "Buff-barred Warbler") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Buff-barred Warbler") %>% .$season , method = "jaccard")


#######
#PERMDISP for each species
perdismig<-betadisper(vegdist(bird_diet %>% filter(mig == "mig") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(mig == "mig") %>% .$season, type = "centroid")
anova(perdismig)

perdisres<-betadisper(vegdist(bird_diet %>% filter(mig == "HE_res") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(mig == "HE_res") %>% .$season, type = "centroid")
anova(perdisres)

perdisBTFU<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Brown-throated Fulvetta") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Brown-throated Fulvetta") %>% .$season, type = "centroid")
anova(perdisBTFU)

perdisRCBA<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Rufous-capped Babbler") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Rufous-capped Babbler") %>% .$season, type = "centroid",)
anova(perdisRCBA)

perdisRVYU<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Rufous-vented Yuhina") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Rufous-vented Yuhina") %>% .$season, type = "centroid")
anova(perdisRVYU)

perdisSBSB<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Streak-breasted Scimitar-Babbler") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Streak-breasted Scimitar-Babbler") %>% .$season, type = "centroid")
anova(perdisSBSB)

perdisSTYU<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Stripe-throated Yuhina") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Stripe-throated Yuhina") %>% .$season, type = "centroid")
anova(perdisSTYU)

perdisCHTE<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Chestnut-headed Tesia") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Chestnut-headed Tesia") %>% .$season, type = "centroid")
anova(perdisCHTE)

perdisRBNI<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Rufous-bellied Niltava") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Rufous-bellied Niltava") %>% .$season, type = "centroid")
anova(perdisRBNI)

perdisRGFL<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Rufous-gorgeted Flycatcher") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Rufous-gorgeted Flycatcher") %>% .$season, type = "centroid")
anova(perdisRGFL)

perdisRWFU<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Rufous-winged Fulvetta") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Rufous-winged Fulvetta") %>% .$season, type = "centroid")
anova(perdisRWFU)

perdisWHWA<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Whistler's Warbler") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Whistler's Warbler") %>% .$season, type = "centroid")
anova(perdisWHWA)

perdisBBWA<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Buff-barred Warbler") %>% select(-c(ID, DATE.x.x, INDID, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Buff-barred Warbler") %>% .$season, type = "centroid")
anova(perdisBBWA)


####### Estimating arthropod order richness across species and season

## Resident species - arthropod order presence and absence in summer and winter
BTFU_Summer<- art_bird %>% filter(SPECIES == "Brown-throated Fulvetta" & season == "Summer") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0)) 
RCBA_Summer<- art_bird %>% filter(SPECIES == "Rufous-capped Babbler" & season == "Summer")%>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RVYU_Summer<- art_bird %>% filter(SPECIES == "Rufous-vented Yuhina" & season == "Summer") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
SBSB_Summer<- art_bird %>% filter(SPECIES == "Streak-breasted Scimitar-Babbler" & season == "Summer") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
STYU_Summer<- art_bird %>% filter(SPECIES == "Stripe-throated Yuhina" & season == "Summer") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
BTFU_Winter<- art_bird %>% filter(SPECIES == "Brown-throated Fulvetta" & season == "Winter") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RCBA_Winter<- art_bird %>% filter(SPECIES == "Rufous-capped Babbler" & season == "Winter")%>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RVYU_Winter<- art_bird %>% filter(SPECIES == "Rufous-vented Yuhina" & season == "Winter") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
SBSB_Winter<- art_bird %>% filter(SPECIES == "Streak-breasted Scimitar-Babbler" & season == "Winter") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
STYU_Winter<- art_bird %>% filter(SPECIES == "Stripe-throated Yuhina" & season == "Winter") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))

# Creating an array of resident species
HE_res_art<-list(BTFU_Summer = BTFU_Summer,RCBA_Summer = RCBA_Summer,RVYU_Summer = RVYU_Summer,SBSB_Summer = SBSB_Summer,STYU_Summer = STYU_Summer,BTFU_Winter = BTFU_Winter,RCBA_Winter = RCBA_Winter,RVYU_Winter=RVYU_Winter,SBSB_Winter=SBSB_Winter,STYU_Winter=STYU_Winter)

# Migrant Species - arthropod order presence and absence in summer and winter
CHTE_Summer<- art_bird %>% filter(SPECIES == "Chestnut-headed Tesia" & season == "Summer") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0)) 
RBNI_Summer<- art_bird %>% filter(SPECIES == "Rufous-bellied Niltava" & season == "Summer")%>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RGFL_Summer<- art_bird %>% filter(SPECIES == "Rufous-gorgeted Flycatcher" & season == "Summer") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RWFU_Summer<- art_bird %>% filter(SPECIES == "Rufous-winged Fulvetta" & season == "Summer") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
WHWA_Summer<- art_bird %>% filter(SPECIES == "Whistler's Warbler" & season == "Summer") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
BBWA_Summer<-art_bird %>% filter(SPECIES == "Buff-barred Warbler" & season == "Summer") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
CHTE_Winter<- art_bird %>% filter(SPECIES == "Chestnut-headed Tesia" & season == "Winter") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RGFL_Winter<- art_bird %>% filter(SPECIES == "Rufous-gorgeted Flycatcher" & season == "Winter")%>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RWFU_Winter<- art_bird %>% filter(SPECIES == "Rufous-winged Fulvetta" & season == "Winter") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
WHWA_Winter<- art_bird %>% filter(SPECIES == "Whistler's Warbler" & season == "Winter") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RBNI_Winter<- art_bird %>% filter(SPECIES == "Rufous-bellied Niltava" & season == "Winter") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
BBWA_Winter<- art_bird %>% filter(SPECIES == "Buff-barred Warbler" & season == "Winter") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))

# Creating an array of migrant species
mig_art<-list(BBWA_Summer = BBWA_Summer, CHTE_Summer = CHTE_Summer,RBNI_Summer = RBNI_Summer,RGFL_Summer = RGFL_Summer,RWFU_Summer = RWFU_Summer,WHWA_Summer = WHWA_Summer, BBWA_Winter = BBWA_Winter, CHTE_Winter = CHTE_Winter,RBNI_Winter = RBNI_Winter,RGFL_Winter=RGFL_Winter,RWFU_Winter=RWFU_Winter,WHWA_Winter=WHWA_Winter)

# Creating an array of all species
all<-list(BTFU_Summer = BTFU_Summer,RCBA_Summer = RCBA_Summer,RVYU_Summer = RVYU_Summer,SBSB_Summer = SBSB_Summer,STYU_Summer = STYU_Summer,BTFU_Winter = BTFU_Winter,RCBA_Winter = RCBA_Winter,RVYU_Winter=RVYU_Winter,SBSB_Winter=SBSB_Winter,STYU_Winter=STYU_Winter,BBWA_Summer = BBWA_Summer,CHTE_Summer = CHTE_Summer,RBNI_Summer = RBNI_Summer,RGFL_Summer = RGFL_Summer,RWFU_Summer = RWFU_Summer,WHWA_Summer = WHWA_Summer,BBWA_Winter = BBWA_Winter,CHTE_Winter = CHTE_Winter,RBNI_Winter = RBNI_Winter,RGFL_Winter=RGFL_Winter,RWFU_Winter=RWFU_Winter,WHWA_Winter=WHWA_Winter)

# Arthropod order presence and absence in summer and winter in High-elevation residents on the whole
HEres_Summer<- art_bird %>% filter(mig == "HE_res" & season == "Summer") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0)) 
HEres_Winter<- art_bird %>% filter(mig == "HE_res" & season == "Winter") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))

# Arthropod order presence and absence in summer and winter in elevational migrants on the whole
mig_Summer<- art_bird %>% filter(mig == "mig" & season == "Summer") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0)) 
mig_Winter<- art_bird %>% filter(mig == "mig" & season == "Winter") %>% mutate(across(c(2:25),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:25)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))

# Creating an array of all migratory strategies
all_merged_art<-list(Resident_Summer=HEres_Summer,Resident_Winter=HEres_Winter,Migrant_Summer=mig_Summer,Migrant_Winter=mig_Winter)

# Estimating richness using hill numbers
# Create a final table with estimated arthropod order richness at the species level and community level
set.seed(47)
DA_res<-estimateD(HE_res_art, datatype = 'incidence_raw', q = 0, base = "size")
DA_mig<-estimateD(mig_art, datatype = 'incidence_raw', q = 0, base = "size")
DA_res$mig<-"Resident"
DA_mig$mig<-"Migrant"
DA_res<- DA_res %>% separate(Assemblage, c('Species', 'Season'))
DA_mig<- DA_mig %>% separate(Assemblage, c('Species', 'Season'))
DA<-estimateD(all_merged_art, datatype = 'incidence_raw', q = 0, base = "size")
DA<- DA %>% separate(Assemblage, c('mig', 'Season'))
DA$Species<-"All Species"
DA_SP<-rbind(DA_res,DA_mig,DA)
DA_SP$mig <- factor(DA_SP$mig, levels=c("Resident", "Migrant"))

# Barplots of estimated arthropod order richness
Order_art<-ggplot(DA_SP, aes(Species, qD, fill = Season))+
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = qD.LCL, ymax = qD.UCL), width=.2, position = position_dodge(0.9))+
  facet_wrap(~mig, scales = "free_x", nrow = 2)+
  labs(title = "Number of Arthropod Orders",
       y = "Estimated Richness") +
  theme_minimal(base_size = 12)+
  scale_fill_grey(start = 0.3,
                  end = 0.7)+
  theme(plot.title = element_text(hjust = 0.5))


##### Number of arthropod orders in each sample - min, max, mean
art_bird$order_no<-rowSums(art_bird[,c(2:25)] != 0)

art_order_avg<-art_bird %>% group_by(SPECIES,season) %>% summarise(order_num = mean(order_no), CI = 1.96*(sd(order_no)/sqrt(n())), n= n())  %>% filter(n>5)


####### Estimating plant family richness across species and season

## Resident species - plant family presence and absence in summer and winter
BTFU_Summer<- pln_bird_family %>% filter(SPECIES == "Brown-throated Fulvetta" & season == "Summer") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0)) 
RCBA_Summer<- pln_bird_family %>% filter(SPECIES == "Rufous-capped Babbler" & season == "Summer")%>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RVYU_Summer<- pln_bird_family %>% filter(SPECIES == "Rufous-vented Yuhina" & season == "Summer") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
SBSB_Summer<- pln_bird_family %>% filter(SPECIES == "Streak-breasted Scimitar-Babbler" & season == "Summer") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
STYU_Summer<- pln_bird_family %>% filter(SPECIES == "Stripe-throated Yuhina" & season == "Summer") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
BTFU_Winter<- pln_bird_family %>% filter(SPECIES == "Brown-throated Fulvetta" & season == "Winter") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RCBA_Winter<- pln_bird_family %>% filter(SPECIES == "Rufous-capped Babbler" & season == "Winter")%>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RVYU_Winter<- pln_bird_family %>% filter(SPECIES == "Rufous-vented Yuhina" & season == "Winter") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
SBSB_Winter<- pln_bird_family %>% filter(SPECIES == "Streak-breasted Scimitar-Babbler" & season == "Winter") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
STYU_Winter<- pln_bird_family %>% filter(SPECIES == "Stripe-throated Yuhina" & season == "Winter") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))

# Creating an array of resident species
HE_res_pln<-list(BTFU_Summer = BTFU_Summer,RCBA_Summer = RCBA_Summer,RVYU_Summer = RVYU_Summer,SBSB_Summer = SBSB_Summer,STYU_Summer = STYU_Summer,BTFU_Winter = BTFU_Winter,RCBA_Winter = RCBA_Winter,RVYU_Winter=RVYU_Winter,SBSB_Winter=SBSB_Winter,STYU_Winter=STYU_Winter)


## Migrant species - plant family presence and absence in summer and winter
CHTE_Summer<- pln_bird_family %>% filter(SPECIES == "Chestnut-headed Tesia" & season == "Summer") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0)) 
RBNI_Summer<- pln_bird_family %>% filter(SPECIES == "Rufous-bellied Niltava" & season == "Summer")%>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RGFL_Summer<- pln_bird_family %>% filter(SPECIES == "Rufous-gorgeted Flycatcher" & season == "Summer") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RWFU_Summer<- pln_bird_family %>% filter(SPECIES == "Rufous-winged Fulvetta" & season == "Summer") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
WHWA_Summer<- pln_bird_family %>% filter(SPECIES == "Whistler's Warbler" & season == "Summer") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
BBWA_Summer<- pln_bird_family %>% filter(SPECIES == "Buff-barred Warbler" & season == "Summer") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
CHTE_Winter<- pln_bird_family %>% filter(SPECIES == "Chestnut-headed Tesia" & season == "Winter") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RGFL_Winter<- pln_bird_family %>% filter(SPECIES == "Rufous-gorgeted Flycatcher" & season == "Winter")%>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RWFU_Winter<- pln_bird_family %>% filter(SPECIES == "Rufous-winged Fulvetta" & season == "Winter") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
WHWA_Winter<- pln_bird_family %>% filter(SPECIES == "Whistler's Warbler" & season == "Winter") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RBNI_Winter<- pln_bird_family %>% filter(SPECIES == "Rufous-bellied Niltava" & season == "Winter") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
BBWA_Winter<- pln_bird_family %>% filter(SPECIES == "Buff-barred Warbler" & season == "Winter") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))

# Creating an array of migrant species
mig_pln<-list(BBWA_Summer = BBWA_Summer, CHTE_Summer = CHTE_Summer,RBNI_Summer = RBNI_Summer,RGFL_Summer = RGFL_Summer,RWFU_Summer = RWFU_Summer,WHWA_Summer = WHWA_Summer, BBWA_Winter = BBWA_Winter, CHTE_Winter = CHTE_Winter,RBNI_Winter = RBNI_Winter,RGFL_Winter=RGFL_Winter,RWFU_Winter=RWFU_Winter,WHWA_Winter=WHWA_Winter)

# Plant family presence and absence in summer and winter in High-elevation residents on the whole
HEres_Summer<- pln_bird_family %>% filter(mig == "HE_res" & season == "Summer") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0)) 
HEres_Winter<- pln_bird_family %>% filter(mig == "HE_res" & season == "Winter") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))

# Plant family presence and absence in summer and winter in elevational migrants on the whole
mig_Summer<- pln_bird_family %>% filter(mig == "mig" & season == "Summer") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0)) 
mig_Winter<- pln_bird_family %>% filter(mig == "mig" & season == "Winter") %>% mutate(across(c(2:103),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:103)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))

# Creating an array of all migratory strategies
all_merged_pln<-list(Resident_Summer=HEres_Summer,Resident_Winter=HEres_Winter,Migrant_Summer=mig_Summer,Migrant_Winter=mig_Winter)

# Estimating richness using hill numbers
# Create a final table with estimated plant family richness at the species level and community level
set.seed(47)
DP_res<-estimateD(HE_res_pln, datatype = 'incidence_raw', q = 0, base = "size")
DP_mig<-estimateD(mig_pln, datatype = 'incidence_raw', q = 0, base = "size")
DP_res$mig<-"Resident"
DP_mig$mig<-"Migrant"
DP_res<- DP_res %>% separate(Assemblage, c('Species', 'Season'))
DP_mig<- DP_mig %>% separate(Assemblage, c('Species', 'Season'))
DP<-estimateD(all_merged_pln, datatype = 'incidence_raw', q = 0, base = "size")
DP<- DP %>% separate(Assemblage, c('mig', 'Season'))
DP$Species<-"All Species"
DP_SP<-rbind(DP_res,DP_mig,DP)
DP_SP$mig <- factor(DP_SP$mig, levels=c("Resident", "Migrant"))

# Barplots of estimated plant family richness
Family_pln<-ggplot(DP_SP, aes(Species, qD, fill = Season))+
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = qD.LCL, ymax = qD.UCL), width=.2, position = position_dodge(0.9))+
  facet_wrap(~mig, scales = "free_x", nrow = 2)+
  labs(title = "Number of Plant Families",
       y = "Estimated Richness") +
  theme_minimal(base_size = 12)+
  scale_fill_grey(start = 0.3,
                  end = 0.7)+
  theme(plot.title = element_text(hjust = 0.5))

##### Number of plant familes in each sample - min, max, mean
pln_bird_family$family_no<-rowSums(pln_bird_family[,c(2:103)] != 0)

pln_fam_avg<-pln_bird_family%>% group_by(SPECIES,season) %>% summarise(family_no = mean(family_no), CI = 1.96*(sd(family_no)/sqrt(n())), n= n())  %>% filter(n>5) 


####### Estimating arthropod ASV richness across species and season

# Create the summary table of counts of each arthropod ASV in each sample
# Fix typos in sample names
result_art_asv <- art_non_zero %>%
  group_by(ID, ASV) %>%
  summarise(Count = sum(Count), .groups = 'drop') %>%
  pivot_wider(names_from = ASV, values_from = Count, values_fill = 0) %>% mutate(ID= recode(ID, "Z2662"="A2662", "A2584" = "Z2584", "A1297"="AB1297", "A2591"="Z2591", "A2594.25422" = "A2690", "A2593"="Z2593",  "A2598"="Z2598",  "A2599"="Z2599", "Z2388" = "A2883", "A1739_10124"="AB1739_10124", "A8224_19124" = "A8224_14124", "A8233_16124"="A8233_15124", "A8317_1224"="Z8317_1224","AB1179_16124" = "AB1179_161223",  "Z5093_19124"="Z5092_19124", "Z5019_22124" = "Z2711_23522",  "A3700_22323"="A3700_22523","A3841_22523" =  "Z3841_22523", "A3865_23523" = "Z3865_23523", "A5119_20124" = "Z5119_20124","Z5152_1224" = "Z2513_30422", "AB0525_18523" = "AB0525_19523", "Z3632_6523" = "A3632_6523", "A3606_18523"="Z3806_18523", "A269807_25422" = "A269087_25422", "AB0373_26422" = "AB0373_26423", "Z2604_5622" = "Z2604_5522", "Z2605_6522"= "Z2605_5522", "Z2612_5622"= "Z2612_5522", "Z3428_26422" = "Z3428_26423", "A3110_141222" = "A3110_151222", "A3165_7123" = "A3165_8123", "A5004_121223" = "A5104_121223","A5200_7124" = "A5200_8124"))

# Merge with ringing and guild data
art_bird_asv<- left_join(result_art_asv,samples_date22) %>% mutate(ID = str_remove(ID, "\\..*"))
art_bird_asv$ID <- ifelse(is.na(art_bird_asv$DATE), art_bird_asv$ID, paste0(art_bird_asv$ID, "_",sub("^0", "", sub("^(\\d+)-0?(\\d+)-(\\d{2})(\\d{2})$", "\\1\\2\\4", art_bird_asv$DATE))))
art_bird_asv<-art_bird_asv%>% left_join(bird_data, by = "ID") %>% left_join(guilds)

## Resident species - Arthropod ASV presence and absence in summer and winter
BTFU_Summer<- art_bird_asv %>% filter(SPECIES == "Brown-throated Fulvetta" & season == "Summer") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0)) 
RCBA_Summer<- art_bird_asv %>% filter(SPECIES == "Rufous-capped Babbler" & season == "Summer")%>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RVYU_Summer<- art_bird_asv %>% filter(SPECIES == "Rufous-vented Yuhina" & season == "Summer") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
SBSB_Summer<- art_bird_asv %>% filter(SPECIES == "Streak-breasted Scimitar-Babbler" & season == "Summer") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
STYU_Summer<- art_bird_asv %>% filter(SPECIES == "Stripe-throated Yuhina" & season == "Summer") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
BTFU_Winter<- art_bird_asv %>% filter(SPECIES == "Brown-throated Fulvetta" & season == "Winter") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RCBA_Winter<- art_bird_asv %>% filter(SPECIES == "Rufous-capped Babbler" & season == "Winter")%>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RVYU_Winter<- art_bird_asv %>% filter(SPECIES == "Rufous-vented Yuhina" & season == "Winter") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
SBSB_Winter<- art_bird_asv %>% filter(SPECIES == "Streak-breasted Scimitar-Babbler" & season == "Winter") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
STYU_Winter<- art_bird_asv %>% filter(SPECIES == "Stripe-throated Yuhina" & season == "Winter") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))

HE_res_art<-list(BTFU_Summer = BTFU_Summer,RCBA_Summer = RCBA_Summer,RVYU_Summer = RVYU_Summer,SBSB_Summer = SBSB_Summer,STYU_Summer = STYU_Summer,BTFU_Winter = BTFU_Winter,RCBA_Winter = RCBA_Winter,RVYU_Winter=RVYU_Winter,SBSB_Winter=SBSB_Winter,STYU_Winter=STYU_Winter)

## Migrant species - Arthropod ASV presence and absence in summer and winter
CHTE_Summer<- art_bird_asv %>% filter(SPECIES == "Chestnut-headed Tesia" & season == "Summer") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0)) 
RBNI_Summer<- art_bird_asv %>% filter(SPECIES == "Rufous-bellied Niltava" & season == "Summer")%>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RGFL_Summer<- art_bird_asv %>% filter(SPECIES == "Rufous-gorgeted Flycatcher" & season == "Summer") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RWFU_Summer<- art_bird_asv %>% filter(SPECIES == "Rufous-winged Fulvetta" & season == "Summer") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
WHWA_Summer<- art_bird_asv %>% filter(SPECIES == "Whistler's Warbler" & season == "Summer") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
BBWA_Summer<- art_bird_asv %>% filter(SPECIES == "Buff-barred Warbler" & season == "Summer") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
CHTE_Winter<- art_bird_asv %>% filter(SPECIES == "Chestnut-headed Tesia" & season == "Winter") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RGFL_Winter<- art_bird_asv %>% filter(SPECIES == "Rufous-gorgeted Flycatcher" & season == "Winter")%>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RWFU_Winter<- art_bird_asv %>% filter(SPECIES == "Rufous-winged Fulvetta" & season == "Winter") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
WHWA_Winter<- art_bird_asv %>% filter(SPECIES == "Whistler's Warbler" & season == "Winter") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RBNI_Winter<- art_bird_asv %>% filter(SPECIES == "Rufous-bellied Niltava" & season == "Winter") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
BBWA_Winter<- art_bird_asv %>% filter(SPECIES == "Buff-barred Warbler" & season == "Winter") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))

mig_art<-list(BBWA_Summer = BBWA_Summer, CHTE_Summer = CHTE_Summer,RBNI_Summer = RBNI_Summer,RGFL_Summer = RGFL_Summer,RWFU_Summer = RWFU_Summer,WHWA_Summer = WHWA_Summer, BBWA_Winter = BBWA_Winter, CHTE_Winter = CHTE_Winter,RBNI_Winter = RBNI_Winter,RGFL_Winter=RGFL_Winter,RWFU_Winter=RWFU_Winter,WHWA_Winter=WHWA_Winter)

# Arthropod ASV presence and absence in summer and winter in High-elevation residents on the whole
HEres_Summer<- art_bird_asv %>% filter(mig == "HE_res" & season == "Summer") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0)) 
HEres_Winter<- art_bird_asv %>% filter(mig == "HE_res" & season == "Winter") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))

# Arthropod ASV presence and absence in summer and winter in elevational migrants on the whole
mig_Summer<- art_bird_asv %>% filter(mig == "mig" & season == "Summer") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0)) 
mig_Winter<- art_bird_asv %>% filter(mig == "mig" & season == "Winter") %>% mutate(across(c(2:3010),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:3010)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))

all_merged_art<-list(Resident_Summer=HEres_Summer,Resident_Winter=HEres_Winter,Migrant_Summer=mig_Summer,Migrant_Winter=mig_Winter)

# Estimating richness using hill numbers
# Create a final table with estimated arthropod ASV richness at the species level and community level
set.seed(47)
DA_res<-estimateD(HE_res_art, datatype = 'incidence_raw', q = 0, base = "size")
DA_mig<-estimateD(mig_art, datatype = 'incidence_raw', q = 0, base = "size")
DA_res$mig<-"Resident"
DA_mig$mig<-"Migrant"
DA_res<- DA_res %>% separate(Assemblage, c('Species', 'Season'))
DA_mig<- DA_mig %>% separate(Assemblage, c('Species', 'Season'))
DA<-estimateD(all_merged_art, datatype = 'incidence_raw', q = 0, base = "size")
DA<- DA %>% separate(Assemblage, c('mig', 'Season'))
DA$Species<-"All Species"
DA_SP<-rbind(DA_res,DA_mig,DA)
DA_SP$mig <- factor(DA_SP$mig, levels=c("Resident", "Migrant"))

# Barplots of estimated arthropod ASV richness
ASV_art<-ggplot(DA_SP, aes(Species, qD, fill = Season))+
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = qD.LCL, ymax = qD.UCL), width=.2, position = position_dodge(0.9))+
  facet_wrap(~mig, scales = "free_x", nrow = 2)+
  labs(title = "Number of Arthropod ASVs",
       y = "Estimated Richness") +
  theme_minimal(base_size = 12)+
  scale_fill_grey(start = 0.3,
                  end = 0.7)+
  theme(plot.title = element_text(hjust = 0.5))

##### Number of arthropod ASVs in each sample - min, max, mean
art_bird_asv$asv<-rowSums(art_bird_asv[,c(2:3010)] != 0)

art_asv_avg<-art_bird_asv %>% group_by(SPECIES,season) %>% summarise(asv_no = mean(asv), CI = 1.96*(sd(asv)/sqrt(n())), n= n())  %>% filter(n>5)


####### Estimating plant ASV richness across species and season

# Create the summary table of counts of each plant ASV in each sample
# Fix typos in sample names
result_pln_asv <- plnF_non_zero %>%
  group_by(ID, ASV) %>%
  summarise(Count = sum(Count), .groups = 'drop') %>%
  pivot_wider(names_from = ASV, values_from = Count, values_fill = 0) %>% mutate(ID= recode(ID, "Z2662"="A2662", "A2584" = "Z2584", "A1297"="AB1297", "A2591"="Z2591", "A2594.25422" = "A2690", "A2593"="Z2593",  "A2598"="Z2598",  "A2599"="Z2599", "Z2388" = "A2883", "A1739_10124"="AB1739_10124", "A8224_19124" = "A8224_14124", "A8233_16124"="A8233_15124", "A8317_1224"="Z8317_1224","AB1179_16124" = "AB1179_161223",  "Z5093_19124"="Z5092_19124", "Z5019_22124" = "Z2711_23522",  "A3700_22323"="A3700_22523","A3841_22523" =  "Z3841_22523", "A3865_23523" = "Z3865_23523", "A5119_20124" = "Z5119_20124","Z5152_1224" = "Z2513_30422", "AB0525_18523" = "AB0525_19523", "Z3632_6523" = "A3632_6523", "A3606_18523"="Z3806_18523", "A269807_25422" = "A269087_25422", "AB0373_26422" = "AB0373_26423", "Z2604_5622" = "Z2604_5522", "Z2605_6522"= "Z2605_5522", "Z2612_5622"= "Z2612_5522", "Z3428_26422" = "Z3428_26423", "A3110_141222" = "A3110_151222", "A3165_7123" = "A3165_8123", "A5004_121223" = "A5104_121223","A5200_7124" = "A5200_8124"))

# Merge with ringing and guild data
pln_bird_asv<- left_join(result_pln_asv,samples_date22) %>% mutate(ID = str_remove(ID, "\\..*"))
pln_bird_asv$ID <- ifelse(is.na(pln_bird_asv$DATE), pln_bird_asv$ID, paste0(pln_bird_asv$ID, "_",sub("^0", "", sub("^(\\d+)-0?(\\d+)-(\\d{2})(\\d{2})$", "\\1\\2\\4", pln_bird_asv$DATE))))
pln_bird_asv<-pln_bird_asv%>% left_join(bird_data, by = "ID") %>% left_join(guilds)

## Resident species - plant ASV presence and absence in summer and winter
BTFU_Summer<- pln_bird_asv %>% filter(SPECIES == "Brown-throated Fulvetta" & season == "Summer") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0)) 
RCBA_Summer<- pln_bird_asv %>% filter(SPECIES == "Rufous-capped Babbler" & season == "Summer")%>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RVYU_Summer<- pln_bird_asv %>% filter(SPECIES == "Rufous-vented Yuhina" & season == "Summer") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
SBSB_Summer<- pln_bird_asv %>% filter(SPECIES == "Streak-breasted Scimitar-Babbler" & season == "Summer") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
STYU_Summer<- pln_bird_asv %>% filter(SPECIES == "Stripe-throated Yuhina" & season == "Summer") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
BTFU_Winter<- pln_bird_asv %>% filter(SPECIES == "Brown-throated Fulvetta" & season == "Winter") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RCBA_Winter<- pln_bird_asv %>% filter(SPECIES == "Rufous-capped Babbler" & season == "Winter")%>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RVYU_Winter<- pln_bird_asv %>% filter(SPECIES == "Rufous-vented Yuhina" & season == "Winter") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
SBSB_Winter<- pln_bird_asv %>% filter(SPECIES == "Streak-breasted Scimitar-Babbler" & season == "Winter") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
STYU_Winter<- pln_bird_asv %>% filter(SPECIES == "Stripe-throated Yuhina" & season == "Winter") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))

HE_res_pln<-list(BTFU_Summer = BTFU_Summer,RCBA_Summer = RCBA_Summer,RVYU_Summer = RVYU_Summer,SBSB_Summer = SBSB_Summer,STYU_Summer = STYU_Summer,BTFU_Winter = BTFU_Winter,RCBA_Winter = RCBA_Winter,RVYU_Winter=RVYU_Winter,SBSB_Winter=SBSB_Winter,STYU_Winter=STYU_Winter)

## Migrant species - plant ASV presence and absence in summer and winter
CHTE_Summer<- pln_bird_asv %>% filter(SPECIES == "Chestnut-headed Tesia" & season == "Summer") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0)) 
RBNI_Summer<- pln_bird_asv %>% filter(SPECIES == "Rufous-bellied Niltava" & season == "Summer")%>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RGFL_Summer<- pln_bird_asv %>% filter(SPECIES == "Rufous-gorgeted Flycatcher" & season == "Summer") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RWFU_Summer<- pln_bird_asv %>% filter(SPECIES == "Rufous-winged Fulvetta" & season == "Summer") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
WHWA_Summer<- pln_bird_asv %>% filter(SPECIES == "Whistler's Warbler" & season == "Summer") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
BBWA_Summer<- pln_bird_asv %>% filter(SPECIES == "Buff-barred Warbler" & season == "Summer") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
CHTE_Winter<- pln_bird_asv %>% filter(SPECIES == "Chestnut-headed Tesia" & season == "Winter") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RGFL_Winter<- pln_bird_asv %>% filter(SPECIES == "Rufous-gorgeted Flycatcher" & season == "Winter")%>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RWFU_Winter<- pln_bird_asv %>% filter(SPECIES == "Rufous-winged Fulvetta" & season == "Winter") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
WHWA_Winter<- pln_bird_asv %>% filter(SPECIES == "Whistler's Warbler" & season == "Winter") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
RBNI_Winter<- pln_bird_asv %>% filter(SPECIES == "Rufous-bellied Niltava" & season == "Winter") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))
BBWA_Winter<- pln_bird_asv %>% filter(SPECIES == "Buff-barred Warbler" & season == "Winter") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))

mig_pln<-list(BBWA_Summer = BBWA_Summer, CHTE_Summer = CHTE_Summer,RBNI_Summer = RBNI_Summer,RGFL_Summer = RGFL_Summer,RWFU_Summer = RWFU_Summer,WHWA_Summer = WHWA_Summer, BBWA_Winter = BBWA_Winter, CHTE_Winter = CHTE_Winter,RBNI_Winter = RBNI_Winter,RGFL_Winter=RGFL_Winter,RWFU_Winter=RWFU_Winter,WHWA_Winter=WHWA_Winter)

# Plant ASV presence and absence in summer and winter in High-elevation residents on the whole
HEres_Summer<- pln_bird_asv %>% filter(mig == "HE_res" & season == "Summer") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0)) 
HEres_Winter<- pln_bird_asv %>% filter(mig == "HE_res" & season == "Winter") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))

# Plant ASV presence and absence in summer and winter in elevational migrants on the whole
mig_Summer<- pln_bird_asv %>% filter(mig == "mig" & season == "Summer") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0)) 
mig_Winter<- pln_bird_asv %>% filter(mig == "mig" & season == "Winter") %>% mutate(across(c(2:421),  ~ ifelse(. > 0, 1, .))) %>%select(c(2:421)) %>% t() %>% as.data.frame() %>% filter_all(any_vars(. != 0))

all_merged_pln<-list(Resident_Summer=HEres_Summer,Resident_Winter=HEres_Winter,Migrant_Summer=mig_Summer,Migrant_Winter=mig_Winter)

# Estimating richness using hill numbers
# Create a final table with estimated plant ASV richness at the species level and community level
set.seed(47)
DP_res<-estimateD(HE_res_pln, datatype = 'incidence_raw', q = 0, base = "size")
DP_mig<-estimateD(mig_pln, datatype = 'incidence_raw', q = 0, base = "size")
DP_res$mig<-"Resident"
DP_mig$mig<-"Migrant"
DP_res<- DP_res %>% separate(Assemblage, c('Species', 'Season'))
DP_mig<- DP_mig %>% separate(Assemblage, c('Species', 'Season'))
DP<-estimateD(all_merged_pln, datatype = 'incidence_raw', q = 0, base = "size")
DP<- DP %>% separate(Assemblage, c('mig', 'Season'))
DP$Species<-"All Species"
DP_SP<-rbind(DP_res,DP_mig,DP)
DP_SP$mig <- factor(DP_SP$mig, levels=c("Resident", "Migrant"))

# Barplots of estimated plant ASV richness
ASV_pln<-ggplot(DP_SP, aes(Species, qD, fill = Season))+
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = qD.LCL, ymax = qD.UCL), width=.2, position = position_dodge(0.9))+
  facet_wrap(~mig, scales = "free_x", nrow = 2)+
  labs(title = "Number of Plant ASVs",
       y = "Estimated Richness") +
  theme_minimal(base_size = 12)+
  scale_fill_grey(start = 0.3,
                  end = 0.7)+
  theme(plot.title = element_text(hjust = 0.5))

#### Number of plant ASVs in each sample - min, max, mean
pln_bird_asv$asv<-rowSums(pln_bird_asv[,c(2:421)] != 0)

pln_asv_avg<-pln_bird_asv %>% group_by(SPECIES,season) %>% summarise(asv_no = mean(asv), CI = 1.96*(sd(asv)/sqrt(n())), n= n())  %>% filter(n>5) 

### Plotting all 4 sets of barplots - Supplementary figure 2
art_div<-ggarrange(ASV_art,Order_art, common.legend = T, ncol = 2, legend =  "right")

pln_div<-ggarrange(ASV_pln,Family_pln, common.legend = T, ncol = 2, legend =  "right")

com_div<-ggarrange(art_div,pln_div, common.legend = T, nrow = 2, legend =  "right", labels = c("A", "B"))

#ggsave(plot = com_div, "FigS2.jpeg", width = 11, height = 9)
