## These scripts are to be run on the output of the dada2 pipeline (per sample ASV counts) and the ecotag script (Taxonomically annotated ASVs) to replicate the faecal DNA metabarcoding related visualisation and analysis in the manuscript

setwd("~/Elev_migration_diet")

library(tidyverse)
library(vegan)
library(ggpubr)
library(gridExtra)
library(ape)
library(brms)
library(marginaleffects)
library(purrr)
library(broom.mixed)
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
  pivot_wider(names_from = order_name, values_from = Count, values_fill = 0) %>% mutate(ID= recode(ID, "Z2662"="A2662", "A2584" = "Z2584", "A1297"="AB1297", "A2591"="Z2591", "A2594.25422" = "A2690", "A2593"="Z2593",  "A2598"="Z2598",  "A2599"="Z2599", "Z2388" = "A2883", "A1739_10124"="AB1739_10124", "A8224_19124" = "A8224_14124", "A8233_16124"="A8233_15124", "A8317_1224"="Z8317_1224","AB1179_16124" = "AB1179_161223",  "Z5093_19124"="Z5092_19124", "Z5019_22124" = "Z2711_23522",  "A3700_22323"="A3700_22523","A3841_22523" =  "Z3841_22523", "A3865_23523" = "Z3865_23523", "A5119_20124" = "Z5119_20124","Z5152_1224" = "Z2513_30422", "AB0525_18523" = "AB0525_19523", "Z3632_6523" = "A3632_6523", "A3606_18523"="Z3806_18523", "A269807_25422" = "A269087_25422", "AB0373_26422" = "AB0373_26423", "Z2604_5622" = "Z2604_5522", "Z2605_6522"= "Z2605_5522", "Z2612_5622"= "Z2612_5522", "Z3428_26422" = "Z3428_26423", "A3110_141222" = "A3110_151222", "A3165_7123" = "A3165_8123", "A5004_121223" = "A5104_121223", "A3627_6523" = "A3673_19523", "A3673_19523" = "A3627_6523","A5200_7124" = "A5200_8124"))


# Read in the data that maps the sample names with information on species, elevation and date of sample collection
bird_data<-read.csv("ring_number_sci.csv")

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
art_bird<-art_bird%>% left_join(bird_data, by = "ID") %>% left_join(guilds) %>% filter(mig == "HE_res" | mig == "mig")

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
  mutate(mig = recode(mig, "mig" = "Migrant")) 

art_all$mig <- factor(art_all$mig, levels=c("Resident", "Migrant"))
  
art_all_plot <-  art_all %>% ggplot(aes(x = mig, y = proportion, fill = season)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~arthropod_order ) +
  labs(x = "Migratory Status",
       y = "Frequency of Occurrence",
       fill = "Season") +
  theme_bw()+
  theme(axis.title = element_text(size =15), legend.title = element_text(size =20), legend.text = element_text(size =15), axis.text = element_text(size =11), strip.text = element_text(size = 15))+
  scale_fill_manual(values = cbbPalette)+
  theme(legend.position = "inside",
        legend.justification = c(0.9,0.2),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  ggtitle("Figure 3")

#ggsave(plot = art_all_plot, "Figure3.tif", device = "tiff", dpi = 600, compression = "lzw", width = 11, height = 7)


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

#######

github.directory <- "https://raw.githubusercontent.com/evolucionario/BigBirdTree/refs/heads/main/"

stage <- "RAGBackbone/"

tree <- "BBtreeC2022.tre"

url <- paste0(github.directory, stage, tree)

BBtree2 <- read.tree(url)
BBtree2$tip.label[BBtree2$tip.label == "Phylloscopus_burkii"] <- "Phylloscopus_whistleri"
BBtree2$tip.label[BBtree2$tip.label == "Phylloscopus_inornatus"] <- "Phylloscopus_pulcher"

Diet_Binary<-art_bird %>%
  mutate(across(Araneae:Trombidiformes, ~ ifelse(. > 0, 1, 0)))%>%
  mutate(Species = str_replace_all(Species, " ", "_"))

Diet_Binary$phylo_id<-Diet_Binary$Species

tree_pruned <- keep.tip(BBtree2, Diet_Binary$Species)
phylo_matrix <- vcv.phylo(tree_pruned, corr = TRUE)
Diet_Binary$mig <- factor(Diet_Binary$mig, 
                          levels = c("mig", "HE_res"), 
                          labels = c("Migrant", "Resident"))


spider<- brm(Araneae ~ season*mig + (1 | gr(phylo_id, cov = phylo)),
             data = Diet_Binary,
             data2 = list(phylo = phylo_matrix),
             family = bernoulli(),
             chains = 4, iter = 4000, cores = 4,
             control = list(adapt_delta = 0.95))

comparisons(spider, variables = "season", by = "mig")

bug<- brm(Hemiptera ~ season*mig + (1 | gr(phylo_id, cov = phylo)),
             data = Diet_Binary,
             data2 = list(phylo = phylo_matrix),
             family = bernoulli(),
             chains = 4, iter = 4000, cores = 4,
             control = list(adapt_delta = 0.95))

comparisons(bug, variables = "season", by = "mig")

beetle<- brm(Coleoptera ~ season*mig + (1 | gr(phylo_id, cov = phylo)),
             data = Diet_Binary,
             data2 = list(phylo = phylo_matrix),
             family = bernoulli(),
             chains = 4, iter = 4000, cores = 4,
             control = list(adapt_delta = 0.95))

comparisons(beetle, variables = "season", by = "mig")

lep<- brm(Lepidoptera ~ season*mig + (1 | gr(phylo_id, cov = phylo)),
             data = Diet_Binary,
             data2 = list(phylo = phylo_matrix),
             family = bernoulli(),
             chains = 4, iter = 4000, cores = 4,
             control = list(adapt_delta = 0.95))

comparisons(lep, variables = "season", by = "mig")

fly<- brm(Diptera ~ season*mig + (1 | gr(phylo_id, cov = phylo)),
             data = Diet_Binary,
             data2 = list(phylo = phylo_matrix),
             family = bernoulli(),
             chains = 4, iter = 4000, cores = 4,
             control = list(adapt_delta = 0.95))

comparisons(fly, variables = "season", by = "mig")

arthropod_models <- list(
  "Araneae" = spider,
  "Hemiptera" = bug,
  "Coleoptera" = beetle,
  "Lepidoptera" = lep,
  "Diptera" = fly
)

arthropod_long_table <- arthropod_models %>%
  map_dfr(function(m) {
    # tidy(m) with effects = c("fixed", "ran_pars") pulls both regression and SD
    tidy(m, effects = c("fixed", "ran_pars"), conf.int = TRUE, rhat = TRUE, ess = TRUE)
  }, .id = "Prey_Order") %>%
  mutate(
    Parameter = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "seasonWinter" ~ "Season_winter",
      term == "migResident" ~ "MigStrat_resident",
      term == "seasonWinter:migResident" ~ "MigStrat*Season",
      term == "sd__(Intercept)" ~ "Phylogenetic SD (sd_Intercept)",
      TRUE ~ term
    ),
    Estimate = round(estimate, 2),
    `95% CI` = paste0("[", round(conf.low, 2), ", ", round(conf.high, 2), "]"),
    Rhat = round(rhat, 3),
    ESS = round(ess, 0)
  ) %>%
  select(Prey_Order, Parameter, Estimate, `95% CI`, Rhat, ESS)

#write.csv(arthropod_long_table, "Table_S5.csv", row.names = FALSE)

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
  pivot_wider(names_from = family_name, values_from = Count, values_fill = 0) %>% mutate(ID= recode(ID, "Z2662"="A2662", "A2584" = "Z2584", "A1297"="AB1297", "A2591"="Z2591", "A2594.25422" = "A2690", "A2593"="Z2593",  "A2598"="Z2598",  "A2599"="Z2599", "Z2388" = "A2883", "A1739_10124"="AB1739_10124", "A8224_19124" = "A8224_14124", "A8233_16124"="A8233_15124", "A8317_1224"="Z8317_1224","AB1179_16124" = "AB1179_161223",  "Z5093_19124"="Z5092_19124", "Z5019_22124" = "Z2711_23522",  "A3700_22323"="A3700_22523","A3841_22523" =  "Z3841_22523", "A3865_23523" = "Z3865_23523", "A5119_20124" = "Z5119_20124","Z5152_1224" = "Z2513_30422", "AB0525_18523" = "AB0525_19523", "Z3632_6523" = "A3632_6523", "A3606_18523"="Z3806_18523", "A269807_25422" = "A269087_25422", "AB0373_26422" = "AB0373_26423", "Z2604_5622" = "Z2604_5522", "Z2605_6522"= "Z2605_5522", "Z2612_5622"= "Z2612_5522", "Z3428_26422" = "Z3428_26423", "A3110_141222" = "A3110_151222", "A3165_7123" = "A3165_8123", "A5004_121223" = "A5104_121223", "A3627_6523" = "A3673_19523", "A3673_19523" = "A3627_6523","A5200_7124" = "A5200_8124"))

## Merging the all the datasets
pln_bird_family<- left_join(result_pln_family,samples_date22) %>% mutate(ID = str_remove(ID, "\\..*"))
pln_bird_family$ID <- ifelse(is.na(pln_bird_family$DATE), pln_bird_family$ID, paste0(pln_bird_family$ID, "_",sub("^0", "", sub("^(\\d+)-0?(\\d+)-(\\d{2})(\\d{2})$", "\\1\\2\\4", pln_bird_family$DATE))))
pln_bird_family<-pln_bird_family%>% left_join(bird_data, by = "ID") %>% left_join(guilds) %>% filter(mig == "HE_res" | mig == "mig")

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

#
plant_Binary<-pln_bird_family %>%
  mutate(across(Acanthaceae:Zingiberaceae, ~ ifelse(. > 0, 1, 0)))%>%
  mutate(Species = str_replace_all(Species, " ", "_"))

plant_Binary$phylo_id<-plant_Binary$Species

tree_pruned <- keep.tip(BBtree2, plant_Binary$Species)
phylo_matrix <- vcv.phylo(tree_pruned, corr = TRUE)
plant_Binary$mig <- factor(plant_Binary$mig, 
                          levels = c("mig", "HE_res"), 
                          labels = c("Migrant", "Resident"))

knotweed<- brm(Polygonaceae ~ season*mig + (1 | gr(phylo_id, cov = phylo)),
             data = plant_Binary,
             data2 = list(phylo = phylo_matrix),
             family = bernoulli(),
             chains = 4, iter = 4000, cores = 4,
             control = list(adapt_delta = 0.98))

comparisons(knotweed, variables = "season", by = "mig")

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
  pivot_wider(names_from = order_name, values_from = Count, values_fill = 0) %>% mutate(ID= recode(ID, "Z2662"="A2662", "A2584" = "Z2584", "A1297"="AB1297", "A2591"="Z2591", "A2594.25422" = "A2690", "A2593"="Z2593",  "A2598"="Z2598",  "A2599"="Z2599", "Z2388" = "A2883", "A1739_10124"="AB1739_10124", "A8224_19124" = "A8224_14124", "A8233_16124"="A8233_15124", "A8317_1224"="Z8317_1224","AB1179_16124" = "AB1179_161223",  "Z5093_19124"="Z5092_19124", "Z5019_22124" = "Z2711_23522",  "A3700_22323"="A3700_22523","A3841_22523" =  "Z3841_22523", "A3865_23523" = "Z3865_23523", "A5119_20124" = "Z5119_20124","Z5152_1224" = "Z2513_30422", "AB0525_18523" = "AB0525_19523", "Z3632_6523" = "A3632_6523", "A3606_18523"="Z3806_18523", "A269807_25422" = "A269087_25422", "AB0373_26422" = "AB0373_26423", "Z2604_5622" = "Z2604_5522", "Z2605_6522"= "Z2605_5522", "Z2612_5622"= "Z2612_5522", "Z3428_26422" = "Z3428_26423", "A3110_141222" = "A3110_151222", "A3165_7123" = "A3165_8123", "A5004_121223" = "A5104_121223", "A3627_6523" = "A3673_19523", "A3673_19523" = "A3627_6523","A5200_7124" = "A5200_8124"))

## Merging the all the datasets
pln_bird_order<- left_join(result_pln_order,samples_date22) %>% mutate(ID = str_remove(ID, "\\..*"))
pln_bird_order$ID <- ifelse(is.na(pln_bird_order$DATE), pln_bird_order$ID, paste0(pln_bird_order$ID, "_",sub("^0", "", sub("^(\\d+)-0?(\\d+)-(\\d{2})(\\d{2})$", "\\1\\2\\4", pln_bird_order$DATE))))
pln_bird_order<-pln_bird_order%>% left_join(bird_data, by = "ID") %>% left_join(guilds) %>% filter(mig == "HE_res" | mig == "mig")



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
x.dist<-vegdist(bird_diet[,-c(1,26:37,78:86)], method = "jaccard", binary = T)
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
PCoA_all<-pcoa_PLOT %>% ggplot(aes(X1, X2)) + geom_point(aes(colour = migsea, shape = season), show.legend = FALSE) + 
  geom_polygon(data = df_ell, aes(x=Dim1, y=Dim2, colour = migsea), fill = NA , linewidth=1)+
  scale_color_manual(values = cbbPalette)+
  xlab("PCoA 1") + ylab("PCoA2")+
  xlim(-0.5,0.5) + ylim(-0.5,0.5)+
  coord_fixed(ratio = 1)+
  theme_bw()+
  theme(axis.title = element_text(size =15), legend.title = element_text(size =15), legend.text = element_text(size =15), axis.text = element_text(size =11), strip.text = element_text(size = 15))+
  guides(color = guide_legend(title = "Migratory Behaviour \nand Season"), shape = "none", linetype = "none", override.aes = list(shape = NA))+
  ggtitle("Figure 4")

#ggsave(plot = PCoA_all, "Figure4.tif",  device = "tiff", dpi = 600, compression = "lzw", width = 11, height = 6)


## PCoA at the species level
# Migratory species
mig_diet<-bird_diet %>% filter(SPECIES %in% c("Whistler's Warbler", "Rufous-gorgeted Flycatcher", "Chestnut-headed Tesia", "Rufous-winged Fulvetta", "Rufous-bellied Niltava", "Buff-barred Warbler"))

mig.dist<-vegdist(mig_diet[,-c(1,26:37,78:86)], method = "jaccard", binary = T)
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
  coord_fixed(ratio = 1)+
  scale_color_manual(values = cbbPalette)+xlim(-0.5,0.5)+
  xlab("PCoA 1") + ylab("PCoA2")+
  guides(color = guide_legend(title = "Migratory Species"), shape = "none", linetype = "none")

## PCoA at the species level
# Resident species
res_diet<-bird_diet %>% filter(SPECIES %in% c("Brown-throated Fulvetta", "Rufous-vented Yuhina", "Streak-breasted Scimitar-Babbler", "Stripe-throated Yuhina", "Rufous-capped Babbler"))

res.dist<-vegdist(res_diet[,-c(1,26:37,78:86)], method = "jaccard", binary = T)
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
  coord_fixed(ratio = 1)+
  scale_color_manual(values = cbbPalette) +xlim(-0.5,0.5)+
  xlab("PCoA 1") + ylab("PCoA2")+
  guides(color = guide_legend(title = "Resident Species"), shape = "none", linetype = "none")

## Plotting the species level PCoA - Supplementary figure 3
PCoA_sp<-ggarrange(resplot, migplot, nrow = 2)

#ggsave(plot = PCoA_sp, "FigS3.jpeg", width = 6, height = 6)

#########

#PERMANOVA for each species
perm_structure_res <- how(blocks = bird_diet %>% filter(mig == "HE_res") %>% .$SPECIES, nperm = 999)

permres<-adonis2(bird_diet %>% filter(mig == "HE_res") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                   bird_diet %>% filter(mig == "HE_res") %>% .$season , method = "jaccard", permutations = perm_structure_res)


perm_structure_mig <- how(blocks = bird_diet %>% filter(mig == "mig") %>% .$SPECIES, nperm = 999)
permmig<-adonis2(bird_diet %>% filter(mig == "mig") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                   bird_diet %>% filter(mig == "mig") %>% .$season , method = "jaccard", permutations = perm_structure_mig)



permBTFU<-adonis2(bird_diet %>% filter(SPECIES == "Brown-throated Fulvetta") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Brown-throated Fulvetta") %>% .$season , method = "jaccard")


permCHTE<-adonis2(bird_diet %>% filter(SPECIES == "Chestnut-headed Tesia") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Chestnut-headed Tesia") %>% .$season , method = "jaccard")

permRCBA<-adonis2(bird_diet %>% filter(SPECIES == "Rufous-capped Babbler") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Rufous-capped Babbler") %>% .$season , method = "jaccard")

permRBNI<-adonis2(bird_diet %>% filter(SPECIES == "Rufous-bellied Niltava") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Rufous-bellied Niltava") %>% .$season , method = "jaccard")

permRVYU<-adonis2(bird_diet %>% filter(SPECIES == "Rufous-vented Yuhina")%>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Rufous-vented Yuhina") %>% .$season , method = "jaccard")

permRGFL<-adonis2(bird_diet %>% filter(SPECIES == "Rufous-gorgeted Flycatcher") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Rufous-gorgeted Flycatcher") %>% .$season , method = "jaccard")

permSBSB<-adonis2(bird_diet %>% filter(SPECIES == "Streak-breasted Scimitar-Babbler")%>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Streak-breasted Scimitar-Babbler") %>% .$season , method = "jaccard")

permRWFU<-adonis2(bird_diet %>% filter(SPECIES == "Rufous-winged Fulvetta") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Rufous-winged Fulvetta") %>% .$season , method = "jaccard")

permSTYU<-adonis2(bird_diet %>% filter(SPECIES == "Stripe-throated Yuhina") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Stripe-throated Yuhina") %>% .$season , method = "jaccard")

permWHWA<-adonis2(bird_diet %>% filter(SPECIES == "Whistler's Warbler") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Whistler's Warbler") %>% .$season , method = "jaccard")

permBBWA<-adonis2(bird_diet %>% filter(SPECIES == "Buff-barred Warbler") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)) ~
                    bird_diet %>% filter(SPECIES == "Buff-barred Warbler") %>% .$season , method = "jaccard")


#######
#PERMDISP for each species

perdismig<-betadisper(vegdist(bird_diet %>% filter(mig == "mig") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(mig == "mig") %>% .$season, type = "centroid")

permutest(perdismig, pairwise = TRUE, permutations = perm_structure_mig)

perdisres<-betadisper(vegdist(bird_diet %>% filter(mig == "HE_res") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(mig == "HE_res") %>% .$season, type = "centroid")

permutest(perdisres, pairwise = TRUE, permutations = perm_structure_res)

perdisBTFU<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Brown-throated Fulvetta") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Brown-throated Fulvetta") %>% .$season, type = "centroid")
anova(perdisBTFU)

perdisRCBA<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Rufous-capped Babbler") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Rufous-capped Babbler") %>% .$season, type = "centroid",)
anova(perdisRCBA)

perdisRVYU<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Rufous-vented Yuhina") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Rufous-vented Yuhina") %>% .$season, type = "centroid")
anova(perdisRVYU)

perdisSBSB<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Streak-breasted Scimitar-Babbler") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Streak-breasted Scimitar-Babbler") %>% .$season, type = "centroid")
anova(perdisSBSB)

perdisSTYU<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Stripe-throated Yuhina")%>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Stripe-throated Yuhina") %>% .$season, type = "centroid")
anova(perdisSTYU)

perdisCHTE<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Chestnut-headed Tesia") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Chestnut-headed Tesia") %>% .$season, type = "centroid")
anova(perdisCHTE)

perdisRBNI<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Rufous-bellied Niltava") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Rufous-bellied Niltava") %>% .$season, type = "centroid")
anova(perdisRBNI)

perdisRGFL<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Rufous-gorgeted Flycatcher") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Rufous-gorgeted Flycatcher") %>% .$season, type = "centroid")
anova(perdisRGFL)

perdisRWFU<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Rufous-winged Fulvetta") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Rufous-winged Fulvetta") %>% .$season, type = "centroid")
anova(perdisRWFU)

perdisWHWA<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Whistler's Warbler") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Whistler's Warbler") %>% .$season, type = "centroid")
anova(perdisWHWA)

perdisBBWA<-betadisper(vegdist(bird_diet %>% filter(SPECIES == "Buff-barred Warbler") %>% select(-c(ID, DATE.x.x, INDID, Time.x, Time.y, Species.x, Species.y, DATE.y.x, DATE.x.y, DATE.y.y, season, mig, SPECIES, migsea, ELEV,specsea, DAY.x, DAY.y, MONTH.x, MONTH.y, YEAR.x, YEAR.y)), method = "jaccard", binary = T), bird_diet %>% filter(SPECIES == "Buff-barred Warbler") %>% .$season, type = "centroid")
anova(perdisBBWA)
