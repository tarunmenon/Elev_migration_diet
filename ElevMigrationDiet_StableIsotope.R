## These scripts are to visualize and analyse the Stable Isotope ratios of Carbon and Nitrogen from whole bird blood which are reported in the manuscript

library(tidyverse)
library(scales)
library(SIBER)
library(nicheROVER)
library(ggpubr)
library(brms)
library(ape)
library(marginaleffects)
library(purrr)
library(broom.mixed)
`%ni%` <- Negate(`%in%`)
cbbPalette <- c("#D55E00", "#56B4E9","#000000", "#CC79A7","#009E73", "#E69F00", "#4B0092", "#0072B2","#DC3220" )

setwd("~/Elev_migration_diet")

# Load bird and baseline plant isotope data
bird<-read.csv("bird_iso.csv")
plant<-read.csv("plant_iso.csv")

## cleaning up bird isotope data, typos in ring numbers and dates
bird<-bird %>% mutate(Samples= recode(Samples, "A3874_121223"="A3874_171223", "A3119_161223" = "A3119_161222", "A8211_14224" = "A8211_14124", "A2862_18522" = "A2682_18522", "Z3626_21523" = "Z3826_21523", "AB1449_5523" = "AB1445_5523", "A1426_27423" = "AB1426_27423", "Z2772_16122" = "Z2772_16123", "A5004_12223"= "A5104_121223", "A0348_8522" = "A0384_8522", "A0417_18422" = "A0417_16422", "A2824_22522" = "A2834_22522", "A8239_22522" = "A2839_22522", "Z2852_23522" = "A2852_23522", "Z2528_17523" = "Z2538_17523")) %>% 
  filter(Samples %ni% c("A8261_20124",  "Z3767_16423"))

# remove repeat samples which were used for testing
bird<-bird %>% separate_wider_delim(Samples, delim = "_", names = c("INDID", "Date","repeat"), too_few = "align_start") %>% filter(is.na(`repeat`))

# converting date so it can be matched with the ring/date data

bird <- bird %>%
  mutate(
    Length_Condition = nchar(Date),
    
    # Extract the first and second parts for numeric comparison
    First_Part_Value = as.numeric(substr(Date, 1, 2)),
    Second_Part_Value = if_else(Length_Condition == 5, as.numeric(substr(Date, 2, 3)), NA_real_),
    
    # Determine if the second part starts with '0'
    Starts_With_Zero = if_else(Length_Condition == 5 & substr(Date, 2, 2) == "0", TRUE, FALSE),
    
    # Determine if the first part value is under 32 and the second part is under 12
    First_Part_Under_32_Second_Under_12 = First_Part_Value < 32 & (Length_Condition == 5 & Second_Part_Value < 12),
    
    DAY = case_when(
      Length_Condition == 4 ~ substr(Date, 1, 1),
      Length_Condition == 6 ~ substr(Date, 1, 2),
      Length_Condition == 5  & (Second_Part_Value > 12 | Starts_With_Zero|First_Part_Under_32_Second_Under_12) ~ substr(Date, 1, 2),
      TRUE ~ substr(Date, 1, 1)
    ),
    
    MONTH = case_when(
      Length_Condition == 4 ~ substr(Date, 2, 2),
      Length_Condition == 6 ~ substr(Date, 3, 4),
      Length_Condition == 5 & (Second_Part_Value > 12 |Starts_With_Zero|First_Part_Under_32_Second_Under_12) ~ substr(Date, 3, 3),
      TRUE ~ substr(Date, 2, 3)
    ),
    
    YR = case_when(
      Length_Condition == 4 ~ substr(Date, 3, 4),
      Length_Condition == 6 ~ substr(Date, 5, 6),
      Length_Condition == 5 & (Second_Part_Value > 12 | Starts_With_Zero|First_Part_Under_32_Second_Under_12) ~ substr(Date, 4, 5),
      TRUE ~ substr(Date, 4, 5)
    )
  )%>%
  select(INDID, DAY, MONTH, YR, d13C, d15N )

bird$YEAR<-ifelse(bird$YR == 24, 2024, ifelse(bird$YR == 23, 2023, 2022))
bird$DAY<-as.numeric(bird$DAY)
bird$MONTH<-as.numeric(bird$MONTH)
bird$Season<-ifelse(bird$MONTH == 4 | bird$MONTH == 5, "Summer", "Winter") #Defining seasons

setwd("E:/IISC/Global_Change_Lab/Data")

#load data linking date and ring number to isotope data
ring_sp<-read.csv("ring_number_sci.csv")

#Loading guild data and merging
guild<-read.csv("guilds.csv")
bird_data<-left_join(ring_sp,guild)

birdiso<-left_join(bird,bird_data)

#renaming variables for ease of plotting
birdiso <- birdiso %>% mutate(mig = recode(mig,  "HE_res" = "Resident"))%>% mutate(mig = recode(mig, "mig" = "Migrant")) %>% filter (mig != "LE_res")

plant<-plant %>% separate_wider_delim(Sample.name, delim = "-", names = c("ELEV", "sample"), too_many = "merge")
plant$ELEV<-as.numeric(plant$ELEV)


### Estimating Trophic Position


#######################
#Estimating trophic positions of all Species in a single model

birdiso$elevcat <- ifelse(birdiso$ELEV > 2200, "high", ifelse(birdiso$ELEV > 1500 & birdiso$ELEV < 2300, "mid", ifelse(birdiso$ELEV > 500 & birdiso$ELEV < 1400, "low", NA)))

plant$elevcat<-ifelse(plant$ELEV > 2200, "high", ifelse(plant$ELEV > 1500 & plant$ELEV < 2300, "mid", ifelse(plant$ELEV > 500 & plant$ELEV < 1400, "low", NA)))

baselines<-plant %>% group_by(elevcat) %>% summarise(d15N_base = mean(d15N), d13C_base = mean(d13C), d13C_min = min(d13C), d13C_max = max(d13C), SE_N = sd(d15N)/sqrt(n()), SE_C = sd(d13C)/sqrt(n()))

birdiso_base<-left_join(birdiso,baselines)
birdiso_base<- birdiso_base %>%
  mutate(SCI_NAME = str_replace_all(Species, " ", "_"))
birdiso_base$TP<-1 + (birdiso_base$d15N - birdiso_base$d15N_base)/2.3

# first-order Taylor series expansion for error propagation
birdiso_base$TP_se <- sqrt( (birdiso_base$SE_N^2 / 2.3^2) + 
                              (((birdiso_base$d15N - birdiso_base$d15N_base)^2 * 0.3^2) / (2.3^4)) )

# get global bird phylogentic tree from Claramunt et al. (2026)
github.directory <- "https://raw.githubusercontent.com/evolucionario/BigBirdTree/refs/heads/main/"

stage <- "RAGBackbone/"

tree <- "BBtreeC2022.tre"

url <- paste0(github.directory, stage, tree)

BBtree2 <- read.tree(url)

# Replacing species in our dataset that are absent in the tree by replacing it with the next closely related species
BBtree2$tip.label[BBtree2$tip.label == "Phylloscopus_burkii"] <- "Phylloscopus_whistleri"
BBtree2$tip.label[BBtree2$tip.label == "Phylloscopus_inornatus"] <- "Phylloscopus_pulcher"

tree_pruned <- keep.tip(BBtree2, birdiso_base$SCI_NAME)
birdiso_base$phylo_id <- birdiso_base$SCI_NAME
phylo_matrix <- vcv.phylo(tree_pruned, corr = TRUE)

# running the bayesian multilevel phylogenetic regression
fit_tp <- brm(
  formula = TP | mi(TP_se) ~ mig*Season + (1 | gr(phylo_id, cov = phylo)),
  data = birdiso_base,
  data2 = list(phylo = phylo_matrix),
  family = gaussian(),
  prior = c(
    prior(normal(3, 1), class = "Intercept"),
    prior(normal(0, 0.5), class = "b")
  ),
  chains = 4, iter = 4000, cores = 4,
  control = list(adapt_delta = 0.99)
)

summary(fit_tp)

comparisons(fit_tp, variables = "Season", by = "mig")


#####################
## Species specific models

sep_species_preds <- birdiso_base %>% filter(SPECIES %in% c("Brown-throated Fulvetta", "Rufous-vented Yuhina", "Streak-breasted Scimitar-Babbler", "Stripe-throated Yuhina", "Rufous-capped Babbler", "Whistler's Warbler", "Rufous-gorgeted Flycatcher", "Chestnut-headed Tesia", "Rufous-winged Fulvetta", "Rufous-bellied Niltava", "Buff-barred Warbler")) %>%
  group_by(SPECIES, mig) %>%
  nest() %>%
  mutate(model = map(data, ~brm(TP ~ Season, data = .x, 
                                family = gaussian(), chains = 4, iter = 4000, cores = 4,
                                prior = c(
                                  prior(normal(3, 1), class = "Intercept"),
                                  prior(normal(0, 0.5), class = "b")),
                                control = list(adapt_delta = 0.99)))) 

species_stats <- sep_species_preds %>%
  mutate(tidy_model = map(model, ~tidy(.x, effects = "fixed", conf.int = TRUE, rhat = TRUE, ess = TRUE))) %>%
  unnest(tidy_model)

# Creating model output table
final_table <- species_stats %>%
  mutate(
    Parameter = case_when(
      term == "(Intercept)" ~ "Summer TP (Intercept)",
      term == "SeasonWinter" ~ "Winter Shift (Beta)",
      TRUE ~ term
    ),
    `95% CI` = paste0("[", round(conf.low, 2), ", ", round(conf.high, 2), "]"),
    Estimate = round(estimate, 2),
    Rhat = round(rhat, 3),
    ESS = round(ess, 0)
  ) %>%
  select(SPECIES, mig, Parameter, Estimate, `95% CI`, Rhat, ESS)

# Export for supplement
#write.csv(final_table, "Table_S3.csv", row.names = FALSE)

############
## Visualization (for Figure 1)

all_species_cond <- predictions(
  fit_tp, 
  newdata = datagrid(mig = unique, Season = unique), 
  re_formula = NULL
) %>%
  as.data.frame() %>%
  mutate(Species = "All Species") %>%
  rename(Mean_TP = estimate, Lower = conf.low, Upper = conf.high) %>%
  select(Species, mig, Season, Mean_TP, Lower, Upper)

sep_species_cond <- sep_species_preds %>%
  mutate(preds = map(model, function(m) {
    # Use datagrid to ensure it is a conditional prediction for that species
    res <- predictions(m, newdata = datagrid(Season = c("Summer", "Winter"))) %>% 
      as.data.frame()
    return(res)
  })) %>%
  unnest(preds) %>%
  mutate(Species = SPECIES) %>%
  rename(Mean_TP = estimate, Lower = conf.low, Upper = conf.high) %>%
  select(Species, SPECIES, mig, Season, Mean_TP, Lower, Upper)

ss<-birdiso_base %>% group_by(SPECIES, Season, mig) %>% summarise(count = n()) %>% rename(Species = SPECIES)
ss_all<-birdiso_base %>% group_by(mig, Season) %>% summarise(count = n())
ss_all$Species<- "All Species"
ss<- rbind(ss,ss_all)

plot_data_final <- bind_rows(all_species_cond, sep_species_cond) %>%
  mutate(
    Species = factor(Species, levels = c("All Species", sort(unique(sep_species_cond$Species)))),
    mig = factor(mig, levels = c("Resident", "Migrant"))
  )%>%
  left_join(ss, by = c("Species", "Season", "mig"))

TP_Bay_all<-ggplot(plot_data_final, aes(x = Species, y = Mean_TP, color = Season)) +
  geom_point(size = 3,position = position_dodge(width = 0.25)) + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge(width = 0.25), width = 0) + 
  ylab("Estimated Trophic Position") +
  geom_text(aes(label = count, group = Season), 
            position = position_dodge(width = 0.75), 
            vjust = -3, size = 3.5, color = "black")+
  theme_bw()+
  theme(axis.title = element_text(size =15), legend.title = element_text(size =20), legend.text = element_text(size =15), axis.text = element_text(size =11), strip.text = element_text(size = 15))+ 
  scale_x_discrete(labels = wrap_format(10)) +
  facet_wrap(~mig,  ncol=1, scales = "free") +
  scale_colour_manual(values = cbbPalette) + 
  ggtitle("Figure 1")
TP_Bay_all

#ggsave(plot = TP_Bay_all, "Figure1.tif", device = "tiff", dpi = 600, compression = "lzw", width = 11, height = 7)


##########
# Visualising Isotopic Niche
# Creating biplots of baseline corrected δ15N and δ13C

# Formatting data for plotting

birdiso_base$d13Ccorr<- (birdiso_base$d13C - birdiso_base$d13C_base)/(birdiso_base$d13C_max-birdiso_base$d13C_min)

birdiso_base$mig<-factor( birdiso_base$mig, levels = c("Resident", "Migrant"))

## Plotting biplots with standard ellipses
p.ell <- pchisq(1,2) #SEAc ellipse parameter

#removing species without enough data
res_ellipse <- birdiso_base %>% filter(mig == "Resident") %>% filter(SPECIES != "Black-faced Laughingthrush") %>% filter(SPECIES != "Black-browed Tit")

# Plot of all residents and individual species
all_res<-res_ellipse
all_res$SPECIES<-"All Species"
res_ellipse<-rbind(res_ellipse,all_res)

res_ellipse_plot<- res_ellipse %>% 
  ggplot(aes(d13Ccorr,TP)) +
  geom_point(aes(color = Season))+
  ylab("Trophic Position")+
  xlab(expression(paste(delta^{13}, "Ccorr"))) + 
  #xlim(c(0.3,1))+ylim(c(2,5.3))+
  theme_bw()+
  theme(text = element_text(size=12)) + 
  facet_wrap(~SPECIES,  ncol=1) +
  ggtitle("Residents") +
  theme(plot.title = element_text(hjust = 0.5, size = 15))+
  scale_fill_manual(values = cbbPalette)+
  stat_ellipse(aes(group = Season, 
                   fill = Season, 
                   color = Season), 
               alpha = 0.25, 
               level = p.ell,
               type = "norm",
               geom = "polygon")+
  theme(axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15))

# Plot of all migrants and individual species
mig_ellipse <- birdiso_base %>% filter(mig == "Migrant")
all_mig<-mig_ellipse
all_mig$SPECIES<-"All Species"
mig_ellipse<-rbind(mig_ellipse,all_mig)

mig_ellipse_plot <- mig_ellipse %>% 
  ggplot(aes(d13Ccorr,TP)) +
  geom_point(aes(color = Season))+
  labs(y = NULL)+
  xlab(expression(paste(delta^{13}, "Ccorr"))) +
  #xlim(c(0.3,1))+ylim(c(2,5.3))+
  theme_bw()+
  theme(text = element_text(size=12)) + 
  facet_wrap(~SPECIES,  ncol=1) +
  ggtitle("Elevation Migrants") +
  theme(plot.title = element_text(hjust = 0.5, size = 15))+
  scale_fill_manual(values = cbbPalette)+
  stat_ellipse(aes(group = Season, 
                   fill = Season, 
                   color = Season), 
               alpha = 0.25, 
               level = p.ell,
               type = "norm",
               geom = "polygon")+
  theme(axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15))

## Putting resident and migrant plots together
std_ellipse_all<-ggarrange(res_ellipse_plot,mig_ellipse_plot, ncol = 2, legend = "bottom", common.legend = T )+
  guides (size = 15)

std_ellipse_all<-annotate_figure(std_ellipse_all,
                top = text_grob("Figure 2",x = 0, hjust = -0.1))

std_ellipse_all<- std_ellipse_all+theme(plot.background = element_rect(fill = "white", color = NA))

library(ragg)
#ggsave(plot = std_ellipse_all, "Figure2.tif", device = ragg::agg_tiff, dpi = 600, compression = "lzw" ,width = 5.5, height = 11)


#### Bayesian analysis of seasonal niche overlap using NicheROVER

# At the community level
# residents
resident<-birdiso_base %>% filter(mig == "Resident")  %>%  select(c(d13Ccorr,TP,Season,mig))
colnames(resident) <- c("iso1", "iso2", "group", "community")
resident <- as.data.frame(resident) 

# Estimating 1000 posterior overlap values where niche size had 95% coverage
res.par<-tapply(1:nrow(resident), resident$group,
                function(ii) niw.post(nsamples = 1000, X = resident[ii,1:2]))

res_over.stat <- overlap(res.par, nreps = 1e3, nprob = 1e3, alpha = .95)

#The median overlap metrics calculated across iterations 
res_over.median <- apply(res_over.stat, c(1, 2), median)*100
round(res_over.median, 2)


# migrants
migrant<-birdiso_base %>% filter(mig == "Migrant") %>%  select(c(d13Ccorr,TP,Season,mig))
colnames(migrant) <- c("iso1", "iso2", "group", "community")
migrant <- as.data.frame(migrant)  
#

# Estimating 1000 posterior overlap values where niche size had 95% coverage
mig.par<-tapply(1:nrow(migrant), migrant$group,
                function(ii) niw.post(nsamples = 1000, X = migrant[ii,1:2]))

mig_over.stat <- overlap(mig.par, nreps = 1e3, nprob = 1e3, alpha = 0.95)

#The median overlap metrics calculated across iterations
mig_over.median <- apply(mig_over.stat, c(1,2), median)*100
round(mig_over.median, 2)

#################

#At the species level
# Same steps as the community level:  Estimating 1000 posterior overlap values where niche size had 95% coverage & The median overlap metrics calculated across iterations

# Brown-throated Fulvetta
BTFU<-birdiso_base %>% filter(SPECIES == "Brown-throated Fulvetta") %>%  select(c(d13Ccorr,TP,Season,mig))
colnames(BTFU) <- c("iso1", "iso2", "group", "community")
BTFU <- as.data.frame(BTFU)  

BTFU.par<-tapply(1:nrow(BTFU), BTFU$group,
                 function(ii) niw.post(nsamples = 1000, X = BTFU[ii,1:2]))

BTFU_over.stat <- overlap(BTFU.par, nreps = 1e3, nprob = 1e3, alpha = 0.95)

BTFU_over.median <- apply(BTFU_over.stat, c(1,2), median)*100
round(BTFU_over.median, 2)
BTFU_over.median<-as.data.frame(BTFU_over.median)

BTFU_over.median$Species<-"Brown-thorated Fulvetta"

# Rufous-capped Babbler
RCBA<-birdiso_base %>% filter(SPECIES == "Rufous-capped Babbler") %>%  select(c(d13Ccorr,TP,Season,mig))
colnames(RCBA) <- c("iso1", "iso2", "group", "community")
RCBA <- as.data.frame(RCBA)  

RCBA.par<-tapply(1:nrow(RCBA), RCBA$group,
                 function(ii) niw.post(nsamples = 1000, X = RCBA[ii,1:2]))

RCBA_over.stat <- overlap(RCBA.par, nreps = 1e3, nprob = 1e3, alpha = 0.95)

#The median overlap metrics calculated across iteratations for both niche 
#region sizes (alpha = .95 and alpha = .99) can be calculated and displayed in an array.
RCBA_over.median <- apply(RCBA_over.stat, c(1,2), median)*100
round(RCBA_over.median, 2)
RCBA_over.median<-as.data.frame(RCBA_over.median)

RCBA_over.median$Species <- "Rufous-capped Babbler"

# Rufous-vented Yuhina
RVYU<-birdiso_base %>% filter(SPECIES == "Rufous-vented Yuhina") %>%  select(c(d13Ccorr,TP,Season,mig))
colnames(RVYU) <- c("iso1", "iso2", "group", "community")
RVYU <- as.data.frame(RVYU)  
#

RVYU.par<-tapply(1:nrow(RVYU), RVYU$group,
                 function(ii) niw.post(nsamples = 1000, X = RVYU[ii,1:2]))

RVYU_over.stat <- overlap(RVYU.par, nreps = 1e3, nprob = 1e3, alpha = 0.95)

RVYU_over.median <- apply(RVYU_over.stat, c(1,2), median)*100
round(RVYU_over.median, 2)
RVYU_over.median<-as.data.frame(RVYU_over.median)

RVYU_over.median$Species<-"Rufous-vented Yuhina"

# Streak-breasted Scimitar-Babbler
SBSB<-birdiso_base %>% filter(SPECIES == "Streak-breasted Scimitar-Babbler") %>%  select(c(d13Ccorr,TP,Season,mig))
colnames(SBSB) <- c("iso1", "iso2", "group", "community")
SBSB <- as.data.frame(SBSB)  

SBSB.par<-tapply(1:nrow(SBSB), SBSB$group,
                 function(ii) niw.post(nsamples = 1000, X = SBSB[ii,1:2]))

SBSB_over.stat <- overlap(SBSB.par, nreps = 1e3, nprob = 1e3, alpha = 0.95)

SBSB_over.median <- apply(SBSB_over.stat, c(1,2), median)*100
round(SBSB_over.median, 2)
SBSB_over.median<-as.data.frame(SBSB_over.median)

SBSB_over.median$Species <- "Streak-breasted Scimitar-Babbler"


# Stripe-throated Yuhina
STYU<-birdiso_base %>% filter(SPECIES == "Stripe-throated Yuhina") %>%  select(c(d13Ccorr,TP,Season,mig))
colnames(STYU) <- c("iso1", "iso2", "group", "community")
STYU <- as.data.frame(STYU)  
#

STYU.par<-tapply(1:nrow(STYU), STYU$group,
                 function(ii) niw.post(nsamples = 1000, X = STYU[ii,1:2]))

STYU_over.stat <- overlap(STYU.par, nreps = 1e3, nprob = 1e3, alpha = 0.95)

STYU_over.median <- apply(STYU_over.stat, c(1,2), median)*100
round(STYU_over.median, 2)
STYU_over.median<-as.data.frame(STYU_over.median)

STYU_over.median$Species<- "Stripe-throated Yuhina"

# Buff-barred Warbler
BBWA<-birdiso_base %>% filter(SPECIES == "Buff-barred Warbler") %>%  select(c(d13Ccorr,TP,Season,mig))
colnames(BBWA) <- c("iso1", "iso2", "group", "community")
BBWA <- as.data.frame(BBWA)  

BBWA.par<-tapply(1:nrow(BBWA), BBWA$group,
                 function(ii) niw.post(nsamples = 1000, X = BBWA[ii,1:2]))

BBWA_over.stat <- overlap(BBWA.par, nreps = 1e3, nprob = 1e3, alpha = 0.95)

BBWA_over.median <- apply(BBWA_over.stat, c(1,2), median)*100
round(BBWA_over.median, 2)
BBWA_over.median<-as.data.frame(BBWA_over.median)

BBWA_over.median$Species<- "Buff-barred Warbler"


# Chestnut-headed Tesia
CHTE<-birdiso_base %>% filter(SPECIES == "Chestnut-headed Tesia") %>%  select(c(d13Ccorr,TP,Season,mig))
colnames(CHTE) <- c("iso1", "iso2", "group", "community")
CHTE <- as.data.frame(CHTE)  

CHTE.par<-tapply(1:nrow(CHTE), CHTE$group,
                 function(ii) niw.post(nsamples = 1000, X = CHTE[ii,1:2]))

CHTE_over.stat <- overlap(CHTE.par, nreps = 1e3, nprob = 1e3, alpha = 0.95)

CHTE_over.median <- apply(CHTE_over.stat, c(1,2), median)*100
round(CHTE_over.median, 2)
CHTE_over.median<-as.data.frame(CHTE_over.median)

CHTE_over.median$Species<- "Chestnut-headed Tesia"

# Rufous-bellied Niltava
RBNI<-birdiso_base %>% filter(SPECIES == "Rufous-bellied Niltava") %>%  select(c(d13Ccorr,TP,Season,mig))
colnames(RBNI) <- c("iso1", "iso2", "group", "community")
RBNI <- as.data.frame(RBNI)  
#

RBNI.par<-tapply(1:nrow(RBNI), RBNI$group,
                 function(ii) niw.post(nsamples = 1000, X = RBNI[ii,1:2]))

RBNI_over.stat <- overlap(RBNI.par, nreps = 1e3, nprob = 1e3, alpha = 0.95)

RBNI_over.median <- apply(RBNI_over.stat, c(1,2), median)*100
round(RBNI_over.median, 2)
RBNI_over.median<- as.data.frame(RBNI_over.median)

RBNI_over.median$Species<- "Rufous-bellied Niltava"

# Rufous-gorgeted Flycatcher
RGFL<-birdiso_base %>% filter(SPECIES == "Rufous-gorgeted Flycatcher") %>%  select(c(d13Ccorr,TP,Season,mig))
colnames(RGFL) <- c("iso1", "iso2", "group", "community")
RGFL <- as.data.frame(RGFL)  
#

RGFL.par<-tapply(1:nrow(RGFL), RGFL$group,
                 function(ii) niw.post(nsamples = 1000, X = RGFL[ii,1:2]))

RGFL_over.stat <- overlap(RGFL.par, nreps = 1e3, nprob = 1e3, alpha = 0.95)

RGFL_over.median <- apply(RGFL_over.stat, c(1,2), median)*100
round(RGFL_over.median, 2)
RGFL_over.median<-as.data.frame(RGFL_over.median)

RGFL_over.median$Species<- "Rufous-gorgeted Flycatcher"


# Rufous-winged Fulvetta
RWFU<-birdiso_base %>% filter(SPECIES == "Rufous-winged Fulvetta") %>%  select(c(d13Ccorr,TP,Season,mig))
colnames(RWFU) <- c("iso1", "iso2", "group", "community")
RWFU <- as.data.frame(RWFU)  
#

RWFU.par<-tapply(1:nrow(RWFU), RWFU$group,
                 function(ii) niw.post(nsamples = 1000, X = RWFU[ii,1:2]))

RWFU_over.stat <- overlap(RWFU.par, nreps = 1e3, nprob = 1e3, alpha = 0.95)

RWFU_over.median <- apply(RWFU_over.stat, c(1,2), median)*100
round(RWFU_over.median, 2)
RWFU_over.median<- as.data.frame(RWFU_over.median)

RWFU_over.median$Species<-"Rufous-winged Fulvetta"


# Whistler's Warbler
WHWA<-birdiso_base %>% filter(SPECIES == "Whistler's Warbler") %>%  select(c(d13Ccorr,TP,Season,mig))
colnames(WHWA) <- c("iso1", "iso2", "group", "community")
WHWA <- as.data.frame(WHWA)  

WHWA.par<-tapply(1:nrow(WHWA), WHWA$group,
                 function(ii) niw.post(nsamples = 1000, X = WHWA[ii,1:2]))

WHWA_over.stat <- overlap(WHWA.par, nreps = 1e3, nprob = 1e3, alpha = 0.95)

WHWA_over.median <- apply(WHWA_over.stat, c(1,2), median)*100
round(WHWA_over.median, 2)
WHWA_over.median<-as.data.frame(WHWA_over.median)

WHWA_over.median$Species<- "Whistler's Warbler"


## Putting all the species overlap values together - Supplementary table 3
overlap_stat<-rbind(BTFU_over.median, RCBA_over.median, RVYU_over.median, SBSB_over.median, STYU_over.median,BBWA_over.median, CHTE_over.median, RBNI_over.median, RGFL_over.median, RWFU_over.median, WHWA_over.median)

