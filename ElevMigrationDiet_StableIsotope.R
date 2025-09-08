## These scripts are to visualize and analyse the Stable Isotope ratios of Carbon and Nitrogen from whole bird blood which are reported in the manuscript

library(tidyverse)
library(scales)
library(SIBER)
library(tRophicPosition)
library(nicheROVER)
library(ggpubr)
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

#load data linking date and ring number to isotope data
ring_sp<-read.csv("ring_numbers.csv") %>% select(-DATE)

#Loading guild data and merging
guild<-read.csv("E:/IISC/Global_Change_Lab/Data/DNA_metabarcoding/guilds.csv")
bird_data<-left_join(ring_sp,guild)

birdiso<-left_join(bird,bird_data, relationship = "many-to-many")

#renaming variables for ease of plotting
birdiso <- birdiso %>% mutate(mig = recode(mig,  "HE_res" = "Resident"))%>% mutate(mig = recode(mig, "mig" = "Migrant")) %>% filter (mig != "LE_res")

#Categorising elevation for Trophic Position calculation
birdiso$elevcat <- ifelse(birdiso$ELEV > 1900, "high", ifelse(birdiso$ELEV > 500 & birdiso$ELEV < 1400, "low", NA))


plant<-plant %>% separate_wider_delim(Sample.name, delim = "-", names = c("ELEV", "sample"), too_many = "merge")
plant$ELEV<-as.numeric(plant$ELEV)

#Categorising elevation for Trophic Position calculation
plant$elevcat<-ifelse(plant$ELEV > 2000, "high", ifelse(plant$ELEV > 500 & plant$ELEV < 1400, "low", NA))


### Estimating Trophic Position

## Putting bird and baseline (plant) isotope values in the same dataframe that can be analysed using the package tRophicPosition
plant_base<-plant %>% select(c(sample,elevcat,d13C,d15N))
plant_base$FG<-"Plant"
birdiso_all<-birdiso %>% mutate(sample = paste(mig, "_",Season, sep="")) %>% select(c(sample,elevcat,d13C,d15N))
birdiso_all$FG<-"Bird"
all_pl<-rbind(birdiso_all,plant_base)
all_pl<-as.data.frame(all_pl)

# Formatting the data for tRophicPosition
# Trophic Discrimination Factor was simulated using a mean of 2.3 and SD of 0.3
all_list<-extractIsotopeData(all_pl, b1 = "Plant", baselineColumn = "FG", consumersColumn = "sample", groupsColumn = "elevcat", d13C = "d13C", d15N = "d15N", deltaN = simulateTDF(nN = length(all_pl$sample), meanN = 2.3, sdN = 0.3, seed = 3))

# Estimating trophic position
all_models <- multiSpeciesTP(all_list, lambda = 1,
                             n.adapt = 10000, n.iter = 10000,
                             burnin = 10000, print = FALSE)

# Putting the results in a dataframe be platted later
plot_all<-all_models$df %>% separate_wider_delim(consumer, delim = "_", names = c("mig", "Season"))
plot_all$mig <- factor(plot_all$mig, levels=c("Resident", "Migrant"))
plot_all$Species<-"All Species"

# Calculating the probability that resident have a higher trophic position in the summer compared to the winter
pres<-compareTwoDistributions(all_models$TPs$high.Resident_Summer.1b, all_models$TPs$high.Resident_Winter.1b, test = ">")

# Calculating the probability that migrant have a higher trophic position in the winter compared to the summer
pmig<-compareTwoDistributions(all_models$TPs$high.Migrant_Summer.1b, all_models$TPs$low.Migrant_Winter.1b, test = "<")

# Calculating the probability that residents had a higher summer trophic position than migrants
pres2<-compareTwoDistributions(all_models$TPs$high.Migrant_Summer.1b, all_models$TPs$high.Resident_Summer.1b, test = "<")


###############
## Estimating trophic positions of Individual Species

## Resident species
## Filter out species with not enough data
birdiso_res<-birdiso %>% filter(mig == "Resident") %>%  filter(SPECIES %ni% c("Black-browed Tit", "Black-faced Laughingthrush")) %>% mutate(sample = paste(SPECIES, "_",Season, sep="")) %>% select(c(sample,elevcat,d13C,d15N))
birdiso_res$FG<-"Bird"

## Putting bird and baseline (plant) isotope values in the same dataframe that can be analysed using the package tRophicPosition
res_pl<-rbind(birdiso_res,plant_base)
res_pl<-as.data.frame(res_pl)

# Formatting the data for tRophicPosition
# Trophic Discrimination Factor was simulated using a mean of 2.3 and SD of 0.3
res_list<-extractIsotopeData(res_pl, b1 = "Plant", baselineColumn = "FG", consumersColumn = "sample", groupsColumn = "elevcat", d13C = "d13C", d15N = "d15N", deltaN = simulateTDF(nN = length(res_pl$sample), meanN = 2.3, sdN = 0.3, seed = 3))

# Estimating Trophic Position
res_models <- multiSpeciesTP(res_list, lambda = 1,
                             n.adapt = 10000, n.iter = 10000,
                             burnin = 10000, print = FALSE)

# Putting the results in a dataframe be platted later
plot_res<-res_models$df %>% separate_wider_delim(consumer, delim = "_", names = c("Species", "Season"))
plot_res$mig<-"Resident"

# Calculating the probability that each resident species have a higher trophic position in the summer compared to the winter
pBTFU<-compareTwoDistributions(res_models$TPs$high.Brown.throated.Fulvetta_Winter.1b, res_models$TPs$high.Brown.throated.Fulvetta_Summer.1b, test = "<")
pRCBA<-compareTwoDistributions(res_models$TPs$high.Rufous.capped.Babbler_Winter.1b,res_models$TPs$high.Rufous.capped.Babbler_Summer.1b, test = "<")
pRVYU<-compareTwoDistributions(res_models$TPs$high.Rufous.vented.Yuhina_Winter.1b,res_models$TPs$high.Rufous.vented.Yuhina_Summer.1b, test = "<")
pSBSB<-compareTwoDistributions(res_models$TPs$high.Streak.breasted.Scimitar.Babbler_Winter.1b,res_models$TPs$high.Streak.breasted.Scimitar.Babbler_Summer.1b, test = "<")
pSTYU<-compareTwoDistributions(res_models$TPs$high.Stripe.throated.Yuhina_Winter.1b,res_models$TPs$high.Stripe.throated.Yuhina_Summer.1b, test = "<")

## migrant species
birdiso_mig<-birdiso %>% filter(mig == "Migrant") %>% mutate(sample = paste(SPECIES, "_",Season, sep="")) %>% select(c(sample,elevcat,d13C,d15N))
birdiso_mig$FG<-"Bird"

## Putting bird and baseline (plant) isotope values in the same dataframe that can be analysed using the package tRophicPosition
mig_pl<-rbind(birdiso_mig,plant_base)
mig_pl<-as.data.frame(mig_pl)

# Formatting the data for tRophicPosition
# Trophic Discrimination Factor was simulated using a mean of 2.3 and SD of 0.3
mig_list<-extractIsotopeData(mig_pl, b1 = "Plant", baselineColumn = "FG", consumersColumn = "sample", groupsColumn = "elevcat", d13C = "d13C", d15N = "d15N", deltaN = simulateTDF(nN = length(mig_pl$sample), meanN = 2.3, sdN = 0.3, seed = 3))

# Estimating Trophic Position
mig_models <- multiSpeciesTP(mig_list, lambda = 1,
                             n.adapt = 10000, n.iter = 10000,
                             burnin = 10000, print = FALSE)

# Putting the results in a dataframe be platted later
plot_mig<-mig_models$df %>% separate_wider_delim(consumer, delim = "_", names = c("Species", "Season"))
plot_mig$mig<-"Migrant"

# Calculating the probability that each migrant species have a higher trophic position in the summer compared to the winter
pBBWA<-compareTwoDistributions(mig_models$TPs$low.Buff.barred.Warbler_Winter.1b, mig_models$TPs$high.Buff.barred.Warbler_Summer.1b, test = ">")
pCHTE<-compareTwoDistributions(mig_models$TPs$low.Chestnut.headed.Tesia_Winter.1b,mig_models$TPs$high.Chestnut.headed.Tesia_Summer.1b, test = ">")
pRBNI<-compareTwoDistributions(mig_models$TPs$low.Rufous.bellied.Niltava_Winter.1b,mig_models$TPs$high.Rufous.bellied.Niltava_Summer.1b, test = ">")
pRGFL<-compareTwoDistributions(mig_models$TPs$low.Rufous.gorgeted.Flycatcher_Winter.1b,mig_models$TPs$high.Rufous.gorgeted.Flycatcher_Summer.1b, test = ">")
pRWFU<-compareTwoDistributions(mig_models$TPs$low.Rufous.winged.Fulvetta_Winter.1b,mig_models$TPs$high.Rufous.winged.Fulvetta_Summer.1b, test = ">")
pWHWA<-compareTwoDistributions(mig_models$TPs$low.Whistler.s.Warbler_Winter.1b,mig_models$TPs$high.Whistler.s.Warbler_Summer.1b, test = ">")

# putting all result dataframes together for plotting
plot_3<-rbind(plot_res,plot_mig, plot_all)
plot_3$mig<-factor(plot_3$mig, levels = c("Resident", "Migrant"))

# Creating seasonal trophic position plot - Figure 2
TP_Bay_all <- plot_3 %>% 
  ggplot(aes(Species, median, colour = Season)) + 
  geom_point(size = 3,position = position_dodge(width = 0.25)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.25), width = 0) + 
  ylab("Posterior Trophic Position") +
  theme(axis.title = element_text(size =20), legend.title = element_text(size =15), legend.text = element_text(size =12), axis.text = element_text(size =12))+ 
  scale_x_discrete(labels = wrap_format(10)) + 
  scale_colour_grey(start = 0,
                    end = 0.7)+
  theme_bw(base_size = 12)+
  facet_wrap(~mig,  ncol=1, scales = "free")

#ggsave(plot = TP_Bay_all, "Figure2.jpeg", width = 10, height = 7)


##########
# Visualising Isotopic Niche
# Creating biplots of baseline corrected δ15N and δ13C

# Formatting data for plotting
birdiso$elevcat <- ifelse(birdiso$ELEV > 2200, "high", ifelse(birdiso$ELEV > 1500 & birdiso$ELEV < 2300, "mid", ifelse(birdiso$ELEV > 500 & birdiso$ELEV < 1400, "low", NA)))

plant$elevcat<-ifelse(plant$ELEV > 2200, "high", ifelse(plant$ELEV > 1500 & plant$ELEV < 2300, "mid", ifelse(plant$ELEV > 500 & plant$ELEV < 1400, "low", NA)))

baselines<-plant %>% group_by(elevcat) %>% summarise(d15N_base = mean(d15N), d13C_base = mean(d13C), d13C_min = min(d13C), d13C_max = max(d13C), SE_N = sd(d15N)/sqrt(n()), SE_C = sd(d13C)/sqrt(n()))

birdiso_base<-left_join(birdiso,baselines) %>% select(c(INDID,SPECIES,ELEV,elevcat,Season,d13C,d15N,mig, d15N_base, d13C_base, d13C_min, d13C_max ))

birdiso_base$TP<-1 + (birdiso_base$d15N - birdiso_base$d15N_base)/2.3
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
  facet_wrap(~SPECIES,  nrow=1) +
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
  ylab("Trophic Position")+
  xlab(expression(paste(delta^{13}, "Ccorr"))) +
  #xlim(c(0.3,1))+ylim(c(2,5.3))+
  theme_bw()+
  theme(text = element_text(size=12)) + 
  facet_wrap(~SPECIES,  nrow=1) +
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
std_ellipse_all<-ggarrange(res_ellipse_plot,mig_ellipse_plot, nrow = 2, legend = "bottom", common.legend = T )+
  guides (size = 15)

#ggsave(plot = std_ellipse_all, "Figure3.jpeg", width = 14, height = 7)


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

