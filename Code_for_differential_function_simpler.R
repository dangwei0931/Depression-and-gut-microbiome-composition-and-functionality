# To analyze the differences of gut microbiome functional profiles between individuals with and without a history of depression
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(Maaslin2)

simpler <- readRDS("/proj/...path.../simpler_study_analysis.RDS") %>%
  select(id, depre_atfecal, depre_bydef, sex, agec, education, bmic, smoke, phy_act, alcohol, probiotics, fiberc, energyc, antibiotics_6mon, crohn_ibs, diabetes, cvd, depre_atfecal_op, dna_plate, dna_read, antid) 

simp2 <- simpler %>%
  filter(antid<3) %>%
  mutate(across(2:21, ~as.factor(.))) 

simp <- simpler %>%
  mutate(across(2:21, ~as.factor(.))) 

# functional data: SIMPLER
# to calculate the relative abundance of GMM and KEGG
simp_gmm <- read.table("/proj/...path.../simpler_metagenomics_gmm_modules_v2.0.tsv", header = TRUE, sep = "\t")
simp_kegg <- read.table("/proj/...path.../simpler_metagenomics_kegg_modules_v2.0.tsv", header = TRUE, sep = "\t")
simp_mgs <- read.table("/proj/...path.../simpler_metagenomics_mgs_relative_abundances_v2.0.tsv", header = TRUE, sep = "\t")

# GMM
mgss <- NULL
mgss <- colnames(simp_mgs)
i <- NULL
gmms <- simp %>% select(id)
for (i in 1:nrow(simp_gmm)) {
  print(i)
  temp <- NULL
  mgs_for_gmm <- NULL
  comm <- NULL
  gmm_name <- NULL
  gmm_name <- simp_gmm$gmm_id[i]
  mgs_for_gmm <- strsplit(simp_gmm$mgs_id[i],",")[[1]]
  comm <- intersect(mgs_for_gmm, mgss)
  temp <- simp_mgs %>%
    rename(id=SIMPKEY) %>%
    select(id, all_of(comm)) %>%
    mutate(sum_abundace=rowSums(across(all_of(comm)), na.rm=TRUE)) %>%
    rename(!!gmm_name :=sum_abundace) %>%
    select(id,!!gmm_name)
  
  gmms <- left_join(gmms, temp, by="id")
}

# KEGG
i <- NULL
keggs <- simp %>% select(id)
for (i in 1:nrow(simp_kegg)) {
  print(i)
  temp <- NULL
  mgs_for_kegg <- NULL
  comm <- NULL
  kegg_name <- NULL
  kegg_name <- simp_kegg$kegg_id[i]
  mgs_for_kegg <- strsplit(simp_kegg$mgs_id[i],",")[[1]]
  comm <- intersect(mgs_for_kegg, mgss)
  temp <- simp_mgs %>%
    rename(id=SIMPKEY) %>%
    select(id, all_of(comm)) %>%
    mutate(sum_abundace=rowSums(across(all_of(comm)), na.rm=TRUE)) %>%
    rename(!!kegg_name :=sum_abundace) %>%
    select(id,!!kegg_name)
  
  keggs <- left_join(keggs, temp, by="id")
}

# sort the tables by id
simp1 <- simp[order(simp$id),]
gmm <- gmms[order(gmms$id),]
kegg <- keggs[order(keggs$id),]

simp1 <- as.data.frame(simp1)
rownames(simp1) <- simp1$id
gmm <- as.data.frame(gmm)
rownames(gmm) <- gmm$id
kegg <- as.data.frame(kegg)
rownames(kegg) <- kegg$id
# remove the column of id
gmm <- gmm[, -c(1)]
kegg <- kegg[, -c(1)]
simp1 <- simp1[, -c(1)]

masterid <- simp2 %>% select(id)
gmm2 <- left_join(masterid, gmms, by="id")
kegg2 <- left_join(masterid, keggs, by="id")
# sort the tables by id
gmm2 <- gmm2[order(gmm2$id),]
kegg2 <- kegg2[order(kegg2$id),]
gmm2 <- as.data.frame(gmm2)
rownames(gmm2) <- gmm2$id
kegg2 <- as.data.frame(kegg2)
rownames(kegg2) <- kegg2$id
# remove the column of id
gmm2 <- gmm2[, -c(1)]
kegg2 <- kegg2[, -c(1)]

# to use the maaslin regression to conduct the differential abundance analysis
maaslin_tra<-function(func, exposure, dir){
  # Set the reference level for covariates
  master$sex <- relevel(master$sex, ref = "0")
  master$agec <- relevel(master$agec, ref = "0")
  master$education <- relevel(master$education, ref = "1")
  master$bmic <- relevel(master$bmic, ref = "1")
  master$smoke <- relevel(master$smoke, ref = "0")
  master$phy_act <- relevel(master$phy_act, ref = "1")
  master$alcohol <- relevel(master$alcohol, ref = "0")
  master$probiotics <- relevel(master$probiotics, ref = "0")
  master$fiberc <- relevel(master$fiberc, ref = "1")
  master$energyc <- relevel(master$energyc, ref = "0")
  master$antibiotics_6mon <- relevel(master$antibiotics_6mon, ref = "0")
  master$crohn_ibs <- relevel(master$crohn_ibs, ref = "0")
  master$diabetes <- relevel(master$diabetes, ref = "0")
  master$cvd <- relevel(master$cvd, ref = "0")
  master$other_psychi <- relevel(master$other_psychi, ref = "0")
  master$dna_plate <- relevel(master$dna_plate, ref = "1")
  master$dna_read <- relevel(master$dna_read, ref = "0")
  
  Maaslin2(func, 
           master, 
           dir,
           # min_abundance=0.0001, 
           min_prevalence=0, # to select all species; PS: those species with all individuals having relative abundance of 0
           random_effects="dna_plate",  #ignored the random effect previously
           fixed_effects = c(exposure, "sex", "agec", "dna_read", "education", "bmic", "smoke", "phy_act", "alcohol", "probiotics", "fiberc", "energyc", "antibiotics_6mon", "crohn_ibs", "diabetes", "cvd", "other_psychi"))
} 

funcs <- c("gmm", "kegg")

fun <- NULL
for (fun in funcs) {
  print(fun)
  master <- NULL
  master <- simp1 %>%
    mutate(other_psychi=depre_atfecal_op)
  
  maaslin_tra(func=get(fun), exposure = "depre_bydef", dir=paste0("/proj/...path.../function_simpler_bydef/", fun, "/"))
}

#--------------------------------------------------------------------------------


