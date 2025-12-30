# Differential abundance analysis: to identify species that are more or less abundance in individuals with and without depression
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(Maaslin2)

# differential abundance analysis
# gut microbiome composition
simp_mgs <- read.table("/proj/...path.../simpler_metagenomics_mgs_relative_abundances_v2.0.tsv", header = TRUE, sep = "\t") 
simpler <- readRDS("/proj/...path.../simpler_study_analysis.RDS") %>%
  select(id, depre_atfecal, depre_bydef, sex, agec, education, bmic, smoke, phy_act, alcohol, probiotics, fiberc, energyc, antibiotics_6mon, crohn_ibs, diabetes, cvd, depre_atfecal_op, dna_plate, dna_read, antid) 

simp2 <- simpler %>%
  filter(antid<3) %>%
  mutate(across(2:21, ~as.factor(.))) 

simp <- simpler %>%
  mutate(across(2:21, ~as.factor(.))) 

sim_id <- simpler %>% select(id)
simp_mgs1 <- simp_mgs %>% 
  rename(id=SIMPKEY) %>% 
  right_join(sim_id, by="id")

# sort the tables by id
simp1 <- simp[order(simp$id),]
simp_mgs1 <- simp_mgs1[order(simp_mgs1$id),]

rownames(simp_mgs1) <- simp_mgs1$id
simp1 <- as.data.frame(simp1)
rownames(simp1) <- simp1$id
# remove the column of id
simp_mgs1 <- simp_mgs1[, -c(1)]
simp1 <- simp1[, -c(1)]

# to use the maaslin regression to conduct the differential abundance analysis
maaslin_tra<-function(exposure, dir){
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
  
  Maaslin2(simp_mgs1, 
           master, 
           dir,
           # min_abundance=0.0001, 
           min_prevalence=0, # to select all species; PS: those species with all individuals having relative abundance of 0 
           random_effects="dna_plate",  #ignored the random effect previously
           fixed_effects = c(exposure, "sex", "agec", "education", "bmic", "smoke", "phy_act", "alcohol", "probiotics", "fiberc", "energyc", "antibiotics_6mon", "crohn_ibs", "diabetes", "cvd", "other_psychi", "dna_read"))
} 

# stratifies analysis by sex for antidepressant use 
maaslin_tras<-function(exposure, dir){
  # Set the reference level for covariates
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
  
  Maaslin2(simp_mgs2, 
           master, 
           dir,
           # min_abundance=0.0001, 
           min_prevalence=0, # to select all species; PS: those species with all individuals having relative abundance of 0 
           random_effects="dna_plate",  #ignored the random effect previously
           fixed_effects = c(exposure, "agec", "education", "bmic", "smoke", "phy_act", "alcohol", "probiotics", "fiberc", "energyc", "antibiotics_6mon", "crohn_ibs", "diabetes", "cvd", "other_psychi", "dna_read"))
} 


# differential abundance analysis
master <- simp1 %>%
  mutate(other_psychi=depre_atfecal_op)
maaslin_tra(exposure = "depre_atfecal", dir="/proj/...path.../composition_simpler_depre_atfecal/")
maaslin_tra(exposure = "depre_bydef", dir="/proj/...path.../composition_simpler_bydef/")

master <- simp2 %>%
  mutate(other_psychi=depre_atfecal_op)
masterid <- master %>% select(id)
simp_mgs1 <- simp_mgs %>% 
  rename(id=SIMPKEY) %>% 
  right_join(masterid, by="id")
master <- master[order(master$id),]
simp_mgs1 <- simp_mgs1[order(simp_mgs1$id),]
rownames(simp_mgs1) <- simp_mgs1$id
master <- as.data.frame(master)
rownames(master) <- master$id
simp_mgs1 <- simp_mgs1[, -c(1)]
master <- master[, -c(1)]
maaslin_tra(exposure = "antid", dir="/proj/...path.../composition_simpler_antid/")

# stratifies analysis by sex
# master$sex <- factor(master$sex, levels = 0:1, labels = c("Women", "Men"))
i <- NULL
for (i in 0:1) {
  master <- NULL
  master <- simp2 %>%
    mutate(other_psychi=depre_atfecal_op) %>%
    filter(sex==i)
  
  sim_id <- NULL
  sim_id <- master %>% select(id)
  simp_mgs2 <- NULL
  simp_mgs2 <- simp_mgs %>% 
    rename(id=SIMPKEY) %>% 
    right_join(sim_id, by="id")
  
  # sort the tables by id
  master <- master[order(master$id),]
  simp_mgs2 <- simp_mgs2[order(simp_mgs2$id),]
  
  rownames(simp_mgs2) <- simp_mgs2$id
  master <- as.data.frame(master)
  rownames(master) <- master$id
  # remove the column of id
  simp_mgs2 <- simp_mgs2[, -c(1)]
  master <- master[, -c(1)]
  
  maaslin_tras(exposure = "depre_atfecal", dir=paste0("/proj/...path.../composition_simpler_depre_atfecal_sex", i, "/"))
}
#--------------------------------------------------------------------------------

# interaction between sex and depression
# remove those with 0 in all individuals
simp_mgs <- simp_mgs %>% rename(id=SIMPKEY)
ids <- simp_mgs %>% select(id)
mgs1 <- simp_mgs %>% select(-id)
col_sums <- colSums(mgs1[,1:4676], na.rm = T)
mgs2 <- mgs1[, col_sums!=0]
simp_mgs <- cbind(ids, mgs2)
mgsid <- print(colnames(simp_mgs[,2:4677]))
mgs <- NULL
results <- NULL
for (mgs in mgsid) {
  print(mgs)
  res <- NULL
  res <- simp_mgs %>% select(id, mgs)
  master <- NULL
  master <- left_join(simpler, res, by="id") %>%
    mutate(outcome=get(mgs),
           other_psychi=depre_atfecal_op,
           exposure=depre_atfecal)
  model <- NULL
  model <- tryCatch({ 
    lmer(outcome ~ exposure*sex + agec + dna_read + education + bmic + smoke + phy_act + alcohol + probiotics + fiberc + energyc + antibiotics_6mon + crohn_ibs + diabetes + cvd + other_psychi + (1|dna_plate), data=master)
  }, error = function(e){
    message("returning NA matrix", e$message)
    return(NULL)
  })
  if (!is.null(model)) { 
    m_exp <- NULL
    m_exp <- data.frame(
      mgs_id=mgs,
      p_interaction = summary(model)$coefficients[, "Pr(>|t|)"]
    )
    m_exp <- m_exp["exposure1:sex1",]
    results <- rbind(results, m_exp)
  } else {
    m_exp <- NULL
    m_exp <- data.frame(
      mgs_id=mgs,
      p_interaction = NA) # error
    results <- rbind(results, m_exp)
  }
}

# adjust for multiple testing
results$q_interaction <- p.adjust(results$p_interaction, method="BH")

write.csv(results, "/proj/...path.../composition_sex_interaction_simpler.csv")
#--------------------------------------------------------------------------------
