# To analyze the association between depression and absence of gut microbial species identified from differential abundance analysis using logistic regression models
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

# SCAPIS
scapis <- readRDS("/proj/...path.../scapis_study_analysis.RDS") %>%
  select(id, depre_atfecal, depre_bydef, sex, agec, education, bmic, smoke, phy_act, alcohol, probiotics, fiberc, energyc, antibiotics_6mon, crohn_ibs, diabetes, cvd, depre_atfecal_op, dna_plate, dna_read, antid) 
scapis_mgs <- read.table("/proj/...path.../scapis_metagenomics_mgs_relative_abundances_v1.0.tsv", header = TRUE, sep = "\t")
scapis_mgs <- scapis_mgs %>% rename(id=scapis_id)
mgsid <- print(colnames(scapis_mgs[,2:4742]))

# overall analysis
results <- NULL
mgs <- NULL
for (mgs in mgsid) {
  print(mgs)
  res <- NULL
  res <- scapis_mgs %>% select(id, mgs)
  master <- NULL
  master <- left_join(scapis, res, by="id") %>%
    mutate(rl_mgs=get(mgs),
           outcome=ifelse(rl_mgs==0, 0, 1)) 
  master$depre_atfecal <- factor(master$depre_atfecal)
  
  model <- NULL
  model <- glm(outcome ~ depre_atfecal + as.factor(sex) + as.factor(agec) + as.factor(dna_read) + as.factor(education) + as.factor(bmic) + as.factor(smoke) + as.factor(phy_act) + as.factor(alcohol) + as.factor(probiotics) + as.factor(fiberc) + as.factor(energyc) + as.factor(antibiotics_6mon) + as.factor(crohn_ibs) + as.factor(diabetes) + as.factor(cvd) + as.factor(depre_atfecal_op), data=master, family = "binomial")
  res1 <- NULL
  res1 <- data.frame(
    coefficient=coef(model)["depre_atfecal1"],
    se=summary(model)$coefficients[, "Std. Error"]["depre_atfecal1"],
    p=summary(model)$coefficients[, "Pr(>|z|)"]["depre_atfecal1"]
  )
  
  mm <- master %>%
    group_by(depre_atfecal) %>%
    summarise(
      n=n(),
      events=sum(outcome==1)
    )
  mm0 <- mm %>% filter(depre_atfecal==0) %>%
    mutate(mgs_id=mgs) %>%
    rename(n0=n, events0=events) %>%
    select(mgs_id, n0, events0)
  mm1 <- mm %>% filter(depre_atfecal==1) %>%
    rename(n1=n, events1=events) %>%
    select(n1, events1)
  res2 <- NULL
  res2 <- cbind(mm0, mm1, res1)
  
  results <- rbind(results, res2)
}

write.csv(results, "/proj/...path.../scapis_absence_depre_atfecal.csv")

#SIMPLER
simpler <- readRDS("/proj/...path.../simpler_study_analysis.RDS") %>%
  select(id, depre_atfecal, depre_bydef, sex, agec, education, bmic, smoke, phy_act, alcohol, probiotics, fiberc, energyc, antibiotics_6mon, crohn_ibs, diabetes, cvd, depre_atfecal_op, dna_plate, dna_read, antid) 
simp_mgs <- read.table("/proj/...path.../simpler_metagenomics_mgs_relative_abundances_v2.0.tsv", header = TRUE, sep = "\t")
simp_mgs <- simp_mgs %>% rename(id=SIMPKEY)
mgsid <- print(colnames(simp_mgs[,2:4677]))

# overall analysis
results <- NULL
mgs <- NULL
for (mgs in mgsid) {
  print(mgs)
  res <- NULL
  res <- simp_mgs %>% select(id, mgs)
  master <- NULL 
  master <- left_join(simpler, res, by="id") %>%
    mutate(rl_mgs=get(mgs),
           outcome=ifelse(rl_mgs==0, 0, 1)) 
  master$depre_atfecal <- factor(master$depre_atfecal)
  
  model <- NULL
  model <- glm(outcome ~ depre_atfecal + as.factor(sex) + as.factor(agec) + as.factor(dna_read) + as.factor(education) + as.factor(bmic) + as.factor(smoke) + as.factor(phy_act) + as.factor(alcohol) + as.factor(probiotics) + as.factor(fiberc) + as.factor(energyc) + as.factor(antibiotics_6mon) + as.factor(crohn_ibs) + as.factor(diabetes) + as.factor(cvd) + as.factor(depre_atfecal_op), data=master, family = "binomial")
  res1 <- NULL
  res1 <- data.frame(
    coefficient=coef(model)["depre_atfecal1"],
    se=summary(model)$coefficients[, "Std. Error"]["depre_atfecal1"],
    p=summary(model)$coefficients[, "Pr(>|z|)"]["depre_atfecal1"]
  )
  
  mm <- master %>%
    group_by(depre_atfecal) %>%
    summarise(
      n=n(),
      events=sum(outcome==1)
    )
  mm0 <- mm %>% filter(depre_atfecal==0) %>%
    mutate(mgs_id=mgs) %>%
    rename(n0=n, events0=events) %>%
    select(mgs_id, n0, events0)
  mm1 <- mm %>% filter(depre_atfecal==1) %>%
    rename(n1=n, events1=events) %>%
    select(n1, events1)
  res2 <- NULL
  res2 <- cbind(mm0, mm1, res1)
  
  results <- rbind(results, res2)
}

write.csv(results, "/proj/...path.../simpler_absence_depre_atfecal.csv")
#--------------------------------------------------------------------------------

