# To analyze the association between depression and alpha diversity of gut microbiome 
# Alpha diversity: Shannon diversity, inverse Simpson diversity, number of observed features (richness), and Pielou's evenness

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(VennDiagram)

rm(list=ls())

# 1. the association between depression and alpha diversity
# load data
simpler <- readRDS("/proj/...path.../simpler_study_analysis.RDS")
simp_tch <- read.table("/proj/...path.../simpler_metagenomics_technical_variables_v2.0.tsv", header = TRUE, sep = "\t")
dna_plate <- simp_tch %>%
  rename(id=SIMPKEY) %>%
  mutate(dna_plate=case_when(
    startsWith(dna_extraction_plate, "1st") ~ "1",
    startsWith(dna_extraction_plate, "Box_") ~ substr(dna_extraction_plate, 5,6),
    dna_extraction_plate=="" ~ "87",
    TRUE ~ dna_extraction_plate),
    dna_read=cut(sequencing_reads,
                 breaks = quantile(sequencing_reads, probs = c(0, 1/5, 2/5, 3/5, 4/5, 1), na.rm = TRUE),
                 labels = c(0:4),
                 include.lowest = TRUE)) %>%
  select(id, dna_plate, dna_read)

# to count the observed features for each study individual
simp_mgs <- read.table("/proj/...path.../simpler_metagenomics_mgs_relative_abundances_v2.0.tsv", header = TRUE, sep = "\t") 
simp_mgs$obs_feature <- rowSums(simp_mgs[,2:4677] !=0)
simp_obs <- simp_mgs %>%
  rename(id=SIMPKEY) %>%
  select(id,obs_feature)

simp_a <- read.table("/proj/...path.../simpler_metagenomics_alpha_diversity_v2.0.tsv", header = TRUE, sep = "\t") %>% rename(id=SIMPKEY)

simpler <- simpler %>%
  left_join(dna_plate, by="id") %>%
  left_join(simp_obs, by="id") %>%
  left_join(simp_a, by="id") %>%
  mutate(agec=case_when(
    age_atfecal<70 ~ 0,
    age_atfecal>69 & age_atfecal<80 ~ 1,
    age_atfecal>79 ~ 2)) %>%
  mutate(shannon=as.numeric(shannon),
         invsimpson=as.numeric(invsimpson),
         pielou=shannon/log(obs_feature)) %>%
  mutate(fiberc=cut(fiber_intake,
                    breaks = quantile(fiber_intake, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                    labels = c(0:2),
                    include.lowest = TRUE),
         energyc=cut(energy_intake,
                     breaks = quantile(energy_intake, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                     labels = c(0:2),
                     include.lowest = TRUE)) 

# SCAPIS
scapis <- readRDS("/proj/...path.../scapis_study_analysis.RDS")
scapis_tch <- read.table("/proj/...path.../scapis_metagenomics_technical_variables_v1.0.tsv", header = TRUE, sep = "\t")
dna_plate <- scapis_tch %>%
  rename(id=scapis_id) %>%
  mutate(dna_plate=substr(extraction_plate, 8, nchar(extraction_plate)),
         dna_read=cut(sequencing_read_pairs,
                      breaks = quantile(sequencing_read_pairs, probs = c(0, 1/5, 2/5, 3/5, 4/5, 1), na.rm = TRUE),
                      labels = c(0:4),
                      include.lowest = TRUE)) %>%
  select(id, dna_plate, dna_read)
scapis_a <- read.table("/proj/...path.../scapis_metagenomics_shannon_diversity_v1.0.tsv", header = TRUE, sep = "\t") %>% rename(id=scapis_id)
scapis_mgs <- read.table("/proj/...path.../scapis_metagenomics_mgs_relative_abundances_v1.0.tsv", header = TRUE, sep = "\t")
scapis_mgs$obs_feature <- rowSums(scapis_mgs[,2:4742] !=0)
scapis_obs <- scapis_mgs %>%
  rename(id=scapis_id) %>%
  select(id,obs_feature)

scapis <- scapis %>%
  inner_join(dna_plate, by="id") %>%
  left_join(scapis_a, by="id") %>%
  left_join(scapis_obs, by="id") %>%
  select(-alcohol) %>%
  mutate(smoke=ifelse(is.na(smoke), 3, smoke),
         alcohol=alcoholc - 1) %>%
  mutate(shannon=as.numeric(shannon),
         invsimpson=as.numeric(invsimpson),
         pielou=shannon/log(obs_feature)) %>%
  mutate(fiberc=cut(fiber_intake,
                    breaks = quantile(fiber_intake, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                    labels = c(0:2),
                    include.lowest = TRUE),
         energyc=cut(energy_intake,
                     breaks = quantile(energy_intake, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                     labels = c(0:2),
                     include.lowest = TRUE)) 

# composition of depression phenotypes
simpler <- simpler %>%
  mutate(diagnosis=ifelse(!is.na(depre_diag_new), 1, 0),
         medication=ifelse(!is.na(depre_anti), 1, 0),
         symptoms=ifelse(!is.na(depre_sym), 1, 0))

diagnosis <- simpler$id[simpler$diagnosis==1]
medication <- simpler$id[simpler$medication==1]
symptoms <- simpler$id[simpler$symptoms==1]

p_simpler <- venn.diagram(
  x=list(
    "Antidepressant use"=medication,
    "Clinical diagnosis"=diagnosis,
    "Depressive symptoms"=symptoms),
  filename = "/proj/...path.../depression_phenotypes_simpler.tiff",
  height = 8,
  width = 8,
  units = "in",
  resolution = 300,
  fill = c("#1A1A1A", "#7A7A7A", "#D9D9D9"),
  alpha=0.5,
  cex=2,
  cat.cex=2,
  print.mode = c("raw", "percent"),
  print.color = "black",
  cat.pos = c(-20, 20, 180),
  cat.dist = c(0.05, 0.05, 0.05),
  main=""
)

scapis <- scapis %>%
  mutate(diagnosis=ifelse(!is.na(depre_diag), 1, 0),
         medication=ifelse(!is.na(depre_anti), 1, 0),
         symptoms=ifelse(!is.na(depre_sym), 1, 0))
diagnosis <- scapis$id[scapis$diagnosis==1]
medication <- scapis$id[scapis$medication==1]
symptoms <- scapis$id[scapis$symptoms==1]

p_scapis <- venn.diagram(
  x=list(
    "Antidepressant use"=medication,
    "Clinical diagnosis"=diagnosis,
    "Depressive symptoms"=symptoms),
  filename = "/proj/...path.../depression_phenotypes_scapis.tiff",
  height = 8,
  width = 8,
  units = "in",
  resolution = 300,
  fill = c("#6699FF", "#0000FF", "#99FFFF"),
  alpha=0.5,
  cex=2,
  cat.cex=2,
  print.mode = c("raw", "percent"),
  print.color = "black",
  cat.pos = c(-20, 20, 180),
  cat.dist = c(0.05, 0.05, 0.05),
  main=""
)
#-------------------------

# to illustrate the distribution of alpha diversity between groups
simpler <- simpler %>% mutate(exposure=ifelse(depre_atfecal==0, "Without depression", "With depression"))
scapis <- scapis %>% mutate(exposure=ifelse(depre_atfecal==0, "Without depression", "With depression"))

# Shannon diversity
p1 <- ggplot(simpler, aes(x=shannon, fill=exposure)) +
  geom_density(alpha=0.6) +
  scale_fill_manual(values = c("Without depression"="skyblue", "With depression"="hotpink")) +
  theme_minimal(base_size=14) +
  theme(
    strip.text = element_text(face="bold", size=14),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  labs(title="SIMPLER", fill="", x="Shannon diversity", y="Density")
ggsave("/proj/...path.../depre_shannon_simpler.pdf", plot=p1, width = 6, height = 6, units = "in", dpi=300)

p2 <- ggplot(scapis, aes(x=shannon, fill=exposure)) +
  geom_density(alpha=0.6) +
  scale_fill_manual(values = c("Without depression"="skyblue", "With depression"="hotpink")) +
  theme_minimal(base_size=14) +
  theme(
    strip.text = element_text(face="bold", size=14),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  labs(title="SCAPIS", fill="", x="Shannon diversity", y="Density")
ggsave("/proj/...path.../depre_shannon_scapis.pdf", plot=p2, width = 6, height = 6, units = "in", dpi=300)

# inverse Simpson diversity
p3 <- ggplot(simpler, aes(x=invsimpson, fill=exposure)) +
  geom_density(alpha=0.6) +
  scale_fill_manual(values = c("Without depression"="skyblue", "With depression"="hotpink")) +
  theme_minimal(base_size=14) +
  theme(
    strip.text = element_text(face="bold", size=14),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  labs(title="SIMPLER", fill="", x="Inverse Simpson diversity", y="Density")
ggsave("/proj/...path.../depre_invsimpson_simpler.pdf", plot=p3, width = 6, height = 6, units = "in", dpi=300)

p4 <- ggplot(scapis, aes(x=invsimpson, fill=exposure)) +
  geom_density(alpha=0.6) +
  scale_fill_manual(values = c("Without depression"="skyblue", "With depression"="hotpink")) +
  theme_minimal(base_size=14) +
  theme(
    strip.text = element_text(face="bold", size=14),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  labs(title="SCAPIS", fill="", x="Inverse Simpson diversity", y="Density")
ggsave("/proj/...path.../depre_invsimpson_scapis.pdf", plot=p4, width = 6, height = 6, units = "in", dpi=300)

# number of observed feature/species
p5 <- ggplot(simpler, aes(x=obs_feature, fill=exposure)) +
  geom_density(alpha=0.6) +
  scale_fill_manual(values = c("Without depression"="skyblue", "With depression"="hotpink")) +
  theme_minimal(base_size=14) +
  theme(
    strip.text = element_text(face="bold", size=14),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  labs(title="SIMPLER", fill="", x="Observed features", y="Density")
ggsave("/proj/...path.../depre_observed_feature_simpler.pdf", plot=p5, width = 6, height = 6, units = "in", dpi=300)

p6 <- ggplot(scapis, aes(x=obs_feature, fill=exposure)) +
  geom_density(alpha=0.6) +
  scale_fill_manual(values = c("Without depression"="skyblue", "With depression"="hotpink")) +
  theme_minimal(base_size=14) +
  theme(
    strip.text = element_text(face="bold", size=14),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  labs(title="SCAPIS", fill="", x="Observed features", y="Density")
ggsave("/proj/...path.../depre_observed_feature_scapis.pdf", plot=p6, width = 6, height = 6, units = "in", dpi=300)

# Pielou's evenness
p7 <- ggplot(simpler, aes(x=pielou, fill=exposure)) +
  geom_density(alpha=0.6) +
  scale_fill_manual(values = c("Without depression"="skyblue", "With depression"="hotpink")) +
  theme_minimal(base_size=14) +
  theme(
    strip.text = element_text(face="bold", size=14),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  labs(title="SIMPLER", fill="", x="Pielou's evenness", y="Density")
ggsave("/proj/...path.../depre_pielou_simpler.pdf", plot=p7, width = 6, height = 6, units = "in", dpi=300)

p8 <- ggplot(scapis, aes(x=pielou, fill=exposure)) +
  geom_density(alpha=0.6) +
  scale_fill_manual(values = c("Without depression"="skyblue", "With depression"="hotpink")) +
  theme_minimal(base_size=14) +
  theme(
    strip.text = element_text(face="bold", size=14),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  labs(title="SCAPIS", fill="", x="Pielou's evenness", y="Density")
ggsave("/proj/...path.../depre_pielou_scapis.pdf", plot=p8, width = 6, height = 6, units = "in", dpi=300)
#-------------------

# alpha diversity in relation to depression
datasets <- c("simpler","scapis")
expos <- c("depre_atfecal", "depre_bydef", "antid")
a_diver <- c("shannon", "invsimpson", "obs_feature", "pielou")
data <- NULL
results <- NULL
for (data in datasets) {
  print(data)
  toana1 <- NULL
  toana1 <- get(data) 
  
  exp <- NULL
  for (exp in expos) {
    print(exp)
    a <- NULL
    
    for (a in a_diver) {
      print(a)
      master <- NULL
      master <- toana1 %>%
        mutate(exposure=get(exp),
               outcome=get(a))
      
      n <- NULL
      n <- n_distinct(master$exposure)
      
      model<- NULL
      model3 <- glm(outcome ~ as.factor(exposure) + as.factor(sex) + age_atfecal + as.factor(dna_read) + as.factor(education) + as.factor(bmic) + as.factor(smoke) + as.factor(phy_act) + as.factor(alcohol) + as.factor(probiotics) + as.factor(fiberc) + as.factor(energyc) + as.factor(antibiotics_6mon) + as.factor(crohn_ibs) + as.factor(diabetes) + as.factor(cvd) + as.factor(other_psychi), family = gaussian(), data=master)
      
      no_ind <- NULL
      no_ind <- master %>%
        group_by(exposure) %>%
        summarise(no.individuals=n()) %>%
        mutate(cohort=data,
               exposures=exp,
               alpha_diversity=a) %>%
        rename(categories=exposure) %>%
        select(cohort, alpha_diversity, exposures, categories, no.individuals)
      means <- NULL
      means <- master %>%
        group_by(exposure) %>%
        summarise(means=mean(outcome, na.rm = TRUE)) %>%
        select(means)
      no_means <- NULL
      no_means <- cbind(no_ind, means)
      
      m_exp <- NULL
      m_exp <- data.frame(
        md = coef(model),
        se = summary(model)$coefficients[, "Std. Error"],
        lci = confint.default(model)[,1],
        uci = confint.default(model)[,2],
        p = summary(model)$coefficients[, "Pr(>|t|)"]
      )
      m_exp <- m_exp %>%
        mutate(md_ci = paste0(round(md, 3), " (", round(lci, 3), " to ", round(uci, 3), ")"))
      m_exp <- m_exp[1:n,]
      m_exp[1,] <- NA
      
      m_exps <- NULL
      m_exps <- cbind(no_means, m_exp)
      
      results <- rbind(results, m_exps)
    }
    
  }
  
}

write.csv(results, paste0("/proj/...path.../depre_alpha_diversity.csv"))  


# interaction between depression and sex
# SIMPLER
master <- simpler %>%
  mutate(exposure=depre_atfecal)
# shannon diversity
model <- NULL
model <- glm(shannon ~ as.factor(exposure) + as.factor(sex) + as.factor(exposure*sex) + agec + as.factor(dna_read) + as.factor(education) + as.factor(bmic) + as.factor(smoke) + as.factor(phy_act) + as.factor(alcohol) + as.factor(probiotics) + as.factor(fiberc) + as.factor(energyc) + as.factor(antibiotics_6mon) + as.factor(crohn_ibs) + as.factor(diabetes) + as.factor(cvd) + as.factor(other_psychi), family = gaussian(), data=master)

# p-value for interaction
summary(model)
# 0.263759

#inverse simpson diversity
model <- glm(invsimpson ~ as.factor(exposure) + as.factor(sex) + as.factor(exposure*sex) + agec + as.factor(dna_read) + as.factor(education) + as.factor(bmic) + as.factor(smoke) + as.factor(phy_act) + as.factor(alcohol) + as.factor(probiotics) + as.factor(fiberc) + as.factor(energyc) + as.factor(antibiotics_6mon) + as.factor(crohn_ibs) + as.factor(diabetes) + as.factor(cvd) + as.factor(other_psychi), family = gaussian(), data=master)

# p-value for interaction
summary(model)
# 0.25541

# number of observed features
model <- glm(obs_feature ~ as.factor(exposure) + as.factor(sex) + as.factor(exposure*sex) + agec + as.factor(dna_read) + as.factor(education) + as.factor(bmic) + as.factor(smoke) + as.factor(phy_act) + as.factor(alcohol) + as.factor(probiotics) + as.factor(fiberc) + as.factor(energyc) + as.factor(antibiotics_6mon) + as.factor(crohn_ibs) + as.factor(diabetes) + as.factor(cvd) + as.factor(other_psychi), family = gaussian(), data=master)

# p-value for interaction
summary(model)
# 0.41589

# Pielou's evenness
model <- glm(pielou ~ as.factor(exposure) + as.factor(sex) + as.factor(exposure*sex) + agec + as.factor(dna_read) + as.factor(education) + as.factor(bmic) + as.factor(smoke) + as.factor(phy_act) + as.factor(alcohol) + as.factor(probiotics) + as.factor(fiberc) + as.factor(energyc) + as.factor(antibiotics_6mon) + as.factor(crohn_ibs) + as.factor(diabetes) + as.factor(cvd) + as.factor(other_psychi), family = gaussian(), data=master)

# p-value for interaction
summary(model)
# 0.285874

# SCAPIS
master <- scapis %>%
  mutate(exposure=depre_atfecal)
# shannon diversity
model <- NULL
model <- glm(shannon ~ as.factor(exposure) + as.factor(sex) + as.factor(exposure*sex) + agec + as.factor(dna_read) + as.factor(education) + as.factor(bmic) + as.factor(smoke) + as.factor(phy_act) + as.factor(alcohol) + as.factor(probiotics) + as.factor(fiberc) + as.factor(energyc) + as.factor(antibiotics_6mon) + as.factor(crohn_ibs) + as.factor(diabetes) + as.factor(cvd) + as.factor(other_psychi), family = gaussian(), data=master)

# p-value for interaction
summary(model)
# 0.205516

#inverse simpson diversity
model <- glm(invsimpson ~ as.factor(exposure) + as.factor(sex) + as.factor(exposure*sex) + agec + as.factor(dna_read) + as.factor(education) + as.factor(bmic) + as.factor(smoke) + as.factor(phy_act) + as.factor(alcohol) + as.factor(probiotics) + as.factor(fiberc) + as.factor(energyc) + as.factor(antibiotics_6mon) + as.factor(crohn_ibs) + as.factor(diabetes) + as.factor(cvd) + as.factor(other_psychi), family = gaussian(), data=master)

# p-value for interaction
summary(model)
# 0.210655

# number of observed features
model <- glm(obs_feature ~ as.factor(exposure) + as.factor(sex) + as.factor(exposure*sex) + agec + as.factor(dna_read) + as.factor(education) + as.factor(bmic) + as.factor(smoke) + as.factor(phy_act) + as.factor(alcohol) + as.factor(probiotics) + as.factor(fiberc) + as.factor(energyc) + as.factor(antibiotics_6mon) + as.factor(crohn_ibs) + as.factor(diabetes) + as.factor(cvd) + as.factor(other_psychi), family = gaussian(), data=master)

# p-value for interaction
summary(model)
# 0.213114

# Pielou's evenness
model <- glm(pielou ~ as.factor(exposure) + as.factor(sex) + as.factor(exposure*sex) + agec + as.factor(dna_read) + as.factor(education) + as.factor(bmic) + as.factor(smoke) + as.factor(phy_act) + as.factor(alcohol) + as.factor(probiotics) + as.factor(fiberc) + as.factor(energyc) + as.factor(antibiotics_6mon) + as.factor(crohn_ibs) + as.factor(diabetes) + as.factor(cvd) + as.factor(other_psychi), family = gaussian(), data=master)

# p-value for interaction
summary(model)
# 0.499739
#-------------------------------------------------------------------------------

