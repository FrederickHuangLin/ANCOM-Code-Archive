library(readr)
library(tidyverse)

otu_data = read_tsv("data/ecam-table-taxa.tsv", skip = 1)
otu_id = otu_data$`feature-id`
meta_data = read_tsv("data/ecam-sample-metadata.tsv")[-1, ]
meta_data = meta_data%>%rename(Sample.ID = `#SampleID`)

# Run ANCOM v1.0
otu_data1 = data.frame(t(otu_data[, -1]), check.names = FALSE)
colnames(otu_data1) = otu_id
otu_data1 = otu_data1%>%rownames_to_column("Sample.ID")%>%arrange(Sample.ID)

meta_data1 = meta_data%>%arrange(Sample.ID)

source("archive/ancom_v1.0.R")

OTUdat = otu_data1; Vardat = meta_data1
adjusted = F; repeated = T; main.var = "delivery"; adj.formula = NULL
repeat.var = "month"; longitudinal = T; random.formula="~ 1 | studyid"
multcorr = 2; sig = 0.05; prev.cut = 0.90

t_start = Sys.time()
out1 = ANCOM.main(OTUdat, Vardat, adjusted, repeated, main.var, adj.formula, repeat.var, 
                  longitudinal, random.formula, multcorr, sig, prev.cut)
t_end = Sys.time()
t_run_old = t_end - t_start # around 30s

res1 = out1$W.taxa

# Run ANCOM v2.1
otu_data2 = data.frame(otu_data[, -1], check.names = FALSE)
rownames(otu_data2) = otu_id

source("scripts/ancom_v2.1.R")

feature_table = otu_data2; sample_var = "Sample.ID"; group_var = NULL
out_cut = 0; zero_cut = 0.90; lib_cut = 0; neg_lb = TRUE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

main_var = "delivery"; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = "~ 1 | studyid"
t_start = Sys.time()
res2 = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
             alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start






