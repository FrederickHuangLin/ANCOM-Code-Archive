library(readr)
library(tidyverse)

otu_data = read_tsv("data/ecam-table-taxa.tsv", skip = 1)
otu_id = otu_data$`feature-id`
otu_data = data.frame(t(otu_data[, -1]), check.names = FALSE)
colnames(otu_data) = otu_id
otu_data = otu_data%>%rownames_to_column("Sample.ID")

meta_data = read_tsv("data/ecam-sample-metadata.tsv")[-1, ]
meta_data = meta_data%>%rename(Sample.ID = `#SampleID`)

otu_data = otu_data%>%arrange(Sample.ID)
meta_data = meta_data%>%arrange(Sample.ID)
all(otu_data$Sample.ID == meta_data$Sample.ID)

main_var = "delivery"; zero_cut = 0.90; p_adjust_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = "~ 1 | studyid"

# Run ANCOM
source("ancom_v2.0.R")

t_start = Sys.time()
res = ANCOM(otu_data, meta_data, main_var,  zero_cut, p_adjust_method, alpha,
            adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start # around 30s

write_csv(res, "outputs/res_ecam.csv")

# Run old ANCOM
source("ancom_v1.0.R")

var_data = meta_data; long = T; rand.formula="~1|studyid"
OTUdat = otu_data; Vardat = meta_data
adjusted=F; repeated=T; main.var="delivery"; adj.formula=NULL
repeat.var="month"; longitudinal=T; random.formula="~1|studyid"
multcorr=2; sig=0.05; prev.cut=0.90

t_start = Sys.time()
longitudinal_comparison_test=ANCOM.main(OTUdat, Vardat, adjusted, repeated,
                                        main.var, adj.formula, repeat.var, longitudinal,
                                        random.formula, multcorr, sig, prev.cut)
t_end = Sys.time()
t_run_old = t_end - t_start # around 30s

res_old = longitudinal_comparison_test$W.taxa




