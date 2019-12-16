library(readr)
library(tidyverse)
library(compositions)

otu_data = read_tsv("data/moving-pics-table.tsv", skip = 1)
otu_id = otu_data$`feature-id`
otu_data = data.frame(otu_data[, -1], check.names = FALSE)
rownames(otu_data) = otu_id

meta_data = read_tsv("data/moving-pics-sample-metadata.tsv")[-1, ]
meta_data = meta_data%>%rename(Sample.ID = SampleID)

source("scripts/ancom_v2.0.R")

# Step 1: Data preprocessing

feature_table = otu_data; sample_var = "Sample.ID"; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM

main_var = "Subject"; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = NULL
t_start = Sys.time()
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start # around 30s

write_csv(res, "outputs/res_moving_pics.csv")

# Step 3 (Optional): Volcano Plot

# Calculate clr
clr_table = apply(feature_table, 2, clr)
# Calculate clr mean difference
eff_size = apply(clr_table, 1, function(y) 
  lm(y ~ x, data = data.frame(y = y, x = meta_data %>% pull(main_var)))$coef[-1])
# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Data frame for the figure
dat_fig = data.frame(taxa_id = res$taxa_id, x = eff_size, y = res$W) %>% 
  mutate(zero_ind = factor(ifelse(is.infinite(y), "Yes", "No"), levels = c("Yes", "No")))
# Replace Inf by (n_taxa - 1) for structural zeros
dat_fig$y = replace(dat_fig$y, is.infinite(dat_fig$y), n_taxa - 1)
# Annotation data
dat_ann = data.frame(x = min(dat_fig$x), y = cut_off["detected_0.7"], label = "W[0.7]")

fig = ggplot(data = dat_fig) + aes(x = x, y = y) + 
  geom_point(aes(color = zero_ind)) + 
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE) + 
  labs(x = "CLR mean difference", y = "W statistic", title = "Volcano Plot") +
  scale_color_discrete(name = "Structural\nzero", drop = FALSE) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
fig  
ggsave("images/moving_pics.jpeg", height=5, width=6.25, units='in', dpi = 300)  
  
  
  
  
  
  