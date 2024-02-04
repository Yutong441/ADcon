library(tidyverse)
library(data.table)
source("utils.R")
source("arrange_figures.R")

root <- "/home/yutong/Data/ADcon/combined/"
save_dir <-  "/home/yutong/Data/ADcon/figures/"
tract_conv <- fread("/home/yutong/Data/tract_data/fiber_abbrev.csv")
covar_list <- list(c("Age", "Sex"),
                   c("Age", "Sex", "WMH_norm"),
                   c("Age", "Sex", "WMH_norm", "Amyloid"),
                   c("Age", "Sex", "Education", "Amyloid"),
                   c("Age", "Sex", "Amyloid")
)

all_df <- fread(paste(root, "/all_tract.csv", sep="")) %>%
    mutate(WMH_norm = log(WMH_vol/TBV)) %>%
    rename(Atlas = Tract) %>%
    rename(Value = Num) %>%
    # mutate(Value = ifelse(Value == NaN, NA, Value)) %>%
    data.frame()

all_df$Atlas <- gsub("_", " ", tract_conv$full_name[
                     match(all_df$Atlas, tract_conv$abbrev)])

eoad <- all_df %>% group_by(Atlas) %>% filter(!duplicated(ID)) %>% ungroup() %>%
    filter(Diagnosis %in% c("YNC", "EOAD"))
Load <- all_df %>% filter(Diagnosis %in% c("ONC", "LOAD"))
eolo <- all_df %>% filter(Diagnosis %in% c("EOAD", "LOAD")) %>%
    mutate(Diagnosis = ifelse(Diagnosis == "LOAD", 1, 0))

out_eoad <- run_logistic_covars(eoad, "AD", covar_list, "logistic", invert=T)
out_load <- run_logistic_covars(Load, "AD", covar_list, "logistic", invert=T)
out_eolo <- run_logistic_covars(eolo, "Diagnosis", covar_list, "logistic", invert=T)
out_eoad$Group <- "EO"
out_load$Group <- "LO"
out_eolo$Group <- "EOLO"

rbind(out_eoad, out_load) %>%
out_eolo %>%
    filter(covar == "Age,Sex") %>%
    select(Module, Group, padj) %>%
    fwrite("/home/yutong/Downloads/tracts_AD.csv")

# identify significant tracts
sig_tracts <- out_load %>%
    filter(padj < 0.05) %>%
    filter(covar == "Age,Sex") %>%
    slice_min(padj, n=10) %>%
    pull(Module)

out_load %>% filter(Module %in% sig_tracts) %>%
    mutate(Val = paste(round(Estimate, 2), gtools::stars.pval(padj))) %>%
    arrange(padj) %>%
    select(Module, covar, Val) %>%
    spread(covar, Val) %>%
    fwrite("/home/yutong/Downloads/eoad_tract.csv")

rbind(out_eoad, out_load) %>%
    filter(padj < 0.05) %>%
    mutate(Val = paste(round(Estimate, 2), gtools::stars.pval(padj))) %>%
    group_by(Group) %>%
    slice_min(padj, n=10) %>%
    ungroup() %>%
    select(Group, Module, Val) %>%
    fwrite("/home/yutong/Downloads/lesioned_tracts.csv")

# ==============================
# PCA
# ==============================
variables <- c("Sex", "Age", "WMH_norm", "AD", "hypertension")
variables <- c("Sex", "Age", "Amyloid", "hypertension")
plot_eoad <- pc_tract(eoad)
eoad_res <- pc_corr_multi(plot_eoad, variables)
plot_load <- pc_tract(Load)
load_res <- pc_corr_multi(plot_load, variables)
plot_eolo <- pc_tract(rbind(eoad, Load))

eoad_res$Group <- "EO"
load_res$Group <- "LO"
rbind(eoad_res, load_res) %>%
    unite(Type, c(Group, PC)) %>%
    spread(Type, padj) %>%
    fwrite("/home/yutong/Downloads/PC.csv")

# plotting variables in PC
var_key <- c("Age"="Age", "Sex"="Sex", "AD"="Alzheimer diagnosis",
             "WMH_norm"="Normalized WMH volume",
             "Hypertension"="Hypertension")

var_key <- c("Age"="Age", "Sex"="Sex", "AD"="Alzheimer diagnosis",
             "WMH_norm"="Normalized WMH volume",
             "Hypertension"="Hypertension")
eoad_list <- plot_eoad %>%
    mutate(Sex = ifelse(Sex==1, "Male", "Female")) %>%
    mutate(AD = ifelse(AD==1, "Alzheimer", "Normal cognition")) %>%
    mutate(Hypertension = ifelse(Hypertension==1, "Hypertension", "No hypertension")) %>%
    get_pc_plots(variables, var_key=var_key)

load_list <- plot_load %>%
    mutate(Sex = ifelse(Sex==1, "Male", "Female")) %>%
    mutate(AD = ifelse(AD==1, "Alzheimer", "Normal cognition")) %>%
    mutate(Hypertension= ifelse(Hypertension==1, "Hypertension", "No hypertension")) %>%
    get_pc_plots(variables, PCs=c("PC1", "PC3"), var_key=var_key)

ggpubr::ggarrange(plotlist=eoad_list)
ggsave("/home/yutong/Downloads/EOAD_PC.jpg", width=10, height=10)
ggpubr::ggarrange(plotlist=load_list)
ggsave("/home/yutong/Downloads/LOAD_PC.jpg", width=10, height=10)

# ======================================================================
# individual cardiovascular risk factors
plot_load$Atlas <- "global"
plot_eoad$Atlas <- "global"
plot_eolo$Atlas <- "global"
all_out <- as.list(unique(plot_eolo$Atlas)) %>% lapply(function(x) {
    sel_df <- plot_eolo %>% 
        filter(Atlas == x)
    if (sel_df %>% pull(PC1) %>% mean() %>% abs() > 0) {
        sel_df <- sel_df %>%
            mutate_if(is.numeric, function(x){(x - mean(x, na.rm=T))/sd(x, na.rm=T)})

        mod <- summary(lm(PC1 ~ Amyloid * Age + Sex + HTN + DM + hyperlipidemia+ smoking + stroke + MI,
                          data=sel_df))
        one_out <- data.frame(mod$coefficients)
        one_out$Atlas <- x
        return(one_out)
    }
})
all_out <- do.call(rbind, all_out)

# global level
all_out %>% rownames_to_column("Feature") %>%
    # filter(grepl("hypercho", Feature)) %>%
    filter(Atlas == "global") %>%
    filter(!Feature %in% c("(Intercept)", "Age", "Sex")) %>%
    mutate(Val = paste(round(Estimate, 2), gtools::stars.pval(Pr...t..))) %>%
    select(Feature, Val) %>%
    arrange(Feature) %>%
    fwrite("/home/yutong/Downloads/corr.csv")

# module level
# most significant tracts
sig_tract <- all_out %>% rownames_to_column("Feature") %>%
    mutate(Feature = gsub("[0-9]+$", "", Feature)) %>%
    filter(Feature == "Amyloid") %>%
    slice_min(Pr...t.., n=10) %>% pull(Atlas)

all_out %>% rownames_to_column("Feature") %>%
    mutate(Feature = gsub("[0-9]+$", "", Feature)) %>%
    filter(!Feature %in% c("(Intercept)", "Age", "Sex")) %>%
    group_by(Feature) %>%
    mutate(padj = p.adjust(Pr...t.., method="fdr")) %>%
    filter(Atlas %in% sig_tract) %>%
    mutate(Val = paste(round(Estimate, 2), gtools::stars.pval(padj))) %>%
    select(Atlas, Feature, Val) %>%
    spread(Feature, Val) %>%
    arrange(Atlas) %>%
    fwrite("/home/yutong/Downloads/corr.csv")

# ======================================================================
all_tracts <- unique(Load$Atlas)
all_IDs <- unique(Load$ID)
tract_mat <- as.list(all_tracts) %>% lapply(function(xx) {
    one_df <- Load %>% filter(Atlas == xx)
    return(one_df[match(all_IDs, one_df$ID), "Value"])
})

tract_mat <- data.frame(do.call(cbind, tract_mat))
colnames(tract_mat) <- all_tracts
rownames(tract_mat) <- all_IDs

pr_mod <- prcomp(tract_mat)
pr_mod$rotation[, "PC3"] %>%
    data.frame() %>%
    magrittr::set_colnames("PC_loading")%>%
    rownames_to_column("Tract") %>%
    slice_max(abs(PC_loading), n=10)
