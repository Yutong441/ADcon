# Graph metrics and submodule analysis
library(tidyverse)
library(data.table)
source("utils.R")

root <- "/home/yutong/Data/ADcon/combined/"
save_dir <-  "/home/yutong/Data/ADcon/figures/"
# specify covariates
covar_list <- list(c("Age", "Sex", "Education"),
                   c("Age", "Sex", "Education", "WMH_norm"))
covar_list2 <- list(c("Age", "Sex"),
                    c("Age", "Sex", "WMH_norm"),
                    c("Age", "Sex", "WMH_norm", "Amyloid"),
                    c("Age", "Sex", "Amyloid", "Education"),
                    c("Age", "Sex", "Amyloid"),
                    c("Age", "Sex", "Dimension3")
)
covar_list3 <- list(c("Age", "Sex", "WMH_norm", "Amyloid"))
covar_list4 <- list(c("Age", "Sex", "Education"),
                    c("Age", "Sex", "Education", "WMH_norm"),
                    c("Age", "Sex", "Education", "WMH_norm", "Amyloid"),
                    c("Age", "Sex", "Education", "Amyloid")
)

graph_df <- fread(paste(root, "/all_conn.csv", sep="/")) %>%
    mutate(WMH_norm = log(WMH_vol/TBV))

eoad <- graph_df %>% filter(Diagnosis %in% c("YNC", "EOAD")) %>% data.frame() %>%
    group_by(Atlas) %>% filter(!duplicated(ID)) %>% ungroup()
Load <- graph_df %>% filter(Diagnosis %in% c("ONC", "LOAD")) %>% data.frame()

# -----logistic regression
out_eoad <- run_logistic_covars(eoad, "AD", covar_list2, method="logistic",
invert=T)
out_load <- run_logistic_covars(Load, "AD", covar_list2, method="logistic",
invert=T)
# now combine the two results and show them in a nice table
# combine_results(out_eoad, out_load, "/home/yutong/Downloads/modules_AD2.csv",
#                 collapse_adj="Age,Sex")
combine_results(out_eoad, out_load, "/home/yutong/Downloads/modules_AD2.csv")

# -----linear regression
variable <- "Amyloid"
out_eoad <- run_logistic_covars(eoad, variable, covar_list2[1:2], method="lingress",
                                invert=T)
out_load <- run_logistic_covars(Load, variable, covar_list2[1:2], method="lingress",
                                invert=T)
combine_results(out_eoad, out_load,
                paste("/home/yutong/Downloads/modules_", variable, ".csv", sep=""))

# -----correlation with cognition
cogn <- c("executive_function", "language", "memory", "visuo_spatial")
eoad_list <- as.list(cogn) %>% lapply(function(xx) {
    one_out <- run_logistic_covars(eoad, xx, covar_list4, method="lingress")
    one_out$cognition <- xx
    return(one_out)
})
eoad_list <- do.call(rbind, eoad_list)

load_list <- as.list(cogn) %>% lapply(function(xx) {
    one_out <- run_logistic_covars(Load, xx, covar_list4, method="lingress")
    one_out$cognition <- xx
    return(one_out)
})
load_list <- do.call(rbind, load_list)

combine_results_cogn(eoad_list, "/home/yutong/Downloads/modules_cogn.csv")

demo_df <- fread(paste(root, "/demo_df.csv", sep="/")) %>%
    mutate(WMH_norm = log(WMH_vol/TBV))

glob_df <- eoad %>% filter(Atlas == "default")
summary(lm(AD ~ Value * Amyloid + Age + Sex, data=glob_df))
summary(lm(AD ~ Value + Amyloid + Age + Sex, data=glob_df))

source("utils.R")
out_eoad <- run_copathology(eoad, "AD", "Amyloid", c("Age", "Sex"))
out_load <- run_copathology(Load, "AD", "Amyloid", c("Age", "Sex"))
combine_results_copath(out_eoad, out_load, "/home/yutong/Downloads/modules_AD2.csv")

# ----------correlation between WMH and amyloid
all_out <- as.list(unique(eoad$Atlas)) %>% lapply(function(x) {
    sel_df <- eoad %>% 
        filter(Atlas == x) %>%
        mutate_if(is.numeric, function(x){(x - mean(x, na.rm=T))/sd(x, na.rm=T)})
    mod <- summary(lm(Value ~ Amyloid + Age + Sex + HTN + DM + hyperlipidemia+ smoking + stroke + MI,
                      data=sel_df))
    one_out <- data.frame(mod$coefficients)
    one_out$Atlas <- x
    return(one_out)
})
all_out <- do.call(rbind, all_out)

# global level
all_out %>% rownames_to_column("Feature") %>%
    # filter(grepl("hypercho", Feature)) %>%
    filter(Atlas == "HCP-MMP") %>%
    filter(!Feature %in% c("(Intercept)", "Age", "Sex")) %>%
    mutate(Val = paste(round(Estimate, 2), gtools::stars.pval(Pr...t..))) %>%
    select(Feature, Val) %>%
    arrange(Feature) %>%
    fwrite("/home/yutong/Downloads/corr.csv")

# module level
all_out %>% rownames_to_column("Feature") %>%
    filter(Atlas != "HCP-MMP") %>%
    mutate(Feature = gsub("[0-9]+$", "", Feature)) %>%
    filter(!Feature %in% c("(Intercept)", "Age", "Sex")) %>%
    group_by(Feature) %>%
    mutate(padj = p.adjust(Pr...t.., method="fdr")) %>%
    ungroup() %>%
    mutate(Val = paste(round(Estimate, 2), gtools::stars.pval(padj))) %>%
    select(Atlas, Feature, Val) %>%
    spread(Feature, Val) %>%
    # mutate(Atlas = module_order[match(Atlas, unique(Atlas))]) %>%
    arrange(Atlas) %>%
    fwrite("/home/yutong/Downloads/corr.csv")
