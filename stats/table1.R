library(tidyverse)
library(data.table)
source("utils.R")

root <- "/home/yutong/Data/ADcon/"
# reformat the apo status for display
all_df <- fread(paste(root, "/combined/demo_df.csv", sep="/")) %>%
    mutate(WMH_vol = WMH_vol/1000) %>%
    mutate(missing_apo = ifelse(is.na(APOE4), 1, 0)) %>%
    mutate(has_apo = APOE4) %>%
    mutate(has_apo = ifelse(is.na(has_apo), 0, has_apo)) %>%
    mutate(no_apo = ifelse(APOE4 == 0, 1, 0)) %>% 
    mutate(no_apo = ifelse(is.na(no_apo), 0, no_apo)) %>% 
    data.frame()

ync <- all_df %>% filter(Diagnosis == "YNC")
eoad <- all_df %>% filter(Diagnosis == "EOAD")
onc <- all_df %>% filter(Diagnosis == "ONC")
Load <- all_df %>% filter(Diagnosis == "LOAD")

sum_list <- list()
sum_list[[1]] <- sum_stat(ync)
sum_list[[2]] <- sum_stat(eoad)
sum_list[[3]] <- sum_stat_pval(list(ync, eoad))
sum_list[[4]] <- sum_stat(onc)
sum_list[[5]] <- sum_stat(Load)
sum_list[[6]] <- sum_stat_pval(list(onc, Load))
sum_list <- do.call(cbind, sum_list)
colnames(sum_list) <- c("YNC", "EOAD", "p", "ONC", "LOAD", "p")
sum_list %>% data.frame() %>% rownames_to_column("Feature") %>%
    fwrite("/home/yutong/Downloads/table1.csv")

# --------------------missing data--------------------
# identify whose data is missing and why
ori_df <- readxl::read_excel(paste(root, "/demographics/EOLOtoYutong.xlsx", sep="/"))
missing_df <- fread(paste(root, "/demographics/missing.csv", sep="/"))
ori_df %>% pull(`Diagnosis(basedonCDR)`) %>% table()
ori_df %>% 
    rename(ID = `Subject ID`) %>%
    merge(missing_df, by="ID") %>%
    select(Description.y, `Diagnosis(basedonCDR)`) %>%
    arrange(`Diagnosis(basedonCDR)`)
