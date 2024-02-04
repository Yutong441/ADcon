# summarize connectome and demographics information
library(tidyverse)
library(data.table)

# prepare the data as: ID, metric, atlas, value, group
root <- "/media/D/ADcon/results/archive/EOAD_conn/"
demo_root <- "/home/yutong/Data/ADcon/demographics/"
save_dir <- "/home/yutong/Data/ADcon/connectivity/"
save_root <- "/home/yutong/Data/ADcon/combined/"

patients <- list.files(root, full.names=T)
patients_brief <- list.files(root, full.names=F)

all_df <- as.list(patients) %>% lapply(function(x) {
    if (file.exists(paste(x, "/graph_metrics.csv", sep=""))) {
        one_df <- fread(paste(x, "/graph_metrics.csv", sep="")) %>%
            filter(measure == "count") %>%
            rename(Metric = V1) %>%
            rename(Value = unnormalized) %>%
            rename(Atlas = atlas) %>%
            select(Metric, Atlas, Value)
        one_df$ID <- basename(x)
        return(one_df)
    }
})

all_df <- do.call(rbind, all_df)
fwrite(all_df, paste(save_dir, "EOAD_graph.csv", sep="/"))

# tract
all_df <- as.list(patients) %>% lapply(function(x) {
    if (file.exists(paste(x, "/metric/tract_num.csv", sep=""))) {
        one_df <- fread(paste(x, "/metric/tract_num.csv", sep="")) %>%
            rename(Tract = V1) %>%
            mutate(Tract = gsub(".trk", "", Tract)) %>%
            rename(Num = num)
        one_df$ID <- basename(x)
        return(one_df)
    }
})

all_df <- do.call(rbind, all_df)
fwrite(all_df, paste(save_dir, "LOAD_tract.csv", sep="/"))

# --------------------summarize demographic information--------------------
sel_fea <- c("ID", "RID", "Age", "Sex", "Education", "Diagnosis", "Amyloid",
             "Tau", "APOE4", "executive_function", "language", "visuo_spatial",
             "memory", "CDRSB", "MMSE", "Archive Date")
demo_df <- readxl::read_excel(paste(demo_root, "/EOLOtoYutong.xlsx", sep="/")) %>%
    rename(ID = `Subject ID`) %>%
    rename(Education = EDU) %>%
    rename(Diagnosis = `Diagnosis(basedonCDR)`) %>%
    rename(Tau = PLASMAPTAU181) %>%
    rename(Amyloid = CL_amyloid) %>%
    rename(executive_function = ADNI_EF) %>%
    rename(language = ADNI_LAN) %>%
    rename(visuo_spatial = ADNI_VS) %>%
    rename(memory = ADNI_MEM) %>%
    select(all_of(sel_fea))

vol_eo <- fread(paste(demo_root, "/WMH_EOAD.csv", sep="/")) %>%
    rename(ID = V1)
vol_lo <- fread(paste(demo_root, "/WMH_LOAD.csv", sep="/")) %>%
    rename(ID = V1) %>%
    filter(!ID %in% vol_eo$ID)
vol <- rbind(vol_eo, vol_lo) %>%
    filter(!is.na(WMH_vol))

# add FLAIR 2D/3D info
flair_dim <- fread(paste(demo_root, "/FLAIR_dim.csv", sep="")) %>%
    filter(grepl("FLAIR", Type)) %>%
    rename(ID = V1) %>%
    mutate(Dimension3 = ifelse(Dim > 2, 1, 0)) %>%
    select(ID, Dimension3)

# vascular risk factor
vasc_df <- fread(paste(save_root, "../medical_history/CV_ADNI123.csv", sep="/"))

demo_df %>% merge(vol, by="ID") %>%
    merge(flair_dim, by="ID") %>%
    mutate(AD = ifelse(grepl("AD", Diagnosis), 1, 0)) %>%
    merge(vasc_df, by="RID",all.x=T) %>%
    mutate(`Archive Date` = as.POSIXct(`Archive Date`, format="%m/%d/%Y")) %>%
    mutate(diff_time = difftime(`Archive Date`, Time, units="days")) %>%
    mutate(diff_time = abs(as.numeric(diff_time))) %>%
    group_by(RID) %>%
    slice_min(diff_time) %>%
    ungroup() %>% 
    select(!one_of(c("diff_time", "Time", "Archive Date"))) %>%
    filter(!is.na(WMH_vol)) %>% 
    filter(!duplicated(RID)) %>%
    fwrite(paste(save_root, "/demo_df.csv", sep="/"))

# ==================================================
# include WMH
# add connectivity (only look at global efficiency)
demo_df <- fread(paste(save_root, "/demo_df.csv", sep="/"))
conn_eoad <- fread(paste(demo_root, "/../connectivity/EOAD_graph.csv", sep="/"))
conn_load <- fread(paste(demo_root, "/../connectivity/LOAD_graph.csv", sep="/"))
conn <- rbind(conn_eoad, conn_load) %>%
    filter(Metric == "global efficiency") %>%
    select(!Metric)

# save connectivity data
conn %>% merge(demo_df, by="ID") %>%
    fwrite(paste(save_root, "all_conn.csv", sep="/"))

# save tract data
tract_eoad <- fread(paste(demo_root, "/../connectivity/EOAD_tract.csv", sep="/"))
tract_load <- fread(paste(demo_root, "/../connectivity/LOAD_tract.csv", sep="/"))
rbind(tract_eoad, tract_load) %>%
    merge(demo_df, by="ID") %>%
    fwrite(paste(save_root, "all_tract.csv", sep="/"))

# --------------------prepare data for VLSM--------------------
root <- "/home/yutong/Data/ADcon/"
demo_df <- fread(paste(root, "/combined/demo_df.csv", sep="/"))
all_files1 <- list.files(paste(root, "/WMH_mask/", "EOAD", sep="/"), full.names=F)
all_files2 <- list.files(paste(root, "/WMH_mask/", "LOAD", sep="/"), full.names=F)
all_files <- gsub(".nii.gz", "", c(all_files1, all_files2))

all_paths1 <- list.files(paste(root, "/WMH_mask/", "EOAD", sep="/"), full.names=T)
all_paths2 <- list.files(paste(root, "/WMH_mask/", "LOAD", sep="/"), full.names=T)
all_paths <- c(all_paths1, all_paths2)
demo_df$paths <- all_paths[match(demo_df$ID, all_files)]

# EOAD patients
demo_df %>% filter(Diagnosis %in% c("EOAD", "YNC")) %>%
    filter(!is.na(paths)) %>%
    mutate(WMH_norm = WMH_vol/TBV) %>%
    fwrite(paste(root, "/VLSM/labels/EOAD.csv", sep="/"))

demo_df %>% filter(Diagnosis %in% c("LOAD", "ONC")) %>%
    filter(!is.na(paths)) %>%
    mutate(WMH_norm = WMH_vol/TBV) %>%
    fwrite(paste(root, "/VLSM/labels/LOAD.csv", sep="/"))
