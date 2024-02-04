library(magrittr)
library(data.table)

dataset <- "EOAD"
root <- "/rds/user/yc441/hpc-work/ADcon/results/"
df_path <- paste(root, dataset, "-T1_6_30_2023.csv", sep="")
mask_path <- paste(root, dataset, "_WMH", sep="")
save_dir <- paste(root, "/labels/", sep="")

# prepare data
xdf <- fread(df_path)
all_imgs <- data.frame(ID = list.files(mask_path),
                       paths = list.files(mask_path, full.names=TRUE)
)
all_imgs <- all_imgs %>% dplyr::mutate(ID = gsub(".nii.gz$", "", ID))

for (i in c(dataset)) {
    sel_df <- xdf %>% 
        dplyr::mutate(new_group = ifelse(
                Group == "CN" | Group == "SMC", 0,
                ifelse(Group == "EMCI" | Group == "LMCI" | Group == "MCI", 1, 2))) %>%
        dplyr::select(all_of(c("Subject", "new_group"))) %>%
        dplyr::filter(!is.na(!!as.symbol("new_group"))) %>%
        merge(all_imgs, by.x="Subject", by.y="ID", all=FALSE)
    sel_df %>% dplyr::select(paths) %>%
        fwrite(paste(save_dir, "/", i, "_img.txt", sep=""), col.names=FALSE)
    sel_df %>% dplyr::select(dplyr::all_of(c(i))) %>%
        fwrite(paste(save_dir, "/", i, "_behavior.txt", sep=""), col.names=FALSE)
}

# df_path <- "/home/yutong/Data/SVDnet/RUNDMC/original_labels/run_long.csv"
# all_df <- fread(df_path)
# grep("depre", colnames(all_df), ignore.case=T, value=T)
# 
# all_df %>% filter(year == 2006) %>% 
#     dplyr::mutate(ID = paste("RUNDMC_", sprintf("%03d", ID), sep="")) %>%
#     pull(Depressive_symptoms_2006) %>% table()
