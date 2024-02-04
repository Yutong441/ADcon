# extract vascular risk factors from ADNI
library(data.table)
library(tidyverse)

#' Identify the presence of certain cardiovascular risk factors in a particular
#' column of a dataframe using regular expression
#'
#' @param one_df a dataframe
#' @param history_col the column containing the cardiovascular risk factor
#' information. It could be a string of text.
#' @param ID_col the column containing the ID of each subject
#' @param time_col the column containing the time at which the cardiovascular
#' risk information was taken.
#'
#' @example An example dataframe would look like:
#' | ID  | Time       | History                        |
#' | ID1 | 2 Sep 2022 | MI; TIA; No evidence of stroke |
#' | ID1 | 3 Sep 2022 | hyperlipidemia                 |
summarize_history <- function(one_df, history_col, ID_col, time_col) {
    one_df <- one_df %>% mutate(PMH = !!as.symbol(history_col))
    one_df$uniqID <- 1:dim(one_df)[1]

    # hypertension
    HTN_ID <- one_df %>%
        filter(grepl("HTN|blood pressure|ertension", PMH, ignore.case=T
                     ) | grepl("HBP", PMH)) %>%
        filter(!grepl("low blood pressure|orthostatic|blood pressure drop",
                      PMH, ignore.case=T)) %>%
        pull(uniqID)
    one_df$HTN <- ifelse(one_df$uniqID %in% HTN_ID, 1, 0)

    # hyperlipidemia
    lipid_ID <- one_df %>%
        filter(grepl("lipid|triglyceride|cholesterol", PMH, ignore.case=T)) %>%
        pull(uniqID)
    one_df$hyperlipidemia <- ifelse(one_df$uniqID %in% lipid_ID, 1, 0)

    # smoking
    smoke_ID <- one_df %>%
        filter(grepl("smoke|smoking|cigar", PMH, ignore.case=T)) %>%
        pull(uniqID)
    one_df$smoking <- ifelse(one_df$uniqID %in% smoke_ID, 1, 0)

    # stroke
    stroke_ID <- one_df %>%
        filter(grepl("stroke", PMH, ignore.case=T)) %>%
        filter(!grepl("negative for stroke|no evidence of stroke", PMH, ignore.case=T)) %>%
        filter(!grepl("no evidence found on mri of a stroke", PMH, ignore.case=T)) %>%
        filter(!grepl("no stroke|ruled out|rule out", PMH, ignore.case=T)) %>%
        filter(!grepl("mini stroke|mini-stroke", PMH, ignore.case=T)) %>%
        pull(uniqID)
    one_df$stroke <- ifelse(one_df$uniqID %in% stroke_ID, 1, 0)

    # MI
    MI_ID <- one_df %>%
        filter(grepl("heart attack|myocardial infarction|coronary artery disease|angina",
                     PMH, ignore.case=T) | grepl("CAD|MI", PMH, ignore.case=F)) %>%
        # the following two lines search for "MI" as a single word and not as a
        # part of other abbreviations, except for the case of STEMI (or NSTEMI)
        filter(!grepl("MI[a-z|A-z]+", PMH)) %>%
        filter(!grepl("[a-z|A-z]+MI", PMH) | grepl("STEMI", PMH)) %>%
        filter(!grepl("no MI", PMH, ignore.case=F)) %>%
        pull(uniqID)
    one_df$MI <- ifelse(one_df$uniqID %in% MI_ID, 1, 0)

    # AF
    AF_ID <- one_df %>% filter(grepl("atrial fibrillation", PMH, ignore.case=T)
                               ) %>% pull(uniqID)
    one_df$AF <- ifelse(one_df$uniqID %in% AF_ID, 1, 0)

    # diabetes
    DM_ID <- one_df %>%
        filter(grepl("diabetes", PMH, ignore.case=T
                     ) | grepl("DM", PMH, ignore.case=F)) %>%
        filter(!grepl("DM[a-z|A-z]+", PMH)) %>%
        filter(!grepl("[a-z|A-z]+DM", PMH) | grepl("NIDDM", PMH)) %>%
        filter(!grepl("gestational|borderline|prediabetes|pre-diabetes", PMH,
                      ignore.case=T)) %>%
        pull(uniqID)
    one_df$DM <- ifelse(one_df$uniqID %in% DM_ID, 1, 0)

    # summarize per time point
    one_df %>%
        mutate(oneTP = !!as.symbol(time_col)) %>%
        mutate(oneID = !!as.symbol(ID_col)) %>%
        group_by(oneTP, oneID) %>%
        # for a particular patient at a particular time point
        # count the number of occurence of particular risk factors
        summarize(AF = sum(AF), hyperlipidemia=sum(hyperlipidemia),
                  MI = sum(MI), stroke=sum(stroke),
                  HTN = sum(HTN), DM=sum(DM), smoking=sum(smoking)
                  ) %>%
        ungroup() %>%
        mutate_at(c("AF", "hyperlipidemia", "MI", "stroke", "HTN", "DM",
                    "smoking"), function(x){ifelse(x > 0, 1, 0)}) %>%
        return()
}


#' Incidence carried forward
#' 
#' @description The input dataframe contains the incidence of an event at
#' different time points. This function outputs how many incidents experimented
#' by a particular patient before a particular time point. 
#'
#' @param one_df input dataframe
#' @param sel_col the incidence of an event
#' @param ID_col ID column
#' @param TP_col time point column
#' @param thres assume that if an event is documented at time A, the event
#' probably happened some time before time A. This parameter specifies the
#' amount of time in days.
assign_time <- function(one_df, sel_col, ID_col, TP_col, thres=365) {
    one_df <- one_df %>%
        mutate(oneID = !!as.symbol(ID_col)) %>%
        mutate(oneTP = !!as.symbol(TP_col)) %>%
        mutate(Feature = !!as.symbol(sel_col))

    # find out when a particular feature becomes 1 (the event occurs)
    index_TP <- one_df %>% filter(Feature == 1) %>%
        group_by(oneID) %>%
        # make sure only one time point is selected (the earlist time point)
        slice_min(oneTP) %>%
        ungroup() %>%
        rename(indexT = oneTP) %>%
        select(oneID, indexT)

    one_df <- one_df %>% merge(index_TP, by="oneID", all.x=T) %>%
        mutate(diff_time = as.numeric(difftime(oneTP, indexT, units="days"))) %>%
        mutate(Feature = ifelse(is.na(diff_time), Feature, 
                                ifelse(diff_time >= -thres, 1, Feature)))
    one_df[, sel_col] <- one_df$Feature
    one_df %>% select(!one_of(c("diff_time", "indexT", "Feature"))) %>% return()
}


extend_window <- function(one_df, var_list, ID_col, time_col, thres) {
    for (i in var_list) {
        one_df <- assign_time(one_df, i, ID_col, time_col, thres)
    }
    return(one_df)
}


root <- "/home/yutong/Data/ADcon/"
adni3 <- fread(paste(root, "/medical_history/ADNI3_PMH.csv", sep="/"))
adni12 <- readxl::read_excel(paste(root, "/medical_history/ADNI12_PMH.xlsx", sep="/"))

adni3_risk <- summarize_history(adni3, "IHDESC", "RID", "USERDATE") %>%
    mutate(oneTP = as.POSIXct(oneTP, format="%Y/%m/%d"))
adni12_risk <- summarize_history(adni12, "MHDESC", "RID", "USERDATE")

# extend timewindow
vasc_vec <- c("HTN", "smoking", "DM", "hyperlipidemia", "AF")
vasc_vec2 <- c("stroke", "MI")
adni12_risk2 <- extend_window(data.frame(adni12_risk), vasc_vec, "oneID", "oneTP", thres=365)
adni12_risk2 <- extend_window(data.frame(adni12_risk2), vasc_vec2, "oneID", "oneTP", thres=0)
adni3_risk2 <- extend_window(data.frame(adni3_risk) %>% filter(!is.na(oneTP)),
                             vasc_vec, "oneID", "oneTP", thres=365)
adni3_risk2 <- extend_window(data.frame(adni3_risk2), vasc_vec2, "oneID", "oneTP", thres=0)

adni3_risk2 %>% 
    rename(RID = oneID) %>% rename(Time = oneTP) %>%
    arrange(RID, Time) %>%
    fwrite(paste(root, "/medical_history/CV_ADNI3.csv", sep="/"))
adni12_risk2 %>%
    rename(RID = oneID) %>% rename(Time = oneTP) %>%
    arrange(RID, Time) %>%
    fwrite(paste(root, "/medical_history/CV_ADNI12.csv", sep="/"))

comb_risk <- rbind(adni3_risk, adni12_risk)
comb_risk2 <- extend_window(data.frame(comb_risk), vasc_vec, "oneID", "oneTP", thres=365)
comb_risk2 <- extend_window(data.frame(comb_risk), vasc_vec2, "oneID", "oneTP", thres=0)
comb_risk2 %>%
    rename(RID = oneID) %>% rename(Time = oneTP) %>%
    arrange(RID, Time) %>%
    fwrite(paste(root, "/medical_history/CV_ADNI123.csv", sep="/"))

# ======================================================================
# compare with the most recent ADNI data release
history1 <- fread(paste(root, "/medical_history/MODHACH_03Feb2024.csv", sep="/")) %>%
    rename(HTN = HMHYPERT) %>%
    rename(stroke = HMSTROKE) %>%
    select(RID, USERDATE, HTN, stroke)

comb_risk2 %>%
    mutate(Feature = stroke) %>%
    # extract baseline feature
    select(oneID, oneTP, Feature) %>%
    group_by(oneID) %>%
    slice_min(oneTP) %>%
    ungroup() %>%
    merge(history1, by.x="oneID", by.y="RID") %>%
    filter(Feature == 1 & stroke == 0) %>% 
    pull(oneID) %>% unique() %>%
    length()
