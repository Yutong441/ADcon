get_conv <- function() {
    return(c(
      Sex = "Sex (Male)",
      Age = "Age",
      Education = "Education (years)",
      has_apo = "APOE4 +",
      no_apo = "APOE4 -",
      missing_apo = "No genetics data",
      Amyloid = "Amyloid",
      Tau = "Plasma tau",
      executive_function = "Executive function",
      language = "Language",
      visuo_spatial = "Visual spatial",
      memory = "Memory",
      CDRSB = "CDRSB",
      WMH_vol = "WMH volume (mL)"
))
}


theme_publication <- function(base_size=15, base_family="arial") {
      (ggthemes::theme_foundation(base_size=base_size, base_family=base_family)
       + ggplot2::theme(
               plot.title = ggplot2::element_text(
                    face = "bold", size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = ggplot2::element_rect(colour = NA),
               plot.background = ggplot2::element_rect(colour = NA),
               panel.border = ggplot2::element_rect(colour = NA),
               axis.title = ggplot2::element_text(face = "bold", size = rel(1)),
               axis.title.y = ggplot2::element_text(angle = 90, vjust = 2),
               axis.title.x = ggplot2::element_text(vjust = -0.2),
               axis.text = ggplot2::element_text(),
               axis.line = ggplot2::element_line(colour = "black"),
               axis.ticks = ggplot2::element_line(),
               panel.grid.major = ggplot2::element_line(colour = "#f0f0f0"),
               panel.grid.minor = ggplot2::element_blank(),
               legend.key = ggplot2::element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.key.size = unit(0.2, "cm"),
               legend.margin = unit(0, "cm"),
               legend.title = ggplot2::element_text(face = "italic"),
               plot.margin = unit(c(10, 5, 5, 5), "mm"),
               strip.background = ggplot2::element_rect(colour="#f0f0f0",
                                                        fill="#f0f0f0"),
               strip.text = ggplot2::element_text(face = "bold")
          ))
}

sprint1 <- function(xx) {sprintf("%.1f", round(xx,1))}
sprint2 <- function(xx) {sprintf("%.2f", round(xx, 2))}
sprint3 <- function(xx) {sprintf("%.3f", round(xx, 3))}

factor2numeric <- function (xx){
        if (is.factor (xx)){
                sum_col <- xx %>% as.character () %>% as.numeric ()
        }else{sum_col <- xx}
        return(sum_col[!is.na (sum_col)])
}


sum_stat <- function(xdf) {
    sum_list <- list()
    keep_col <- names(get_conv())
    for (i in keep_col){
        proceed <- sum(c("Date", "POSIXt", "character") %in% class(xdf[, i]))
        if (proceed == 0){
            # convert factors to numeric
            sum_col <- factor2numeric(xdf[, i])
            # percentages for 2-level features
            within_01 <- mean(sum_col %in% c(0, 1))
            if (length(unique(sum_col)) == 2 & within_01 == 1){
                sum_levels <- sum_col %>% as.factor() %>% levels()
                num_lev <- sum(sum_col == sum_levels[2])
                # perc_lev <- sprint1(mean(sum_col == sum_levels[2])*100)
                perc_lev <- sprint1(num_lev/sum(!is.na(sum_col))*100)
                sum_text <- paste(num_lev, " (", perc_lev, ")", sep="")
            } else if (length(unique(sum_col)) == 1) {
                if (mean(sum_col) == 0) {
                    sum_text <- paste(0, " (", 0, ")", sep="")
                } else {
                    sum_text <- paste(length(sum_col), " (", 100.0, ")", sep="")
                }
            } else if (length(sum_col) == 0) {
                sum_text <- "NA"
            }else {
                normality<-shapiro.test (sum_col)$p.value
                normality <- 0
                if (normality > 0.05){
                    sum_text <- paste (
                    sprint1(mean (sum_col, na.rm=T)), " \u00B1", 
                    sprint1(sd (sum_col, na.rm=T)), sep="")
                }else{
                    sum_text <- paste(
                    sprint1(median(sum_col, na.rm=T)), " (",
                    sprint1(quantile(sum_col, 0.25, na.rm=T)), "-",
                    sprint1(quantile(sum_col, 0.75, na.rm=T)), 
                    ")", sep="")
                }
            }
            sum_list [[i]] <- sum_text
        }
    }
    sum_df <- do.call(c, sum_list) %>% data.frame() %>% 
        rownames_to_column("Variable") %>%
        magrittr::set_colnames(c("Variable", "Value")) %>%
        mutate(Variable = as.character(Variable))
    conversion <- get_conv()
    sum_df <- sum_df [match(names(conversion), sum_df$Variable), ]
    sum_df$Variable  <- as.character(conversion)
    for (i in 1:dim(sum_df)[1]) {
        if (grepl("\\(", sum_df$Value[i]) & !grepl("-", sum_df$Value[i])){
            sum_df$Variable[i] <- paste(sum_df$Variable[i], ", n (%)", sep="")
        }else if(grepl("\\(", sum_df$Value[i]) & grepl("-", sum_df$Value[i])){
            sum_df$Variable[i] <- paste(sum_df$Variable[i], ", median (IQR)", sep="")
        }else if(!grepl("\\(", sum_df$Value[i])){
            sum_df$Variable[i] <- paste(sum_df$Variable[i], ", mean \u00B1 std",
                                        sep="")
        }
    }
    rownames (sum_df) <- sum_df$Variable
    return(sum_df %>% select(!Variable))
}

chi2_test <- function(vec_list) {
    all_vec <- do.call(c, vec_list)
    pos_class <- unique(all_vec)[1]
    pos_neg <- function(vec) {
        vec_pos <- sum(vec == pos_class)
        vec_neg <- length(vec) - vec_pos
        return(data.frame("pos" = vec_pos, "neg" = vec_neg))
    }
    vec_list <- vec_list %>% lapply(pos_neg)
    vec_df <- do.call(rbind, vec_list)
    return(chisq.test(t(vec_df)))
}


continuous_test <- function(vec_list) {
    all_vec <- do.call(c, vec_list)
    normality <- shapiro.test(all_vec[1:5000])$p.value
    N <- length(vec_list)

    if (N == 2) {
        if (normality > 0.05) {
            pval <- t.test(vec_list[[1]], vec_list[[2]])$p.value
        }else {
            pval <- wilcox.test(vec_list[[1]], vec_list[[2]])$p.value
        }
    } else {
        samples <- lapply(as.list(1:N), function(x) {
                        M <- length(vec_list[[x]])
                        paste("sample", x) %>% rep(M) %>% return()
                        })
        vec_df <- data.frame(values=all_vec, samples=do.call(c, samples))
        if (normality > 0.5) {
            pval <- oneway.test(values ~ samples, vec_df)$p.value
        } else {
            pval <- kruskal.test(values ~ samples, vec_df)$p.value
        }
    }
    return(pval)
}


sum_stat_pval <- function(df_list_ori) {
    sum_list <- list()
    conv_names <- names(get_conv())
    df_list <- df_list_ori %>% lapply(
        function(x) {
            x %>% dplyr::select(dplyr::all_of(conv_names))
        })
    for (i in colnames(df_list[[1]])){
        if (i %in% colnames(df_list_ori[[1]])) {
            if (!is.character(df_list[[1]][, i])) {
                # percentages for 2-level features
                if (mean(is.na(df_list[[1]][, i])) < 1) {
                    # convert factors to numeric
                    sum_cols <- df_list %>% lapply(
                        function(x) {x[, i] %>% factor2numeric()})
                    all_sum <- do.call(c, sum_cols)
                    if (length(unique(all_sum)) == 2) {
                        pval <- chi2_test(sum_cols)$p.value
                    }else {
                        pval <- continuous_test(sum_cols)
                    }
                } else {pval <- NA}
                sum_list [[i]] <- pval
            }
        } else {sum_list[[i]] <- NA}
    }
    sum_list <- do.call(c, sum_list) %>% round(digit = 3) %>% 
        format(nsmall = 3) %>% as.character()
    for (i in 1:length(sum_list)) {
        sum_list[i] <- ifelse(sum_list[i] == "0", "<0.001", sum_list[i])
    }
    sum_list <- data.frame(pval = sum_list)
    # rownames(sum_list) <- get_conv()
    return(sum_list)
}


scale2 <- function(x) (x - mean(x, na.rm=T)) / sd(x, na.rm=T)

run_one_logistic <- function(one_df, yvar, covars, method="logistic", invert=FALSE,
                             estimate="lingress") {
    if (length(covars) != 0) {
        formu <- paste(covars, collapse="+")
        formu <- paste("+", formu)
    } else { formu <- ""}
    # Colname "Value": represent global efficiency or tract number
    if (invert) {
        formu <- as.formula(paste("Value ~", yvar, formu))
    } else {formu <- as.formula(paste(yvar, " ~ Value ", formu))}
    include_var <- c(covars, "Value")
    # scale the y variable only in linear regression, but not logistic
    # regression
    if (method != "logistic") {
        include_var <- c(include_var, yvar)
    } else {
        one_df[, yvar] <- as.factor(one_df[, yvar])
    }
    one_df <- one_df %>% mutate_at(include_var, scale2)

    if (method == "logistic" && invert==FALSE) {
        model <- glm(formu, data=one_df, family="binomial")
    } else {
        model <- lm(formu, data=one_df)
    }
    sum_tab <- summary(model)$coefficients

    if ("Pr(>|t|)" %in% colnames(sum_tab)) {
        colnames(sum_tab)[colnames(sum_tab) == "Pr(>|t|)"] <- "Pr(>|z|)"
    }

    # save the slope and p value
    out_tab <- data.frame(sum_tab) %>% 
        select(Pr...z.., Estimate) %>%
        rename(p=Pr...z..) %>%
        rownames_to_column("Variable") %>%
        filter(Variable != "(Intercept)")
        # mutate(padj = p.adjust(p, method="fdr")) %>%

    if (invert == TRUE) {
        out_tab <- out_tab %>% filter(Variable == paste(yvar, "1", sep=""))
    } else {out_tab <- out_tab %>% filter(Variable == "Value")}

    if (estimate == "correlation") {
        cor_val <- cor(one_df$Value, one_df %>% pull(!!as.symbol(yvar)))
        out_tab$Estimate <- cor_val
    }
    out_tab %>% select(!Variable) %>% return()
}

#' @param thres at least how many percentage of participants without missing
#' value would the variable be considered for regression 
run_logistic_modules <- function(one_df, yvar, covars, method="logistic",
                                 invert=FALSE, thres=0.5) {
    all_modules <- one_df %>% pull(Atlas) %>% unique()
    all_out <- as.list(all_modules) %>% lapply(function(xx) {
        # Submodules have the colname "Atlas"
        one_df_sub <- one_df %>% filter(Atlas == xx)
        all_na <- mean(is.na(one_df_sub$Value))
        all_zero <- sum(one_df_sub$Value, na.rm=T)
        uniq_len <- length(unique(one_df_sub[!is.na(one_df_sub$Value), yvar]))
        if (all_na < thres && uniq_len > 1 && all_zero > 0) {
            one_out <- run_one_logistic(one_df_sub, yvar, covars, method=method,
                                        invert=invert)
            one_out$Module <- xx
            return(one_out)
        }
    })
    all_out <- do.call(rbind, all_out)
    # Perform FDR correction across the submodules (plus global graph metrics)
    all_out$padj <- p.adjust(all_out$p, method="fdr")
    return(all_out)
}

#' Regression analysis with covariate adjustment
#'
#' @description This function uses a particular metric of interest, either
#' global efficiency of networks, or the number of tracts. This metric can be
#' regional or global. If regional, multiple testing correction will be
#' performed across all the regions using FDR. All input variables are scaled so
#' normalized coefficients will be obtained.
#'
#' @param one_df dataframe with all the data. The metric of interest must have
#' the name "Value". Each region where the metric is obtained must have the name
#' "Atlas".
#' @param yvar which y variable to associate to.
#' @param covar_list an R list of covariate. Each list contains a set of
#' covariates to adjust.
#' @param method either "logistic" or "lingress".
#'
#' @example
#' covar_list <- list(c("Age", "Sex"), c("Age", "Sex", "Education"))
#' out_table <- run_logistic_covars(one_df, "Amyloid", covar_list, "lingress")
run_logistic_covars <- function(one_df, yvar, covar_list, method, invert=FALSE) {
    all_out <- covar_list %>% lapply(function(yy) {
        print(yy)
        one_out <- run_logistic_modules(one_df, yvar, yy, method=method,
                                        invert=invert)
        one_out$covar <- paste(yy, collapse=",")
        return(one_out)
    })
    all_out <- do.call(rbind, all_out)
    return(all_out)
}


#' Perform PCA on tracts
pc_tract <- function(all_df) {
    all_tracts <- unique(all_df$Atlas)
    all_IDs <- unique(all_df$ID)
    tract_mat <- as.list(all_tracts) %>% lapply(function(xx) {
        one_df <- all_df %>% filter(Atlas == xx)
        return(one_df[match(all_IDs, one_df$ID), "Value"])
    })

    tract_mat <- data.frame(do.call(cbind, tract_mat))
    colnames(tract_mat) <- all_tracts
    rownames(tract_mat) <- all_IDs

    pr_mod <- prcomp(tract_mat)
    prs <- pr_mod$x %>% data.frame() %>% rownames_to_column("ID") %>%
        select(ID, PC1, PC2, PC3)

    all_df %>% filter(Atlas == "Anterior Commissure") %>%
        merge(prs, by="ID") %>%
        select(!one_of(c("Atlas", "Value"))) %>%
        return()
}


#' Correlation of variables with each PC
pc_corr <- function(plot_df, variables) {
    all_pc <- as.list(c("PC1", "PC2", "PC3")) %>% lapply(function(xx) {
        one_pc <- as.list(variables) %>% lapply(function(yy) {
            if (yy %in% c("Diagnosis", "Sex", "AD")) {
                formu <- as.formula(paste(xx, " ~ ", yy, sep=""))
                p <- wilcox.test(formu, data=plot_df)$p.value
            } else {
                formu <- as.formula(paste(" ~ ", yy, " + ", xx, sep=""))
                p <- cor.test(formu, data=plot_df)$p.value
            }
            return(data.frame(PC=xx, Variable=yy, p=p))
        })
        one_pc <- do.call(rbind, one_pc)
        return(one_pc)
    })
    do.call(rbind, all_pc) %>% data.frame() %>%
        # group_by(PC) %>%
        mutate(padj = p.adjust(p, method="fdr")) %>%
        # ungroup() %>%
        return()
}


pc_corr_multi <- function(plot_df, variables, discrete=c("Sex", "AD")) {
    all_pc <- as.list(c("PC1", "PC2", "PC3")) %>% lapply(function(xx) {
        formu <- paste(variables, collapse="+")
        formu <- as.formula(paste(xx, "~", formu))
        if (!is.null(discrete)) {
            plot_df <- plot_df %>% mutate_at(discrete, as.factor)
        }
        mod <- summary(lm(formu, data=plot_df))
        mod$coefficients %>% data.frame() %>%
            rename(p=Pr...t..) %>%
            select(p, Estimate) %>%
            rownames_to_column("Variable") %>%
            filter(Variable != "(Intercept)") %>%
            mutate(PC = xx) %>%
            return()
    })
    do.call(rbind, all_pc) %>% data.frame() %>%
        group_by(Variable) %>%
        mutate(padj = p.adjust(p, method="fdr")) %>%
        ungroup() %>%
        select(Variable, PC, padj) %>%
        mutate(padj = round(padj, 3)) %>%
        mutate(Variable = gsub("1$", "", Variable)) %>%
        return()
}


get_pc_plots <- function(plot_df, variables,
                         discrete=c("Sex", "Diagnosis", "AD", "Hypertension"),
                         PCs=c("PC1", "PC2"),
                         var_key=NULL) {
    all_plots <- as.list(variables) %>% lapply(function(xx) {
        p <- plot_df %>%
            ggplot(aes_string(x = PCs[1], y= PCs[2], color = xx)) +
            geom_point(size = 2) +
            labs(color = "")

        if (!xx %in% discrete) {
            p <- p + scale_color_viridis_c()
        } else {
            p <- p + scale_color_viridis_d()
        }

        if (!is.null(var_key)) {
            p <- p + ggtitle(var_key[xx])
        } else {
            p < p + ggtitle(xx)
        }
        p <- p + theme_publication() + 
            theme(legend.key.width = unit(1.5, "cm"),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.ticks.y = element_blank(),
            )
        return (p)
    })
}

num2symbol <- function(one_int) {
    return(paste(rep("+", as.numeric(one_int)), collapse=""))
}

combine_results <- function(out_eoad, out_load, save_path, collapse_adj=NULL) {
    # order of modules in the final table
    module_order <- c("Global", "Default", "Frontoparietal", "Ventral attention",
                      "Dorsal attention", "Somatomotor", "Visual", "Limbic")
    out_eoad$Group <- "EO"
    out_load$Group <- "LO"
    all_out <- rbind(out_eoad, out_load) %>%
        # Capitulize first letters of the module names, remove underscores,
        # order the modules
        mutate(Module = gsub("HCP-MMP", "global", Module)) %>%
        mutate(Module = DescTools::StrCap(Module)) %>%
        mutate(Module = gsub("_", " ", Module)) %>%
        mutate(Module = factor(Module, levels=module_order))

    if (is.null(collapse_adj)) {
        all_out %>% 
            # print normalized slopes in 2 decimal places, append the star symbols
            # for p values
            mutate(Value = paste(sprint2(Estimate), gtools::stars.pval(padj))) %>%
            # mutate(Covar = ifelse(grepl("WMH_norm", covar), "WMH", "no_WMH")) %>%
            mutate(Covar = covar) %>%
            mutate(Type = paste(Group, Covar)) %>%
            select(Module, Type, Value) %>%
            spread(Type, Value) %>%
            arrange(Module) %>%
            fwrite(save_path)
    } else {
        estimate <- all_out %>% 
            filter(covar == collapse_adj) %>%
            mutate(Value = sprint2(Estimate)) %>%
            mutate(Type = paste(Group, Module, sep="__")) %>%
            select(Type, Value)

        sig_level <- all_out %>%
            mutate(Type = paste(Group, Module, sep="__")) %>%
            group_by(Type) %>%
            summarize(sig = as.character(sum(padj < 0.05))) %>%
            ungroup() %>%
            select(Type, sig)

        for (i in 1:dim(sig_level)[1]) {
            sig_level[i, "sig"] <- num2symbol(sig_level[i, "sig"])
        }
        
        estimate %>% merge(sig_level, by="Type") %>%
            separate(Type, c("Group", "Module"), sep="__") %>%
            mutate(Val = paste(Value, sig, sep=" ")) %>%
            # mutate(Val = gsub("_0", "", Val)) %>%
            select(Group, Module, Val) %>%
            spread(Group, Val) %>%
            mutate(Module = factor(Module, levels=module_order)) %>%
            arrange(Module) %>%
            fwrite(save_path)
    }
}

combine_results_cogn <- function(all_out, save_path) {
    # order of modules in the final table
    module_order <- c("Global", "Default", "Frontoparietal", "Ventral attention",
                      "Dorsal attention", "Somatomotor", "Visual", "Limbic")
    all_out %>% 
        # print normalized slopes in 2 decimal places, append the star symbols
        # for p values
        mutate(Value = paste(sprint2(Estimate), gtools::stars.pval(padj))) %>%
        mutate(Module = gsub("HCP-MMP", "global", Module)) %>%
        mutate(Type = paste(cognition, covar)) %>%
        select(Module, Type, Value) %>%
        spread(Type, Value) %>%
        # Capitulize first letters of the module names, remove underscores,
        # order the modules
        mutate(Module = DescTools::StrCap(Module)) %>%
        mutate(Module = gsub("_", " ", Module)) %>%
        mutate(Module = factor(Module, levels=module_order)) %>%
        arrange(Module) %>%
        fwrite(save_path)
}


run_one_colinear <- function(one_df, yvar, inter_val, covars) {
    if (length(covars) != 0) {
        formu <- paste(covars, collapse="+")
        formu <- paste("+", formu)
    } else { formu <- ""}
    formu <- as.formula(paste(yvar, " ~ Value *", inter_val, "+", formu))
    include_var <- c(covars, "Value")
    # scale the y variable only in linear regression, but not logistic
    # regression
    include_var <- c(include_var, yvar)
    one_df <- one_df %>% mutate_at(include_var, scale2)

    model <- lm(formu, data=one_df)
    sum_tab <- summary(model)$coefficients

    if ("Pr(>|t|)" %in% colnames(sum_tab)) {
        colnames(sum_tab)[colnames(sum_tab) == "Pr(>|t|)"] <- "Pr(>|z|)"
    }

    # save the slope and p value
    out_tab <- data.frame(sum_tab) %>% 
        select(Pr...z.., Estimate) %>%
        rename(p=Pr...z..) %>%
        rownames_to_column("Variable") %>%
        filter(Variable != "(Intercept)") %>%
        filter(!Variable %in% covars)
    out_tab %>% return()
}


run_copathology <- function(one_df, yvar, inter_val, covars, thres=0.5) {
    all_modules <- one_df %>% pull(Atlas) %>% unique()
    all_out <- as.list(all_modules) %>% lapply(function(xx) {
        # Submodules have the colname "Atlas"
        one_df_sub <- one_df %>% filter(Atlas == xx)
        all_na <- mean(is.na(one_df_sub$Value))
        all_zero <- sum(one_df_sub$Value, na.rm=T)
        uniq_len <- length(unique(one_df_sub[!is.na(one_df_sub$Value), yvar]))
        if (all_na < thres && uniq_len > 1 && all_zero > 0) {
            one_out <- run_one_colinear(one_df_sub, yvar, inter_val, covars)
            one_out$Module <- xx
            return(one_out)
        }
    })
    do.call(rbind, all_out) %>%
        group_by(Variable) %>%
        mutate(padj = p.adjust(p, method="fdr")) %>%
        ungroup() %>%
        return()
}


combine_results_copath <- function(out_eoad, out_load, save_path) {
    module_order <- c("Global", "Default", "Frontoparietal", "Ventral attention",
                      "Dorsal attention", "Somatomotor", "Visual", "Limbic")
    out_eoad$Group <- "EO"
    out_load$Group <- "LO"
    all_out <- rbind(out_eoad, out_load) %>%
        # Capitulize first letters of the module names, remove underscores,
        # order the modules
        mutate(Module = gsub("HCP-MMP", "global", Module)) %>%
        mutate(Module = DescTools::StrCap(Module)) %>%
        mutate(Module = gsub("_", " ", Module)) %>%
        mutate(Module = factor(Module, levels=module_order))

    all_out %>% 
        mutate(Value = paste(sprint2(Estimate), gtools::stars.pval(padj))) %>%
        mutate(Type = paste(Group, Variable)) %>%
        select(Module, Type, Value) %>%
        spread(Type, Value) %>%
        arrange(Module) %>%
        fwrite(save_path)
}
