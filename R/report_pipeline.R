# report_pipeline.R
library(aucC)
library(dplyr)
library(doFuture)
plan(multicore, workers=2L)
registerDoFuture()
options(future.globals.maxSize=3*(1024^3))
start <- Sys.time()

rerun_boot <- F

### Getting Features ###
main_data_path <- "./data/IBI-1-28-22.csv"

raw <- read.csv(main_data_path,
                header=T)

## Pre-processing
library(dplyr)
data_raw <- raw %>%
  mutate(sex= factor(sex, labels = c("male", "female")),
         cm_condition = factor(cmc, labels = c("no", "yes")),
         age = as.numeric(age),
         gestationalage = as.numeric(gestationalage),
         illappearing = factor(illappearing, labels = c("well", "ill","unknown")),
         temp = as.numeric(temp),
         numdaysill = as.numeric(numberdaysillness),
         numberdaysfever = as.numeric(numberdaysfever),
         coughpresent = factor(coughpresent, labels = c("yes", "no","unknown" )),
         urinarytractinflamm = factor(urinarytractinflamm, 
                                      labels = c("no", "yes", "unknown")),
         leukesterase = factor(leukesterase, 
                               labels = c("no", "yes", "unknown")),
         IBI = factor(ibi, labels = c("no", "yes")),
         IBI2 = factor(ibi2, labels = c("no", "yes")),
         BI = factor(bacterialinfection, labels = c("no", "yes"))
  ) %>% 
  select(recordid, 
         sex,
         cmc,
         age,
         gestationalage,
         illappearing,
         temp,
         numberdaysfever,
         numberdaysillness,
         coughpresent,
         urinarytractinflamm,
         leukesterase,
         IBI, IBI2, BI)

data_raw <- data_raw %>% mutate(
  greater30 = (31 <= age & age <= 90 ),
  preterm=gestationalage <= 36,
  sens1=(!preterm) & ( 8 <= age & age <= 60 ),
  sens2=(!preterm) & greater30,
  wellappearing=illappearing == "well"
)
record_ids <- data_raw$recordid

analysis_1 <- which(data_raw$sens1 & data_raw$wellappearing)
analysis_2 <- which(data_raw$greater30)
analysis_3 <- which(data_raw$wellappearing)
########################

lr_fun <- function(sens, spec){
  if(1-spec < .Machine$double.eps) {
    lr_pos <- Inf
  } else {
    lr_pos <- sens/(1-spec)
  }
  
  if(spec < .Machine$double.eps) {
    lr_neg <- Inf
  } else {
    lr_neg <- (1-sens)/spec
  }
  
  return(
    c(LRp=lr_pos, LRn=lr_neg)
  )
}

config.selection <- function(boot.mat, max=TRUE){
  # accepts matrix of aucs, which we maximize
  summ.mat <- apply(boot.mat, 2, function(col){
    c(median(col, na.rm=T), min(col, na.rm=T))
  })
  
  if(max){
    best.med <- max(summ.mat[1,])
  } else {
    best.med <- min(summ.mat[1,])
  }
  best.idxs <- which(summ.mat[1,] == best.med)
  
  if(length(best.idxs) > 1){
    if(max){
      final.best.idx <- best.idxs[which.max(summ.mat[2, best.idxs])]
    } else {
      final.best.idx <- best.idxs[which.min(summ.mat[2, best.idxs])]
    }
  } else {
    final.best.idx <- best.idxs
  }
  return(final.best.idx)
}
B <- 5e4
analysis <- c("Main", "Sensitivity", "LE")
conf_stat_names <- c("Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value")
cut.pts <- seq(.005, .1, by=.005)

for(classifier_idx in 1:2){ 
  is_sl <- (classifier_idx == 1) # 1: SL, 2: GLM
  for(j in 1:3){
    le_B <- j==3
    sens_outcome_B <- j==2
    
    # set output files
    result_dir <- sprintf("./results/%s/", analysis[j])
    results_file <- paste0(result_dir, "cv_predictions.rds")
    if(!is_sl)
      result_dir <- sprintf("./results/%s/glm-", analysis[j])
    
    boot_preds_file <- paste0(result_dir, "boot_preds.RData")
    report_table_file <- paste0(result_dir, "validation_results.csv")
    sens1_table_file <- paste0(result_dir, "8to60_validation_results.csv")
    sens2_table_file <- paste0(result_dir, "30plus_validation_results.csv")
    sens5_table_file <- paste0(result_dir, "wellappearing_validation_results.csv")
    rec_table_file <- paste0(result_dir, "rec_validation_results.csv")
    
    roc_file <- paste0(result_dir, "other_auc_results.txt")
    
    
    missed_records_file <- paste0(result_dir, "false_negatives.csv")
    rec_missed_records_file <- paste0(result_dir, "rec_false_negatives.csv")
    low_risk_file <- paste0(result_dir, "low_risk.csv")
    
    auc_plot_file <- paste0(result_dir, "validation_roc.png")
    
    pred_list <- readRDS(results_file)
    num_w_conf <- length(pred_list)
    Y <- pred_list[[1]]$Y
    n <- length(Y)
    set.seed(2022, kind="L'Ecuyer-CMRG")
    bootIdxs <- replicate(B, sample(1:n, size=n, replace=T), simplify = F)
    outBootIdxs <- lapply(bootIdxs, function(boot_idx) setdiff(1:n, boot_idx))
    outBootIdxs <- lapply(outBootIdxs, function(boot_idx) 
      sample(boot_idx, size=length(boot_idx), replace=T))
    
    pred_arr_sl <- matrix(NA, nrow=nrow(pred_list[[1]]), ncol=num_w_conf)
    w_names <- rep("", num_w_conf)
    for(w_idx in 1:num_w_conf) {
      w_names[w_idx] <- pred_list[[w_idx]]$Weights[1]
      pred_arr_sl[,w_idx] <- if(is_sl){
        pred_list[[w_idx]]$SL
      } else {
        pred_list[[w_idx]]$GLM
      }
    }
    
    
    Y_fac <- factor(Y, levels=c(0,1), labels=c("No", "Yes"))
    obs_aucs_config <- apply(pred_arr_sl, 2, function(x) aucC_wrap(Y_fac, x))
    obs_best_idx <- config.selection(matrix(obs_aucs_config, nrow=1), max=T)
    best_preds <- pred_arr_sl[,obs_best_idx]
    obs_auc <- obs_aucs_config[obs_best_idx]
    
    if(rerun_boot) {
    tmp_output <- foreach(b=1:B) %dopar% {
      boot_idx <- bootIdxs[[b]]
      out_boot_idx <- outBootIdxs[[b]]
      
      boot_chars <- apply(pred_arr_sl[boot_idx,], 2, 
                          function(x) aucC_wrap(Y_fac[boot_idx], scores=x))
      best_config_idx <- config.selection(matrix(boot_chars, nrow=1), max=TRUE)
      auc <- boot_chars[best_config_idx]
      
      out_sample_preds <- pred_arr_sl[out_boot_idx,best_config_idx]
      # Full Data
      out_sample_auc <- aucC_wrap(Y_fac[out_boot_idx], out_sample_preds)
      these_boot_stats <- sapply(cut.pts, function(cut.pt){
        confusion_mat <- factor(out_sample_preds < cut.pt, levels=c(T, F), labels=c("No", "Yes")) %>% 
          table(Y_fac[out_boot_idx], exclude = F) %>%
          caret::confusionMatrix(positive="Yes")
        
        op_char <- confusion_mat$byClass[1:4]
        c(auc, op_char, lr_fun(sens=op_char[1], spec=op_char[2]))
      })
      
      # Sub-pop 1
      analysis_1_idx <- out_boot_idx[which(out_boot_idx %in% analysis_1)]
      out_sample_preds <- pred_arr_sl[analysis_1_idx,best_config_idx]
      analysis_1_auc <- aucC_wrap(Y_fac[analysis_1_idx], out_sample_preds)
      analysis_1_boot_stats <- sapply(cut.pts, function(cut.pt){
        confusion_mat <- factor(out_sample_preds < cut.pt, levels=c(T, F), labels=c("No", "Yes")) %>% 
          table(Y_fac[analysis_1_idx], exclude = F) %>%
          caret::confusionMatrix(positive="Yes")
        
        op_char <- confusion_mat$byClass[1:4]
        c(auc, op_char, lr_fun(sens=op_char[1], spec=op_char[2]))
      })
      
      # Sub-pop 2: Use Risk calculator; everyone else: high risk
      analysis_2_idx <- out_boot_idx
      out_sample_preds <- pred_arr_sl[,best_config_idx]
      out_sample_preds[-analysis_2] <- 1 # hard-code infants outside 31-90 as high-risk
      out_sample_preds <- out_sample_preds[analysis_2_idx]
      analysis_2_auc <- aucC_wrap(Y_fac[analysis_2_idx], out_sample_preds)
      analysis_2_boot_stats <- sapply(cut.pts, function(cut.pt){
        confusion_mat <- factor(out_sample_preds < cut.pt, levels=c(T, F), labels=c("No", "Yes")) %>% 
          table(Y_fac[analysis_2_idx], exclude = F) %>%
          caret::confusionMatrix(positive="Yes")
        
        op_char <- confusion_mat$byClass[1:4]
        c(auc, op_char, lr_fun(sens=op_char[1], spec=op_char[2]))
      })
      
      # Sub-pop 3
      analysis_3_idx <- out_boot_idx[which(out_boot_idx %in% analysis_3)]
      out_sample_preds <- pred_arr_sl[analysis_3_idx,best_config_idx]
      analysis_3_auc <- aucC_wrap(Y_fac[analysis_3_idx], out_sample_preds)
      analysis_3_boot_stats <- sapply(cut.pts, function(cut.pt){
        confusion_mat <- factor(out_sample_preds < cut.pt, levels=c(T, F), labels=c("No", "Yes")) %>%
          table(Y_fac[analysis_3_idx], exclude = F) %>%
          caret::confusionMatrix(positive="Yes")
        
        op_char <- confusion_mat$byClass[1:4]
        c(auc, op_char, lr_fun(sens=op_char[1], spec=op_char[2]))
      })
      
      these_best_configs <- best_config_idx
      return(list(these_boot_stats, these_best_configs, 
                  analysis_1_boot_stats, analysis_2_boot_stats,
                  analysis_3_boot_stats))
    }
    save.image(file=boot_preds_file)
    } else {
    load(boot_preds_file)
    }

    create_table <- function(boot_loc=1, analysis_idxs, output_file) {
      if(boot_loc == 4) {# analysis 2, use all data
	the_preds <- best_preds
	the_preds[-analysis_idxs] <- 1
	the_Y <- Y_fac
      } else {
        the_preds <- best_preds[analysis_idxs]
        the_Y <- Y_fac[analysis_idxs]
      }
      obs_stats <- sapply(cut.pts, function(cut.pt){
        confusion_mat <- factor(the_preds < cut.pt,
                                levels=c(T, F), labels=c("No", "Yes")) %>% 
          table(the_Y, exclude = F) %>%
          caret::confusionMatrix(positive="Yes")
        
        op_char <- confusion_mat$byClass[1:4]
        c(obs_auc, op_char, lr_fun(sens=op_char[1], spec=op_char[2]))
      })
      dimnames(obs_stats) <- list(c("AUC", conf_stat_names, "LR+", "LR-"), cut.pts)
      
      boot_stats <- array(dim=c(7, length(cut.pts), B),
                          dimnames = list(c("AUC", conf_stat_names, "LR+", "LR-"),
                                          cut.pts,
                                          paste0("B", 1:B)))
      for(b in 1:B) {
        boot_stats[,,b] <- tmp_output[[b]][[boot_loc]]
      }
      cis <- matrix(NA, ncol=length(cut.pts), nrow=dim(obs_stats)[1])
      for(i in 1:nrow(cis)) {
        for(j in 1:length(cut.pts)) {
          ci <- quantile(boot_stats[i,j,], probs=c(.025, .975), na.rm=T)
          cis[i,j] = sprintf("%.3f (%.3f, %.3f)", obs_stats[i,j], ci[1], ci[2])
        }
      }
      dimnames(cis) <- dimnames(obs_stats)
      write.csv(cis, file=output_file)
      invisible()
    }
    create_table(1, 1:n, report_table_file)
    create_table(3, analysis_1, sens1_table_file)
    create_table(4, analysis_2, sens2_table_file)
    create_table(5, analysis_3, sens5_table_file)
    
  }
  
  stop <- Sys.time()
  message(sprintf("Completed bootstrap %i in %.2f minutes.", classifier_idx,
                  difftime(stop, start, units="mins")))
  
}
