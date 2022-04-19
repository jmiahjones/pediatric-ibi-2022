# additional output - AUC plot plus predictions for xlsx filtering sheet

library(aucC)
library(dplyr)
library(doFuture)
plan(multicore, workers=2L)
registerDoFuture()
options(future.globals.maxSize=3*(1024^3))
start <- Sys.time()

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

j <- 1
le_B <- j==3
sens_outcome_B <- j==2

# set output files
result_dir <- sprintf("./results/%s/", analysis[j])
results_file <- paste0(result_dir, "cv_predictions.rds")

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

pred_arr_sl <- matrix(NA, nrow=nrow(pred_list[[1]]), ncol=num_w_conf)
w_names <- rep("", num_w_conf)
for(w_idx in 1:num_w_conf) {
  w_names[w_idx] <- pred_list[[w_idx]]$Weights[1]
  pred_arr_sl[,w_idx] <- pred_list[[w_idx]]$SL
}

Y_fac <- factor(Y, levels=c(0,1), labels=c("No", "Yes"))
obs_aucs_config <- apply(pred_arr_sl, 2, function(x) aucC_wrap(Y_fac, x))
obs_best_idx <- config.selection(matrix(obs_aucs_config, nrow=1), max=T)
best_preds <- pred_arr_sl[,obs_best_idx]
obs_auc <- obs_aucs_config[obs_best_idx]

roc <- pROC::roc(Y_fac, best_preds)


png(auc_plot_file, width=1200, height=1200, pointsize = 36)
plot(roc, col="red",
     xlab="Specificity (1 - False Positive %)",
     ylab="Sensitivity (True Positive %)", lwd=6)

load(boot_preds_file)
boot_stats <- array(dim=c(7, length(cut.pts), B),
                          dimnames = list(c("AUC", conf_stat_names, "LR+", "LR-"),
                                          cut.pts,
                                          paste0("B", 1:B)))
for(b in 1:B) {
  boot_stats[,,b] <- tmp_output[[b]][[1L]]
}

lims <- quantile(boot_stats[1,1,], probs=c(.025, .975))
legend(x=0.55, y=0.15,
       legend=sprintf(paste0("AUC: %0.2f CI: (%0.2f, %0.2f)"), obs_auc, lims[1], lims[2]),
       col="red",
       lty=1, lwd=6,
       cex=0.9)
dev.off()



pred_df <- data.frame(ID=record_ids, IBI=Y, SL=best_preds)
write.csv(pred_df, "./results/Main/tmp.csv", row.names = F)
