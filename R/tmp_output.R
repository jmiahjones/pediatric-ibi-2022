glm <- T
# set output files
if(glm) {
  result_dir <- sprintf("./results/%s/glm-", "Main")
} else { 
  result_dir <- sprintf("./results/%s/", "Main")
}
boot_preds_file <- paste0(result_dir, "boot_preds.RData")
load(boot_preds_file)
pred_df <- data.frame(ID=record_ids, IBI=Y, pred=best_preds)
write.csv(pred_df, "./results/Main/tmp.csv", row.names = F)
