# Output Functions

create_report_table <- function(
  the_df, cut.pts, table_file,
  clNum=20,
  B=2000L, # number of bootstraps
  alpha=0.05, # 1-confidence level
  digit_format="%.3f", # number of digits
  conf_stat_names=c('Sensitivity', 'Specificity', 'Pos Pred Value', 'Neg Pred Value')
){
  
  sprintf_ci_format <- paste0("(", digit_format, ", ", digit_format, ")")
  sprintf_est_ci_format <- paste0(digit_format, " %s")
  
  
  if(clNum > 1){
    cl <- parallel::makeCluster(clNum, "FORK")
    registerDoParallel(cl)
  } else {
    registerDoSEQ()
  }
  
  tbl <- foreach(cut.point=cut.pts, .combine=cbind) %do% {
    
    sl_conf <- the_df %>%
      select(SL, BI) %>%
      mutate(SL=factor(1*(SL > cut.point), levels=0:1)) %>%
      table %>%
      caret::confusionMatrix(positive="1")
    sl_conf <- sl_conf$byClass[conf_stat_names]
    
    glm_conf <- the_df %>%
      select(GLM, BI) %>%
      mutate(GLM=factor(1*(GLM > cut.point), levels=0:1)) %>%
      table %>%
      caret::confusionMatrix(positive="1")
    glm_conf <- glm_conf$byClass[conf_stat_names]
    
    ret <- c(sl_conf, glm_conf)
    names(ret) <- lapply(c("SL.", "GLM."),
                         paste0, conf_stat_names) %>% do.call(c, .)
    
    
    sl.lrp <- ret["SL.Sensitivity"]/(1-ret["SL.Specificity"])
    glm.lrp <- ret["GLM.Sensitivity"]/(1-ret["GLM.Specificity"])
    sl.lrn <- (1-ret["SL.Sensitivity"])/ret["SL.Specificity"]
    glm.lrn <- (1-ret["GLM.Sensitivity"])/ret["GLM.Specificity"]
    lr <- c(sl.lrp, glm.lrp, sl.lrn, glm.lrn)
    names(lr) <- c("SL.LR+", "GLM.LR+",
                   "SL.LR-", "GLM.LR-")
    
    ret <- c(ret, lr)
    
    return(ret)
  }
  
  ci_tbl <- foreach(cut.point=cut.pts, .combine=cbind) %do% {
    foreach(b=1:B, .combine=cbind) %dopar% {
      
      set.seed(34598+b)
      boot <- sample(nrow(the_df), replace=T)
      boot_df <- the_df[boot,]
      sl_conf <- boot_df %>%
        select(SL, BI) %>%
        mutate(SL=factor(1*(SL > cut.point), levels=0:1)) %>%
        table %>%
        caret::confusionMatrix(positive="1")
      sl_conf <- sl_conf$byClass[conf_stat_names]
      
      glm_conf <- boot_df %>%
        select(GLM, BI) %>%
        mutate(GLM=factor(1*(GLM > cut.point), levels=0:1)) %>%
        table %>%
        caret::confusionMatrix(positive="1")
      glm_conf <- glm_conf$byClass[conf_stat_names]
      
      ret <- c(sl_conf, glm_conf)
      names(ret) <- lapply(c("SL.", "GLM."),
                           paste0, conf_stat_names) %>% do.call(c, .)
      
      
      sl.lrp <- ret["SL.Sensitivity"]/(1-ret["SL.Specificity"])
      glm.lrp <- ret["GLM.Sensitivity"]/(1-ret["GLM.Specificity"])
      sl.lrn <- (1-ret["SL.Sensitivity"])/ret["SL.Specificity"]
      glm.lrn <- (1-ret["GLM.Sensitivity"])/ret["GLM.Specificity"]
      lr <- c(sl.lrp, glm.lrp, sl.lrn, glm.lrn)
      names(lr) <- c("SL.LR+", "GLM.LR+",
                     "SL.LR-", "GLM.LR-")
      
      ret <- c(ret, lr)
      ret[is.nan(ret)] <- 1e3
      
      return(ret)
    } %>% apply(1, quantile, probs=c(alpha/2, 1-(alpha/2))) %>% t %>% 
      apply(1, function(x) sprintf(sprintf_ci_format, x[1], x[2]))
  }
  
  stopifnot(all.equal(dim(tbl), dim(ci_tbl)))
  final_tbl <- foreach(a=1:nrow(ci_tbl), .combine=rbind) %:% 
    foreach(b=1:ncol(ci_tbl), .combine=cbind) %do% 
    {
      sprintf(sprintf_est_ci_format, tbl[a,b], ci_tbl[a,b])
    }
  colnames(final_tbl) <- cut.pts
  rownames(final_tbl) <- rownames(tbl)
  write.csv(final_tbl, file=table_file)
  if(clNum > 1)
    stopCluster(cl)
  
  Sys.sleep(2)
  
}

print_roc_results <- function(
  the_df, ci_format="%.2f, (%.2f, %.2f)"
) {
  GLM.roc <- pROC::roc(BI ~ GLM, data=the_df)
  SL.roc <- pROC::roc(BI ~ SL, data=the_df)
  GLM.ci <- (pROC::ci.auc(GLM.roc))
  SL.ci <- (pROC::ci.auc(SL.roc))
  
  paste0(
    "GLM: ",
    sprintf(ci_format,
            GLM.ci[2], GLM.ci[1], GLM.ci[3]),
    "\nSL: ",
    sprintf(ci_format,
            SL.ci[2], SL.ci[1], SL.ci[3])
  )
  
}
