
rm(list=ls())

library(dplyr)
library(foreach)
library(SuperLearner)
library(doFuture)
plan(multicore, workers=3L)
registerDoFuture()

start <- Sys.time()


main_data_path <- "./data/IBI-1-28-22.csv"

raw <- read.csv(main_data_path,
                header=T)

source("R/output_funtions.R")

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

only_features = data_raw %>% 
  select(sex,
         cmc,
         age,
         gestationalage,
         illappearing,
         temp,
         numberdaysfever,
         numberdaysillness,
         coughpresent,
         urinarytractinflamm,
         leukesterase)

dim(only_features)

##Explore data set
missingness <- (apply(data_raw, 2, function(x) sum(is.na(x)) > 0))
if(any(missingness > 0))
  stop("Missing data found!")


# Encode factors
features_encoded <- data.frame(model.matrix(~., data = only_features)[,-1])
record_ids <- data_raw$recordid

IBI <- as.numeric(data_raw$IBI) - 1
IBI2 <- as.numeric(data_raw$IBI2) - 1
BI <- BI_full <- as.numeric(data_raw$BI) - 1

# Y <- IBI
BI <- as.numeric(data_raw$BI) - 1

weight_function <- function(w, Y, BI) {
  w[1]*BI*(1-Y) + w[2]*(1-BI) + w[3]*Y
}

W_configs <- matrix(c( 1, 1, 4,
                       1, 1, 2,
                       1, 1, 8,
                       2, 1, 4,
                       2, 1, 8,
                       2, 1, 16,
                       1, 1, 16), ncol=3, byrow=T)

num_w_conf <- nrow(W_configs)
source("R/SuperLearner_config.R")


analysis <- c("Main", "Sensitivity", "LE")
conf_stat_names <- c("Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value")
cut.pts <- c(.001, .005, .01, .03, .05)

void <- foreach(j=1:3) %dopar% {
# for(j in 1:3) {
  le_B <- j==3
  
  if(le_B){
    features_encoded <- data.frame(model.matrix(~., data = 
      only_features %>% select(-urinarytractinflamm))[,-1])
  } else {
    features_encoded <- data.frame(model.matrix(~., data = 
      only_features %>% select(-leukesterase))[,-1])
  }
  
  sens_outcome_B <- j==2
  
  # set output files
  result_dir <- sprintf("./results/%s/", analysis[j])
  results_file <- paste0(result_dir, "cv_predictions.rds")
  
  Y <- if(sens_outcome_B) { IBI2 } else { IBI }
  
  ##########################
  ## Fit SuperLearner ####
  ##
  pred_list <- vector("list", num_w_conf)
  for(w_idx in 1:num_w_conf) {
    this_w_config <- W_configs[w_idx,]
    W <- weight_function(this_w_config, Y, BI)
    # W <- 20*W/mean(W)
    
    trained_cv_sl = CV.SuperLearner(
      Y = BI,
      X = features_encoded,
      obsWeights = W,
      family = binomial(),
      method = "method.AUC",
      SL.library = SL.lib,
      cvControl = outer_control,
      innerCvControl = inner_controls,
      parallel = "multicore"
    )
    
    trained_cv_glm = CV.SuperLearner(
      Y = BI,
      X = features_encoded,
      family = binomial(),
      method = "method.AUC",
      SL.library = c("SL.glm"),
      cvControl = outer_control,
      innerCvControl = inner_controls
    )
    
    
    pred_list[[w_idx]] <- data.frame(
      Y=Y,
      SL=trained_cv_sl$SL.predict,
      GLM=trained_cv_glm$SL.predict,
      Weights = paste(this_w_config, collapse=",")
    )
    
  }
  
  saveRDS(pred_list, results_file)
  
  return(NULL)
}

stop <- Sys.time()
message(sprintf(
  "Training completed in %.1f minutes.",
  as.numeric(difftime(stop, start, units = "mins"))
))
