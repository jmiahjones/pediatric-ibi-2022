# Define the SuperLearner library and other configuration.

#--------------- SVM ---------------#
svmlrn <- create.Learner(
  "SL.ksvm",
  params = list(kernel="rbfdot", kpar="automatic", 
                cache = 1024
                # class.weights=c("0"=1, "1"=prop_wt)
  ),
  tune = list(C=2^seq(-2,2, length.out = 5))
)

#-------------- xgbTree --------------#

xgblrn <- create.Learner(
  "SL.xgboost", 
  # params=list(scale_pos_weight = prop_wt),
  tune = list(
    nrounds = c(100, 500, 1000),
    max_depth = 1:2
    # ,
    # eta = c(0.01,0.1),
    # gamma = c(1, 2, 3),
    # colsample_bytree = c(0.5, 1.0),
    # min_child_weight = c(0.5, 1),
    # subsample = c(0.8, 1.0)
  )
)

rflrn=create.Learner("SL.randomForest", 
                     tune = list(mtry = seq(6,11,2),
                                 ntrees=c(1000)))
SL.lib <- c(
  "SL.earth",
  "SL.gam",
  "SL.glm",
  "SL.glmnet",
  rflrn$names,
  xgblrn$names,
  svmlrn$names
)
outer_cv_num <- 10L
inner_cv_num <- 5L


# use the same folds to compare against
folds <- caret::createFolds(IBI2, k=outer_cv_num)
names(folds) <- paste0("Outer", names(folds))
inner_folds <- lapply(folds, function(outer_idx) {
  inner_Y <- IBI2[-outer_idx]
  inner_folds <- caret::createFolds(inner_Y, k=inner_cv_num)
  names(inner_folds) <- paste0("Inner", names(inner_folds))
  return(inner_folds)
})

outer_control <- SuperLearner.CV.control(V=outer_cv_num, stratifyCV=T, validRows = folds)
inner_controls <- lapply(1:outer_cv_num, function(idx){
  SuperLearner.CV.control(V=inner_cv_num, stratifyCV=T, validRows = inner_folds[[idx]])
})
