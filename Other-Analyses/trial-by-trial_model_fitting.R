
#--------------------------------------------------------------------------#
#                              FIT DELTA MODEL                             #
#--------------------------------------------------------------------------#

movingwindow <- function(model_param, obs, distort = FALSE){
  window_size <- ceiling(abs(model_param['window_size']))
  get_window_ends <- function(end_i){
    start_i <- end_i - window_size + 1
    if(start_i < 1) start_i <- 1
    return(c(start_i, end_i))
  }
  windows <- lapply(1:length(obs), get_window_ends)
  
  estimates <- sapply(windows, function(window_ends){
    window_indices <- window_ends[1]:window_ends[2]
    mean(obs[window_indices])
  })
  return(estimates)
}

deltamodel = function(model_param, obvs, distort = FALSE) {
  learning_rate <- model_param['learning_rate']
  starting_est <- model_param['starting_est']
  new = c()
  for(obv in 1:length(obvs)) {
    starting_anew <- obv==1 | obv==334 | obv==667
    if (starting_anew) {
      new[obv] <- starting_est + learning_rate*(obvs[obv]-starting_est)
    } else {
      new[obv] <- new[obv-1] + learning_rate*(obvs[obv] - new[obv-1])
      if (distort){new[obv] <- rnorm(1,new[obv-1],0.25) + learning_rate*(obvs[obv] - new[obv-1])}
    }
  }
  return(new)
} 

kullback_leibler <- function(p_o, phat_g){
  if(p_o <= 0) p_o <- 0.0000000001
  else if (p_o >= 1) p_o <- 0.999999999
  if(phat_g <= 0) phat_g <- 0.0000000001
  else if (phat_g >= 1) phat_g <- 0.999999999
  p_o * log(p_o/phat_g) + (1 - p_o)*log((1-p_o)/(1-phat_g))
}

threshold_model <- function(obsv, learning_method, threshold_method, model_param,
                            distort = FALSE){
  threshold <- abs(model_param['threshold'])
  #model_param['learning_rate'] <- abs(model_param['learning_rate'])
  #assume the mental model is updated every trial according to the delta learning model
  #generate delta estimates and store as list "model_ests"
  model_ests <- switch (learning_method,
    "delta" = deltamodel(model_param, obsv, distort),
    "movingwindow" = movingwindow(model_param, obsv, distort)
  )
  #create a second list which keeps track of actual estimates, which may hold steady
  #rather than be updated every trial
  thresholded_ests <- model_ests
  #run through every trial, changing a response to the last trial's respon
  #responses remain steady when the threshold is not met
  #if threshold is met, change running_est to that trial's delta est
  running_est <- model_ests[1]
  threshold_exceeded <- FALSE
  #switch statement that, according to the threshold criteria used, sets the 
  #procedures (e.g. what expression evaluates to true if threshold is exceeded)
  #unusual actions are stored as expressions and called only if they exist
  switch(threshold_method,
         "none" = {  
           threshld_met <- TRUE
         },
         "drift" = {
           #drift diffusion has to sum up the error over multiple trials
           #so i need to make a variable called "running_error" and accumulate
           #it on every trial 
           #and return it to zero on trials where the threshold is exceeded
           running_error <- 0
           trial_proc <- expression(running_error <- running_error + trial_error)
           threshld_met <- expression(abs(running_error) >= threshold)
           threshld_met_action <- expression(running_error <- 0)
         },
         "meandiff" = {
           threshld_met <- expression(abs(running_est - thresholded_ests[est]) >= threshold)
           },
         "kl" = {
           num_red <- 0
           num_obs <- 0
           trial_proc <- expression({num_red = num_red + obsv[est]; num_obs = num_obs + 1})
           threshld_met <- expression(if (num_obs > 0) {kullback_leibler(running_est, thresholded_ests[est]) > threshold} else {FALSE})
           threshld_met_action <- expression({num_red <- 0;num_obs <- 0})
          })
  for (est in 1:length(thresholded_ests)) {
    trial_error <- running_est - thresholded_ests[est]
    if(exists("trial_proc")){eval(trial_proc)}
    threshold_exceeded <- eval(threshld_met)
    if(!is.logical(threshold_exceeded) | is.na(threshold_exceeded)){
      print(running_est)
      print(model_param)
      }
    if(threshold_exceeded){
      running_est <- thresholded_ests[est]
      if(exists("threshld_met_action")){eval(threshld_met_action)}
    }
    thresholded_ests[est] <- running_est
  }
  return(thresholded_ests)  
}

get_ppt_fit_error_function <- function(pptid, learning_method, threshold_method, fix_params = c()) {
  ppt_rows <- task_df$PptId==pptid
  pptests <- task_df$bernoulliEstimate[ppt_rows]
  obsv <- task_df$TrialRingBoolean[ppt_rows]
  force(fix_params)
  return(
    function(param) {
      fixed_p <- names(fix_params)
      param[fixed_p] <- fix_params
      
      modelests <- threshold_model(obsv, learning_method, threshold_method, param)
      error <- modelests - pptests
      return(sum(error**2))
    }
  )
}

#--------------------------------------------------------------------------#
#                                 FITTING                                  #
#--------------------------------------------------------------------------#
param_startval <- function(param_name, threshold_method){
  switch(param_name,
         "learning_rate" = abs(rnorm(1,0.087,0.05)),
         "starting_est" = 0.5,
         "window_size" = round(runif(1,1,100)),
         "threshold" = 0#switch(threshold_method,
                              #"none" = 0,
                              #"drift" = abs(rnorm(1,0,0.5)),
                              #"meandiff" = runif(1),
                              #"kl" = runif(1,0,2))
  )
}

get_param_start_val <- function(learning_method, threshold_method, param_component){
  switch(learning_method,
         "delta" = {
           switch (param_component,
             "param_val" = param_start_vals <- c(learning_rate = param_startval("learning_rate",threshold_method),
                                                 starting_est = param_startval("starting_est",threshold_method),
                                                 threshold = param_startval("threshold",threshold_method)),
             "whether_fixed" = param_start_vals <- c(learning_rate = FALSE,
                                                     starting_est = FALSE,
                                                     threshold = switch(threshold_method,
                                                                        "none" = TRUE,
                                                                        FALSE))
           )
         }, 
         "movingwindow"= {
           switch (param_component,
                   "param_val" = param_start_vals <- c(window_size = param_startval("window_size",threshold_method),
                                                       threshold = param_startval("threshold",threshold_method)),
                   "whether_fixed" = param_start_vals <- c(window_size = FALSE,
                                                           threshold = switch(threshold_method,
                                                                              "none" = TRUE,
                                                                              FALSE)))
         })
  return(param_start_vals)
}

start_ppt_params <- function(learning_method, threshold_method){
  ppt_params <- list()
  ppt_params[['values']] <- sapply(ppt_ids, function(pptid){get_param_start_val(learning_method, threshold_method, "param_val")})
  ppt_params[['fixed']]<- sapply(ppt_ids, function(pptid){get_param_start_val(learning_method, threshold_method, "whether_fixed")})
  return(ppt_params)
}

update_ppt_params <- function(params, new_vals, fit_all_param = FALSE){
  uniform_reset <- function(data_structure, replace_with){
    sapply(data_structure, function(item_in_structure){
      item_in_structure <- replace_with
    })
  }
  reset_to_false <- function(pptid){
    ppt_params <- params[['fixed']][,pptid]
    ppt_params <- uniform_reset(ppt_params,FALSE)
  }
  #if we want to set ALL 3 parameters, then we have to set "fixed" to false for everyone
  if (fit_all_param){
    params[['fixed']] <- sapply(ppt_ids, reset_to_false)
    }
  #if we are not to reset every "fixed" value to FALSE then invert values
  #what was previously taken as fixed (threshold being fixed at zero) is now free to vary
  #what was previously taken as free to vary (the two other parameters to be fit) is now fixed
  else {params[['fixed']] <- !params[['fixed']]}
  params[['values']] <- sapply(ppt_ids, function(pptid){new_vals[[pptid]]})
  return(params)
}

ppt_fit_model <- function(pptid, fixed_par, unfixed_par, learning_method, threshold_method, optim_method){
  error <- get_ppt_fit_error_function(pptid, learning_method, threshold_method, fixed_par)
  upper_bounds <- c(learning_rate=0.5, starting_est = 1, window_size=100, threshold=10)
  lower_bounds <- c(learning_rate=0, starting_est = 0, window_size=0, threshold=0)
  upper_bounds <- upper_bounds[names(upper_bounds) %in% names(unfixed_par)]
  lower_bounds <- lower_bounds[names(lower_bounds) %in% names(unfixed_par)]
  if (optim_method == "GenSA"){
    optim <- GenSA(par=unfixed_par, fn=error, lower=lower_bounds, upper_bounds)
  } else {
    optim <- optim(error, par=unfixed_par, method=optim_method)
  }
  return(c(optim$par, fixed_par))
}

ppt_model_ests <- function(pptid, ppt_params){
  ppt_rows <- task_df$PptId==pptid
  obsv <- task_df$TrialRingBoolean[ppt_rows]
  ests <- threshold_model(obsv, learning_method, threshold_method, ppt_params)
  return(ests)
}

get_num_changepoints <- function(ppt_ests){
  #ppt_ests <- round(ppt_ests,2)
  cp <- diff(ppt_ests) != 0
  return(sum(cp))
}

get_sum_sqrd_error <- function(ppt_ests, model_ests){
  stand_dev <- sd(ppt_ests - model_ests)
  model_error <- sum((model_ests - ppt_ests)**2)
  return(log(model_error))
}

get_loglik <- function(ppt_ests,model_ests){
  error <- ppt_ests - model_ests
  stand_dev <- sqrt(sum(error**2))
  prob_of_error <- dnorm(error, 0, stand_dev)
  log_of_error <- log(prob_of_error)
  return(sum(log_of_error))
}

model_fit_variables <- function(ppt_fit, ppt_ests, ppt_obsv, learning_method, threshold_method){
  model_ests <- threshold_model(ppt_obsv, learning_method, threshold_method, ppt_fit)
  num_cp <- get_num_changepoints(model_ests)
  loglik <- get_loglik(ppt_ests, model_ests)
  sqrd_error <- get_sum_sqrd_error(ppt_ests, model_ests)
  return(list(
    params = ppt_fit,
    model_ests = model_ests,
    num_cp = num_cp,
    loglik = loglik,
    sqrd_err = sqrd_error
  ))
}

get_ppt_fit_vals <- function(pptid, learning_method, threshold_method, ppt_params, optim_method, csv_filename = NA){
  #print(count_up())
  print(pptid)
  ppt_ests <- task_df$bernoulliEstimate[task_df$PptId==pptid]
  ppt_obsv <- task_df$TrialRingBoolean[task_df$PptId==pptid]
  
  whether_fixed <- unlist(ppt_params$fixed[,pptid])
  fixed_par <- unlist(ppt_params$values[,pptid])[whether_fixed] #ppt_params$values[whether_fixed,pptid]
  unfixed_par <- unlist(ppt_params$values[,pptid])[!whether_fixed]
  
  ppt_fit_params <- ppt_fit_model(pptid, fixed_par, unfixed_par, learning_method, threshold_method, optim_method)
  ppt_fit <- model_fit_variables(ppt_fit_params, ppt_ests, ppt_obsv, learning_method, threshold_method)
  
  if(!is.na(csv_filename)){
    temp <- as.matrix(ppt_fit)
    colnames(temp) <- pptid
    temp <- as.list(temp)
    append_fits_csv(csv_filename, temp, learning_method, threshold_method)
  }
  return(ppt_fit)
}

append_fits_csv <- function(csv_filename,ppt_fit,learning_method, threshold_method){
  newcsv <- turn_fits_to_csv(ppt_fit,learning_method, threshold_method)
  if(file.exists(csv_filename)){
    oldcsv <- read.csv(csv_filename)
    newcsv <- rbind(oldcsv,newcsv)
  }
  write.csv(newcsv,csv_filename,row.names = FALSE)
}

fit_many_times <- function(learning_method, threshold_method, times, optim_method, pptids,
                           ppt_params = start_ppt_params(learning_method, threshold_method)){
  timer <- start_timer()
  fits <- c()
  for (i in 1:times) {
    print(count_up())
    fits[[i]] <- sapply(pptids, get_ppt_fit_vals, learning_method, threshold_method, ppt_params, optim_method)
  }
  print(timer())
  return(fits)
}

find_best_fit <- function(learning_method, threshold_method, n_times, optim_method, pptids=ppt_ids,  ...){
  many_fits <- fit_many_times(learning_method,threshold_method,n_times,optim_method, pptids,  ...)
  logliks <- sapply(many_fits, function(fit){fit['loglik',]})
  m_logliks <- matrix(as.numeric(logliks), ncol=n_times)
  indices_of_max_loglik <- apply(X=m_logliks,MARGIN=1,FUN=function(ppt_logliks){which(ppt_logliks==max(ppt_logliks))})
  names(indices_of_max_loglik) <- colnames(many_fits[[1]])
  best_fits <- sapply(pptids,function(pptid){
    best_fit_index <- indices_of_max_loglik[[pptid]][1]
    fit_vals <- many_fits[[best_fit_index]][,pptid]
    fit_vals
  })
  return(best_fits)
}

turn_fits_to_csv <- function(fits, learning_method, threshold_method){
  csv_df <- t(as.data.frame(fits['params',]))
  csv_df <- data.frame(csv_df)
  #print(csv_df)
  csv_df$PptId <- colnames(fits)#sapply(rownames(csv_df),substr,2,50)
  csv_df$learning_method <- learning_method
  csv_df$threshold_method <- threshold_method
  switch(learning_method,
         "delta" = {csv_df$window_size <- NA},
         "movingwindow" = {csv_df$learning_rate <- NA; csv_df$starting_est <- NA})
  csv_df <- add_loglik_to_csv(csv_df)
  return(csv_df)
}

turn_csv_to_fits_obj <- function(csv){
  lms <- unique(csv$learning_method)
  tms <- unique(csv$threshold_method)
  fits <- list()
  for(learning_method in lms){
    lm_df <- csv[csv$learning_method==learning_method,]
    switch(learning_method,
           "delta" = {param_names <- c("starting_est","learning_rate","threshold")},
           "movingwindow" <- {param_names <- c("window_size","threshold")})
    for (threshold_method in tms) {
      fits[[learning_method]][[threshold_method]] <- sapply(ppt_ids, function(pptid){
        row <- lm_df$PptId==pptid & lm_df$learning_method==learning_method & lm_df$threshold_method==threshold_method
        param <- unlist(lm_df[row,param_names])
        ppt_ests <- task_df$bernoulliEstimate[task_df$PptId==pptid]
        ppt_obsv <- task_df$TrialRingBoolean[task_df$PptId==pptid]
        ppt_fit <- model_fit_variables(param, ppt_ests, ppt_obsv, learning_method, threshold_method)
        return(ppt_fit)
      })
    }
  }
  return(fits)
}

fix_loglik <- function(csv){
  csv <- csv[,colnames(csv)[colnames(csv)!="loglik"]]
  csv <- add_loglik_to_csv(csv)
  return(csv)
}
add_loglik_to_csv <- function(csv){
  csv$loglik <- NA
  for (row in 1:nrow(csv)){
    param <- unlist(c(csv[row,c("threshold","starting_est","window_size","learning_rate")]))
    param <- param[!is.na(param)]
    pptid <- csv[row,"PptId"]
    threshold_method <- csv[row,"threshold_method"]
    learning_method <- csv[row,"learning_method"]
    ppt_ests <- task_df$bernoulliEstimate[task_df$PptId==pptid]
    ppt_obsv <- task_df$TrialRingBoolean[task_df$PptId==pptid]
    ppt_modelests <- threshold_model(ppt_obsv, learning_method, threshold_method, param)
    csv$loglik[row] <- get_loglik(ppt_ests,ppt_modelests)
  }
  return(csv)
}

reduce_csv_to_bestest_fits <- function(csv){
  is_max <- rep(FALSE,nrow(csv))
  lms <- unique(csv$learning_method)
  tms <- unique(csv$threshold_method)
  for(learning_method in lms){
    lm_bool <- csv$learning_method==learning_method
    for(threshold_method in tms){
      lm_tm_bool <- lm_bool & csv$threshold_method==threshold_method
      ppt_ids <- unique(csv$PptId)
      for (pptid in ppt_ids){
        lm_tm_ppt_bool <- lm_tm_bool & csv$PptId==pptid
        logliks <- csv$loglik[lm_tm_ppt_bool]
        max_loglik <- max(logliks)
        max_loglik_row_index <- which(lm_tm_ppt_bool & csv$loglik == max_loglik)[[1]]
        is_max[max_loglik_row_index] <- TRUE
      }
    }
  }
  csv <- csv[is_max,]
  return(csv)
}

#or what if i want the model fits to start where the without-threshold best fits were
old_param <- function(pptid, csv_df, learning_method, threshold_method){
  rows <- csv_df$PptId==pptid & csv_df$learning_method==learning_method & csv_df$threshold_method==threshold_method
  param_as_df <- csv_df[rows,c('learning_rate','starting_est','window_size','threshold')]
  param_as_list <- as.numeric(param_as_df)
  names(param_as_list) <- colnames(param_as_df)
  param_as_list <- param_as_list[!is.na(param_as_list)]
  param_as_list
}

use_old_param_to_fit <- function(csv_df,
                                 lms = c("delta","movingwindow"),
                                 tms = c("drift","meandiff","kl"),
                                 n_times = 1,
                                 optim_method = "Nelder-Mead",
                                 csv_filename = NA){
  csv <- data.frame(csv_df[FALSE,]) #creates an empty data frame with the right column headings
  if(n_times > 1 | is.na(csv_filename)){
    for (learning_method in lms){
      for(threshold_method in tms){ #,"meandiff"
        fits <- find_best_fit(learning_method,threshold_method,times=n_times,optim_method=optim_method,ppt_unthresh_param)
        fits_as_csv <- turn_fits_to_csv(fits, learning_method, threshold_method)
        fits_as_csv$optim_method <- optim_method
        fits_as_csv$start_est_old_param <- TRUE
        fits_as_csv$startthresh0 <- NA
        fits_as_csv$num_fits <- n_times
        csv <- rbind(csv,fits_as_csv)
      }
    }
  } else {
    #learning_method <- lms[[1]]
    #threshold_method <- tms[[1]]
    for(learning_method in lms){
      for(threshold_method in tms){
        print(paste(learning_method,threshold_method))
        ppt_unthresh_param <- start_ppt_params(learning_method,"drift")
        ppt_unthresh_param$values <- sapply(ppt_ids,old_param,csv_df, learning_method, "none")
        
        if(file.exists(csv_filename)){
          oldcsv <- read.csv(csv_filename)
          already_fit_ids <- unique(oldcsv$PptId[oldcsv$learning_method==learning_method & oldcsv$threshold_method==threshold_method])
          ids_left_to_fit <- ppt_ids[!(ppt_ids %in% already_fit_ids)]
        } else {
          oldcsv <- NULL
          ids_left_to_fit <- ppt_ids
        }
        fits <- sapply(ids_left_to_fit, get_ppt_fit_vals, learning_method, threshold_method, ppt_unthresh_param, optim_method, csv_filename = csv_filename)
      }
    }
  }
}

my_AIC <- function(loglik,k){2*loglik + 2*k}

get_param_names <- function(learning_method){
  switch(learning_method,
         "delta" = c("starting_est","learning_rate","threshold"),
         "movingwindow" = c("window_size","threshold"))
}
add_AIC_to_df <- function(csv_df){
  num_param <- function(row){
    lm <- csv_df$learning_method[row]
    param_names <- get_param_names(lm)
    tm <- csv_df$threshold_method[row]
    if(tm=="none") param_names <- param_names[param_names!="threshold"]
    return(length(param_names))
  }
  num_p <- sapply(1:nrow(csv_df),num_param)
  csv_df$AIC <- my_AIC(csv_df$loglik,num_p)
  return(csv_df)
}

reduce <- function(csv_df,lm,tm){
  return(csv_df[csv_df$learning_method==lm & csv_df$threshold_method==tm,])
}

append_missing_ppts <- function(){
  csv_df <- read.csv("modelfitting/allfits.csv")
  csv_df <- csv_df[,colnames(csv_df)[!(colnames(csv_df) %in% c("X","X.1"))]]
  missing_ids <- ppt_ids[!(ppt_ids %in% csv_df$PptId)]
  
  #1. get starting Nelder-Mead estimates by fitting 60 times
  csv_df <- append_missing_ppts_Neldermead60(csv_df,missing_ids)
  write.csv(csv_df,"modelfitting/allfits.csv",row.names=FALSE)
  
  #2. now fit Nelder-Mead threshold models to best-fit no-threshod models
  use_old_param_to_fit(csv_df,
                       lms = c("delta","moving_window"),
                       tms = c("drift","meandiff","kl"),
                       n_times = 1,
                       optim_method = "Nelder-Mead",
                       csv_filename = "modelfitting/new.csv")
  
  neldermead_df <- read.csv("modelfitting/new.csv")
  meta_info <- data.frame(start_est_old_param = TRUE,
                          num_fits = 1,
                          startthresh0 = FALSE,
                          optim_method = "Nelder-Mead")
  neldermead_df <- addcol(neldermead_df,meta_info)
  csv_df <- rbind(csv_df,neldermead_df)
  write.csv(csv_df,"modelfitting/neldermead.csv")
  
  #3. GenSA - none fits
  csv_df <- reduce_csv_to_bestest_fits(csv_df)
  use_old_param_to_fit(csv_df,
                       lms = c("delta","movingwindow"),
                       tms = c("none"),
                       n_times = 1,
                       optim_method = "GenSA",
                       csv_filename = "modelfitting/gensa-copy.csv")
  
  #4. GenSA - delta drift fits
  use_old_param_to_fit(csv_df,
                       lms = c("delta","movingwindow"),
                       tms = c("none"),
                       n_times = 1,
                       optim_method = "GenSA",
                       csv_filename = "modelfitting/gensa-copy.csv")
  csv_df <- read.csv("modelfitting/gensa-copy.csv")
  use_old_param_to_fit(csv_df,
                       lms = c("delta"),
                       tms = c("drift"),
                       n_times = 1,
                       optim_method = "GenSA",
                       csv_filename = "modelfitting/gensa-copy.csv")
}

append_missing_ppts_Neldermead60 <- function(csv_df, missing_ids){
  #missing_ids <- as.character(c(369466,382882,413191,413704,416563,428839,430738,430759,431011,432931,433078))
  lms = c("delta","movingwindow")
  tms = c("none","drift","meandiff","kl")
  for(lm in lms){
    for(tm in tms){
      #count_up <- make_count_function()
      with(parent.frame(),count_up <- make_count_function())
      print(paste(lm,tm))
      optim_method <- "Nelder-Mead"
      learning_method <- lm
      threshold_method <- tm
      n_times <- 30
      best_fit <- find_best_fit(learning_method, threshold_method, n_times, optim_method, missing_ids)
      csv_df1 <- turn_fits_to_csv(best_fit, learning_method, threshold_method)
      meta_info <- data.frame(start_est_old_param = FALSE,
                              num_fits = n_times,
                              startthresh0 = FALSE,
                              optim_method = "Nelder-Mead")
      csv_df1 <- addcol(csv_df1,meta_info)
      print(colnames(csv_df)[!(colnames(csv_df) %in% colnames(csv_df1))])
      print(colnames(csv_df1)[!(colnames(csv_df1) %in% colnames(csv_df))])
      csv_df <- rbind(csv_df, csv_df1)
    }
  }
  return(csv_df)
}

learning_methods <- c("delta","movingwindow")
threshold_methods <- c("none","drift","meandiff","kl")
csv_df <- read.csv("modelfitting/allfits.csv")
csv_df <- csv_df[csv_df$PptId %in% ppt_ids,]
csv_df <- reduce_csv_to_bestest_fits(csv_df)
fits <- turn_csv_to_fits_obj(csv_df)

model_df <- task_df[,c('PptId','Condition','AccCorr','TrialNum','bernoulliEstimate','TrueBernoulliParam')]
#prod(unique(model_df$PptId) == colnames(delta_fits[['meandiff']]))

for (learning_method in c("delta","movingwindow")){
  switch(learning_method,
         "delta" = {lm_col <- "Delta"; lm_fits <- fits[['delta']]},
         "movingwindow" = {lm_col <- "Window"; lm_fits <- fits[['movingwindow']]})
  for (threshold_method in c("none","meandiff","drift","kl")){
    switch(threshold_method,
           "none" = tm_col <- "_N_",
           "meandiff" = tm_col <- "_T_M_",
           "drift" = tm_col <- "_T_D_",
           "kl" = tm_col <- "_T_K_")
    col_id <- paste0(lm_col, tm_col)
    model_df[,paste0(col_id,"ModelEsts")] <- unlist(lm_fits[[threshold_method]]['model_ests',])
    model_df[,paste0(col_id,"NumCP")] <- as.numeric(rep(lm_fits[[threshold_method]]['num_cp',],each=999))
    model_df[,paste0(col_id,"LogLik")] <- rep(as.numeric(lm_fits[[threshold_method]]['loglik',]),each=999)
    thresholds <- sapply(ppt_ids,function(pptid){lm_fits[[threshold_method]][['params',pptid]]['threshold']})
    model_df[,paste0(col_id,"Threshold")] <- rep(as.numeric(thresholds),each=999)
  }
}

#display model fits; participant responses versus best-fit threshold model

model_df_key <- function(learning_method, threshold_method){
  lm <- switch(learning_method,
               "delta" = "Delta",
               "movingwindow"="Window")
  tm <- switch(threshold_method,
               "none" = "N",
               "drift" = "T_D",
               "meandiff" = "T_M",
               "kl" = "T_K")
  paste(lm,tm,sep = "_")
}
