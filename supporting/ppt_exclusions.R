#----------------------------------------------------------------------------------#
#                                   EXCLUSIONS                                     #
#----------------------------------------------------------------------------------#

ppt_exclusions <- function(f_df){
  flag_df <- data.frame(PptId = unique(f_df$PptId))
  flag_df$BlockedAccCorr <- flag_df$PptId %in% exclude_on_acc(f_df)
  flag_df$TimeOnTask <- flag_df$PptId %in% exclude_on_taskdur(f_df)
  #flag_df$SurveyResponses <- flag_df$PptId %in% exclude_on_surveyQ()
  return(flag_df)
}

collapse_ppt_exclusions <- function(flag_df){
  temp_df <- flag_df[,-1]
  temp_df[temp_df!=FALSE]<-TRUE
  temp_df <- apply(temp_df,2,as.logical)
  
  is_excluded <- apply(temp_df,1,FUN=sum) >= 1
  excluded_ppts <- flag_df$PptId[is_excluded]
  return(excluded_ppts)
}

get_excluded_ppts <- function(f_df){
  flag_df <- ppt_exclusions(f_df)
  if(ncol(flag_df) == 1){
    excluded_ppts <- flag_df
  } else {
    excluded_ppts <- collapse_ppt_exclusions(flag_df)
  }
  return(excluded_ppts)
}

#--------------------------------------------------------------------------#
#                           EXCLUSIONS                                     #
#--------------------------------------------------------------------------#

          #----------------------------------------------#
          #               TASK PERFORMANCE               #
          #----------------------------------------------#
# BECAUSE YOU CAN SEE PERFORMANCE BIFURICATE INTO FAIR/GOOD AND BAD/AT-CHANCE BY THE 
# FINAL BLOCK; THIS CODE EXCLUDES ON THE CRITERIA OF CONSISTENT ACCURACY OVER TASK 

#given a ppt id, returns their accuracy (correlation to true bernoulli parameter)
#separated by the three blocks

## ---- ppt-exclusions

# this section is looking at how accuracy changes over the course of the task
make_accorr_plot <- function(f_df){
  p <- ggplot(gg_df)+
    geom_histogram(bins=40,fill="grey75",color="black")+
    facet_wrap(~Block,scale='free')+
    labs(
      x = "Accuracy (Beta coefficient predicting true bernoulli with bernoulli estimate)",
      y = "Number of participants"
    )
  return(p)
}


## ----
exclude_on_acc <- function(f_df, critical_r = critical.r(333)){
  blocknames <- c("ACC_Block1","ACC_Block2","ACC_Block3")
  acc_df <- unique(f_df[,c("PptId",blocknames)])
  exclude <- apply(acc_df[,blocknames],1,function(accs){sum(accs<critical_r,rm.na=T)>1})
  flag_ppts <- acc_df$PptId[exclude | is.na(exclude)]
  print(length(flag_ppts))
  return(flag_ppts)
}
critical.r <- function( n, alpha = .05 ) {
  df <- n - 2
  critical.t <- qt(alpha/2, df, lower.tail = F)
  critical.r <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
  return(critical.r)
}

  
          #----------------------------------------------#
          #               TIME SPENT ON TASK             #
          #----------------------------------------------#

#------------- FULL TASK DURATION -------------#

exclude_on_taskdur <- function(f_df = task_df,
                               full_ids = unique(f_df$PptId)){
  #functions:
  get_time_on_task <- function(pptid){
    starttime <- f_df$TimeElapsed[f_df$PptId==pptid & f_df$TrialNum==1]
    stoptime <- f_df$TimeElapsed[f_df$PptId==pptid & f_df$TrialNum==999]
    time_in_ms <- stoptime-starttime
    time_in_min <- ms_in_min(time_in_ms)
    return(time_in_min)
  }
  get_max_trial_length <- function(pptid){
    trial_lengths_ms <- f_df$rt[f_df$PptId==pptid]
    trial_lengths_min <- ms_in_min(trial_lengths_ms)
    return(max(trial_lengths_min))
  }
  
  #1. get participant task durations
  tasktime_list <- sapply(full_ids, get_time_on_task)
  maxtrial_list <- sapply(full_ids, get_max_trial_length)
  tasktime_df <- data.frame(PptId=full_ids, 
                            TimeOnTask_min=tasktime_list, 
                            MaxTrialLength_min=maxtrial_list)
  
  #2. what is an unusually long time to take on the task
  #here defined: +/- 1.5 IQR from above 75% quartile
  outlier_timeontask <- outlierCriticalValues(tasktime_list)
  #91.70802 min
  
  #3. flag outlier participants
  flag_ppts <- tasktime_df$PptId[tasktime_df$TimeOnTask_min>91.7]
  tasktime_df[tasktime_df$PptId %in% flag_ppts,]
  
  return(flag_ppts)
}

