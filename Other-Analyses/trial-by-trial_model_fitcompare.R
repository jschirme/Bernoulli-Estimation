curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
subdir <- paste0(curdir, "/supporting/")
setwd(subdir)
source("bernEst_startup.R")
setwd(curdir)
source("trial-by-trial_model_GenSA.R")
csv_df <- add_AIC_to_df(csv_df)

#---------------------------------------------------------#
#     WHICH LEARNING METHOD (LM) PROVIDES BETTER FITS:    #
#---------------------------------------------------------#
ppt_modelfit_goodness <- function(pptid, 
                            lms = c("delta","movingwindow"),
                            tms = c("none","drift","meandiff","kl"),
                            goodness_var="loglik"){
  rows <- csv_df$PptId==pptid
  if (length(lms)==1) {
    lm <- lms[[1]]
    rows <- rows & csv_df$learning_method==lm
    sapply(tms,function(tm){
      rows <- rows & csv_df$threshold_method==tm
      csv_df[[goodness_var]][rows]
    })
  } else {
    sapply(lms,function(lm){
      rows <- rows & csv_df$learning_method==lm
      sapply(tms,function(tm){
        rows <- rows & csv_df$threshold_method==tm
        csv_df[[goodness_var]][rows]
      })
    })
  }
  
}
compare_ppt_lm_fit <- function(pptid, ...){
  apply(X=ppt_modelfit_goodness(pptid,...),MARGIN=1,FUN=diff)
}

#ALMOST UNIVERSALLY, DELTA PROVIDES A BETTER FIT THAN MOVING WINDOW, EVEN USING AIC
d_to_mw <- sapply(ppt_ids,compare_ppt_lm_fit,goodness_var="AIC")
d_to_mw <- as.data.frame(t(d_to_mw))
gg_df <- melt(d_to_mw)
ggplot(gg_df,aes(x=value))+
  geom_histogram(bins=50)+
  facet_wrap(~variable)+
  themeAPA()


#----------------------------------------------------------#
#     WHICH THRESHOLD METHOD (TM) PROVIDES BETTER FITS:    #
#----------------------------------------------------------#
none_to_thresh_improvement <- function(pptid, 
                                       lms = c("delta","movingwindow"),
                                       tms = c("drift","meandiff","kl"),
                                       goodness_var="loglik"){
  fit_goodness <- ppt_modelfit_goodness(pptid,goodness_var=goodness_var)
  sapply(lms,function(lm){
    sapply(tms,function(tm){
      fit_goodness[tm,lm] - fit_goodness['none',lm]
    })
  })
}
#if we try to plot each participant's loglik/AIC to compare the fit of tms
#we see that individual variation far exceeds differences in model fits
# delta_fit_goodness <- sapply(ppt_ids,ppt_modelfit_goodness,"delta")
# gg_df <- melt(delta_fit_goodness)
# ggplot(gg_df,aes(x=Var1,y=value))+
#   geom_boxplot()+
#   themeAPA()

#so lets normalize the model fits by subtracting out the model-fit for a no-threshold model
#this is like a repeated-measures correction where a participant value is estimated and adjusted for
#another way of thinking about this is we are looking at relative improvements to the no-threshold model
delta_improvements <- sapply(ppt_ids,none_to_thresh_improvement,lms=c("delta"))
row.names(delta_improvements) <- c("drift","meandiff","kl")
gg_df <- melt(delta_improvements)
ggplot(gg_df,aes(x=Var1,y=value))+
  geom_boxplot()+
  themeAPA()

  #wow. its certainly not a clear cut case which threshold method is best
#for each type of threshold mechanism, we get some outliers with better improvements to AIC
#lets try to find out if these are the same or different participants
lm_improvements <- sapply(ppt_ids,none_to_thresh_improvement)
row.names(lm_improvements) <- paste(rep(c("delta","movingwindow"),each=3),c("drift","meandiff","kl"),sep = ".")
improves_with_thresh <- function(pptid,
                                 critval = 50){
  names(which(lm_improvements[,pptid]>critval))
}
thresh_ppts <- sapply(ppt_ids,improves_with_thresh,critval=50)
thresh_ppts <- thresh_ppts[sapply(thresh_ppts,function(improvements){length(improvements)>0})]
thresh_ppt_ids <- names(thresh_ppts)
thresh_ppts 
c('331711','347989','348394','372892','384880','385213','386290','398974','401923','407335','418369','427924','428860','42940')
#it seems many of those participants who got better results with thresholds added
#did so with multiple threshold methods (sometimes various learning methods)


#------------------------------------------------------------#
#     INVESTIGATING WHY SOME THRESHOLD-MODELS FIT BETTER:    #
#------------------------------------------------------------#

ppt_best_model <- function(pptid, 
                           lms = learning_methods,
                           tms = threshold_methods,
                           goodness_measure = "AIC"){
  logliks <- sapply(lms, ppt_modelfit_goodness, pptid = pptid, goodness_var=goodness_measure)
  logliks <- logliks[row.names(logliks) %in% tms,]
  if(is.null(dim(logliks))){ #code above messes up array formatting if lms or tms is only length 1
    logliks <- t(array(logliks, dim = c(length(lms),length(tms)), dimnames = list(lms,tms)))
  }
  #which(is.max(logliks),arr.ind = T)
  max <- which_names(is.max(logliks)) 
  return(c(learning_method=max[2],threshold_method=max[1]))
}

plot_ppt_best_fit <- function(ppts_to_display,
                              num_at_once=4){
  #count <- make_count_function(length(ppts_to_display))
  num_pages <- ceiling(length(ppts_to_display)/num_at_once)
  best_models <- sapply(ppts_to_display,ppt_best_model)
  plot_ppt_fit(ppts_to_display, best_models)
}

make_labeler <- function(colname,replace_dict){
  function(labels, multi_line = TRUE, sep = ": ")
  {
    labels[[colname]] <- as.character(labels[[colname]])
    labels[[colname]] <- listreplace(labels[[colname]],names(replace_dict),replace_dict)
    list(unname(unlist(labels[[colname]])))
  }
}

plot_ppt_fit <- function(ppts_to_display, model_type,
                         num_at_once=4, block="1", 
                         order_by = NULL,
                         labeller = label_both){
  #this section is just the contigency plans for being passed different types of objects into "model_type"
  #i know it is going to have a threshold_method and a learning_method way to index the data
  #but i dont know whether it will be of length one (and thus uniform across participants)
  #or an array (matrix) where i index both by ppt id AND model type
  #what i need that info for, is to know which columsn of model_df to get the model estimates from
  if(class(model_type)[[1]]=="list"){model_type<-list_to_array(model_type)}
  else {model_type <- as.array(model_type)}
  indv_model_assignment <- length(dim(model_type)) == 2
  if(indv_model_assignment){
    model_ests <- sapply(ppts_to_display,function(pptid){
      #find the corresponfing column of model df to get the model estimates from
      colname <- paste(model_df_key(model_type['learning_method',pptid],model_type['threshold_method',pptid]),"ModelEsts",sep = "_")
      model_df[model_df$PptId==pptid,colname]
    })
  } else {
    model_ests <- sapply(ppts_to_display,function(pptid){
      learning_method <- model_type['learning_method']
      threshold_method <- model_type['threshold_method']
      colname <- paste(model_df_key(learning_method,threshold_method),"ModelEsts",sep = "_")
      model_df[model_df$PptId==pptid,colname]
    })
  }
  #now that i have the model estimates, i need to make the df that ggplot2 will use
  gg_df <- task_df[task_df$PptId %in% ppts_to_display,c("PptId","Condition","TrialNum","bernoulliEstimate","TrueBernoulliParam","BlockNum","AdjSign")]
  #make sure its in the order of ppts_to_display
  gg_df$PptId <- factor(gg_df$PptId, levels=ppts_to_display)
  gg_df <- gg_df[order(gg_df$PptId),]
  
  #in order to colour the trial-to-trial movement by correct vs incorrect direction
  #i need the adjsign to be shifted up one
  gg_df$AdjSign <- shifter(gg_df$AdjSign,1)
  #to get the one column of participant estimates and the one column of model estimates
  #i am making a loong data-frame, that has a column that assigns whetehr the model or 
  #the participant made the estimates
  #then i can facet by this variable
  gg_df <- rbind(gg_df,gg_df)
  gg_df$estimater <- rep(c('model','ppt'),each=nrow(gg_df)/2)
  gg_df$bernoulliEstimate[1:(nrow(gg_df)/2)] <- as.numeric(model_ests)
  #for only the model estimates, i am setting the adjsign (whether the movement was 
  #right or wrong to NA - which will make the line grey)
  gg_df$AdjSign[1:(nrow(gg_df)/2)] <- NA
  #let me add AIC values for these fits
  model_AICs <- ppt_model_AIC(model_type)
  nothresh_model <- sapply(ppts_to_display, ppt_best_model, tms = "none")
  nothresh_AIC <- ppt_model_AIC(nothresh_model)
  gg_df$AIC <- rep(model_AICs[ppts_to_display],each=999)
  gg_df$AIC_Improvement <- gg_df$AIC - rep(nothresh_AIC[ppts_to_display],each=999)
  
  #now i order the participants if so instructed
  if(!is.null(order_by)){
    gg_df$PptId <- factor(gg_df$PptId,levels = unique(gg_df$PptId[order(gg_df[,order_by])]))
  }
  num_pages <- ceiling(length(ppts_to_display)/num_at_once)
  for(page_num in 1:num_pages){
    p <- ggplot(gg_df[gg_df$PptId %in% ppts_to_display & grepl(as.character(block),gg_df$BlockNum),], aes(x=TrialNum))+
      geom_line(aes(y=TrueBernoulliParam),color="black",size=0.5,linetype="dashed")+
      #geom_line(aes(y=modelEstimate),color="lightblue",size=1)+
      geom_line(aes(y=bernoulliEstimate,color=AdjSign),size=1)+
      facet_wrap_paginate(PptId~estimater,scales = 'free',
                          labeller = labeller,
                          ncol=2,nrow=num_at_once,page=page_num)+
      coord_cartesian(ylim = c(0,1))+
      themeAPA()+
      labs(x = "Trial Number", y = "Bernoulli Parameter")+
      #theme(strip.text.x = element_blank())+
      scale_colour_gradient(low = "#ff0000", high = "#00FF28")
    print(make_APA(p,margin = c(x=0,y=0.05)))
    #img_name <- paste("images/",order_by,"_block",block,"_page",page_num,".png")
    #ggsave(img_name,p)
  }
}

ppt_model_AIC <- function(model_type){
  #contingency plan if passed a vector for all participants 
  #make it into a matrix with identical vales for all participants
  model_type <- as.array(model_type)
  if(length(dim(model_type))==1){
    model_type <- sapply(ppt_ids,function(pptid){model_type})
  }
  #now functions
  get_k_for_AIC <- function(pptid){
    lm <- model_type['learning_method',pptid]
    tm <- model_type['threshold_method',pptid]
    k <- 1
    if(lm=="delta")k <- k+1
    if(tm!="none")k<-k+1
    return(k)
  }
  get_loglik <- function(pptid){
    lm <- model_type['learning_method',pptid]
    tm <- model_type['threshold_method',pptid]
    csv_df$loglik[csv_df$PptId==pptid & csv_df$learning_method==lm & csv_df$threshold_method==tm]
  }
  #now procedure
  k_per_ppt_model <- sapply(colnames(model_type),get_k_for_AIC)
  loglik_per_ppt_model <- sapply(colnames(model_type),get_loglik)
  mapply(my_AIC,loglik_per_ppt_model,k_per_ppt_model)
}

gg_df <- data.frame(PptId = ppt_ids,Condition=cond_dict[ppt_ids])
no_thresh <- sapply(ppt_ids,ppt_best_model,tms="none")
gg_df$NT_lm <- no_thresh['learning_method',]
gg_df$NT_tm <- no_thresh['threshold_method',]
thresh <- sapply(ppt_ids,ppt_best_model,tms=c("drift","meandiff","kl"))
gg_df$T_lm <- thresh['learning_method',]
gg_df$T_tm <- thresh['threshold_method',]
gg_df$NT_AIC <- ppt_model_AIC(no_thresh)
gg_df$T_AIC <- ppt_model_AIC(thresh)
gg_df$AIC_diff <- gg_df$T_AIC - gg_df$NT_AIC
gg_df$AIC_diff[gg_df$AIC_diff<0] <- 0

#lets see if the groups are different in how much a threshold improved the
#fit of their models
p <- ggplot(gg_df,aes(x=AIC_diff,fill=Condition))+
  geom_histogram(color="black")+
  facet_wrap(~Condition)
gg_bernEst(p)

men_AICdiff <- gg_df$AIC_diff[gg_df$Condition=="Manual"]
auto_AICdiff <- gg_df$AIC_diff[gg_df$Condition=="Automatic"]
ks.test(men_AICdiff, auto_AICdiff)

#now much much can we trust the AIC improvement to tell us whether the ppt showed
#thresholding behavior
#lets sort participants by how much adding the threshold improved their model
gg_df$PptId <- factor(gg_df$PptId,levels=gg_df$PptId[order(gg_df$AIC_diff)])
# now lets plot every participants first block best threshold fit model
# sorted by how much the threshold improved their model

#make a labeller
cond_and_AICdiff_dict <- sapply(ppt_ids,function(pptid){
  paste(cond_dict[pptid],round(gg_df$AIC_diff[gg_df$PptId==pptid]))
})
labeller <- make_labeler("PptId",cond_and_AICdiff_dict)
plot_ppt_fit(sample(ppt_ids,4),thresh,order_by = "AIC_Improvement",labeller = labeller)
