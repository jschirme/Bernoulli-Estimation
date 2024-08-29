## ---- plot-unusual-ppts

plot_ppts <- function(ppts_to_display, 
                      blocknum=c(1),
                      n_col=2,
                      n_row=2,
                      n_trials=333){
  gg_df <- full_df[full_df$PptId %in% ppts_to_display,c("PptId","Condition","AdjSize","TrueBernoulliParam","bernoulliEstimate","BlockNum","TrialNum")]
  #1. set-up adjsign for the colouring of the participant response line (green, red, or neutral)
  gg_df$AdjSign <- NA
  gg_df$AdjSign <- shifter(sign(gg_df$AdjSize),1)
  
  #2. order participants 
  gg_df$PptId <- factor(gg_df$PptId,levels = ppts_to_display,ordered=T)
  
  #3. reduce the df to only have the blocks i want to display
  gg_df <- gg_df[gg_df$BlockNum %in% blocknum & (gg_df$TrialNum-1)%%333 <= n_trials,]
  
  make_labeler <- function(colname,replace_dict){
    function(labels, multi_line = TRUE, sep = ": ")
    {
      labels[[colname]] <- as.character(labels[[colname]])
      labels[[colname]] <- listreplace(labels[[colname]],names(replace_dict),replace_dict)
      list(unname(unlist(labels[[colname]])))
    }
  }
  f <- make_labeler('PptId',cond_dict)

    
  p <- ggplot(gg_df, aes(x=TrialNum))+
    geom_line(aes(y=TrueBernoulliParam),color="black",size=0.5,linetype="dashed")+
    geom_line(aes(y=bernoulliEstimate,color=AdjSign),size=1.0,show.legend = F)+
    facet_wrap(~PptId,scales = 'free',
                        labeller = f,
                        ncol=n_col,nrow=n_row)+
    coord_cartesian(ylim = c(0,1))+
    themeAPA()+
    labs(x = "Trial Number", y = "Bernoulli Parameter")+
    scale_colour_gradient(low = "#ff0000", high = "#00FF28")
  p <- make_APA(p)
  p <- p + theme(axis.title.x = element_text(size = 15),
                 axis.title.y = element_text(size = 15),
                 strip.text.x = element_text(size = 15))
  p
}

#------------------------------------------------------------#
#                   MANUSCRIPT PLOT                          #
#------------------------------------------------------------#

#ppts_to_display <- c(c("361756","347989"),c("410113","430963"))
#plot_ppts(ppts_to_display) + theme(aspect.ratio=1.5/3)

#------------------------------------------------------------#
#                  SUPPLEMENTARY PLOT                        #
#------------------------------------------------------------#

## ---- supplementary-plot-ppts

plot_all_ppts <- function(ppts_to_display, ppt_labels,
                          blocknum=c(1,3),
                          n_col=2,
                          n_row=4,
                          n_trials=333){
  gg_df <- full_df[full_df$PptId %in% ppts_to_display,c("PptId","Condition","AdjSize","TrueBernoulliParam","bernoulliEstimate","BlockNum","TrialNum")]
  #1. set-up adjsign for the colouring of the participant response line (green, red, or neutral)
  gg_df$AdjSign <- NA
  gg_df$AdjSign <- shifter(sign(gg_df$AdjSize),1)
  
  #2. order participants 
  gg_df$PptId <- factor(gg_df$PptId,levels = ppts_to_display,ordered=T)
  
  #3. reduce the df to only have the blocks i want to display
  gg_df <- gg_df[gg_df$BlockNum %in% blocknum & (gg_df$TrialNum-1)%%333 <= n_trials,]
  
  make_labeler <- function(colname,replace_dict){
    function(labels, multi_line = TRUE, sep = ": ")
    {
      labels[[colname]] <- as.character(labels[[colname]])
      labels[[colname]] <- listreplace(labels[[colname]],names(replace_dict),replace_dict)
      list(unname(unlist(labels[[colname]])))
    }
  }
  f <- make_labeler('PptId',ppt_labels) #cond_dict
  
  num_pages <- ceiling(length(ppts_to_display)/n_row)
  num_ppt_per_page <- n_row
  for (pagenum in 1:num_pages){
    start_ppt_index <- ((pagenum-1)*num_ppt_per_page) + 1
    pptids <- ppts_to_display[start_ppt_index:(start_ppt_index+num_ppt_per_page-1)]
    p <- ggplot(gg_df, aes(x=TrialNum))+
      geom_line(aes(y=TrueBernoulliParam),color="black",size=0.5,linetype="dashed")+
      geom_line(aes(y=bernoulliEstimate,color=AdjSign),size=1.0,show.legend = F)+
      facet_wrap_paginate(PptId~BlockNum,scales = 'free',
                          labeller = f,
                          ncol=n_col,nrow=n_row,page=pagenum)+
      coord_cartesian(ylim = c(0,1))+
      themeAPA()+
      labs(x = "Trial Number", y = "Bernoulli Parameter")+
      scale_colour_gradient(low = "#ff0000", high = "#00FF28")
    p <- make_APA(p)
    p <- p + theme(axis.title.x = element_text(size = 15),
                   axis.title.y = element_text(size = 15),
                   strip.text.x = element_text(size = 15))
    #print(p)
    ggsave(filename = paste0(pagenum,"PPT_data_visual.pdf"),plot = p,width=21,height=20,unit="cm")
  }
}

full_df$RespMain <- full_df$AdjSize==0
a_df <- aggregate_rf(full_df,"RespMain",c("PptId","Condition"),function(x){sum(x,na.rm=T)})
a_df <- a_df[order(a_df$RespMain),]
ids <- a_df$PptId
plot_ids <- ids[ids %in% ppt_ids]
plot_labels <- sapply(plot_ids,function(pptid){
  num_adj <- 998 - a_df$RespMain[a_df$PptId==pptid]
  cond <- a_df$Condition[a_df$PptId==pptid]
  paste0("Condition: ",cond,"\n Number Adjustments:", num_adj)
})
plot_all_ppts(plot_ids, plot_labels)


