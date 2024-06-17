makeAdjDF <- function(f_df){
  
  ppt_ids <- unique(f_df$PptId)
  
  #quality assurance
  x <- unique(f_df[,c('PptId','BlockNum')])
  excluded_by_block <- nrow(x)/3 != length(ppt_ids)

  
  #to start, lets make a data frame that has a row for every possible adjustment size
  #and, for each participant, gives how many adjustments were made of each size
  
  #-----------------------------------------------------#
  #### AGGREGATE NUMBER OF ADJUSTMENTS: BY DIRECTION ####
  #-----------------------------------------------------#
  
    get_numadj_byadjsize <- function(df, ppt_rows, adjsize){
      #given the participant's id and the adjustment size returns 1 value: num of adjustments
      sum(df$AdjSize[ppt_rows]==adjsize,na.rm = TRUE)
    }
    adj_df <- sapply_to_df(ppt_ids, function(pptid){
      ppt_rows <- f_df$PptId==pptid
      sapply(-100:100,function(adjsize){
        get_numadj_byadjsize(f_df, ppt_rows, adjsize)
      })
    },"PptId","NumAdjs")
    
    #-----------------------------------------------------#
    ####             ADD SOME MORE COLUMNS             ####
    #-----------------------------------------------------#
    
    adj_df$AdjSize <- -100:100
    adj_df <- merge(adj_df, unique(f_df[,c('PptId','Condition','AccCorr')]))
    adj_df$CorrectDirection <- adj_df$AdjSize > 0
    adj_df$CorrectDirection[adj_df$AdjSize==0] <- NA
    adj_df$AdjSize_Abs <- abs(adj_df$AdjSize)
    
    #-----------------------------------------------------#
    ####                SANITY CHECK                   ####
    #-----------------------------------------------------#
    
    num_trials <- length(unique(f_df$TrialNum))
    ppts_sum_num_adj <- aggregate(adj_df$NumAdjs,by=list(adj_df$PptId),sum)$x
    
    check1 <- prod(ppts_sum_num_adj==num_trials | ppts_sum_num_adj == (num_trials - 1))
    check2 <- prod(adj_df$AdjSize == rep(-100:100,length(unique(adj_df$PptId))))
    if (!check1){
      stop(paste("the number of adjustments (for at least one participant)",
                 "made is less than the total number of trials. Since a zero-size",
                 "adjustment (no movement) is still counted as an adjustment,",
                 "something has gone wrong since not all trials are accounted for"))
    }
    if(!check2){
      stop(paste("for some reason, your adjustment size column is out of order. it",
                 "is not simply the list -100:100 repeated as many times as there are",
                 "participants"))
    }
    
    #---------------------------------------------------------#
    #### AGGREGATE NUMBER OF ADJUSTMENTS: BY LEFT VS RIGHT ####
    #---------------------------------------------------------#
    
    agg_adjs <- sapply(ppt_ids,function(pptid){
      temp_df <- f_df[f_df$PptId==pptid,]
      sapply(-100:100,function(adjsize){
        sum(temp_df$AdjSize_LeftRight==adjsize,na.rm = TRUE)
      })
    })
    temp_df <- as.data.frame(agg_adjs)
    colnames(temp_df) <- ppt_ids
    mdf <- melt(temp_df)
    colnames(mdf) <- c("PptId","NumAdjs_LeftRight")
    mdf$AdjSize <- -100:100
    
    mdf <- mdf[order(mdf$PptId, mdf$AdjSize),]
    adj_df <- adj_df[order(adj_df$PptId,adj_df$AdjSize),]
    
    check1 <- prod(mdf$PptId==adj_df$PptId)
    check2 <- prod(mdf$AdjSize == rep(-100:100,length(unique(mdf$PptId))))
    if(!check1 | !check2){
      #stop("Something went wrong at checkpoint 2")
    }
    
    adj_df$NumAdjs_LeftRight <- mdf$NumAdjs_LeftRight
    
    return(adj_df)
  
}
