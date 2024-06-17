curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
subdir <- paste0(curdir, "/supporting/")
setwd(subdir)
source("bernEst_startup.R")

#------------------------------------------------------------#
# DO THE TWO CONDITIONS MAKE DIFFERNT NUMBERS OF ADJUSTMENTS #
#------------------------------------------------------------#

## ---- numadj-groupcompare

# 1. T-TEST: do the two groups make different number of adjustments
# do Manual participants make more adjustments?

#a_df <- adj_df_full[adj_df_full$AdjSize!=0,]
a_df <- adj_df[adj_df$AdjSize!=0,]

adf <- aggregate(a_df$NumAdjs,by=list(Condition=a_df$Condition,PptId=a_df$PptId),sum)

p <- ggplot(adf,aes(x=Condition,y=x,fill=Condition,color=Condition)) +
  geom_violin(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE,alpha=0.5)+
  geom_point()+
  ylab("Number of Adjustments")
p <- gg_bernEst(p,scales = c(x="discrete",y="continuous"),margin = c(x=.25,y=0))

t_test <- t.test(adf$x[adf$Condition=="Manual"],adf$x[adf$Condition=="Automatic"],var.equal = T)
t_test <- APA_format_vals(t_test)

#------------------------------------------------------------#
#   HOW/WHERE DO MANUAL PARTICIPANTS MAKE MORE ADJUSTMENTS   #
#------------------------------------------------------------#

# 2. NOISE hypothesis: Higher proportion negative adjustments
num_adj_by_sign <- function(pptid){
  out <- sapply(c(-1,1),function(c_sign){sum(task_df$AdjSign[task_df$PptId==pptid]==c_sign,na.rm = T)})
  out <- as.table(out)
  row.names(out) <- c("-1","+1")
  out
}
manual_ids <- unique(task_df$PptId[task_df$Condition=="Manual"])
auto_ids <- unique(task_df$PptId[task_df$Condition=="Automatic"])
manual_avg_adj <- apply(sapply(manual_ids,num_adj_by_sign),1,mean)
automatic_avg_adj <- apply(sapply(auto_ids,num_adj_by_sign),1,mean)


# T TESTS
#more negative adjustments?
ttest_more_neg <- t.test(sapply(manual_ids,num_adj_by_sign)["-1",],sapply(auto_ids,num_adj_by_sign)["-1",])
#higher portion negative adjustments?
portion_neg <- function(ids){
  sapply(ids,num_adj_by_sign)["-1",]/(sapply(ids,num_adj_by_sign)["+1",])
}
ttest_higher_proportion_neg <- t.test(portion_neg(manual_ids),portion_neg(auto_ids))
ttest_prop <- APA_format_vals(ttest_higher_proportion_neg)

#3. ADDED ADJUSTMENTS hypothesis: 
#are added adjustments positively skewed?

M <- as.table(rbind(manual_avg_adj, automatic_avg_adj))
dimnames(M) <- list(sign = c("Manual", "Automatic"),
                    condition = c("-1","+1"))

pos_bias_x1 <- sapply(manual_ids,num_adj_by_sign)["-1",] - M["Automatic","-1"]
pos_bias_x2 <- sapply(manual_ids,num_adj_by_sign)["+1",] - M["Automatic","+1"]
ttest_remaining_positive_bias <- t.test(pos_bias_x1,pos_bias_x2)
ttest_pos_bias <- APA_format_vals(ttest_remaining_positive_bias)
