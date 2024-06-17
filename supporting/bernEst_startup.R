
source("misc_general_use_functions.R")
source("ppt_exclusions.R")
source("setup_adjsize_df.R")

importlibraries()
print("get df")
r_df <- read.csv("bernEstData.csv")

excluded_ppts <- get_excluded_ppts(r_df)
full_ids <- unique(r_df$PptId)
ppt_ids <- full_ids[!(full_ids %in% excluded_ppts)]
task_df <- r_df[r_df$PptId %in% ppt_ids,]

cond_dict <- sapply(full_ids,function(pptid){unique(r_df$Condition[r_df$PptId==pptid])})
names(cond_dict) <- full_ids

adj_df <- makeAdjDF(task_df)
adj_df_full <- makeAdjDF(r_df) 


