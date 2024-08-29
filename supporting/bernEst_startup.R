
source("misc_general_use_functions.R")
source("ppt_exclusions.R")
source("setup_adjsize_df.R")

importlibraries()
print("get df")
full_df <- read.csv("bernEstData.csv")
pract_df <- read.csv("bernEstPract.csv")

excluded_ppts <- get_excluded_ppts(full_df)
full_ids <- unique(full_df$PptId)
ppt_ids <- full_ids[!(full_ids %in% excluded_ppts)]
task_df <- full_df[full_df$PptId %in% ppt_ids,]

cond_dict <- sapply(full_ids,function(pptid){unique(full_df$Condition[full_df$PptId==pptid])})
names(cond_dict) <- full_ids

adj_df <- makeAdjDF(task_df)
adj_df_full <- makeAdjDF(full_df) 


