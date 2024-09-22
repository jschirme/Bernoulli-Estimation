#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------- GENERAL --------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
{
  #DESCRIPTIONS FOR FUNCTIONS
  DESCR <- function(function_name = DESCR){
    function_name <- deparse(substitute(function_name))
    typeof(function_name)
    switch(function_name,
           DESCR = description <- "Tells you about each function in this library,<br>
            to use, this function will ignore all the whitespace in the string you pass<br>
            description that is more than 1 space. <br>If you want to add a new line, you<br>
            have to use the newline keyword for html 'br' ",
           themeAPA = description <- "returns a theme that can be passed to the 
            ggplot2 parameter theme.",
           outlierBool = description <- "Accepts a vector of numeric valies. Returns 
            a list of booleans, the length of the input, determining whether each 
            value is an outlier or not (+/- 1.5 IQR)",
           getData = description <- "If not passed a directory, assumes the data is in 
            a subdirectory (called \"Data\") from the location 
            that the current script is being run. <br><br> If not 
            passed a list of filenames, it will use all the csv
            files in taskdir"
    )
    description <- gsub("[\n]","",description)
    description <- gsub("  ","",description)
    description <- gsub("<br>","\n",description)
    writeLines(paste0("\n",description))
  }
  
  #GENERALLY NICE TO HAVE LIBRARIES AT DISPOSAL
  importlibraries <- function(){
    library(reshape2)
    library(ggplot2)
    library(rebus)
    library(stringr)
    library(plyr)
    #library(lme4)
    #library(lmerTest)
    library(ggforce)
    library(ggthemes)
  }
  
  #FORMATING
  themeAPA <- function(){
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
           panel.background = element_blank(), plot.margin=unit(c(0,0,0,0),"null"),
           strip.background = element_rect(fill="white"),
           #text=element_text(size=16,  family="TT Arial")
          )
  }
  
  make_APA <- function(gg_plot, 
                       scales = c(x="continuous",y="continuous"),
                       margin = c(x=0,y=0)){
    gg_plot <- gg_plot + themeAPA()
    #if(rm.margin){
      gg_plot <- switch(scales[['x']],
                        "continuous" =  gg_plot + scale_x_continuous(expand=c(margin[['x']],margin[['x']])),
                        "discrete" =  gg_plot + scale_x_discrete(expand=c(margin[['x']],margin[['x']])))
      gg_plot <- switch(scales[['y']],
                        "continuous" =  gg_plot + scale_y_continuous(expand=c(margin[['y']],margin[['y']])),
                        "discrete" =  gg_plot + scale_y_discrete(expand=c(margin[['y']],margin[['y']])))
    #}
    return(gg_plot)
  }
  
  gg_fixfacet <- function(gg_plot,
                          lims = get_xylims(gg_plot)){
    get_xylims <- function(gg_plot){
      x_col <- as.character(gg_plot$mapping$x)[[2]]
      y_col <- as.character(gg_plot$mapping$y)[[2]]
      x <- gg_plot$data[,x_col]
      y <- gg_plot$data[,y_col]
      xlims <- c(min(x,na.rm=T),max(x,na.rm=T))
      ylims <- c(min(y,na.rm=T),max(y,na.rm=T))
      return(list(x=xlims,y=ylims))
    }
    gg_plot <- gg_plot + theme_tufte() + theme(axis.line=element_line())+
      #scale_x_continuous(limits=lims[['x']], expand=c(0,0))+ 
      #scale_y_continuous(limits=lims[['y']], expand=c(0,0))
      coord_cartesian(xlim = lims[['x']], ylim = lims[['y']])
    gg_plot
  }
  APA_result_string <- function(stat_test, df_value,f_value, p_value){
    p_report <- paste0("p = ",round(p_value,3))
    if (p_value > 0.001) {p_report <- "p < 0.001"}
    paste0(stat_test, "(", round(df_value,2), ") = ", round(f_value,2), ", ",
           p_report)
  }
  APA_result_ttest <- function(model){
    f <- model$statistic
    p <- model$p.value
    df <- model$parameter
    APA_result_string('t',df,f,p)
  }
  APA_results.summary.lm <- function(summary_obj){
    f <- summary_obj$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    paste0("F(",f['numdf'],", ", f['dendf'], ") = ", round(f['value'],2), ", ",
           "p = ",round(p,3))
  }
  APA_results.anova <- function(anova_obj){
    df_value <- anova_obj$Df
    f_value <- anova_obj$F
  }
  APA_results.htest <- function(x){
    p_value <- x$p.value
    if(p_value < 0.001){p_value = " > 0.001"}
    else{p_value = paste0(" > ", round(p_value,3))}
    f_value <- x$statistic
    df_value <- x$parameter
    APA_result_string("t",df_value,f_value,p_value)
  }
  APA_results.chisq <- function(model){
    APA_result_string('x',model$parameter,model$statistic,model$p.value)
  }
  
  APA_format_vals <- function(model){
    model_class <- class(model)[[1]]
    switch(model_class,
           "htest" = {vals <- c(f = unname(round(model$statistic,2)),
                                p = unname(round(model$p.value,3)),
                                df = unname(round(model$parameter,2)))},
           "anova" = {vals <- list(f = model$F,
                                   p = model$`Pr(>F)`,
                                   df = model$Df,
                                   res.df = model$Res.Df)
                      vals <- sapply(vals,function(val){val})},
           "lm" = {sum_obj <- summary(model)
                   f_obj <- sum_obj$fstatistic
                   p_val <- pf(f_obj[1],f_obj[2],f_obj[3],lower.tail=F)
                   vals <- c(f =  unname(round(f_obj[['value']],3)),
                                p =  unname(round(p_val,3)),
                                df = f_obj[['dendf']],
                                r_sq = unname(round(sum_obj$r.squared,3)))}
    )
    return(vals)
  }
  
  ANOVA_APA <- function(aov_output){
    aov_df <- summary(aov_output)[[1]]
    rownames <- unname(sapply(row.names(aov_df),function(rname){unlist(strsplit(rname," "))[1]}))
    #row_index <- which(rowname==rownames)
    out <- data.frame(F = aov_df$`F value`,
                      Df1 = aov_df$`Df`,
                      Df2 = aov_df$`Df`[nrow(aov_df)],
                      p = aov_df$`Pr(>F)`)
    row.names(out) <- rownames
    out <- round(out,3)
    out$F <- round(out$F,2)
    out
  }
  
  TTEST_APA <- function(ttest_output){
    data.frame(t = unname(round(ttest_output$statistic,2)),
    p = unname(round(ttest_output$p.value,3)),
    Df = unname(round(ttest_output$parameter,2)))
  }
  
  APA_report <- function(stat_vals){
    if("F" %in% colnames(stat_vals)){
      out <- paste0("F(", stat_vals$Df1, ",", stat_vals$Df2, ") = ",stat_vals$F)
    } else if ("t" %in% colnames(stat_vals)) {
      out <- paste0("t(", stat_vals$Df, ") = ",stat_vals$t)
    }
    if(stat_vals$p==0){
      out <- paste0(out,", p < .001")
    } else {
      out <- paste0(out,", p = ",stat_vals$p)
    }
    out
  }
  
  
  
  #SOME GENERAL USE FUNCTIONS
  
  #math functions
  
  std_err <- function(x) sd(x)/sqrt(length(x))
  ci <- function(x){1.96*std_err(x)}
  specify_decimal <- function(x, k) as.numeric(trimws(format(round(x, k), nsmall=k)))
  ms_in_min <- function(ms){
    (ms/1000)/60
  }
  ms_in_sec <- function(ms){
    (ms/1000)
  }
  
  outlierBool <- function(variable){
    maxval <- quantile(variable)["75%"] + IQR(variable)
    minval <- quantile(variable)["25%"] - IQR(variable)
    return( variable > maxval | variable < minval)
  }
  
  outlierCriticalValues <- function(variable, quant_range = 50){
    upper_quant <- 50 + quant_range/2
    lower_quant <- 50 - quant_range/2
    upper_limit <- quantile(variable, upper_quant/100)
    lower_limit <- quantile(variable, lower_quant/100)
    range <- upper_limit - lower_limit
    minval <- upper_limit + range*1.5
    maxval <- lower_limit - range*1.5
    return(c(minval,maxval))
  }
  
  standardise <- function(variable){
    if(equality_list(variable)){
      return(rep(0,length(variable)))
    } else {
      (variable-mean(variable,na.rm = T)) / sd(variable,na.rm = T) 
    }
  }
  
  center <- function(variable,
                     mean = mean(variable)){
    variable-mean
  }
  
  invlogit <- function(x){
    1 / (1 + exp(-x))
  }
  
  #convenient formatting functions
  
  matrix_to_list <- function(matrix,dim){
    if(dim==1){
      lapply(seq_len(nrow(matrix)),function(i) matrix[i,])
    } else {
      lapply(seq_len(ncol(matrix)),function(i) matrix[,i])
    }
  }

  ifexists <- function(x){
    out <- NA
    varname <- as.character(enexpr(x))
    if(exists(varname)){
      out <- x
    }
    out
  }
  
  count_by <- function(df,count_var,by_var){
    if(length(by_var)>0){
      colname <- by_var[[1]]
      colvar <- unique(df[,colname])
      count_dfs <- sapply(colvar,function(var){
        var_df <- df[df[,colname]==var,]
        count_by(var_df,count_var,by_var[-1]) 
      })
      count_dfs <- apply(count_dfs,2,as.data.frame)
      for(i in 1:length(colvar)){
        count_dfs[[i]][,colname] <- colvar[[i]]
      }
      count_df <- rbind_list(count_dfs)
      #count_df[,colname] <- rep(colvar,each=nrow(count_dfs[[1]]))
      
    } else {
      count_df <- count(df,count_var)
    }
    return(count_df)
  }
  

  shifter <- function(guy,amount=1){
    prev_order <- 0:(length(guy)-1)
    new_order <- (prev_order + amount) %% length(guy)
    new_order <- new_order + 1
    return(guy[new_order])
  }
  
  aggregate_rf <- function(df, variable_colname, groupby_colnames, func = mean){
    #rt stands for retain format; instead of returning a list, returns a df
    groupby_list <- lapply(groupby_colnames, function(colname){df[,colname]})
    df <- aggregate(df[,variable_colname],by=groupby_list,func)
    colnames(df) <- c(groupby_colnames,variable_colname)
    return(df)
  }
  
  sapply_to_df <- function(sapplyto, func, colname.id = "name", colname.var = "output", ...){
    #instead of a list with names, returns a df    
    sapplyto <- as.character(sapplyto)
    x <- lapply(sapplyto, func, ...)
    names_rep <- sapply(x,length)
    names <- rep(sapplyto, names_rep)
    df <- data.frame(names, unlist(x), row.names = 1:length(unlist(x)))
    colnames(df) <- c(colname.id, colname.var)
    return(df)
  }
    
  listreplace <- function(list, values, replacements){
    #e.g. c(1,2,3,4,4) replacing 4s with 7s to make c(1,2,3,7,7)
    #would be (c(1,2,3,4,4), 4, 7)
    newlist <- list
    for (i in 1:length(values)) {
      oldval <- values[i]
      newval <- replacements[i]
      newlist[list==oldval] <- replacements[i]
    }
    return(newlist)
  }
  
  factor_and_reorder <- function(df, colname, reordercol){
    df[,colname] <- factor(df[,colname])
    df[,colname] <- reorder(df[,colname], df[,reordercol])
    return(df)
  }
  
  #misc
  
  make_count_function <- function(up_to=Inf){
    force(up_to)
    x <- -1
    return(function(count = TRUE){
      if(count) x <<- (x + 1) %% up_to
      return(x+1)
    })
  }
  
  start_timer <- function(){
    starttime <- Sys.time()
    return(function(){
      timepassed <- round(as.double(difftime(Sys.time(),starttime,u='mins')),4)
      return(timepassed)
    })
  }
  
  OR <- function(a_list){
    for(el in a_list){
      if(el) return(TRUE)
    }
    return(FALSE)
  }
  
  is.max <- function(a_list){
    out <- a_list==max(a_list)
    #names(out) <- names(a_list)
    out
  }
  
  which_names <- function(an_array){
    out <- which(an_array, arr.ind = T)
    num_dim <- length(dim(an_array))
    all_dimnames <- dimnames(an_array)
    x <- c()
    for(i in 1:num_dim){
      current_dimnames <- all_dimnames[[i]]
      current_index <- out[1,i]
      dimname <- current_dimnames[current_index] 
      x <- append(x,dimname)
    }
    #first_row <- index_array(out,1)
    #row_name <- names(an_array[1,])[first_row[['row']]]
    #col_name <- names(an_array[,1])[first_row[['col']]]
    #return(c(row_name,col_name))
    return(x)
  }
  
  equality <- function(...){
    #pass it as many arguments as you like, and it will tell you whether theyre equal
    equal <- TRUE
    for (i in 1:...length()){
      if (!(prod(...elt(1) == ...elt(i)))) equal <- FALSE
    }
    equal
  }
  
  equality_list <- function(alist){
    prod(mapply(`==`,alist[1:(length(alist)-1)],alist[-1]))
  }
  
  make_recursive <- function(func){
    re_func <- function(alist, ...){
      if (length(alist) > 1){
        func(alist[[1]],re_func(alist[2:length(alist)],...),...)
      }
      else {
        alist[[1]]
      }
    }
    
  }

  
  rbind_list <- make_recursive(rbind)
  paste_list <- make_recursive(paste)
  
  list_to_array <- function(nested_list,
                            array_names = lapply(recursive_names(nested_list),c))
    {
      #x <- array(c(1,1,1,1,1,1,1,1),dim = list(2,2,2),dimnames = list(c('one','two'),c('one','two'),c('one','two')))
      unraveled_list <- unlist(nested_list)
      num_dims <- lapply(array_names,length)
      return(array(unraveled_list,dim=num_dims,dimnames=array_names))
  }
  recursive_names <- function(nested_list){
    if(typeof(nested_list[[1]])=="list"){
      #this section is just to catch potential problems
      if(length(nested_list) > 1){
        names_at_sub_level <- lapply(nested_list, names2)
        names_are_equal <- prod(mapply(`==`,names_at_sub_level[1:(length(names_at_sub_level)-1)],names_at_sub_level[-1]))
        if(!names_are_equal){
          stop(paste0("your list indices do not match at nest level,", nest_counter, " a list of consistent subindices will not work for this level across its elements"))
        }
      }
      #this is the actual work thats done
      sublist_names <- recursive_names(nested_list[[1]])
      new_names <- list(names2(nested_list))
      names <- append(sublist_names,new_names)
      nest_counter <<- nest_counter + 1
    }else{
      names <- list(names2(nested_list))
      nest_counter <<- 1
    }
    return(names)
  }
  
  names2 <- function(a_list){
    names <- names(a_list)
    if(is.null(names)){names <- as.character(1:length(a_list))}
    return(names)
  }
  
  array_to_list <- function(an_array){
    df <- as.data.frame(an_array)
    df_to_list(df)
  }
  
  df_to_list <- function(df, colsep = '\\.'){
    colnames <- colnames(df)
    varnames <- lapply(colnames,strsplit,"\\.")
    varnames <- lapply(varnames,`[[`,1)
    current_var <- unlist(lapply(varnames,`[[`,1))
    num_var_left <- length(varnames[[1]]) - 1
    #num_var <- length(varnames[[1]])
    df_list <- list()
    for(var in unique(current_var)){
      var_col <- colnames[grepl(var,colnames)]
      var_df <- df[,var_col]
      if(num_var_left==0){
        var_names <- row.names(df)
        var_values <- unlist(var_df)
        names(var_values) <- var_names
        df_list[[var]] <- var_values
      } else {
        colnames(var_df) <- gsub(paste0(var,"\\."),"",colnames(var_df))
        df_list[[var]] <- df_to_list(var_df, colsep)
      }
    }
    return(df_list)
  }
  
  array_to_df <- function(an_array){
    df<- as.data.frame(an_array)
    library(reshape2)
    mdf <- melt(df)
    vars <- sapply(as.character(mdf$variable),strsplit,"\\.")
    
    #num_vars <- stringr::str_count(as.character(mdf$variable[1]),"\\.") + 1
    #colnames <- paste0("Var",1:num_vars)
    #mdf[,colnames] <- NA
    # for(row_i in 1:nrow(mdf)){
    #   var_string <- as.character(mdf$variable[row_i]) 
    #   vars <- strsplit(var_string,"\\.")
    #   for(col_i in 1:length(colnames)){
    #     colname <- colnames[col_i]
    #     mdf[row_i,colname] <- vars[col_i]
    #   }
    # }
  }
  
  index_list <- function(a_list,indices){
    indices <- unlist(indices)
    dim = length(indices)
    current_index <- indices[1]
    remaining_indices <- indices[2:dim]
    if(dim>1){
      if(is.na(current_index)){
        sublists <- sapply(a_list,index,remaining_indices)
        return(sublists)
      } else {
        sublist <- a_list[[current_index]]
        return(index(sublist,remaining_indices))
        }
    } else {
      if(is.na(current_index)){
        return(unlist(a_list))
      } else {
        return(a_list[[current_index]])
      }
    }
  }
  #x <- list( list(list(1,2,3),list(4,5,6),list(7,8,9)), list(list(1,2,3),list(4,5,6),list(7,8,9)))
  #index(x,c(1,NA,1))
  
  index_array <- function(an_array,
                          dim,
                          loc,
                          replacement=NULL){
    pre_dim_index <- paste_list(rep(",",dim-1))
    post_dim_index <- paste_list(rep(",",length(dim(an_array)) - dim))
    out_string <- paste0("an_array[",pre_dim_index,loc,post_dim_index,"]")
    if(!is.null(replacement)){
      out_string <- paste0(out_string," <- ",replacement)
      out_expr <- str2expression(out_string)
      eval(out_expr)
      an_array
    } else {
      out_expr <- str2expression(out_string)
      eval(out_expr)
    }
  }
  

  
  append_array <- function(list_of_arrays){
    #check first
    arrays_are_uniform <- prod(apply(X=sapply(array_list,dimnames),MARGIN = 1,FUN = equality_list))
    print(arrays_are_uniform)
    if(!arrays_are_uniform){stop("arrays are not the same structure")}
    #now make appended array
    list_names <- list(names(list_of_arrays))
    array_names <- dimnames(list_of_arrays[[1]])
    array_names <- append(list_names,array_names,0)
    num_dim <- dim(list_of_arrays[[1]])
    num_dim <- append(num_dim,length(list_of_arrays),0)
    return(array(unlist(list_of_arrays),dim=num_dim,dimnames = array_names))
  }
  
  dimnames_i <- function(an_array, dim_index){
    return(dimnames(an_array)[[dim_index]])
  }
  # nested_list1 <-list(yes=list(hello=1,goodbye=2),no=list(hello=1,goodbye=2))
  # nested_list2 <- list(yes=list(mer=1,ci=2,fer=3),no=list(mer=1,ci=2,fer=3))
  # as_array1 <- list_to_array(nested_list1)
  # as_array2 <- list_to_array(nested_list2)
  # array_list <- list(x=as_array2,y=as_array2)
  # x <- append_array(array_list)
  
  addcol <- function(df_to, 
                     df_from, 
                     id.col = colnames(df_from)[colnames(df_from) %in% colnames(df_to)], #id col is assumed to be common to both dfs (unless specified otherwise)
                     val.col = colnames(df_from)[!(colnames(df_from) %in% colnames(df_to))] #new col is assumed not to be in the col it is moved into (unless specified otherwise)
  ){
    #if dfs have different row lengths, need to expand df_from
    #assumes the longer df is df_to
    nrow_to <- nrow(df_to)
    nrow_from <- nrow(df_to)
    num_rep <- nrow(df_to) / nrow(df_from)
    
    requires_extending <- num_rep > 1
    if(requires_extending){
      if (nrow(df_to) %% nrow(df_from) != 0){
        error("your two dfs have differing numbers of rows")
      } else {
        store_df_from <- df_from
        for (rep in 1:(num_rep-1)) {
          df_from <- rbind(df_from, store_df_from)
        }
      }
    }
    #if already all the id values align, no need to reorder
    requires_reordering <- !(prod(df_from[,id.col] ==  df_to[,id.col]))
    if(requires_reordering){
      if (typeof(df_from[,id.col]) != typeof(df_to[,id.col])){
        error("your id vars dont match in type")
      } else {
        #lazy solution to the id values not aligning is to reoder both dfs by that column
        df_from <- df_from[order(df_from[,id.col]),]
        df_to <- df_to[order(df_to[,id.col]),]
      }
    }
    df_to[,val.col] <- df_from[,val.col]
    return(df_to)
  }

  
  #GENERAL USE GET DATA FUNCTION
  getData <- function(taskdir = paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/Data"),
                      filenames = list.files(taskdir)
                      ){
    df <- data.frame()
    for (filename in filenames) {
      csv_ext_check <- substr(filename, nchar(filename)-3, nchar(filename))
      if(csv_ext_check == ".csv"){
        filename <- paste0(taskdir,'/',filename)
        df <- rbind.fill(df, as.data.frame(read.csv(filename)))
      }
    }
    return(df)
  }
  
  #note to self
  #i did a few things to permanenetly change te data to make it more convinient to 
  #use, and i no longer have the code that does this
  #one of those things is to correct the trial_index that is automated by jsPsych
  #and so it restarts the numbering whenver a new timeline (block) is introduced
  #to fix that numbering, run 
  #df <- fix_BlockCSVRestartedNumbering(df, "trial_index")
  
  #participant  "386290" seems to be one of the first participants i ran
  #and so they are missing some columns like "BlockNum" and "SlideClass/SlideType"
  #I manually filled those in using excel
  
}
