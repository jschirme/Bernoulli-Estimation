curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
subdir <- paste0(curdir, "/supporting/")
setwd(subdir)
source("bernEst_startup.R")

#------------------------------------------------------------#
#                     FITTING THE MODELS                     #
#------------------------------------------------------------#

make_lower <- function(param_names){
  sapply(param_names,function(name){if(grepl("mean",name)){return(-5)}else{return(0)}})
}
make_start <- function(p_names, extensive=FALSE){
  if(extensive){
    steps <- c(100,1,0.25)
    bounds <- matrix(c(100,0,0.25, #100,1,0.25,600,3,1.5
                       600,3,1.5),ncol=2)
  } else {
    steps <- c(200,2,0.5)
    bounds <- matrix(c(200,1,0.25,
                       600,3,1.25),ncol=2)
  }
  param_name_to_start_bounds <- function(p_name){
    if(grepl("scale",p_name)){
      seq(bounds[1,1],bounds[1,2],steps[1])
    } else if (grepl("mean",p_name)){
      seq(bounds[2,1],bounds[2,2],steps[2])
    } else {
      seq(bounds[3,1],bounds[3,2],steps[3])
    }
  }
  out_bounds = lapply(p_names,param_name_to_start_bounds)
  out <- expand.grid(out_bounds)
  colnames(out) <- p_names
  out
}

nls_fit <- function(start, a_df, m_expr){
  nls(m_expr,
      data=a_df,
      start = start,
      lower = make_lower(names(start)),
      algorithm="port")
}

find_nls_fit <- function(a_df,model_name,
                         find_best=TRUE,
                         model_expr = NULL){
  search_for_nls <- function(){
    search_while_true <-  expr((is.null(out_nls) | find_best) & i < nrow(starts))
    out_nls <- NULL
    while(eval(search_while_true)){
      i = i +  1
      start = starts[i,]
      if(length(start)==1){names(start) <- dimnames(starts)[[2]]}
      tryCatch({
        new_nls <- nls_fit(start,a_df,model_expr)
        if(is.null(out_nls)){
          out_nls <- new_nls
        } else if(logLik(new_nls) > logLik(out_nls)){
          out_nls <- new_nls
        }
      },error = function(e){})
    }
    return(out_nls)
  }
  if(is.null(model_expr)){
    g <- quote(NumAdjs)
    h <- model_exprs[[model_name]]
    model_expr <- substitute(g~h, list(g=g, h=h))
  }
  params <- expr_to_param(model_expr)
  starts <- make_start(params)
  lower <- make_lower(params)

  out_nls <- NULL
  i = 0
  
  out_nls <- search_for_nls()
  if(is.null(out_nls)){
    print("second try")
    starts <- make_start(params,extensive=TRUE)
    find_best <- FALSE
    out_nls <- search_for_nls()
    if(is.null(out_nls)){
      print("double failure")
    }
  }
  return(out_nls)
}

find_ppt_fit <- function(pptid,model_name,find_best=TRUE,a_df = adj_df){
  print(pptid)
  ppt_df <- a_df[a_df$PptId==pptid,]
  find_nls_fit(ppt_df,model_name,find_best)
}

## ---- fit-models

# takes a long time to run this code below, so i've commented it out and saved the rds
 # nls_auto_avg = find_nls_fit(adj_df[!manual_boolean,],"auto2")
 # auto2_param <-  nls_auto_avg$m$getAllPars()   
 # individual_nls_fits <- list(auto1 = lapply(automatic_ids,find_ppt_fit,"auto1"),
 #                            auto2 = lapply(automatic_ids,find_ppt_fit,"auto2"),
 #                            auto_avg = nls_auto_avg,
 #                            man2 = lapply(manual_ids,find_ppt_fit,"man2_man"),
 #                            man3 = lapply(manual_ids,find_ppt_fit,"man3_man"),
 #                            man2_auto = lapply(automatic_ids,find_ppt_fit,"man2"),
 #                            man3_auto = lapply(automatic_ids,find_ppt_fit,"man3"))
#saveRDS(individual_nls_fits,"ind_ppt_gauss_series.rds")

# nls_auto_avg = find_nls_fit(adj_df[!manual_boolean,],"auto2")
# auto2_param <-  nls_auto_avg$m$getAllPars()   
# group_nls_fits <- list(auto1 = find_nls_fit(adj_df[!manual_boolean,],"auto1"),
#                       auto2 =  nls_auto_avg,
#                       man2_man = find_nls_fit(adj_df[manual_boolean,],"man2"),
#                       man3_man = find_nls_fit(adj_df[manual_boolean,],"man3"),
#                       man2_auto = find_nls_fit(adj_df[!manual_boolean,],"man2"),
#                       man3_auto = find_nls_fit(adj_df[!manual_boolean,],"man3"))
#saveRDS(group_nls_fits,"group_gauss_series.rds")

individual_nls_fits <- readRDS(file = "ind_ppt_gauss_series.rds")
group_nls_fits <- readRDS(file = "group_gauss_series.rds")

#------------------------------------------------------------#
#         FUNCTIONS FOR WORKING WITH THE MODELS              #
#------------------------------------------------------------#

## ---- gaussian-series

expr_to_param <- function(model_expr){
  ws <- unlist(str_extract_all(deparse(model_expr), "\\w+"))
  is_param_name <- sapply(ws,function(w){grepl(".+\\d",w)})
  if("auto2_param" %in% ws){
    not_params <- c(which(ws=="auto2_param"), which(ws=="auto2_param")+1)
    is_param_name[not_params] <- FALSE
    #is_param_name <- sapply(ws,function(w){grepl(".+\\d",w) & !grepl("2",w)})
  }
  ws[is_param_name]
}

manual_boolean <- adj_df$Condition=="Manual"
manual_ids <- unique(adj_df$PptId[manual_boolean])
automatic_ids <- unique(adj_df$PptId[!manual_boolean])

model_exprs <- list(auto1 = expr(scale1*dnorm(AdjSize,mean1,stdev1)),
                    auto2 = expr(scale1*dnorm(AdjSize,0,stdev1) + scale2*dnorm(AdjSize,mean2,stdev2)),
                    man2 = expr(scale1*dnorm(AdjSize,0,stdev1) + auto2_param[["scale2"]]*dnorm(AdjSize,auto2_param[["mean2"]],auto2_param[["stdev2"]])),
                    man3 = expr(scale1*dnorm(AdjSize,0,stdev1) + auto2_param[["scale2"]]*dnorm(AdjSize,auto2_param[["mean2"]],auto2_param[["stdev2"]]) + scale3*dnorm(AdjSize,mean3,stdev3)),
                    man3_alt  = expr(scale1*dnorm(AdjSize,0,stdev1) + scale2*dnorm(AdjSize,mean2,stdev2) + scale3*dnorm(AdjSize,mean3,stdev3)))
param_names <- sapply(model_exprs,expr_to_param)
model_names <- c("auto1","auto2","man2_man","man3_man")
model_descriptions <- c(auto1 = "A. Automatic\n One Gaussian",auto2 = "B. Automatic\n Two Gaussian",man2 = "C. Manual\n Two Gaussian",man3 = "D. Manual\n Three Gaussian")
model_conditions <- c("Automatic","Automatic","Manual","Manual")

auto2_param <-  individual_nls_fits[["auto_avg"]]$m$getAllPars()   

unlist_nls <- function(nls_fits,
                       fit_names = names(nls_fits)){
  unlist(nls_fits[fit_names],recursive=F)
}

AIC_diffs <- function(nls_fits,
                      model_names = list(c("man2_man","man3_man"),c("man2_auto","man3_auto")),
                      f = AIC){
  lapply(model_names,function(model_name){
    mod1 <- nls_fits[[model_name[1]]]
    mod2 <- nls_fits[[model_name[2]]]
    valid_nls <- !sapply(mod1,is.null) & !sapply(mod2,is.null)
    sapply(mod1[valid_nls],f) - sapply(mod2[valid_nls],f)
  })
}

nls_to_function <- function(nls_obj,
                            rescale_nls=FALSE){
  rescale <- function(params){
    is_scale <- sapply(names(params),grepl,pattern="scale")
    total_scale <- sum(params[is_scale])
    params[is_scale] <- params[is_scale]/total_scale
    params
  }
  set_params <- function(){
    params <- nls_obj$m$getAllPars()
    if(rescale_nls){params <- rescale(params)}
    for(param_name in names(params)){
      ex <- str2expression(paste0(param_name,"<-",params[param_name]))
      env <- parent.env(environment())
      eval(ex,env)
    }
  }
  str_model_expr <- paste_list(deparse(nls_obj$m$formula()))
  str_model_expr <- strsplit(str_model_expr," ~ ")[[1]][2]
  model_expr <- str2lang(str_model_expr)
  set_params()
  return(function(AdjSize){
    eval(model_expr)
  })
}

nls_to_model_name <- function(nls_obj){
  str_model_expr <- paste_list(deparse(nls_obj$m$formula()))
  str_model_expr <- strsplit(str_model_expr," ~ ")[[1]][2]
  model_expr <- str2lang(str_model_expr)
  str_model_exprs <- sapply(model_exprs,function(mexpr){paste_list(deparse(mexpr))})
  names(which(sapply(str_model_exprs,`==`,str_model_expr)))
}

est_G_numadj <- function(model_name,G_id){
  if(grepl("man",model_name) & G_id==2){ 
    #manual participants are not fit to gaussian 2; they take it from automatic participants (their distribution of known updating events)
    f <- function(x){auto2_param["scale2"]*dnorm(x,auto2_param["mean2"],auto2_param["mean2"])}
    est_numadj <- integrate(f,-100,100)$value
  } else {
    mk_f <- function(id_num){
      return(function(AdjSize){
        G_scale <- mod_param[paste0("scale",G_id),id_num]
        G_sd <- mod_param[paste0("stdev",G_id),id_num]
        G_scale*dnorm(AdjSize,mean=0,sd=G_sd)
      })
    }
    mod_param <- get_nls_param(model_name)
    fs <- sapply(1:dim(mod_param)[2],mk_f)
    est_numadj <- sapply(fs,function(f){
      integrate(f,-100,100)$value
    })
  }
  round(est_numadj,2)
}

total_est_numadj <- function(model_name){
  mod_param_names <- dimnames(get_nls_param(model_name))[[1]]
  mod_Gs <- unique(sapply(mod_param_names,function(param_name){gsub("[a-z]","",param_name)}))
  num_G <- max(as.numeric(mod_Gs))
  G_estadj <- lapply(1:num_G,est_G_numadj,model_name=model_name)
  sum(sapply(G_estadj,mean))
}

quick_plot <- function(nls_mods, mod_descr,
                       window_bounds = c(-15,15)){
  window <- window_bounds[1]:window_bounds[2]
  in_window <-  c(-100:100) >= window_bounds[1] & c(-100:100) <= window_bounds[2]
  
  nls_mods <- nls_mods[!sapply(nls_mods,is.null)]
  mod_fs <- lapply(nls_mods,nls_to_function)
  mod_names <- sapply(nls_mods,nls_to_model_name)
  mod_names <- mod_descr[mod_names]
  nls_adjsize <- seq(window_bounds[1],window_bounds[2],0.25)
  mod_estimates <- lapply(mod_fs,function(f){f(nls_adjsize)})
  mod_data <- lapply(nls_mods,function(nls_obj){nls_obj$m$getEnv()$NumAdjs[nls_obj$m$getEnv()$AdjSize >= window_bounds[1] & nls_obj$m$getEnv()$AdjSize <= window_bounds[2]]})
  
  nls_df <- data.frame(PptIndex = rep(1:length(nls_mods),each=length(nls_adjsize)),
                       AdjSize = rep(nls_adjsize,length(nls_mods)),
                       Estimates = round(unlist(mod_estimates),4),
                       Model = rep(mod_names,each=length(nls_adjsize)))
  
  gg_df <- data.frame(PptIndex = rep(1:length(nls_mods),each=length(window)),
                      AdjSize=rep(window,length(nls_mods)),
                      NumAdj=unlist(mod_data),
                      Model=rep(mod_names,each=length(window)))
  p <- ggplot(gg_df,aes(x=AdjSize,y=NumAdj))+
    stat_summary(geom="line",fun=mean,color="red",size=1,alpha=0.5)+
    stat_summary(geom="line",fun=mean,data=nls_df,aes(x=AdjSize,y=Estimates))+  
    xlab("Adjustment Size")+
    ylab("Number of Adjustments")
  p <- p   +  facet_wrap(~Model,ncol = length(unique(mod_names)))
  p <- make_APA(p)
  p <- p + theme(axis.title.x = element_text(size = 13),
                 axis.title.y = element_text(size = 13),
                 strip.text.x = element_text(size = 11),
                 panel.spacing = unit(1.5, "lines"))
  return(p)
}

get_nls_param <- function(model_name){
  nls_fits <- individual_nls_fits[[model_name]]
  nls_fits <- nls_fits[!sapply(nls_fits,is.null)]
  sapply(nls_fits,function(mod){mod$m$getAllPars()})
}

#------------------------------------------------------------#
#                     MODEL ANALYSES                         #
#------------------------------------------------------------#

## ---- manuscript-analyses

AIC_diff_auto12 <- unlist(AIC_diffs(individual_nls_fits,list(c("auto1","auto2"))))
AIC_diff_man23 <- unlist(AIC_diffs(individual_nls_fits,list(c("man2_man","man3_man"))))
AIC_diff_auto23 <- unlist(AIC_diffs(individual_nls_fits,list(c("man2_auto","man3_auto"))))

p <- quick_plot(unlist_nls(individual_nls_fits,c("auto1","auto2","man2_man","man3_man")),model_descriptions)

ttest_stdev1 <- APA_format_vals(t.test(get_nls_param("auto2")["stdev1",],get_nls_param("man2_man")["stdev1",]))

