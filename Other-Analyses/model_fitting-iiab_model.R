#--------------------------------------------------------------------------#
#               GET PPT DATA AND FUNCTIONS FROM GETDATA.R                  #
#--------------------------------------------------------------------------#

curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
subdir <- paste0(curdir, "/supporting/")
setwd(subdir)
source("bernEst_startup.R")
setwd(curdir)

#-------------------------------------------#
#### GET LOG LIKLIHOOD OF DATA ####
#-------------------------------------------#
##### FUNCTIONS #####
BernCP <- function(Data,AlphaB,BetaB,AlphaA,BetaA,pc) {
  
  L <- 1:length(Data)
  Sf = cumsum(Data)
  Ff = L - Sf
  Sr = sum(Data) - Sf
  Lr = seq((length(Data)-1),0,-1)
  Fr = Lr - Sr
  
  LgPstLkf = log(beta(AlphaB + Sf, BetaB + Ff))
  LgPstLkr = log(beta(AlphaA + Sr, BetaA + Fr))
  LgLkFun = LgPstLkf + LgPstLkr
  LgLkFunN = LgLkFun + abs(max(LgLkFun))
  LkFun = exp(LgLkFunN)
  CPmean = sum((1:length(Data))*LkFun)/sum(LkFun)
  
  RelLkhd = (sum(LkFun[1:length(LkFun)])/(length(Data)))/LkFun[length(LkFun)]
  
  pNC = 1-pc
  PriorOdds = L[length(L)]*pc/pNC
  Odds = PriorOdds*RelLkhd
  BF = RelLkhd
  CP = round(CPmean)
  
  return(list(CP,Odds,LkFun,CPmean))
  
}

berndkl <- function(p1,p2){
  q1 = 1 - p1
  q2 = 1 - p2
  
  LVp1_0 = p1==0
  LVp1_1 = p1==1
  
  LVp2 = (p2==0) | (p2==1)
  
  LVeq = p1==p2
  
  Dkl = rep(NaN,length(p1))
  
  Dkl[LVeq] = 0
  Dkl[!LVeq & LVp2] = Inf
  Dkl[!LVeq & LVp1_0] = -log(q2[!LVeq&LVp1_0])
  Dkl[!LVeq&LVp1_1] = -log(p2[!LVeq&LVp1_1])
  
  Dkl[!LVeq&!LVp1_0&!LVp1_1&!LVp2] = p1[!LVeq&!LVp1_0&!LVp1_1&!LVp2] * 
    log(p1[!LVeq&!LVp1_0&!LVp1_1&!LVp2] / p2[!LVeq&!LVp1_0&!LVp1_1&!LVp2]) + 
    q1[!LVeq&!LVp1_0&!LVp1_1&!LVp2]*log(q1[!LVeq&!LVp1_0&!LVp1_1&!LVp2] / 
                                           q2[!LVeq&!LVp1_0&!LVp1_1&!LVp2])
  return(Dkl)
}
BernCPKLfun <- function(Data, alpha_p, beta_p, alpha_c, beta_c, KLcrit, BFcrit){
  p_hat = alpha_p/(alpha_p + beta_p)
  alpha_c_h = alpha_c
  beta_c_h = beta_c
  pc_hat = alpha_c/(alpha_c + beta_c)
  hyper_c_rec = list(c(alpha_c, beta_c))
  hyper_p_rec = list(c(alpha_p, beta_p))
  
  alpha_a = alpha_p
  beta_a = beta_p
  
  Nc = alpha_c
  NcInit = Nc
  
  cp = c(0)
  DP = c(0)
  RTcpRprt = c() #real time change point report
  RTphatRprt = rep(p_hat,length(Data))
  D = Data
  Det = 0
  
  while( length(D) > 20) { #20 ){
    N = 1:length(D)
    p_o = cumsum(D)/N
    KL = berndkl(p_o, matrix(rep(p_hat,length(N)),length(N))) #repmat(p_hat,N[length(N)],1) )
    E = N*KL
    
    if((Det >= length(E))  |  ( !any( E[(Det+1):length(E)] > KLcrit) )){ #E[(Det+1):length(E)] > KLcrit) ){ #for some reason this gave me a lot of trouble and i am not even sure why Det is being chancged so like whatever
      break #break out of while loop
    } else {
      Det = which(E[(Det+1):length(E)] > KLcrit)[1] + Det
    }
    
    DP = append(DP, (cp[length(cp)]) + Det, length(DP))
    alpha_c = Nc + alpha_c_h
    beta_c = DP[length(DP)]-Nc + beta_c_h
    pc_hat = alpha_c/(alpha_c + beta_c)
    bern_CP_list <- BernCP(D[1:Det],alpha_p,beta_p,alpha_a,beta_a,pc_hat)
    CP <- unlist(bern_CP_list[1])
    cp_local = CP
    cp_global = CP + cp[length(cp)] #- 1
    PstO <- unlist(bern_CP_list[2])
    LF <- unlist(bern_CP_list[3])
    
    #print(PstO)
    
    if(PstO > BFcrit){
      RTcpRprt = append(RTcpRprt,DP[length(DP)],length(RTcpRprt))
      cp = append(cp,cp_global,length(cp))
      #cp = append(cp,CP+cp[length(cp)],length(cp))
      #CL = append(CL,Lmts,length(CL))
      alpha_p = 1 + sum(D[(cp_local+1):Det])
      beta_p = 1 + length(D[(cp_local+1):Det])-sum(D[(cp_local+1):Det])
      hyper_p_rec = append(hyper_p_rec, list(c(alpha_p,beta_p)), length(hyper_p_rec))
      hyper_c_rec = append(hyper_c_rec, list(c(alpha_c,beta_c)), length(hyper_c_rec))
      D <- D[(cp_local+1):length(D)]
      Det <- 1 + DP[length(DP)] - cp[length(cp)]
      p_hat = alpha_p/(alpha_p + beta_p)
      Nc = Nc+1
    } else if (PstO < BFcrit && length(cp)>1) { #uhhh idk what the line here (580 says)
      TestPrevCP_data = Data[ (cp[length(cp)-1]+1) :DP[length(DP)]]
      
      bern_CP_list <- BernCP(TestPrevCP_data,alpha_p,beta_p,alpha_a,beta_a,pc_hat)
      CPr <- unlist(bern_CP_list[1])
      cpr_local <- CPr
      cpr_global <- (cp[length(cp)-1]) + CPr #- 1
      Or <- unlist(bern_CP_list[2])
      LF <- unlist(bern_CP_list[3])
      
      #previous change point was in retrospect not justified
      if(Or < BFcrit){
        cp = cp[1:(length(cp)-1)]
        #CL = CL[1:length(CL)-1]
        hyper_p_rec= hyper_p_rec[1:(length(hyper_p_rec)-1)]
        hyper_c_rec= hyper_c_rec[1:(length(hyper_c_rec)-1)]
        Nc = Nc-1
        
        alpha_p = unlist(hyper_p_rec[[length(hyper_p_rec)]][1]) + sum(TestPrevCP_data)
        beta_p = unlist(hyper_p_rec[[length(hyper_p_rec)]][2]) + length(TestPrevCP_data)- sum(TestPrevCP_data)
        p_hat = alpha_p/(alpha_p + beta_p)
        
        alpha_c = Nc + alpha_c_h
        beta_c = DP[length(DP)] - Nc + beta_c_h
        D = Data[(cp[length(cp)]+1):length(Data)]
        Det = 1 + DP[length(DP)] - cp[length(cp)]
        
      } else { #previous estimate of p_g was wrong
        cp[length(cp)] = cpr_global
        #CL(end,:) = Lmtsr
        alpha_p = 1 + sum(Data[ (cp[length(cp)]+1) :DP[length(DP)]]) #sum(Data[ (cp[length(cp)]+1) :DP[length(DP)]])
        beta_p = 1 + length(Data[ (cp[length(cp)]+1) :DP[length(DP)]]) - sum(Data[ (cp[length(cp)]+1) :DP[length(DP)]])
        p_hat = alpha_p/(alpha_p + beta_p)
        D = Data[ (cp[length(cp)]+1) :length(Data)]# i quote "Is this necessary?"
        Det = 1 + DP[length(DP)] - cp[length(cp)]
      }
    } else {
      #there was no apparent change point and no change point has so far
      # been found, in which case, it must be that the initial estimate of
      # p_g was wrong
      alpha_p = 1 + sum(D[1:Det])
      beta_p = 1 + length(D[1:Det]) - sum(D[1:Det])
      p_hat = alpha_p/(alpha_p + beta_p)
    }
    RTphatRprt[(DP[length(DP)]+1):length(Data)] <- p_hat
  } #end of while loop
  
  CmRec = cumsum(Data)
  CmTrialNum = 1:length(Data)
  cp = append(cp, length(CmRec), length(cp)) #makes final trial a nominal change point
  
  Succ = unlist(sapply(2:length(cp),function(thiscp){CmRec[cp[thiscp]] - CmRec[cp[thiscp-1]]}))
    #CmRec[cp[2] : diff(CmRec([cp(2:end);length(Data)])) #number of successes from one change point to the next
  
  NumTrls = unlist(sapply(2:length(cp),function(thiscp){CmTrialNum[cp[thiscp]] - CmTrialNum[cp[thiscp-1]]}))
  
  ps = Succ/NumTrls
  
  return(list(cp,RTphatRprt,DP))
  
}

##### variable assignment #####

alpha_p = .5
beta_p = .5
alpha_c = .04 
beta_c = .005

KLcrit = 1.92 #corresponds of a pvalue of 0.05
BFcrit = 1 #corresponds to even odds


#### fit to participant data ####

#TEST
#trialrings = df$TrialRingBoolean[1:(999+51)]
#trialrings = trialrings[!is.na(trialrings)]
#BernCPKLfun(trialrings, 0.5,0.5,0.04,0.005,1.92,1)
#BernCP(trialrings, alpha_p, beta_p,alpha_p, beta_p,0.005)

create_sqrderror_function <- function(pptests, obsv) {
  force(pptests)
  force(obsv)
  return(
    function(param) {
      KLcrit <- param[1]
      BLcrit <- param[2]
      alpha_p = .5
      beta_p = .5
      alpha_c = .04 
      beta_c = .005
      iiab_resp <- BernCPKLfun(obsv, alpha_p, beta_p, alpha_c, beta_c, KLcrit, BFcrit)
      iiab_resp <- iiab_resp[[2]]
      return(sum((iiab_resp - pptests)**2,na.rm = TRUE))
    }
  )
}

##now we pull it all together##

#example ppt: 
# pptid <- '348394'
# error <- create_sqrderror_function(pptdf$bernoulliEstimate, pptdf$TrialRingBoolean)
# opt <- optim(error,par=c(KLcrit=1.92, BFcrit=1))
# params <- opt$par
# KLcrit <- 0#as.numeric(params['KLcrit'])
# BFcrit <- as.numeric(params['BFcrit'])
#trialrings = task_df$TrialRingBoolean[task_df$PptId==pptid]
#CPKL_ests <- BernCPKLfun(trialrings, alpha_p, beta_p, alpha_c, beta_c, KLcrit, BFcrit)
#CPKL_ests <- CPKL_ests[[2]][1:999]

split_df <- split(reduced_df, reduced_df$PptId)

ppt_error_fncts <- sapply(split_df, function(x){
  create_sqrderror_function(x$bernoulliEstimate, x$TrialRingBoolean)
})

tmp <- sapply(ppt_error_fncts,optim,par=c(KLcrit=1.92, BFcrit=1))
optim_KLcrits <- sapply(colnames(tmp), function(x) {
  params <- tmp['par',x][[1]]
  return(as.numeric(params['KLcrit']))
})
optim_BFcrits <- sapply(colnames(tmp), function(x) {
  params <- tmp['par',x][[1]]
  return(as.numeric(params['BFcrit']))
})

optim_berncpkl_ests <- sapply(names(optim_KLcrits), function(x){
  print(x)
  KLcrit <- optim_KLcrits[[x]]
  BFcrit <- optim_BFcrits[[x]]
  observations <- split_df[[x]][['TrialRingBoolean']]
  #print(c(length(observations),x))
  berncpkl = BernCPKLfun(observations, alpha_p, beta_p, alpha_c, beta_c, KLcrit, BFcrit)
  berncpkl[[2]][1:999]
})

#names(optim_berncpkl_ests) <- unique(df$PptId)
iiab_df <- melt(data.frame(optim_berncpkl_ests))
iiab_df$variable <- sub('X','',iiab_df$variable)
#missing_ppt <- unique(iiab_df$variable)[!(unique(iiab_df$variable) %in% unique(reduced_df$PptId))]
#iiab_df <- iiab_df[iiab_df$variable != missing_ppt,]
if (prod(iiab_df$variable == reduced_df$PptId)) {
  reduced_df$IiabEsts <- iiab_df$value
}

ggplot(data=reduced_df) +
  geom_line(aes(x=TrialNum,y=bernoulliEstimate,color=Condition)) +
  geom_line(aes(x=TrialNum,y=IiabEsts)) +
  facet_wrap(~Condition+PptId,labeller = labeller(PptId = optim_BFcrits))

#-------------------------------------------#
#### GET LOG LIKLIHOOD OF DATA ####
#-------------------------------------------#

get_loglik <- function(ppt_ests,model_ests){
  stand_dev <- sd(ppt_ests - model_ests)
  model_error <- pnorm(model_ests - ppt_ests, 0, stand_dev, log = TRUE)
  #print(model_error)
  #mod_err_log <- log(model_error)
  return(sum(model_error))
}

ppt_loglik_iiab <- sapply(names(optim_KLcrits), function(x){
  KLcrit <- optim_KLcrits[[x]]
  BFcrit <- optim_BFcrits[[x]]
  observations <- split_df[[x]][['TrialRingBoolean']]
  berncpkl = BernCPKLfun(observations, alpha_p, beta_p, alpha_c, beta_c, KLcrit, BFcrit)
  model_ests <- berncpkl[[2]][1:999]
  ppt_ests <- split_df[[x]][['bernoulliEstimate']]
  get_loglik(model_ests, ppt_ests)
})

ppt_numUpdates_iiab <- sapply(names(optim_KLcrits), function(pptid){
  ppt_df <- reduced_df[reduced_df$PptId==pptid,]
  #View(data.frame(ppt_df$ThresholdedDeltaEsts[1:(length(ppt_df$TrialNum)-1)],ppt_df$ThresholdedDeltaEsts[2:(length(ppt_df$TrialNum))]))
  update_event <- !(ppt_df$IiabEsts[1:(length(ppt_df$TrialNum)-1)] == ppt_df$IiabEsts[2:(length(ppt_df$TrialNum))])
  return(sum(update_event))
})

temp_df <- data.frame(PptId = names(ppt_loglik_iiab),LogLikIIAB = unlist(ppt_loglik_iiab))
summarystats_df <- merge(summarystats_df, temp_df)
temp_df <- data.frame(PptId = names(ppt_numUpdates_iiab),NumUpdatesIIABModel = unlist(ppt_numUpdates_iiab))
summarystats_df <- merge(summarystats_df, temp_df)

ggplot(summarystats_df,aes(x=LogLikIIAB-LogLikDelta))+
  geom_histogram()+
  facet_wrap(~Condition)
# 
# ggplot()+
#   geom_histogram(data=summarystats_df[summarystats_df$Condition=="Manual",],aes(x=LogLikIIAB-LogLikDelta),fill="coral2",alpha=0.5,bins=50)+
#   geom_histogram(data=summarystats_df[summarystats_df$Condition=="Automatic",],aes(x=LogLikIIAB-LogLikDelta),fill="cyan3",alpha=0.5,bins=50)
# 

ggplot(summarystats_df)+
  geom_histogram(aes(x=LogLikIIAB), fill="cyan3",alpha=0.5,bins=50)+
  geom_histogram(aes(x=LogLikDelta), fill="coral2",alpha=0.5,bins=50)

# ggplot(summarystats_df)+
#   geom_histogram(aes(x=LogLikIIAB), fill="cyan3",alpha=0.5,bins=50)+
#   geom_histogram(aes(x=LogLikDelta), fill="coral2",alpha=0.5,bins=50)


