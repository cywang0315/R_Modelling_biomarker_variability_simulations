###Simulation study: case 4
###Program in a parallel way, all functions required are stored in "source_function_20.R".

library(snowfall)
library(parallel)

################
sfInit(parallel = TRUE, cpus = detectCores())
sfLibrary(survival)
sfLibrary(JM)
sfLibrary(joineR)
sfLibrary(dplyr) # near
sfLibrary(statmod)
sfLibrary(progress)
sfLibrary(MASS)
sfLibrary(mvtnorm)
sfLibrary(tensor)
sfSource(here::here("source_function_20.R"))

SIMULATE=function(s){
  ##Generate simulation data
  Data=gendat(s,obsmax=6,gammatrue=-1,alpha1true=0.3,alpha2true=0.3,sigmatrue=0.4)
  data=Data$data
  
  ##Determine the number and locations of interior knots:
  ##(1)specify candidate patterns (possible knot sequences);
  ##(2)knot.sel: calculate AIC and BIC under each pattern 
  ##   by fitting the corresponding longitudinal model to the simulated data;
  ##(3)select the best one among the candidate patterns based on AIC and BIC.
  
  pattern_set=list()
  pattern_candi=c(pi/2,3*pi/2)
  pattern_set[[1]]=pattern_candi
  knot_1=knot.sel(data,pattern_candi,ctl=lmeControl (msMaxIter=100))

  pattern_candi=c(pi)
  pattern_set[[2]]=pattern_candi
  knot_2=knot.sel(data,pattern_candi,ctl=lmeControl (msMaxIter=100))

  pattern_candi=c(pi/2,pi,3*pi/2)
  pattern_set[[3]]=pattern_candi
  knot_3=knot.sel(data,pattern_candi,ctl=lmeControl (msMaxIter=100))

  pattern_candi=quantile(data$time,0.5)
  pattern_set[[4]]=pattern_candi
  knot_4=knot.sel(data,pattern_candi,ctl=lmeControl (msMaxIter=100))

  pattern_candi=quantile(data$time,0.75)
  pattern_set[[5]]=pattern_candi
  knot_5=knot.sel(data,pattern_candi,ctl=lmeControl (msMaxIter=100))

  pattern_candi=quantile(data$time,c(0.25,0.5,0.75))
  pattern_set[[6]]=pattern_candi
  knot_6=knot.sel(data,pattern_candi,ctl=lmeControl (msMaxIter=100))

  pattern_candi=c(3)
  pattern_set[[7]]=pattern_candi
  knot_7=knot.sel(data,pattern_candi,ctl=lmeControl (msMaxIter=100))

  knot_aic_set=c(knot_1$aic,knot_2$aic,knot_3$aic,knot_4$aic,knot_5$aic,knot_6$aic,knot_7$aic)
  knot_bic_set=c(knot_1$bic,knot_2$bic,knot_3$bic,knot_4$bic,knot_5$bic,knot_6$bic,knot_7$bic)
  if(which.min(knot_aic_set)==which.min(knot_bic_set)){
    internal_knots=pattern_set[[which.min(knot_aic_set)]]
  }else{
    return(NULL) # be careful when AIC and BIC produce different selection result.
  }
  
  
  q=length(internal_knots)+4
  M=as.data.frame(mycubicbs(data$time,internal_knots,boundary_knots=c(0,obsmax))$mat)
  names(M)=paste0("time",c(1:q))
  data=cbind(data,M)
  data.id=data[!duplicated(data$id),]

  
  knots=mycubicbs(0,internal_knots,boundary_knots)$knots
  Q.2=matrix(0,nrow=q,ncol=q-2)
  for(l in 1:(q-2)){
    Q.2[l,l]=6/((knots[l+4]-knots[l+2])*(knots[l+4]-knots[l+1]))
    Q.2[l+2,l]=6/((knots[l+4]-knots[l+2])*(knots[l+5]-knots[l+2]))
    Q.2[l+1,l]=-(Q.2[l,l]+Q.2[l+2,l])
  }
  
  ######Get initial values: estimation results from TS and â€œNo Variabilityâ€?.
  initialvalue=try(inival(M,data,data.id,ctl=lmeControl (msMaxIter=100),internal_knots,Q.2))
  if('try-error' %in% class(initialvalue)){
    return(NULL)
  }
  
  
  ##Estimate under the true model: fit a cox regression model with true covariates
  cox.sin=cox.nonpar(data,data.id,Data$sin_ran,Data$sin_int)
  print(cox.sin$res_par)
  
  ###Estimate by JM
  riskset=function(t){# individuals in the risk set
    unique(data$id[(data$Time>=t)])## "=" is important
  }
  
  cc=c()
  for(t in sort(unique(data$Time[data$delta==1]))){
    rs.id=riskset(t)
    s=sum(exp(0*data.id$W[rs.id]))

    cc=c(cc,length(data.id$id[(data.id$Time==t)&(data.id$delta==1)])/s)
  }


  cumbase=data.frame(hazard=0,time=c(0,sort(unique(data$Time[data$delta==1]))))
  cumbase[,1]=cumsum(c(0,cc))

  res.jm1=try(est(data,data.id,0,0,0,
              cumbase,initialvalue$beta,initialvalue$sigma2,initialvalue$D,internal_knots,Q.2))
 
  if('try-error' %in% class(res.jm1)){
    return(NULL)
  }
  
  if(is.null(res.jm1)){
    return(NULL)
  }
  
  res_jm.par=res.jm1$par
  print(res_jm.par)

  
  res_ts.par=c(initialvalue$coxts_coxph,initialvalue$coxts_my,initialvalue$beta,initialvalue$sigma2,initialvalue$coxts_coxph_novar)

  res_non.par=c(cox.sin$res_par$par)

  return(c(which.min(knot_aic_set),res_jm.par,res_ts.par,res_non.par))
}

RES=sfLapply(1:100,SIMULATE)
sfStop()



