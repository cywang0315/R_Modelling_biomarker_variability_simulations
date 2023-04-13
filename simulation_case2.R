### Simulation study: case 2.
### Program in a parallel way, all functions required are stored in "source_function_18.R".


library(snowfall)
library(parallel)

sfInit(parallel = TRUE, cpus = detectCores()-4)
sfLibrary(MCMCglmm) #rtnorm
sfLibrary(distr)# mixture distribution
sfLibrary(survival)
sfLibrary(JM)
sfLibrary(joineR)
sfLibrary(dplyr) # near
sfLibrary(statmod)
sfLibrary(progress)
sfLibrary(MASS)
sfLibrary(mvtnorm)
sfLibrary(tensor)
sfSource(here::here("source_function_18.R"))
####################

SIMULATE=function(s){
  knots=mycubicbs(0,internal_knots,boundary_knots)$knots
  Q.2=matrix(0,nrow=q,ncol=q-2)
  for(l in 1:(q-2)){
    Q.2[l,l]=6/((knots[l+4]-knots[l+2])*(knots[l+4]-knots[l+1]))
    Q.2[l+2,l]=6/((knots[l+4]-knots[l+2])*(knots[l+5]-knots[l+2]))
    Q.2[l+1,l]=-(Q.2[l,l]+Q.2[l+2,l])
  }
  
  ## Generate simulation data.
  data=gendat(s,obstim=c(0,0.041,0.082,0.166,0.25,0.5,0.75,seq(1,5.8,by=0.25)),obsmax=5.8,gammatrue=c(-0.4,-0.1,0.71,0.22,0.36),
              alpha1true=0.37,alpha2true=0.24,betatrue=c(6.2,4.5,5,4.7,5,4.5),beta1true=c(6.8,4.2,4.5,4.1,4.5,4.2),beta2true=c(6.4,4.5,4.6,4.1,4.5,4.3),
              etatrue=c(-0.05,0.03,0.09),xitrue=0.08,
              D0=diag(c(0.2,0.2,0.4,0.5,0.3)),sigmatrue=0.42,knots,Q.2)
  
  M=as.data.frame(mycubicbs(data$time,internal_knots,boundary_knots)$mat)
  names(M)=paste0("time",c(1:q))
  data=cbind(data,M)
  data.id=data[!duplicated(data$id),]
  
  ### Initial values of parameters
  initialvalue=inival(data,data.id,ctl=lmeControl (msMaxIter=100),knots,Q.2)

  gamma=initialvalue$gamma 
  beta=initialvalue$beta
  beta1=initialvalue$beta1
  beta2=initialvalue$beta2
  eta1=initialvalue$eta1
  eta2=initialvalue$eta2
  sigma2=initialvalue$sigma2
  D=initialvalue$D
  alpha1=initialvalue$alpha1
  alpha2=initialvalue$alpha2
  cumbase=initialvalue$cumbase
  
  ## Result from the two-stage method
  res.ts=c(gamma,alpha1,alpha2,beta,beta1,beta2,eta1,eta2,sigma2,diag(D))
  
  # Evaluate the log-likelihood function
  logLikmc.or=logLik(data,data.id,gamma,alpha1,alpha2,beta,beta1,beta2,eta1,eta2,sigma2,D,cumbase,knots,Q.2)
  
  ## Estimate by the joint modelling method
  res.jm=est(data,data.id,gamma,alpha1,alpha2,beta,beta1,beta2,eta1,eta2,sigma2,D,cumbase,knots,Q.2)
  return(c(c(res.jm),c(res.ts)))
}

RES=sfLapply(1:3,SIMULATE)

sfStop()

