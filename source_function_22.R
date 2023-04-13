
m=1000
obsmax=6
boundary_knots=c(0,obsmax)
####### B-spline basis #########
## my cubic B-spline basis functions ####
mycubicbs=function(s,internal_knots,boundary_knots){
  df=length(internal_knots)+3+1
  r=numeric(df)
  out_knots=c(seq(from=min(boundary_knots),by=-0.1,length=4),seq(from=max(boundary_knots),by=0.1,length=4))
  knots=sort(c(internal_knots,out_knots))
  temp=function(x){
    my0bs=function(x,k){
      a=ifelse(x<knots[k+1]&x>=knots[k],1,0)
      return(a)
    }
    my1bs=function(x,k){
      (x-knots[k])/(knots[k+1]-knots[k])*my0bs(x,k)+
        (knots[k+2]-x)/(knots[k+2]-knots[k+1])*my0bs(x,k+1)
    }
    
    my2bs=function(x,k){
      (x-knots[k])/(knots[k+2]-knots[k])*my1bs(x,k)+
        (knots[k+3]-x)/(knots[k+3]-knots[k+1])*my1bs(x,k+1)
    }
    my3bs=function(x,k){
      (x-knots[k])/(knots[k+3]-knots[k])*my2bs(x,k)+
        (knots[k+4]-x)/(knots[k+4]-knots[k+1])*my2bs(x,k+1)
    }
    for(k in 1:df){
      r[k]=my3bs(x,k)
    }
    return(r)}
  return=list(mat=t(sapply(s,temp)),knots=knots)# a matrix with ncol=df
}

R.2=function(t,knots){
  q=length(knots)-4
  R.2=matrix(0,ncol=q-2,nrow=q-2)
  for (l in 1:(q-2)){
    if((t>knots[l+2])&(t<=knots[l+3])){
      R.2[l,l]=(t-knots[l+2])^3/(3*(knots[l+3]-knots[l+2])^2)
    }else{
      if((t>knots[l+3])&(t<knots[l+4])){
        R.2[l,l]=(knots[l+4]-knots[l+2])/3-(knots[l+4]-t)^3/(3*(knots[l+4]-knots[l+3])^2)
        R.2[l,l+1]=(-t^3/3+(knots[l+4]+knots[l+3])/2*t^2-knots[l+4]*knots[l+3]*t-(knots[l+3])^3/6+
                      (knots[l+3])^2*knots[l+4]/2)/(knots[l+4]-knots[l+3])^2
        R.2[l+1,l]= R.2[l,l+1]
      }else{
        if(t>=knots[l+4]){
          R.2[l,l]=(knots[l+4]-knots[l+2])/3
          R.2[l,l+1]=(knots[l+4]-knots[l+3])/6
          R.2[l+1,l]=(knots[l+4]-knots[l+3])/6
          
        }
      }
    }
  }
  R.2[1,1]=ifelse(t>=knots[5],(knots[5]-knots[4])/3,(knots[5]-knots[4])/3-(knots[5]-t)^3/(3*(knots[5]-knots[4])^2))
  return(R.2)
}

#The true trajectory behind longitudinal observations is a polynomial function
## with random coefficients (ran_par and ran_int).
longit.true=function(t,beta_cub,beta_slo,ran_int,ran_par){
  ran_par*(beta_cub*(t-3)^3+beta_slo*(t-3))+ran_int
}




adv=function(x){
  grad=diag(1,length(x))
  grad2=matrix(0,ncol=length(x),nrow = length(x))
  attr(x,"grad")=grad
  attr(x,"grad2")=grad2
  class(x)="adv"
  x
}


quad.adv=function(A,b){
  grad.b=attr(b,"grad")
  grad2.b=attr(b,"grad2")
  b=as.numeric(b)
  d=c(t(b)%*%A%*%b) # scalar
  attr(d,"grad")=grad.b%*%(2*A%*%b)
  attr(d,"grad2")=grad.b%*%(2*A)
  class(d)="adv"
  d
}


"^.adv"=function(b,alpha){ # the input is a scalar function of b
  grad.b=attr(b,"grad")
  grad2.b=attr(b,"grad2")
  b=as.numeric(b)
  d=b^alpha
  attr(d,"grad")=alpha*b^(alpha-1)*grad.b
  attr(d,"grad2")=alpha*(alpha-1)*b^(alpha-2)*tcrossprod(grad.b)+alpha*b^(alpha-1)*grad2.b
  class(d)="adv"
  d
}


linear_sca.adv=function(A,b,j){ ### derivative of <Ab>j with respect to b, A is matrix independent of b
  grad.b=attr(b,"grad")
  b=as.numeric(b) ## transform to a vector
  d=c(A%*%b)[j]
  K=numeric(nrow(A))
  K[j]=1
  attr(d,"grad")=t(A)%*%K
  attr(d,"grad2")=matrix(0,ncol=length(b),nrow=length(b))
  class(d)="adv"
  d
}

linear.adv=function(A,b){## derivative of Ab w.r.t b, A is a matrix or a column vector
  grad.b=attr(b,"grad")
  b=as.numeric(b)
  check=length(A)==length(b)
  if(check){
    d=t(A)%*%b
    attr(d,"grad")=A
    attr(d,"grad2")=matrix(0,ncol=length(b),nrow = length(b))
  }else{
    d=A%*%b
    attr(d,"grad")=A
  }
  
  class(d)="adv"
  d
}

exp.adv=function(b){ ## the input is a scalar function of b
  grad.b=attr(b,"grad")
  grad2.b=attr(b,"grad2")
  b=as.numeric(b)
  d=exp(b)
  attr(d,"grad")=exp(b)*grad.b
  attr(d,"grad2")=exp(b)*tcrossprod(grad.b)+exp(b)*grad2.b
  class(d)="adv"
  d
}



"*.adv"=function(a1,a2){#derivative of a1(b)*a2(b) w.r.t. b, where a1 is a scalar function of b and a2 is a scalar or vector function of b.
  
  if(is.null(attr(a1,"grad"))){
    grad.a2=attr(a2,"grad")
    grad2.a2=attr(a2,"grad2")
    a2=as.numeric(a2)
    d=a1*a2
    attr(d,"grad")=a1*grad.a2
    attr(d,"grad2")=a1*grad2.a2
    
  }else{
    grad.a1=attr(a1,"grad")
    grad2.a1=attr(a1,"grad2")
    grad.a2=attr(a2,"grad")
    grad2.a2=attr(a2,"grad2")
    a1=as.numeric(a1)
    a2=as.numeric(a2)
    d=a1*a2
    check=length(a2)==1
    if(check){
      attr(d,"grad")=grad.a1*a2+a1*grad.a2
      attr(d,"grad2")=grad.a1%*%t(grad.a2)+grad.a2%*%t(grad.a1)+a2*grad2.a1+a1*grad2.a2
    }else{
      attr(d,"grad")=a2%*%t(grad.a1)+a1*grad.a2
    }
  }
  
  class(d)="adv"
  d
}


"+.adv"=function(a,b){ #derivative of a(b)+b(b) where a and b are vector or scalar functions of b, and a(.) can be not function of b
  if(is.null(attr(a,"grad"))){
    grad.b=attr(b,"grad")
    grad2.b=attr(b,"grad2")
    b=as.numeric(b)
    d=a+b
    attr(d,"grad")=grad.b
    attr(d,"grad2")=grad2.b
    
  }else{
    grad.a=attr(a,"grad")
    grad2.a=attr(a,"grad2")
    
    grad.b=attr(b,"grad")
    grad2.b=attr(b,"grad2")
    
    a=as.numeric(a)
    b=as.numeric(b)
    d=a+b
    attr(d,"grad")=grad.a+grad.b
    attr(d,"grad2")=grad2.a+grad2.b
    
  }
  
  class(d)="adv"
  d
}

##################Generate data 
gendat=function(s,obsmax,gammatrue,alpha1true,alpha2true,beta_cub,beta_slo,sigmatrue){
  set.seed(s)
  basetrue=function(t){ifelse(t<2,0,exp(-3))} ## 0 when t<2
  #b=mvrnorm(m,mu=rep(0,q),Sigma=D0) # random effects in longitudinal submodel
  W=sample(c(0,1),size=m,replace = TRUE)
  ran_par=runif(m,0.5,2)
  ran_int=runif(m,4,8)
  #plot(x=seq(0,6,by=0.02),y=2*sin(seq(0,6,by=0.02))+4)
  #plot(x=seq(0,6,by=0.02),y=2^2*(seq(0,6,by=0.02)/2-sin(2*seq(0,6,by=0.02))/4))

  TD.fun=function(i,t){# i is a scale, t can be a vector
    tt=function(t){ #ran_par*(beta_cub*t^3+beta_slo*t)+beta_int
      c(basetrue(t)*exp(gammatrue*W[i]+alpha1true*(ran_par[i]*(beta_cub*(t-3)^3+beta_slo*(t-3))+ran_int[i])+
                          alpha2true*sqrt(12*beta_cub^2*ran_par[i]^2*((t-3)^3+27))))
    }
    #hazard=sapply(t,tt)
    #p=1-exp(-hazard*0.01)
    r=rbinom(n=length(t),size=1,prob=1-exp(-sapply(t,tt)*0.002)) ##prob: event occurs
    Time1=obsmax
    delta1=0
    if(max(r)>0){
      Time1=t[min(which(r==1))]
      delta1=1
    }
    return(list(Time1=Time1,delta1=delta1))
  }
  
  TD=sapply(c(1:m),TD.fun,t=seq(2,obsmax,by=0.002))
  Time1=as.numeric(TD[1,])
  delta1=as.numeric(TD[2,])
  sum(delta1)
  hist(Time1)
  ####generate censoring 
  censtim=sample(seq(3,10,by=0.001),size=m,replace = TRUE)
  delta=ifelse(Time1<=censtim,delta1,0)
  sum(delta)
  Time=ifelse(Time1<=censtim,Time1,censtim)
  hist(Time)
  #############generate longitudinal data
  Y=c()
  l=numeric(m)
  time=c()
  for (i in 1:m){
    l[i]=length(c(seq(0,1.6,by=0.4),seq(2,Time[i],by=0.5)))
    for(j in c(seq(0,1.6,by=0.4),seq(2,Time[i],by=0.5))){
      YY=longit.true(j,beta_cub,beta_slo,ran_int[i],ran_par[i])+rnorm(1,sd=sigmatrue)
      Y=rbind(Y,YY)
      time=rbind(time,j)
    }
  }
  
  Y=as.vector(Y)
  time=as.vector(time)
  ### construct data frame. Note here we use "id" not "sub"
  data=data.frame(id=rep(c(1:m),l),Time=rep(Time,l),delta=rep(delta,l),Y=Y,W=rep(W,l),
                  time=time, start=time,stop=time+1/2,event=numeric(length(Y))) ##should adjust "stop" and "event"
  
  data$stop=ifelse(data$stop<=data$Time,data$stop,data$Time)
  data$event=ifelse((data$stop==data$Time)&(data$delta==1),1,data$event)
  return(list(data=data,ran_int=ran_int,ran_slo=ran_par))
}

## Knot selection
knot.sel=function(data,pattern_candi,ctl=lmeControl (msMaxIter=100)){
  internal_knots=pattern_candi
  M=as.data.frame(mycubicbs(data$time,internal_knots,boundary_knots)$mat)
  q=length(internal_knots)+4
  names(M)=paste0("time",c(1:q))
  data.M=cbind(data,M)
  fix_for=as.formula(paste("Y ~ ", "0","+",paste(names(M), collapse= "+")))
  ran_for=as.formula(paste("~ ", "0","+",paste(names(M), collapse= "+")))
  lme.data=lme(fixed=fix_for,random=list(id=pdDiag(form= ran_for)),data = data.M,control=ctl)
  
  return(list(aic=AIC(lme.data),bic=BIC(lme.data)))
  
}

############################# Get initial value ##################
inival=function(M,data,data.id,ctl,internal_knots,Q.2){
  knots=mycubicbs(0,internal_knots,boundary_knots)$knots
  q=dim(M)[2]
  fix_for=as.formula(paste("Y ~ ", "0","+",paste(names(M), collapse= "+")))
  ran_for=as.formula(paste("~ ", "0","+",paste(names(M), collapse= "+")))
  lme.data=lme(fixed=fix_for,random=list(id=pdDiag(form= ran_for)),data = data,control=ctl)

  beta=lme.data$coefficients[[1]]
  sigma2=lme.data$sigma^2
  D=diag(c(as.numeric(VarCorr(lme.data)[1:q])))
  
  #######################
  #tt=seq(0,6,by=0.02)
  #fit_tra=numeric(length(tt))
  
  # for(i in 1:length(tt)){
  #   fit_tra[i]=c(mycubicbs(tt[i],internal_knots,boundary_knots=c(0,obsmax))$mat%*%as.numeric(beta+random.effects(lme.data)[100,]))
  # }
  # 
  # points(x=tt,y=fit_tra,col=4)
  # plot(x=tt,y=fit_tra)
  
  #### two stage, using coxph function ####
  data$stop[data$start==data$stop]=data$stop[data$start==data$stop]+0.01
  data.ts=data
  l=as.data.frame(table(data$id))$Freq
  raneff=apply(random.effects(lme.data),2,rep,times=l)
  data.ts$preY=rowSums(data[,10:(10+q-1)]*t(beta+t(raneff)))
  data.ts$prevar=sapply(c(1:(dim(data)[1])),function(x) sqrt(c(t(as.numeric(beta+raneff[x,]))%*%(Q.2%*%R.2(data$time[x],knots)%*%t(Q.2))%*%(as.numeric(beta+raneff[x,])))))
  coxts.data=coxph(Surv(start, stop, event)~W+preY+prevar,data=data.ts)
  coxts_coxph=coxts.data$coefficients
  cumbase_coxph=rbind(c(0,0),basehaz(coxts.data, centered=FALSE))
  
  #######No variability
  coxts.data_novar=coxph(Surv(start, stop, event)~W+preY,data=data.ts)
  coxts_coxph_novar=coxts.data_novar$coefficients
  #### two stage, self-coding, exact two-stage ####
  conloglik_3=function(par){
    gamma=par[1]
    alpha1=par[2]
    alpha2=par[3]
    conloglik=0
    for(i in data.id$id[data.id$delta==1]){
      Ti=data.id$Time[data.id$id==i]
      rs.id=data.id$id[data.id$Time>=Ti]
      conloglik=conloglik+gamma*data.id$W[i]+alpha1*sum(c(mycubicbs(Ti,internal_knots,boundary_knots)$mat)*(beta+t(random.effects(lme.data)[i,])))+
        alpha2*sqrt(sum((beta+t(random.effects(lme.data)[i,]))*(Q.2%*%R.2(Ti,knots)%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[i,])))))-
        log(sum(exp(gamma*data.id$W[rs.id]+alpha1*colSums(c(mycubicbs(Ti,internal_knots,boundary_knots)$mat)*(beta+t(random.effects(lme.data)[rs.id,])))+
                      alpha2*sqrt(colSums((beta+t(random.effects(lme.data)[rs.id,]))*(Q.2%*%R.2(Ti,knots)%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[rs.id,]))))))))
    }
    return(-conloglik)
    
  }
  
  #optim(c(0,0,0),conloglik_3,method = "BFGS")
  
  
  gg_conloglik_3=function(par){
    gamma=par[1]
    alpha1=par[2]
    alpha2=par[3]
    gg_gamma=0
    gg_alpha1=0
    gg_alpha2=0
    i_alpha1=0
    i_alpha2=0
    for(i in data.id$id[data.id$delta==1]){
      Ti=data.id$Time[data.id$id==i]
      rs.id=data.id$id[data.id$Time>=Ti]
      s=sum(exp(gamma*data.id$W[rs.id]+alpha1*colSums(c(mycubicbs(Ti,internal_knots,boundary_knots)$mat)*(beta+t(random.effects(lme.data)[rs.id,])))+
                  alpha2*sqrt(colSums((beta+t(random.effects(lme.data)[rs.id,]))*(Q.2%*%R.2(Ti,knots)%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[rs.id,])))))))
      s1=sum(exp(gamma*data.id$W[rs.id]+alpha1*colSums(c(mycubicbs(Ti,internal_knots,boundary_knots)$mat)*(beta+t(random.effects(lme.data)[rs.id,])))+
                   alpha2*sqrt(colSums((beta+t(random.effects(lme.data)[rs.id,]))*(Q.2%*%R.2(Ti,knots)%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[rs.id,]))))))*data.id$W[rs.id])
      
      s2=sum(exp(gamma*data.id$W[rs.id]+alpha1*colSums(c(mycubicbs(Ti,internal_knots,boundary_knots)$mat)*(beta+t(random.effects(lme.data)[rs.id,])))+
                   alpha2*sqrt(colSums((beta+t(random.effects(lme.data)[rs.id,]))*(Q.2%*%R.2(Ti,knots)%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[rs.id,]))))))*colSums(c(mycubicbs(Ti,internal_knots,boundary_knots)$mat)*(beta+t(random.effects(lme.data)[rs.id,]))))
      
      i2=sum(exp(gamma*data.id$W[rs.id]+alpha1*colSums(c(mycubicbs(Ti,internal_knots,boundary_knots)$mat)*(beta+t(random.effects(lme.data)[rs.id,])))+
                   alpha2*sqrt(colSums((beta+t(random.effects(lme.data)[rs.id,]))*(Q.2%*%R.2(Ti,knots)%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[rs.id,]))))))*
               (colSums(c(mycubicbs(Ti,internal_knots,boundary_knots)$mat)*(beta+t(random.effects(lme.data)[rs.id,]))))^2)
      
      s3=sum(exp(gamma*data.id$W[rs.id]+alpha1*colSums(c(mycubicbs(Ti,internal_knots,boundary_knots)$mat)*(beta+t(random.effects(lme.data)[rs.id,])))+
                   alpha2*sqrt(colSums((beta+t(random.effects(lme.data)[rs.id,]))*(Q.2%*%R.2(Ti,knots)%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[rs.id,]))))))*sqrt(colSums((beta+t(random.effects(lme.data)[rs.id,]))*(Q.2%*%R.2(Ti,knots)%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[rs.id,]))))))
      
      i3=sum(exp(gamma*data.id$W[rs.id]+alpha1*colSums(c(mycubicbs(Ti,internal_knots,boundary_knots)$mat)*(beta+t(random.effects(lme.data)[rs.id,])))+
                   alpha2*sqrt(colSums((beta+t(random.effects(lme.data)[rs.id,]))*(Q.2%*%R.2(Ti,knots)%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[rs.id,]))))))*
               (colSums((beta+t(random.effects(lme.data)[rs.id,]))*(Q.2%*%R.2(Ti,knots)%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[rs.id,]))))))
      
      gg_gamma=gg_gamma+data.id$W[i]-s1/s
      gg_alpha1=gg_alpha1+sum(c(mycubicbs(Ti,internal_knots,boundary_knots)$mat)*(beta+t(random.effects(lme.data)[i,])))-s2/s
      i_alpha1=i_alpha1+i2/s-(s2/s)^2
      gg_alpha2=gg_alpha2+sqrt(sum((beta+t(random.effects(lme.data)[i,]))*(Q.2%*%R.2(Ti,knots)%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[i,])))))-s3/s
      i_alpha2=i_alpha2+i3/s-(s3/s)^2
      
    }
    return(c(-gg_gamma,-gg_alpha1,-gg_alpha2))
    #return(c(-gg_alpha1,-gg_alpha2))
    
  }
  #
  res=optim(c(0,0,0),conloglik_3,gg_conloglik_3,method = "BFGS")
  coxts_my=res$par
  
  gamma=res$par[1]
  alpha1=res$par[2]
  alpha2=res$par[3]
  #
  # ### CALCULATE CUMULATIVE HAZARD
  riskset=function(t){# individuals in the risk set
    unique(data$id[(data$Time>=t)])## "=" is important
  }
  
  cc=c()
  for(t in sort(unique(data$Time[data$delta==1]))){
    rs.id=riskset(t)
    s=sum(exp(gamma*data.id$W[rs.id]+alpha1*colSums(c(mycubicbs(t,internal_knots,boundary_knots)$mat)*(beta+t(random.effects(lme.data)[rs.id,])))+
                alpha2*sqrt(colSums((beta+t(random.effects(lme.data)[rs.id,]))*(Q.2%*%R.2(t,knots)%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[rs.id,])))))))
    
    cc=c(cc,length(data.id$id[(data.id$Time==t)&(data.id$delta==1)])/s)
  }
  
  cumbase_my=data.frame(hazard=0,time=c(0,sort(unique(data$Time[data$delta==1]))))
  cumbase_my[,1]=cumsum(c(0,cc))
  
  out=list(beta=beta,sigma2=sigma2,D=D,coxts_coxph=coxts_coxph,cumbase_coxph=cumbase_coxph,coxts_coxph_novar=coxts_coxph_novar,coxts_my=coxts_my,cumbase_my=cumbase_my,b=random.effects(lme.data))
  
}


### Estimate parameters in Cox regression model if we know the true trajectory and its second derivative
cox.nonpar=function(data,data.id,beta_cub, beta_slo,ran_int,ran_par){
  
  conloglik_nonpar=function(par){
    gamma=par[1]
    alpha1=par[2]
    alpha2=par[3]
    conloglik=0
    for(i in data.id$id[data.id$delta==1]){
      Ti=data.id$Time[data.id$id==i]
      rs.id=data.id$id[data.id$Time>=Ti]
      conloglik=conloglik+gamma*data.id$W[i]+alpha1*(ran_par[i]*(beta_cub*(Ti-3)^3+beta_slo*(Ti-3))+ran_int[i])+
        alpha2*sqrt(12*beta_cub^2*ran_par[i]^2*((Ti-3)^3+27))-
        log(sum(exp(gamma*data.id$W[rs.id]+alpha1*(ran_par[rs.id]*(beta_cub*(Ti-3)^3+beta_slo*(Ti-3))+ran_int[rs.id])+
                      alpha2*sqrt(12*beta_cub^2*ran_par[rs.id]^2*((Ti-3)^3+27)))))
    }
    return(-conloglik)
  }
  
  res=optim(c(0,0,0),conloglik_nonpar,method = "BFGS")
  
  gamma=res$par[1]
  alpha1=res$par[2]
  alpha2=res$par[3]
  
  riskset=function(t){# individuals in the risk set
    unique(data$id[(data$Time>=t)])## "=" is important
  }
  
  cc=c()
  for(t in sort(unique(data$Time[data$delta==1]))){
    rs.id=riskset(t)
    s=sum(exp(gamma*data.id$W[rs.id]+alpha1*(ran_par[rs.id]*(beta_cub*(t-3)^3+beta_slo*(t-3))+ran_int[rs.id])+
                alpha2*sqrt(12*beta_cub^2*ran_par[rs.id]^2*((t-3)^3+27))))
    
    cc=c(cc,length(data.id$id[(data.id$Time==t)&(data.id$delta==1)])/s)
  }
  
  
  cumbase=data.frame(hazard=0,time=c(0,sort(unique(data$Time[data$delta==1]))))
  cumbase[,1]=cumsum(c(0,cc))
  
  
  return(list(res_par=res,cumbase=cumbase))
  
}


##############Likelihood
logLik=function(data,data.id,gamma,alpha1,alpha2,cumbase,beta,sigma2,D,internal_knots,Q.2,L){
  
  lambda0=function(t){
    if(t %in% unique(data$Time[data$delta==1])){
      x=cumbase[which(near(cumbase[,2],t)),1]-cumbase[which(near(cumbase[,2],t))-1,1] ##warn: options(digits = 22) default 7
      if(length(x)>1){
        x=x[1]
      }
      
    }else{
      x=0}
    return(x)
  }
  
  q=length(internal_knots)+4
  knots=mycubicbs(0,internal_knots,boundary_knots)$knots
  Time=data$Time[!duplicated(data$id)]
  delta=data$delta[!duplicated(data$id)]
  W=data$W[!duplicated(data$id)]
  des.Y=as.matrix(data[,c(paste0("time",c(1:q)))])
  des.T=mycubicbs(data.id$Time[data.id$delta==1],internal_knots,boundary_knots)$mat
  
  
  ### collection of K(t)=Q.2%*%R.2(t)%*%t(Q.2) for t=data.id$Time[data.id$delta==1]
  K.2=array(c(sapply(data.id$Time[data.id$delta==1],function(t) Q.2%*%R.2(t,knots)%*%t(Q.2))),
            dim=c(q,q,length(data.id$Time[data.id$delta==1]))) 
  
  log.p.Y=numeric(m)
  #log.p.Y2=numeric(m)
  log.s_set=matrix(0,nrow=m,ncol=L)
  ZTimeb=matrix(0,nrow=sum(delta==1),ncol=L) ##Z_Time%*%b
  bKb=matrix(0,nrow=sum(delta==1),ncol=L)
  rmc=rmvnorm(L,mean=rep(0,q))
  for(i in 1:m){
    Zi=des.Y[data$id==i,]
    Xi=Zi
    Yi=data$Y[data$id==i]
    Sigma=solve(solve(D)+c(sigma2^(-1))*t(Zi)%*%Zi)
    mu=c(Sigma%*%(t(Zi)/c(sigma2))%*%(Yi-Xi%*%beta))
    log.p.Y[i]=dmvnorm(Yi,mean=c(Xi%*%beta),sigma=sigma2*diag(length(Yi))+Zi%*%D%*%t(Zi),log=TRUE)
    #log.p.Y2[i]=dmvnorm(Yi,mean=c(Xi%*%betanew),sigma=sigma2new*diag(length(Yi))+Zi%*%diag(diag(Dnew))%*%t(Zi),log=TRUE)
    
    rmci=t(mu+t(chol(Sigma))%*%t(rmc))## matrix with ncol=q,nrow=L, MC sample from bi|yi
    
    if(delta[i]==1){
      ZTimeb[match(i,which(delta==1)),]=des.T[match(i,which(delta==1)),]%*%t(rmci)
      bKb[match(i,which(delta==1)),]=colSums((t(rmci)+beta)*(K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(t(rmci)+beta)))
    }
    
    a=unique(data.id$Time[(data.id$Time<=Time[i])&(data.id$delta==1)])
    if(length(a)>0){
      if(length(a)==1){
        btKb=colSums((t(rmci)+beta)*(K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(t(rmci)+beta)))
      }else{
        btKb=apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) colSums((t(rmci)+beta)*(y%*%(t(rmci)+beta))))
      }
      log.s_set[i,]=apply(sapply(a,lambda0)*exp(gamma*W[i])*
                            exp(alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(t(rmci)+beta)+
                                  alpha2*t(sqrt(btKb))),2,sum)
      
    }
    
  }
  
  
  log.hazard=matrix(0,ncol=L,nrow=m)
  log.hazard[which(delta==1),]=log(sapply(Time[delta==1],lambda0))+gamma*W[delta==1]+alpha1*(c(des.T%*%beta)+ZTimeb)+
    alpha2*sqrt(bKb)
  log.p.Tb=log.hazard-log.s_set
  p.Tb=exp(log.p.Tb)
  logLikmc=sum(log.p.Y+log(rowMeans(p.Tb)))
  return(logLikmc)
}


###############Estimation##########

est=function(data,data.id,gamma,alpha1,alpha2,cumbase,beta,sigma2,D,internal_knots,Q.2){
  q=length(internal_knots)+4
  knots=mycubicbs(0,internal_knots,boundary_knots)$knots
  Time=data$Time[!duplicated(data$id)]
  delta=data$delta[!duplicated(data$id)]
  W=data$W[!duplicated(data$id)]
  des.Y=as.matrix(data[,c(paste0("time",c(1:q)))])
  des.T=mycubicbs(data.id$Time[data.id$delta==1],internal_knots,boundary_knots)$mat
  N=nrow(des.Y)
  l=as.data.frame(table(data$id))$Freq
  diXZ=matrix(0,ncol=q*m,nrow=N) ##block diagnal version of XZ
  for(i in 1:m){
    diXZ[c((cumsum(l)[i]-l[i]+1):(cumsum(l)[i])),c((q*(i-1)+1):(q*i))]=des.Y[c((cumsum(l)[i]-l[i]+1):(cumsum(l)[i])),]
  }
  ### collection of K(t)=Q.2%*%R.2(t)%*%t(Q.2) for t=data.id$Time[data.id$delta==1]
  K.2=array(c(sapply(data.id$Time[data.id$delta==1],function(t) Q.2%*%R.2(t,knots)%*%t(Q.2))),
            dim=c(q,q,length(data.id$Time[data.id$delta==1]))) 
  ### EM algorithm  ###
  K=10
  gamma_set=numeric(K)
  alpha1_set=numeric(K)
  alpha2_set=numeric(K)
  beta_set=matrix(0,ncol=q,nrow=K)
  sigma2_set=numeric(K)
  D_set=array(0,dim=c(q,q,K))
  Q.fun_set=numeric(K)
  logLik_set=numeric(K)
  b_set=matrix(0,ncol=q,nrow=m)
  inv.Fish_set=array(0,dim=c(q,q,m)) ## V(bi)~=inv.Fishi
  ## Estimation
  lambda0=function(t){
    if(t %in% unique(data$Time[data$delta==1])){
      x=cumbase[which(near(cumbase[,2],t)),1]-cumbase[which(near(cumbase[,2],t))-1,1] ##warn: options(digits = 22) default 7
      if(length(x)>1){
        x=x[1]
      }
      
    }else{
      x=0}
    return(x)
  }
  
  riskset=function(t){# individuals in the risk set
    unique(data$id[(data$Time>=t)])## "=" is important
  }
  
  for (k in 1:K){
    
    for(i in 1:m){
      Zi=des.Y[data$id==i,]
      Xi=Zi
      Yi=data$Y[data$id==i]
      chol_invSig=chol(solve(D)+c(sigma2^(-1))*t(Zi)%*%Zi)
      Q=solve(chol_invSig)
      bi=c(Q%*%t(Q)%*%(t(Zi)/c(sigma2))%*%(Yi-Xi%*%beta)) #E(bi|Yi,theta),initial value of bi
      a=unique(data.id$Time[(data.id$Time<=Time[i])&(data.id$delta==1)])
      if(length(a)>0){
        if(length(a)>1){
          ### initial value of p(Ti,deltai,Yi,bi;theta)
          log.p.b=dmvnorm(bi,mean=rep(0,q),sigma=D,log=TRUE)
          log.p.Yb=sum(dnorm(Yi,mean=c(Xi%*%beta+Zi%*%bi),sd=sqrt(sigma2),log=TRUE)) # p(Y|b)
          log.s=exp(gamma*W[i])*sum(sapply(a,lambda0)*exp(alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(beta+bi)+
                                                            alpha2*apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) sqrt(t(beta+bi)%*%y%*%(beta+bi)))))
          
          log.hazard=ifelse(delta[i]==0,0,log(lambda0(Time[i]))+gamma*W[i]+alpha1*des.T[match(Time[i],data.id$Time[data.id$delta==1]),]%*%(beta+bi)+
                              alpha2*sqrt(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi)))
          
          log.p.YTb=log.hazard-log.s+log.p.Yb+log.p.b
          
          for(kk in 1:10){
            wei=c(exp(gamma*W[i])*sapply(a,lambda0)*exp(alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(beta+bi)+
                                                          alpha2*apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) sqrt(t(beta+bi)%*%y%*%(beta+bi)))))
            ZK=alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]+alpha2*t(apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) (c(t(beta+bi)%*%y%*%(beta+bi)))^(-1/2)*y%*%(beta+bi)))
            Sbi=-colSums(wei*ZK)+t(Zi)%*%(Yi-Xi%*%beta-Zi%*%bi)/sigma2-solve(D)%*%bi
            
            Fishi=tensor(array(apply(ZK,1,tcrossprod),dim=c(q,q,length(a)))+
                           alpha2*array(apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) (c(t(beta+bi)%*%y%*%(beta+bi)))^(-1/2)*y),dim=c(q,q,length(a)))-
                           alpha2*array(apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) (c(t(beta+bi)%*%y%*%(beta+bi)))^(-3/2)*tcrossprod(y%*%(beta+bi))),dim=c(q,q,length(a))),
                         wei,3,1)+t(Zi)%*%Zi/sigma2+solve(D)
            if(delta[i]==1){
              Sbi=alpha1*des.T[match(Time[i],data.id$Time[data.id$delta==1]),]+alpha2*c(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-1/2)*K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi)+Sbi
              Fishi=-alpha2*(-c(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-3/2)*tcrossprod(K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))+
                               c(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-1/2)*K.2[,,match(Time[i],data.id$Time[data.id$delta==1])])+Fishi
            }
            ## update bi using Newton-Raphson: bi=bi+solve(I)S
            binew=c(bi+solve(Fishi)%*%Sbi)
            if(all(abs((binew-bi)/bi)[which(is.nan(abs((binew-bi)/bi))!=TRUE)]<10^(-3))){
              break
            }
            ## Check
            log.p.b.new=dmvnorm(binew,mean=rep(0,q),sigma=D,log=TRUE)
            log.p.Yb.new=sum(dnorm(Yi,mean=c(Xi%*%beta+Zi%*%binew),sd=sqrt(sigma2),log=TRUE)) # p(Y|b)
            log.s.new=exp(gamma*W[i])*sum(sapply(a,lambda0)*exp(alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(beta+binew)+
                                                                  alpha2*apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) sqrt(t(beta+binew)%*%y%*%(beta+binew)))))
            log.hazard.new=ifelse(delta[i]==0,0,log(lambda0(Time[i]))+gamma*W[i]+alpha1*des.T[match(Time[i],data.id$Time[data.id$delta==1]),]%*%(beta+binew)+
                                    alpha2*sqrt(t(beta+binew)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+binew)))
            
            log.p.YTb.new=log.hazard.new-log.s.new+log.p.Yb.new+log.p.b.new
            if(log.p.YTb.new<=log.p.YTb){
              #print("oh")
              break
            }
            bi=binew
            log.p.YTb=log.p.YTb.new
            #print(log.p.YTb)
            
          }
        }else{
          ### initial value of p(Ti,deltai,Yi,bi;theta)
          log.p.b=dmvnorm(bi,mean=rep(0,q),sigma=D,log=TRUE)
          log.p.Yb=sum(dnorm(Yi,mean=c(Xi%*%beta+Zi%*%bi),sd=sqrt(sigma2),log=TRUE)) # p(Y|b)
          
          log.s=exp(gamma*W[i])*sum(sapply(a,lambda0)*exp(alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(beta+bi)+
                                                            alpha2*sqrt(t(beta+bi)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi))))
          
          log.hazard=ifelse(delta[i]==0,0,log(lambda0(Time[i]))+gamma*W[i]+alpha1*des.T[match(Time[i],data.id$Time[data.id$delta==1]),]%*%(beta+bi)+
                              alpha2*sqrt(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi)))
          
          log.p.YTb=log.hazard-log.s+log.p.Yb+log.p.b
          
          for(kk in 1:10){
            wei=c(exp(gamma*W[i])*sapply(a,lambda0)*exp(alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(beta+bi)+
                                                          alpha2*sqrt(t(beta+bi)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi))))
            ZK=alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]+alpha2*c(t(beta+bi)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-1/2)*K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi)
            Sbi=-colSums(wei*ZK)+t(Zi)%*%(Yi-Xi%*%beta-Zi%*%bi)/sigma2-solve(D)%*%bi
            
            Fishi=wei*(tcrossprod(ZK)-alpha2*c(t(beta+bi)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-3/2)*tcrossprod(K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi))+
                         alpha2*c(t(beta+bi)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-1/2)*K.2[,,match(a,data.id$Time[data.id$delta==1])])+t(Zi)%*%Zi/sigma2+solve(D)
            if(delta[i]==1){
              Sbi=alpha1*des.T[match(Time[i],data.id$Time[data.id$delta==1]),]+alpha2*c(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-1/2)*K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi)+Sbi
              Fishi=-alpha2*(-c(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-3/2)*tcrossprod(K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))+
                               c(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-1/2)*K.2[,,match(Time[i],data.id$Time[data.id$delta==1])])+Fishi
            }
            ## update bi using Newton-Raphson: bi=bi+solve(I)S
            binew=c(bi+solve(Fishi)%*%Sbi)
            if(all(abs((binew-bi)/bi)[which(is.nan(abs((binew-bi)/bi))!=TRUE)]<10^(-3))){
              break
            }
            ## Check
            log.p.b.new=dmvnorm(binew,mean=rep(0,q),sigma=D,log=TRUE)
            log.p.Yb.new=sum(dnorm(Yi,mean=c(Xi%*%beta+Zi%*%binew),sd=sqrt(sigma2),log=TRUE)) # p(Y|b)
            log.s.new=exp(gamma*W[i])*sum(sapply(a,lambda0)*exp(alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(beta+binew)+
                                                                  alpha2*sqrt(t(beta+binew)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+binew))))
            log.hazard.new=ifelse(delta[i]==0,0,log(lambda0(Time[i]))+gamma*W[i]+alpha1*des.T[match(Time[i],data.id$Time[data.id$delta==1]),]%*%(beta+binew)+
                                    alpha2*sqrt(t(beta+binew)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+binew)))
            
            log.p.YTb.new=log.hazard.new-log.s.new+log.p.Yb.new+log.p.b.new
            if(log.p.YTb.new<=log.p.YTb){
              #print("oh")
              break
            }
            bi=binew
            log.p.YTb=log.p.YTb.new
            #print(log.p.YTb)
            
          }
          
        }
        
      }else{
        Fishi=t(Zi)%*%Zi/sigma2+solve(D)
      }
      
      b_set[i,]=bi
      # Fish_set[,,i]=Fishi
      inv.Fish_set[,,i]=solve(Fishi)
    }
    
    tr=0 #sum of trace
    for(i in 1:m){
      tra=sum(diag(t(des.Y[data$id==i,])%*%des.Y[data$id==i,]
                   %*%inv.Fish_set[,,i]))
      tr=tr+tra
    }
    
    ###################### Updata parameters using b_set
    sgamma=0
    salpha1=0
    salpha2=0
    sbeta=numeric(q)
    Igamma=0
    Ialpha1=0
    Ialpha2=0
    Ibeta=matrix(0,ncol=q,nrow=q)
    Igamalp1=0
    Igamalp2=0
    Ialp12=0
    
    
    for(i in data.id$id[data.id$delta==1]){
      
      rs=riskset(Time[i])
      Exp.f=numeric(length(rs))
      Exp.fBb=numeric(length(rs))
      Exp.fBb2=numeric(length(rs))
      Exp.fbetaKb=numeric(length(rs))
      sqrbetaKb=numeric(length(rs))
      Exp.sqrfbetaKb=numeric(length(rs))
      Exp.fBbbetaKb=numeric(length(rs))
      Exp.scorebeta=matrix(0,ncol=q,nrow=length(rs))
      
      
      K2=K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]
      desT=des.T[match(Time[i],data.id$Time[data.id$delta==1]),]
      
      betabi=adv(c(beta+b_set[i,]))
      misqrbetaKbKb=c(c(t(beta+b_set[i,])%*%K2%*%(beta+b_set[i,]))^(-1/2)*K2%*%(beta+b_set[i,]))+
        sapply(c(1:q),function(s) sum(diag(inv.Fish_set[,,i]%*%attr((quad.adv(K2,betabi))^(-1/2)*linear_sca.adv(K2,betabi,s),"grad2")))/2)
      
      
      
      for(j in rs){
        f=c(exp(alpha1*desT%*%b_set[j,]+
                  alpha2*sqrt(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))))
        inv.Fishj=inv.Fish_set[,,j]
        
        
        fg1=c(f*(alpha1*desT+alpha2*c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-1/2)*K2%*%(beta+b_set[j,])))
        
        fg2=f*(tcrossprod(alpha1*desT+alpha2*c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-1/2)*K2%*%(beta+b_set[j,]))-
                 alpha2*c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-3/2)*tcrossprod(K2%*%(beta+b_set[j,]))+
                 alpha2*c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-1/2)*K2)
        
        h=c(sqrt(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,])))
        
        hg1=c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-1/2)*K2%*%(beta+b_set[j,])
        
        hg2=-c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-3/2)*tcrossprod(K2%*%(beta+b_set[j,]))+
          c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-1/2)*K2
        
        Exp.f[which(rs==j)]=f+sum(diag(inv.Fishj%*%fg2))/2
        
        Bb=c(t(desT)%*%b_set[j,])
        
        Exp.fBb[which(rs==j)]=f*Bb+
          sum(diag(inv.Fishj%*%(fg2*Bb+fg1%*%t(desT)+desT%*%t(fg1))))/2
        
        Exp.fBb2[which(rs==j)]=f*Bb^2+
          sum(diag(inv.Fishj%*%(Bb^2*fg2+2*Bb*(fg1%*%t(desT)+desT%*%t(fg1))+2*f*desT%*%t(desT))))/2
        
        betaKb=c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))
        
        Exp.fbetaKb[which(rs==j)]=f*betaKb+
          sum(diag(inv.Fishj%*%(betaKb*fg2+2*fg1%*%t(K2%*%(beta+b_set[j,]))+2*(K2%*%(beta+b_set[j,]))%*%t(fg1)+2*f*K2)))/2
        
        sqrbetaKb[which(rs==j)]=h+sum(diag(inv.Fishj%*%hg2))/2
        
        Exp.sqrfbetaKb[which(rs==j)]=f*h+sum(diag(inv.Fishj%*%(h*fg2+fg1%*%t(hg1)+hg1%*%t(fg1)+f*hg2)))/2
        
        Exp.fBbbetaKb[which(rs==j)]=f*h*Bb+sum(diag(inv.Fishj%*%(Bb*(h*fg2+fg1%*%t(hg1)+hg1%*%t(fg1)+f*hg2)+(h*fg1+f*hg1)%*%t(desT)+desT%*%t((h*fg1+f*hg1)))))/2
        
        betabj=adv(c(beta+b_set[j,]))
        Exp.scorebeta[which(rs==j),]=c(f*alpha2*c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-1/2)*K2%*%(beta+b_set[j,]))+
          sapply(c(1:q),function(s) sum(diag(inv.Fishj%*%attr(exp(alpha1*linear.adv(desT,adv(b_set[j,])))*exp(alpha2*(quad.adv(K2,betabj))^(1/2))*(alpha2*(quad.adv(K2,betabj))^(-1/2))*linear_sca.adv(K2,betabj,s),
                                                              "grad2")))/2)
        
      }
      
      ### score of gamma; Fisher of gamma
      d=exp(gamma*W[rs])#*(exp(alpha*xbeta))# a vector, can be cancelled 
      Deno=c(t(d)%*%Exp.f)
      ssgamma=W[i]-t(W[rs]*d)%*%Exp.f/Deno
      sgamma=sgamma+ssgamma
      igamma=t(W[rs]^2*d)%*%Exp.f/Deno-(t(W[rs]*d)%*%Exp.f/Deno)^2
      Igamma=Igamma+igamma
      
      ### score of alpha1; Fisher of alpha1
      xbeta=matrix(rep(c(t(desT)%*%beta),length(rs)),ncol=1)#collection of Xj(Ti)%*%beta
      ssalpha1=desT%*%(beta+b_set[i,])-t(d)%*%(xbeta*Exp.f+Exp.fBb)/Deno
      salpha1=salpha1+ssalpha1
      ialpha1=t(d)%*%(xbeta*Exp.fBb+Exp.fBb2)/Deno-(t(d)%*%(xbeta*Exp.f+Exp.fBb))*(t(d)%*%Exp.fBb)/(Deno^2)
      Ialpha1=Ialpha1+ialpha1
      
      ### score of alpha2; Fisher of alpha2
      ssalpha2=sqrbetaKb[which(rs==i)]-t(d)%*%Exp.sqrfbetaKb/Deno
      salpha2=salpha2+ssalpha2
      ialpha2=t(d)%*%Exp.fbetaKb/Deno-(t(d)%*%Exp.sqrfbetaKb/Deno)^2
      Ialpha2=Ialpha2+ialpha2
      
      ##############
      igamalp1=t(d)%*%(Exp.fBb*W[rs])/Deno-(t(d)%*%Exp.fBb)*(t(W[rs]*d)%*%Exp.f)/(Deno^2)
      Igamalp1=Igamalp1+igamalp1
      
      igamalp2=t(d)%*%(Exp.sqrfbetaKb*W[rs])/Deno-(t(d)%*%Exp.sqrfbetaKb)*(t(W[rs]*d)%*%Exp.f)/(Deno^2)
      Igamalp2=Igamalp2+igamalp2
      
      ialp12=t(d)%*%(Exp.fBbbetaKb+c(desT%*%beta)*Exp.sqrfbetaKb)/Deno-(t(d)%*%Exp.sqrfbetaKb)*(t(d)%*%(Exp.fBb+c(desT%*%beta)*Exp.f))/(Deno^2)
      Ialp12=Ialp12+ialp12
      
      ### part of score of beta; part of Fisher of beta
      
      ssbeta=alpha1*c(desT)+alpha2*misqrbetaKbKb-(Deno*alpha1*desT+colSums(d*Exp.scorebeta))/Deno
      sbeta=sbeta+ssbeta
      ibeta=tcrossprod(ssbeta)
      Ibeta=Ibeta+ibeta
      
    }
    
    
    res=c(data$Y-apply(des.Y*t(beta+t(b_set[rep(c(1:m),l),])),1,sum))
    sigma.2=1/N*(t(data$Y-des.Y%*%beta)%*%(data$Y-des.Y%*%beta-2*diXZ%*%c(t(b_set)))+
                   tr+t(diXZ%*%c(t(b_set)))%*%diXZ%*%c(t(b_set)))
    sigma.2=c(sigma.2)
    sbeta=c(sbeta)+1/sigma.2*t(des.Y)%*%res
    Ibeta=Ibeta+(-(sigma.2)^(-2)*2/N)*t(des.Y)%*%res%*%t(res)%*%des.Y+(sigma.2)^(-1)*t(des.Y)%*%des.Y
    
    Igamalp12=matrix(0,ncol=3,nrow=3)
    Igamalp12[upper.tri(Igamalp12)]=c(Igamalp1,Igamalp2,Ialp12)
    Igamalp12=Igamalp12+t(Igamalp12)
    diag(Igamalp12)=c(Igamma,Ialpha1,Ialpha2)
    
    if(k>1){
      loglik=logLik_set[k-1]
    }
    
    for(v in c(10,0:9)){
      
      stepsize=2^(-v)
      par.sur=c(gamma,alpha1,alpha2)+stepsize*solve(Igamalp12,c(sgamma,salpha1,salpha2))
      gammanew=par.sur[1]
      alpha1new=par.sur[2]
      alpha2new=par.sur[3]
      betanew=c(beta+stepsize*solve(Ibeta)%*%sbeta)
      
      sigma2new=1/N*(t(data$Y-des.Y%*%betanew)%*%(data$Y-des.Y%*%betanew-2*apply(des.Y*b_set[rep(c(1:m),l),],1,sum))+tr+sum(apply(des.Y*b_set[rep(c(1:m),l),],1,sum)^2))
      sigma2new=c(sigma2new)
      Dnew=1/m*(apply(inv.Fish_set,c(1,2),sum)+t(b_set)%*%b_set)
      Dnew=diag(diag(Dnew))
      
      cc=c()
      for(t in sort(unique(data$Time[data$delta==1]))){
        rs=riskset(t)
        Exp.f=numeric(length(rs))
        for(j in rs){
          f=c(exp(alpha1new*des.T[match(t,data.id$Time[data.id$delta==1]),]%*%b_set[j,]+
                    alpha2new*sqrt(t(betanew+b_set[j,])%*%K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(betanew+b_set[j,]))))
          
          inv.Fishj=inv.Fish_set[,,j]
          
          fg2=f*(tcrossprod(alpha1new*des.T[match(t,data.id$Time[data.id$delta==1]),]+alpha2new*c(t(betanew+b_set[j,])%*%K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(betanew+b_set[j,]))^(-1/2)*K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(betanew+b_set[j,]))-
                   alpha2new*c(t(betanew+b_set[j,])%*%K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(betanew+b_set[j,]))^(-3/2)*tcrossprod(K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(betanew+b_set[j,]))+
                   alpha2new*c(t(betanew+b_set[j,])%*%K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(betanew+b_set[j,]))^(-1/2)*K.2[,,match(t,data.id$Time[data.id$delta==1])])
          
          
          Exp.f[which(rs==j)]=f+sum(diag(inv.Fishj%*%fg2))/2
        }
        
        xbeta=matrix(rep(des.T[match(t,data.id$Time[data.id$delta==1]),],length(rs)),ncol=q,byrow=TRUE)%*%betanew #collection of Xj(t)%*%beta1
        cc=c(cc,length(data.id$id[(data.id$Time==t)&(data.id$delta==1)])/(c(t(exp(gammanew*W[rs]+alpha1new*xbeta))%*%Exp.f)))
      }
      
      cumbase=data.frame(hazard=0,time=c(0,sort(unique(data$Time[data$delta==1]))))
      cumbase[,1]=cumsum(c(0,cc))
      ### Log-likelihood 
      logLikmcnew=logLik(data,data.id,gammanew,alpha1new,alpha2new,cumbase,betanew,sigma2new,Dnew,internal_knots,Q.2,L=2000)
      if(v==10){
        logLikcom=logLikmcnew
        print(logLikcom)
      }else{
        if(logLikmcnew>logLikcom){
          print(v)
          print(logLikmcnew)
          break
        }
      }
      
    }
    
    
    
    # if(v==9){
    #   break
    # 
    # }
    
    if(k>1){
      if(((abs(logLikmcnew-logLik_set[k-1])/abs(logLik_set[k-1]))<10^(-4))){
        if(logLikmcnew>logLik_set[k-1]){
          iter=c(k,logLikmcnew)
          par=c(gammanew,alpha1new,alpha2new,betanew,sigma2new,diag(Dnew)) 
          return(list(iter=iter,par=par,cumbase=cumbase,b_set=b_set))
          
        }else{
          iter=c(k-1,logLik_set[k-1])
          par=c(gamma,alpha1,alpha2,beta,sigma2,diag(D))
          return(list(iter=iter,par=par,cumbase=cumbase,b_set=b_set))
          
        }
        
      }
    }
    
    print(k)
    logLik_set[k]=logLikmcnew
    gamma_set[k]=gammanew
    print(gamma_set[k])
    alpha1_set[k]=alpha1new
    print(alpha1_set[k])
    alpha2_set[k]=alpha2new
    print(alpha2_set[k])
    beta_set[k,]=betanew
    sigma2_set[k]=sigma2new
    D_set[,,k]=Dnew
    
    gamma=gammanew
    alpha1=alpha1new
    alpha2=alpha2new
    beta=betanew
    sigma2=sigma2new
    D=Dnew
    #Q.fun_set[k]=Q.fun.up
  }
  return(NULL)
  
}


