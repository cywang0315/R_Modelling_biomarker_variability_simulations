m=3700
q=6
internal_knots=c(0.25,1.5)
boundary_knots=c(0,5.8)

####### B-spline basis #########
## my cubic B-spline basis functions ####
mycubicbs=function(s,internal_knots,boundary_knots){
  df=length(internal_knots)+3+1
  r=numeric(df)
  
  out_knots=c(seq(from=min(boundary_knots),by=-0.1,length=4),seq(from=max(boundary_knots),by=0.1,length=4))
  #out_knots=c(seq(from=min(boundary_knots),by=-2,length=4),seq(from=max(boundary_knots),by=2,length=4))
  
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

# true value of longitudinal biomarker (no measurement error)
longit.true=function(t,x1,x2,fix_eff,ran_eff,eta,xi){
  desmat=mycubicbs(t,internal_knots=internal_knots,boundary_knots=boundary_knots)$mat# design matrix
  x1%*%eta+desmat%*%(fix_eff+ran_eff)+x2*xi
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
# cumtrue=function(t){ifelse(t<1.5,0,0.003*(t-1.5))}
# cumwei=function(t,b,k){ifelse(t<1.5,0,b/k*(t-1.5)^{k})}
# cumgom=function(t,b,theta){ifelse(t<1.5,0,b/theta*(exp(theta*(t-1.5))-1))}

gendat=function(s,obstim,obsmax,gammatrue,alpha1true,alpha2true,betatrue,beta1true,beta2true,etatrue,xitrue,D0,sigmatrue,knots,Q.2){
  set.seed(s)
  basetrue=function(t){ifelse(t<1.5,0,0.0001*exp(0.7*(t-1.5)))} ## 0 when t<1.5
  b=cbind(numeric(m),mvrnorm(m,mu=rep(0,q-1),Sigma=D0))# random effects in longitudinal submodel
  W=matrix(0,ncol=length(gammatrue),nrow=m)
  tmtind=runif(m)
  ###W=(tmt1, tmt2,sex,smoke,age/10)
  W=cbind(ifelse(tmtind<=1/4,1,0),ifelse((tmtind<=1/2)&(tmtind>1/4),1,0),
          rbinom(m,1,prob = 0.42),rbinom(m,1,prob = 0.164),rtnorm(m,mean=7,sd=sqrt(0.29),lower = 6.5,upper = 7.4))
  desmat=mycubicbs(seq(1.5,obsmax,by=0.005),internal_knots=internal_knots,boundary_knots=boundary_knots)$mat# design matrix  
  #set.seed(1)
  TD.fun=function(i,t){# i is a scale, t can be a vector
    betai=W[i,1]*beta1true+W[i,2]*beta2true+(1-W[i,1]-W[i,2])*betatrue
    tt=function(t){
      c(basetrue(t)*exp(W[i,]%*%gammatrue+alpha1true*desmat[match(t,seq(1.5,obsmax,by=0.005)),]%*%(betai+b[i,])+
                          alpha2true*sqrt(t(betai+b[i,])%*%Q.2%*%(R.2(t,knots)-R.2(0.25,knots))%*%t(Q.2)%*%(betai+b[i,]))))
    }
    #hazard=sapply(t,tt)
    #p=1-exp(-hazard*0.01)
    r=rbinom(n=length(t),size=1,prob=1-exp(-sapply(t,tt)*0.005)) ##prob: event occurs
    Time1=obsmax
    delta1=0
    if(max(r)>0){
      Time1=t[min(which(r==1))]
      delta1=1
    }
    return(list(Time1=Time1,delta1=delta1))
  }
  
  TD=sapply(c(1:m),TD.fun,t=seq(1.5,obsmax,by=0.005))
  Time1=as.numeric(TD[1,])
  delta1=as.numeric(TD[2,])
  #hist(Time1)
  sum(delta1)
  ####generate censoring 
  #censtim=sample(seq(12,20,by=0.001),size=m,replace = TRUE)
  cendis=UnivarMixingDistribution(Truncate(Norm(3,1), upper=5.8,lower=1.5), 
                                   Truncate(Norm(5.2,0.3), lower=1.5,upper=5.8), 
                                   mixCoeff=c(0.6, 0.4))
  rcenstim=r(cendis)
  censtim=rcenstim(m)
  #hist(censtim)
  delta=ifelse(Time1<=censtim,delta1,0)
  sum(delta)
  Time=ifelse(Time1<=censtim,Time1,censtim)
  #hist(Time)
  #hist(Time[delta==1])
  #hist(Time[delta==0])
  #############generate longitudinal data
  #set.seed(1)
  X1=W[,3:5]
  Y=c()
  l=numeric(m)
  time=c()
  for (i in 1:m){
    l[i]=length(obstim[obstim<=Time[i]])
    betai=W[i,1]*beta1true+W[i,2]*beta2true+(1-W[i,1]-W[i,2])*betatrue
    for(j in obstim[obstim<=Time[i]]){
      YY=longit.true(j,X1[i,],ifelse(j %in% c(0,1,2,3,4,5),1,0),fix_eff=betai,ran_eff=b[i,],etatrue,xitrue)+rnorm(1,sd=sigmatrue)
      Y=rbind(Y,YY)
      time=rbind(time,j)
    }
  }
  
  Y=as.vector(Y)
  time=as.vector(time)
  ### construct data frame. Note here we use "id" not "sub"
  data=data.frame(id=rep(c(1:m),l),Time=rep(Time,l),delta=rep(delta,l),Y=Y,W=apply(W,2,rep,times=l),
                  time=time, start=time,event=numeric(length(Y))) ##should adjust "stop" and "event"
  data$start[data$start!=0]=data$start[data$start!=0]-0.01
  data$stop=c(data$start[-1],0)
  data$stop=ifelse(data$stop<data$start,data$Time,data$stop)
  data$event=numeric(dim(data)[1])
  data$event=ifelse((data$stop==data$Time)&(data$delta==1),1,data$event)
  #data$event[which((data$start==data$stop)&(data$delta==1))-1]=0
  sum(data$event)
  #data$stop[data$start==data$stop]=data$stop[data$start==data$stop]+0.01
  data$WC=ifelse(data$time %in% c(0,1,2,3,4,5),1,0)
  data$tmt=numeric(dim(data)[1])+3
  data$tmt=ifelse(data$W.1==1,1,data$tmt)
  data$tmt=ifelse(data$W.2==1,2,data$tmt)
  return(data)
  
}

############################# Get initial value ##################
inival=function(data,data.id,ctl,knots,Q.2){
  lme.data=lme(fixed=Y~-1+time1+time2+time3+time4+time5+time6+W.3+W.4+W.5+WC+
                 W.1:(time1+time2+time3+time4+time5+time6)+
                 W.2:(time1+time2+time3+time4+time5+time6),
               random=list(id=pdDiag(form=~-1+time2+time3+time4+time5+time6)),data = data,control=ctl)
  
  
  beta=lme.data$coefficients[[1]][1:6]
  beta1=lme.data$coefficients[[1]][1:6]+lme.data$coefficients[[1]][11:16]
  beta2=lme.data$coefficients[[1]][1:6]+lme.data$coefficients[[1]][17:22]
  eta1=lme.data$coefficients[[1]][c(7,8,9)]
  eta2=lme.data$coefficients[[1]][10]
  sigma2=lme.data$sigma^2
  D=diag(c(as.numeric(VarCorr(lme.data)[1:5])))
  
  cox.data=coxph(Surv(Time,delta)~W.1+W.2+W.3+W.4+W.5,data=data.id,x=TRUE,method="breslow")

  #### two stage ####
  data.ts=data
  l=as.data.frame(table(data$id))$Freq
  raneff=apply(random.effects(lme.data),2,rep,times=l)
  data.ts$preY=numeric(dim(data)[1])
  ### Baseline covariates in longitudinal analysis
  X=cbind(data.id$W.3,data.id$W.4,data.id$W.5)
  for(i in 1:m){
    coe=lme.data$coefficients[[1]][1:6]+data.id$W.1[i]*lme.data$coefficients[[1]][11:16]+
       data.id$W.2[i]*lme.data$coefficients[[1]][17:22]+c(0,as.numeric(random.effects(lme.data)[i,]))
    
    data.ts$preY[data.ts$id==i]=as.numeric(mycubicbs(data.ts$time[data.ts$id==i],internal_knots = c(0.25,1.5),boundary_knots = c(0,5.8))$mat%*%coe)+
      c(X[i,]%*%lme.data$coefficients[[1]][7:9])
  }
  ##### prevar: variability 
  
  data.ts$prevar0.25c=numeric(dim(data.ts)[1])
  for(i in 1:dim(data.ts)[1]){
    coe=lme.data$coefficients[[1]][1:6]+data.ts$W.1[i]*lme.data$coefficients[[1]][11:16]+
      data.ts$W.2[i]*lme.data$coefficients[[1]][17:22]+c(0,as.numeric(random.effects(lme.data)[data.id$id==data.ts$id[i],]))

    if(data.ts$time[i]>0.25){
      data.ts$prevar0.25c[i]=sqrt(t(coe)%*%(Q.2%*%(R.2(data.ts$time[i],knots)-R.2(0.25,knots))%*%t(Q.2))%*%(coe))
    }
  }
  
  coxts.data=coxph(Surv(start, stop, event)~W.1+W.2+W.3+W.4+W.5+preY+prevar0.25c,data=data.ts,method="breslow")
  gamma=coxts.data$coefficients[1:5]
  alpha1=coxts.data$coefficients[6]
  alpha2=coxts.data$coefficients[7]
  cumbase=rbind(c(0,0),basehaz(coxts.data, centered=FALSE))
  out=list(beta=beta,beta1=beta1,beta2=beta2,eta1=eta1,eta2=eta2,sigma2=sigma2,D=D,
           gamma=gamma,alpha1=alpha1,alpha2=alpha2,cumbase=cumbase)

}

###calculate log likelihood 

logLik=function(data,data.id,gamma,alpha1,alpha2,beta,beta1,beta2,eta1,eta2,sigma2,D,cumbase,knots,Q.2){
  
  lambda0=function(t){
    if(t %in% unique(data.id$Time[data.id$delta==1])){
      cumbase[which(cumbase[,2]==t),1]-cumbase[which(cumbase[,2]==t)-1,1]} ##warn: options(digits = 22) default 7
    else{0}
  }
  
  des.Y=mycubicbs(data$time,internal_knots,boundary_knots)$mat
  des.T=mycubicbs(data.id$Time[data.id$delta==1],internal_knots,boundary_knots)$mat
  Time=data.id$Time
  delta=data.id$delta
  W=cbind(data.id$W.1,data.id$W.2,data.id$W.3,data.id$W.4,data.id$W.5)
  X1=matrix(c(data.id$W.3,data.id$W.4,data.id$W.5),ncol=3)
  K.2=array(c(sapply(data.id$Time[data.id$delta==1],function(t) Q.2%*%(R.2(t,knots)-R.2(0.25,knots))%*%t(Q.2))),
            dim=c(q,q,length(data.id$Time[data.id$delta==1]))) 
  l=as.data.frame(table(data$id))$Freq
  m=dim(data.id)[1]
  Beta=matrix(0,ncol=q,nrow=m)
  for (i in 1:m){
    if(data.id$tmt[i]==1){
      Beta[i,]=beta1
    }else{
      if(data.id$tmt[i]==2){
        Beta[i,]=beta2
      }else{
        Beta[i,]=beta
      }
    }
  }
  set.seed(1)
  L=2000
  log.p.Y=numeric(m)
  log.s_set=matrix(0,nrow=m,ncol=L)
  ZTimeb=matrix(0,nrow=sum(delta==1),ncol=L) ##Z_Time%*%b
  bKb=matrix(0,nrow=sum(delta==1),ncol=L)
  rmc=rmvnorm(L,mean=rep(0,q-1))
  for(i in 1:m){
    Zi=des.Y[data$id==i,-1]
    Xi=cbind(matrix(rep(X1[data.id$id==i,],l[data.id$id==i]),ncol=dim(X1)[2],byrow=T),
             des.Y[data$id==i,],data$WC[data$id==i])
    Yi=data$Y[data$id==i]
    
    chol_invSig=chol(solve(D)+c(sigma2^(-1))*t(Zi)%*%Zi)
    Q=solve(chol_invSig)
    if(data.id$tmt[data.id$id==i]==1){
      betai=beta1
    }else{
      if(data.id$tmt[data.id$id==i]==2){
        betai=beta2
      }else{
        betai=beta
      }
    }
    
    mu=c(Q%*%t(Q)%*%(t(Zi)/c(sigma2))%*%(Yi-Xi%*%c(eta1,betai,eta2))) #E(bi|Yi,theta),initial value of bi
    
    Sigma=solve(solve(D)+c(sigma2^(-1))*t(Zi)%*%Zi)
    
    log.p.Y[data.id$id==i]=dmvnorm(Yi,mean=c(Xi%*%c(eta1,betai,eta2)),sigma=sigma2*diag(length(Yi))+Zi%*%D%*%t(Zi),log=TRUE)
    rmci=t(mu+t(chol(Sigma))%*%t(rmc))## matrix with ncol=q,nrow=L, MC sample from bi|yi
    if(delta[data.id$id==i]){
      ZTimeb[match(i,data.id$id[data.id$delta==1]),]=des.T[match(i,data.id$id[data.id$delta==1]),-1]%*%t(rmci)
      bKb[match(i,data.id$id[data.id$delta==1]),]=colSums((t(cbind(numeric(L),rmci))+betai)*(K.2[,,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%(t(cbind(numeric(L),rmci))+betai)))
    }
    
    a=unique(data.id$Time[(data.id$Time<=Time[data.id$id==i])&(data.id$delta==1)])
    if(length(a)>0){
      if(length(a)==1){
        btKb=colSums((t(cbind(numeric(L),rmci))+betai)*(K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(t(cbind(numeric(L),rmci))+betai)))
      }else{
        btKb=apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) colSums((t(cbind(numeric(L),rmci))+betai)*(y%*%(t(cbind(numeric(L),rmci))+betai))))
      }
      
      log.s_set[data.id$id==i,]=apply(sapply(a,lambda0)*exp(c(gamma%*%W[data.id$id==i,]))*
                                       exp(alpha1*(c(X1[data.id$id==i,]%*%eta1)+des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(t(cbind(numeric(L),rmci))+betai))+
                                             alpha2*t(sqrt(btKb))),2,sum)
      
    }
    
  }
  
  log.hazard=matrix(0,ncol=L,nrow=m)
  log.hazard[which(delta==1),]=log(sapply(Time[delta==1],lambda0))+c(W[delta==1,]%*%gamma)+alpha1*(c(X1[delta==1,]%*%eta1+apply(des.T*Beta[delta==1,],1,sum))+ZTimeb)+
    alpha2*sqrt(bKb)
  log.p.Tb=log.hazard-log.s_set
  p.Tb=exp(log.p.Tb)
  logLikmc.or=sum(log.p.Y+log(rowMeans(p.Tb)))
  return(logLikmc.or)              
  
}

################ Estimation ###############
est=function(data,data.id,logLikmc.or,gamma,alpha1,alpha2,beta,beta1,beta2,eta1,eta2,sigma2,D,cumbase,knots,Q.2){
  lambda0=function(t){
    if(t %in% unique(data.id$Time[data.id$delta==1])){
      cumbase[which(cumbase[,2]==t),1]-cumbase[which(cumbase[,2]==t)-1,1]} ##warn: options(digits = 22) default 7
    else{0}
  }
  
  riskset=function(t){# individuals in the risk set
    data.id$id[data.id$Time>=t]## "=" is important
  }
  
  ########################
  des.Y=mycubicbs(data$time,internal_knots,boundary_knots)$mat
  des.T=mycubicbs(data.id$Time[data.id$delta==1],internal_knots,boundary_knots)$mat
  Time=data.id$Time
  delta=data.id$delta
  W=cbind(data.id$W.1,data.id$W.2,data.id$W.3,data.id$W.4,data.id$W.5)
  X1=matrix(c(data.id$W.3,data.id$W.4,data.id$W.5),ncol=3)
  X2=data$WC
  K.2=array(c(sapply(data.id$Time[data.id$delta==1],function(t) Q.2%*%(R.2(t,knots)-R.2(0.25,knots))%*%t(Q.2))),
            dim=c(q,q,length(data.id$Time[data.id$delta==1]))) 
  l=as.data.frame(table(data$id))$Freq
  m=dim(data.id)[1]
  N=dim(data)[1]
  ### EM algorithm  ###
  K=10
  gamma_set=matrix(0,nrow=K,ncol=length(gamma))
  alpha1_set=numeric(K)
  alpha2_set=numeric(K)
  eta1_set=matrix(0,nrow=K,ncol=length(eta1))
  eta2_set=numeric(K)
  beta_set=matrix(0,ncol=q,nrow=K)
  beta1_set=matrix(0,ncol=q,nrow=K)
  beta2_set=matrix(0,ncol=q,nrow=K)
  sigma2_set=numeric(K)
  D_set=array(0,dim=c(q-1,q-1,K)) 
  Q.fun_set=numeric(K)
  logLik_set=numeric(K)
  b_set=matrix(0,ncol=q-1,nrow=m)
  inv.Fish_set=array(0,dim=c(q-1,q-1,m)) ## V(bi)~=inv.Fishi
  ########################
  logLik_set[1]=logLikmc.or
  gamma_set[1,]=gamma
  alpha1_set[1]=alpha1
  alpha2_set[1]=alpha2
  eta1_set[1,]=eta1
  eta2_set[1]=eta2
  beta_set[1,]=beta
  beta1_set[1,]=beta1
  beta2_set[1,]=beta2
  sigma2_set[1]=sigma2
  D_set[,,1]=D
  
  for(k in 2:K){
    ### Get bi=argmax p(Ti,deltai,Yi,bi;theta) using Newton-Raphson
    # Score of bi: dlog(p(Ti,deltai,Yi,bi;theta))/dbi
    for(i in data.id$id){
      Zi=des.Y[data$id==i,-1]
      Xi=cbind(matrix(rep(X1[data.id$id==i,],l[data.id$id==i]),ncol=dim(X1)[2],byrow=T),
               des.Y[data$id==i,],data$WC[data$id==i])
      Yi=data$Y[data$id==i]
      chol_invSig=chol(solve(D)+c(sigma2^(-1))*t(Zi)%*%Zi)
      Q=solve(chol_invSig)
      if(data.id$tmt[data.id$id==i]==1){
        betai=beta1
      }else{
        if(data.id$tmt[data.id$id==i]==2){
          betai=beta2
        }else{
          betai=beta
        }
      }
      
      bi=c(Q%*%t(Q)%*%(t(Zi)/c(sigma2))%*%(Yi-Xi%*%c(eta1,betai,eta2))) #E(bi|Yi,theta),initial value of bi
      a=unique(data.id$Time[(data.id$Time<=Time[data.id$id==i])&(data.id$delta==1)])
      
      if(length(a)>0){
        if(length(a)>1){
          ### initial value of p(Ti,deltai,Yi,bi;theta)
          log.p.b=dmvnorm(bi,mean=rep(0,q-1),sigma=D,log=TRUE)
          log.p.Yb=sum(dnorm(Yi,mean=c(Xi%*%c(eta1,betai,eta2)+Zi%*%bi),sd=sqrt(sigma2),log=TRUE)) # p(Y|b)
          log.s=exp(c(gamma%*%W[data.id$id==i,]))*sum(sapply(a,lambda0)*exp(alpha1*(c(X1[data.id$id==i,]%*%eta1)+des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(betai+c(0,bi)))+
                                                                             alpha2*apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) sqrt(t(betai+c(0,bi))%*%y%*%(betai+c(0,bi))))))
          
          log.hazard=ifelse(delta[data.id$id==i]==0,0,log(lambda0(Time[data.id$id==i]))+c(gamma%*%W[data.id$id==i,])+alpha1*(X1[data.id$id==i,]%*%eta1+des.T[match(Time[data.id$id==i],data.id$Time[data.id$delta==1]),]%*%(betai+c(0,bi)))+
                              alpha2*sqrt(t(betai+c(0,bi))%*%K.2[,,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%(betai+c(0,bi))))
          
          log.p.YTb=log.hazard-log.s+log.p.Yb+log.p.b
          
          for(kk in 1:10){
            wei=c(exp(c(gamma%*%W[data.id$id==i,]))*sapply(a,lambda0)*exp(alpha1*(c(X1[data.id$id==i,]%*%eta1)+des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(betai+c(0,bi)))+
                                                                           alpha2*apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) sqrt(t(betai+c(0,bi))%*%y%*%(betai+c(0,bi))))))
            ZK=alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),-1]+alpha2*apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) (c(t(betai+c(0,bi))%*%y%*%(betai+c(0,bi))))^(-1/2))*
              t(apply(K.2[,-1,match(a,data.id$Time[data.id$delta==1])],3,function(y) t(y)%*%betai)+apply(K.2[-1,-1,match(a,data.id$Time[data.id$delta==1])],3,function(y) y%*%bi))
            Sbi=-colSums(wei*ZK)+t(Zi)%*%(Yi-Xi%*%c(eta1,betai,eta2)-Zi%*%bi)/sigma2-solve(D)%*%bi
            
            Fishi=tensor(array(apply(ZK,1,tcrossprod),dim=c(q-1,q-1,length(a)))+
                           alpha2*array(apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) (c(t(betai+c(0,bi))%*%y%*%(betai+c(0,bi))))^(-1/2)*y[-1,-1]),dim=c(q-1,q-1,length(a)))-
                           alpha2*array(apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) (c(t(betai+c(0,bi))%*%y%*%(betai+c(0,bi))))^(-3/2)*tcrossprod(t(y[,-1])%*%betai+y[-1,-1]%*%bi)),dim=c(q-1,q-1,length(a))),
                         wei,3,1)+t(Zi)%*%Zi/sigma2+solve(D)
            if(delta[data.id$id==i]==1){
              Sbi=alpha1*des.T[match(Time[data.id$id==i],data.id$Time[data.id$delta==1]),-1]+alpha2*(c(t(betai+c(0,bi))%*%K.2[,,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%(betai+c(0,bi))))^(-1/2)*
                (t(K.2[,-1,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])])%*%betai+K.2[-1,-1,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%bi)+Sbi
              Fishi=-alpha2*(-(c(t(betai+c(0,bi))%*%K.2[,,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%(betai+c(0,bi))))^(-3/2)*tcrossprod(t(K.2[,-1,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])])%*%betai+K.2[-1,-1,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%bi)+
                               +(c(t(betai+c(0,bi))%*%K.2[,,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%(betai+c(0,bi))))^(-1/2)*K.2[-1,-1,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])])+Fishi
            }
            ## update bi using Newton-Raphson: bi=bi+solve(I)S
            
            for (v in 0:10){
              
              stepsize=2^(-v)
              binew=c(bi+stepsize*solve(Fishi)%*%Sbi)
              ## Check
              log.p.b.new=dmvnorm(binew,mean=rep(0,q-1),sigma=D,log=TRUE)
              log.p.Yb.new=sum(dnorm(Yi,mean=c(Xi%*%c(eta1,betai,eta2)+Zi%*%binew),sd=sqrt(sigma2),log=TRUE)) # p(Y|b)
              log.s.new=exp(c(gamma%*%W[data.id$id==i,]))*sum(sapply(a,lambda0)*exp(alpha1*(c(X1[data.id$id==i,]%*%eta1)+des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(betai+c(0,binew)))+
                                                                                     alpha2*apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) sqrt(t(betai+c(0,binew))%*%y%*%(betai+c(0,binew))))))
              log.hazard.new=ifelse(delta[data.id$id==i]==0,0,log(lambda0(Time[data.id$id==i]))+c(gamma%*%W[data.id$id==i,])+alpha1*(X1[data.id$id==i,]%*%eta1+des.T[match(Time[data.id$id==i],data.id$Time[data.id$delta==1]),]%*%(betai+c(0,binew)))+
                                      alpha2*sqrt(t(betai+c(0,binew))%*%K.2[,,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%(betai+c(0,binew))))
              
              log.p.YTb.new=log.hazard.new-log.s.new+log.p.Yb.new+log.p.b.new
              
              if(log.p.YTb.new>=log.p.YTb){
                print(v)
                break
              }
            }
            
            
            
            if(all(abs((binew-bi)/bi)[which(is.nan(abs((binew-bi)/bi))!=TRUE)]<10^(-3))){
              cat("kk",kk,"\n")
              break
            }
            
            bi=binew
            log.p.YTb=log.p.YTb.new
            
          }
        }else{
          ### initial value of p(Ti,deltai,Yi,bi;theta)
          log.p.b=dmvnorm(bi,mean=rep(0,q-1),sigma=D,log=TRUE)
          log.p.Yb=sum(dnorm(Yi,mean=c(Xi%*%c(eta1,betai,eta2)+Zi%*%bi),sd=sqrt(sigma2),log=TRUE)) # p(Y|b)
          
          log.s=exp(c(gamma%*%W[data.id$id==i,]))*sum(sapply(a,lambda0)*exp(alpha1*(c(X1[data.id$id==i,]%*%eta1)+des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(betai+c(0,bi)))+
                                                                             alpha2*sqrt(t(betai+c(0,bi))%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(betai+c(0,bi)))))
          
          log.hazard=ifelse(delta[data.id$id==i]==0,0,log(lambda0(Time[data.id$id==i]))+c(gamma%*%W[data.id$id==i,])+alpha1*(X1[data.id$id==i,]%*%eta1+des.T[match(Time[data.id$id==i],data.id$Time[data.id$delta==1]),]%*%(betai+c(0,bi)))+
                              alpha2*sqrt(t(betai+c(0,bi))%*%K.2[,,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%(betai+c(0,bi))))
          
          log.p.YTb=log.hazard-log.s+log.p.Yb+log.p.b
          
          for(kk in 1:10){
            wei=c(exp(c(gamma%*%W[data.id$id==i,]))*sapply(a,lambda0)*exp(alpha1*(c(X1[data.id$id==i,]%*%eta1)+des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(betai+c(0,bi)))+
                                                                           alpha2*sqrt(t(betai+c(0,bi))%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(betai+c(0,bi)))))
            ZK=alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),-1]+
              alpha2*((c(t(betai+c(0,bi))%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(betai+c(0,bi))))^(-1/2)*(t(K.2[,-1,match(a,data.id$Time[data.id$delta==1])])%*%betai+K.2[-1,-1,match(a,data.id$Time[data.id$delta==1])]%*%bi))
            
            Sbi=-colSums(wei*ZK)+t(Zi)%*%(Yi-Xi%*%c(eta1,betai,eta2)-Zi%*%bi)/sigma2-solve(D)%*%bi
            
            Fishi=wei*(tcrossprod(ZK)-alpha2*(c(t(betai+c(0,bi))%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(betai+c(0,bi))))^(-3/2)*tcrossprod(t(K.2[,-1,match(a,data.id$Time[data.id$delta==1])])%*%betai+K.2[-1,-1,match(a,data.id$Time[data.id$delta==1])]%*%bi)+
                         alpha2*(c(t(betai+c(0,bi))%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(betai+c(0,bi))))^(-1/2)*K.2[-1,-1,match(a,data.id$Time[data.id$delta==1])])+t(Zi)%*%Zi/sigma2+solve(D)
            
            if(delta[data.id$id==i]==1){
              Sbi=alpha1*des.T[match(Time[data.id$id==i],data.id$Time[data.id$delta==1]),-1]+alpha2*(c(t(betai+c(0,bi))%*%K.2[,,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%(betai+c(0,bi))))^(-1/2)*
                (t(K.2[,-1,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])])%*%betai+K.2[-1,-1,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%bi)+Sbi
              Fishi=-alpha2*(-(c(t(betai+c(0,bi))%*%K.2[,,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%(betai+c(0,bi))))^(-3/2)*tcrossprod(t(K.2[,-1,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])])%*%betai+K.2[-1,-1,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%bi)+
                               +(c(t(betai+c(0,bi))%*%K.2[,,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%(betai+c(0,bi))))^(-1/2)*K.2[-1,-1,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])])+Fishi
              
            }
            ## update bi using Newton-Raphson: bi=bi+solve(I)S
            for (v in 0:10){
              
              stepsize=2^(-v)
              binew=c(bi+stepsize*solve(Fishi)%*%Sbi)
              
              ## Check
              log.p.b.new=dmvnorm(binew,mean=rep(0,q-1),sigma=D,log=TRUE)
              log.p.Yb.new=sum(dnorm(Yi,mean=c(Xi%*%c(eta1,betai,eta2)+Zi%*%bi),sd=sqrt(sigma2),log=TRUE)) # p(Y|b)
              log.s.new=exp(c(gamma%*%W[data.id$id==i,]))*sum(sapply(a,lambda0)*exp(alpha1*(c(X1[data.id$id==i,]%*%eta1)+des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(betai+c(0,binew)))+
                                                                                     alpha2*sqrt(t(betai+c(0,binew))%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(betai+c(0,binew)))))
              log.hazard.new=ifelse(delta[data.id$id==i]==0,0,log(lambda0(Time[data.id$id==i]))+c(gamma%*%W[data.id$id==i,])+alpha1*(X1[data.id$id==i,]%*%eta1+des.T[match(Time[data.id$id==i],data.id$Time[data.id$delta==1]),]%*%(betai+c(0,binew)))+
                                      alpha2*sqrt(t(betai+c(0,binew))%*%K.2[,,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%(betai+c(0,binew))))
              
              log.p.YTb.new=log.hazard.new-log.s.new+log.p.Yb.new+log.p.b.new
              
              if(log.p.YTb.new>=log.p.YTb){
                print(v)
                break
              }
              
            }
            
            if(all(abs((binew-bi)/bi)[which(is.nan(abs((binew-bi)/bi))!=TRUE)]<10^(-3))){
              cat("kk",kk,"\n")
              break
            }
            
            bi=binew
            log.p.YTb=log.p.YTb.new
            
          }
          
        }
        
      }else{
        
        Fishi=t(Zi)%*%Zi/sigma2+solve(D)
      }
      
      b_set[data.id$id==i,]=bi
      # Fish_set[,,i]=Fishi
      inv.Fish_set[,,data.id$id==i]=solve(Fishi)
      
    }
    
    
    tr=0 #sum of trace
    for(i in data.id$id){
      tra=sum(diag(t(des.Y[data$id==i,-1])%*%des.Y[data$id==i,-1]
                   %*%inv.Fish_set[,,data.id$id==i]))
      tr=tr+tra
      
    }
    
    ###################### Updata parameters using b_set
    sgamma=numeric(dim(W)[2])
    seta1t=numeric(dim(X1)[2])
    salpha1=0
    salpha2=0
    sbetat=numeric(q)
    sbeta1t=numeric(q)
    sbeta2t=numeric(q)
    seta2=0
    Igamma=matrix(0,ncol=dim(W)[2],nrow=dim(W)[2])
    Ieta1t=matrix(0,ncol=dim(X1)[2],nrow=dim(X1)[2])
    Ialpha1=0
    Ialpha2=0
    Ieta2=0
    Igamalp1=0
    Igamalp2=0
    Ialp12=0
    Ibetat=matrix(0,ncol=q,nrow=q)
    Ibeta1t=matrix(0,ncol=q,nrow=q)
    Ibeta2t=matrix(0,ncol=q,nrow=q)
    
    for (i in data.id$id[data.id$delta==1]){
      rs=riskset(Time[data.id$id==i]) ##id
      Exp.f=numeric(length(rs))
      Exp.fBb=numeric(length(rs))
      Exp.fBb2=numeric(length(rs))
      Exp.fbetaKb=numeric(length(rs))
      # Exp.fbetaKb2=numeric(length(rs))
      # Exp.fb=matrix(0,nrow=length(rs),ncol=q-1)
      # Exp.fbb=array(0,dim=c(q-1,q-1,length(rs)))
      K2=K.2[,,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])]
      desT=des.T[match(Time[data.id$id==i],data.id$Time[data.id$delta==1]),]
      d1=numeric(length(rs))
      Bbeta=numeric(length(rs))
      # betaKbeta=numeric(length(rs))
      # Kbeta=matrix(0,ncol=q,nrow=length(rs))
      
      
      sqrbetaKb=numeric(length(rs))
      Exp.sqrfbetaKb=numeric(length(rs))
      Exp.fBbbetaKb=numeric(length(rs))
      Exp.scorebeta=matrix(0,ncol=q,nrow=length(rs))
      
      
      if(data.id$tmt[data.id$id==i]==1){
        betai=beta1
      }else{
        if(data.id$tmt[data.id$id==i]==2){
          betai=beta2
        }else{
          betai=beta
        }
      }
      betabi=adv(c(betai+c(0,b_set[data.id$id==i,])))
      misqrbetaKbKb=c(c(t(betai+c(0,b_set[data.id$id==i,]))%*%K2%*%(betai+c(0,b_set[data.id$id==i,])))^(-1/2)*K2%*%(betai+c(0,b_set[data.id$id==i,])))+
        sapply(c(1:q),function(s) sum(diag(inv.Fish_set[,,data.id$id==i]%*%attr((quad.adv(K2,betabi))^(-1/2)*linear_sca.adv(K2,betabi,s),"grad2")[-1,-1]))/2)
      
      for(j in rs){
        
        if(data.id$tmt[data.id$id==j]==1){
          betaj=beta1
        }else{
          if(data.id$tmt[data.id$id==j]==2){
            betaj=beta2
          }else{
            betaj=beta
          }
        }
        
        d1[which(rs==j)]=exp(c(alpha1*(desT%*%betaj)))
        Bbeta[which(rs==j)]=c(t(desT)%*%betaj)
        # betaKbeta[which(rs==j)]=t(betaj)%*%K2%*%betaj
        # Kbeta[which(rs==j),]=K2%*%betaj
        
        
        f=c(exp(alpha1*desT[-1]%*%b_set[data.id$id==j,]+
                  alpha2*sqrt(t(betaj+c(0,b_set[data.id$id==j,]))%*%K2%*%(betaj+c(0,b_set[data.id$id==j,])))))
        inv.Fishj=inv.Fish_set[,,data.id$id==j]
        
        
        fg1=c(f*(alpha1*desT[-1]+alpha2*c(t(betaj+c(0,b_set[data.id$id==j,]))%*%K2%*%(betaj+c(0,b_set[data.id$id==j,])))^(-1/2)*(t(K2[,-1])%*%betaj+K2[-1,-1]%*%b_set[data.id$id==j,])))
        
        fg2=f*(tcrossprod(alpha1*desT[-1]+alpha2*c(t(betaj+c(0,b_set[data.id$id==j,]))%*%K2%*%(betaj+c(0,b_set[data.id$id==j,])))^(-1/2)*(t(K2[,-1])%*%betaj+K2[-1,-1]%*%b_set[data.id$id==j,]))-
                 alpha2*c(t(betaj+c(0,b_set[data.id$id==j,]))%*%K2%*%(betaj+c(0,b_set[data.id$id==j,])))^(-3/2)*tcrossprod(t(K2[,-1])%*%betaj+K2[-1,-1]%*%b_set[data.id$id==j,])+
                 alpha2*c(t(betaj+c(0,b_set[data.id$id==j,]))%*%K2%*%(betaj+c(0,b_set[data.id$id==j,])))^(-1/2)*K2[-1,-1])
        
        h=c(sqrt(t(betaj+c(0,b_set[data.id$id==j,]))%*%K2%*%(betaj+c(0,b_set[data.id$id==j,]))))
        
        hg1=c(t(betaj+c(0,b_set[data.id$id==j,]))%*%K2%*%(betaj+c(0,b_set[data.id$id==j,])))^(-1/2)*(t(K2[,-1])%*%betaj+K2[-1,-1]%*%b_set[data.id$id==j,])
        
        hg2=-c(t(betaj+c(0,b_set[data.id$id==j,]))%*%K2%*%(betaj+c(0,b_set[data.id$id==j,])))^(-3/2)*tcrossprod(t(K2[,-1])%*%betaj+K2[-1,-1]%*%b_set[data.id$id==j,])+
          c(t(betaj+c(0,b_set[data.id$id==j,]))%*%K2%*%(betaj+c(0,b_set[data.id$id==j,])))^(-1/2)*K2[-1,-1]
        
        Exp.f[which(rs==j)]=f+sum(diag(inv.Fishj%*%fg2))/2
        
        Bb=c(t(desT[-1])%*%b_set[data.id$id==j,])
        
        Exp.fBb[which(rs==j)]=f*Bb+
          sum(diag(inv.Fishj%*%(fg2*Bb+fg1%*%t(desT[-1])+desT[-1]%*%t(fg1))))/2
        
        Exp.fBb2[which(rs==j)]=f*Bb^2+
          sum(diag(inv.Fishj%*%(Bb^2*fg2+2*Bb*(fg1%*%t(desT[-1])+desT[-1]%*%t(fg1))+2*f*desT[-1]%*%t(desT[-1]))))/2
        
        betaKb=c(t(betaj+c(0,b_set[data.id$id==j,]))%*%K2%*%(betaj+c(0,b_set[data.id$id==j,])))
        
        Exp.fbetaKb[which(rs==j)]=f*betaKb+
          sum(diag(inv.Fishj%*%(betaKb*fg2+2*fg1%*%t(t(K2[,-1])%*%betaj+K2[-1,-1]%*%b_set[data.id$id==j,])+2*(t(K2[,-1])%*%betaj+K2[-1,-1]%*%b_set[data.id$id==j,])%*%t(fg1)+2*f*K2[-1,-1])))/2
        
        sqrbetaKb[which(rs==j)]=h+sum(diag(inv.Fishj%*%hg2))/2
        
        Exp.sqrfbetaKb[which(rs==j)]=f*h+sum(diag(inv.Fishj%*%(h*fg2+fg1%*%t(hg1)+hg1%*%t(fg1)+f*hg2)))/2
        
        betabj=adv(c(betaj+c(0,b_set[data.id$id==j,])))
        
        Exp.scorebeta[which(rs==j),]=c(f*alpha2*c(t(betaj+c(0,b_set[data.id$id==j,]))%*%K2%*%(betaj+c(0,b_set[data.id$id==j,])))^(-1/2)*K2%*%(betaj+c(0,b_set[data.id$id==j,])))+
          sapply(c(1:q),function(s) sum(diag(inv.Fishj%*%attr(exp(alpha1*linear.adv(desT,adv(c(0,b_set[data.id$id==j,]))))*exp(alpha2*(quad.adv(K2,betabj))^(1/2))*(alpha2*(quad.adv(K2,betabj))^(-1/2))*linear_sca.adv(K2,betabj,s),
                                                              "grad2")[-1,-1]))/2)
        
        Exp.fBbbetaKb[which(rs==j)]=f*h*Bb+sum(diag(inv.Fishj%*%(Bb*(h*fg2+fg1%*%t(hg1)+hg1%*%t(fg1)+f*hg2)+(h*fg1+f*hg1)%*%t(desT[-1])+desT[-1]%*%t((h*fg1+f*hg1)))))/2
        
        # fb=array(c(rep(fg2,q-1)),dim=c(q-1,q-1,q-1))
        # fb1=array(c(sapply(c(1:(q-1)),function(x){fb[,,x]=fb[,,x]*b_set[mrc1.id$id==j,x]
        # fb[x,,x]=fb[x,,x]+2*fg1
        # return(fb[,,x])})),dim=c(q-1,q-1,q-1))
        # 
        # 
        # Exp.fb[which(rs==j),]=f*b_set[mrc1.id$id==j,]+apply(fb1,3,function(y) sum(diag(inv.Fishj%*%y))/2)
        # 
        # Exp.fbb[,,which(rs==j)]=f*b_set[mrc1.id$id==j,]%*%t(b_set[mrc1.id$id==j,])
        
        
      }
      
      
      
      ### score of gamma; Fisher of gamma
      d=c(exp(W[data.id$id%in%rs,]%*%gamma+alpha1*X1[data.id$id%in%rs,]%*%eta1))*d1 
      Deno=c(t(d)%*%Exp.f)
      ssgamma=W[data.id$id==i,]-colSums((d*Exp.f)*W[data.id$id%in%rs,])/Deno
      sgamma=sgamma+ssgamma
      igamma=tensor(array(apply(W[data.id$id%in%rs,],1,tcrossprod),dim=c(dim(W)[2],dim(W)[2],length(rs))),d*Exp.f,3,1)/Deno-tcrossprod(colSums((d*Exp.f)*W[data.id$id%in%rs,]))/(Deno^2)
      Igamma=Igamma+igamma
      
      #######Partial score of eta1;Partial Fisher of eta1
      sseta1=alpha1*X1[data.id$id==i,]-colSums((alpha1*d*Exp.f)*X1[data.id$id%in%rs,])/Deno
      seta1t=seta1t+sseta1
      ieta1=tensor(array(apply(X1[data.id$id%in%rs,],1,tcrossprod),dim=c(dim(X1)[2],dim(X1)[2],length(rs))),alpha1^2*d*Exp.f,3,1)/Deno-tcrossprod(colSums((alpha1*d*Exp.f)*X1[data.id$id%in%rs,]))/(Deno^2)
      Ieta1t=Ieta1t+ieta1
      
      
      ### score of alpha1; Fisher of alpha1
      ssalpha1=c(X1[data.id$id==i,]%*%eta1+desT%*%(betai+c(0,b_set[data.id$id==i,])))-t(d)%*%((Bbeta+c(X1[data.id$id%in%rs,]%*%eta1))*Exp.f+Exp.fBb)/Deno
      salpha1=salpha1+ssalpha1
      ialpha1=c(t(d)%*%((Bbeta+c(X1[data.id$id%in%rs,]%*%eta1))^2*Exp.f+Exp.fBb2+2*(Bbeta+c(X1[data.id$id%in%rs,]%*%eta1))*Exp.fBb)/Deno)-
        c(t(d)%*%((Bbeta+c(X1[data.id$id%in%rs,]%*%eta1))*Exp.f+Exp.fBb)/Deno)^2
      Ialpha1=Ialpha1+ialpha1
      
      ### score of alpha2; Fisher of alpha2
      ssalpha2=sqrbetaKb[which(rs==i)]-t(d)%*%Exp.sqrfbetaKb/Deno
      salpha2=salpha2+ssalpha2
      ialpha2=t(d)%*%Exp.fbetaKb/Deno-(t(d)%*%Exp.sqrfbetaKb/Deno)^2
      Ialpha2=Ialpha2+ialpha2
      
      ### part of score of beta; part of Fisher of beta
      ssbeta=ifelse(data.id$tmt[data.id$id==i]<3,numeric(q),alpha1*desT+alpha2*misqrbetaKbKb)-
        (c(t(d*ifelse(data.id$tmt[data.id$id==i]<3,0,1))%*%Exp.f)*alpha1*desT+colSums(d*ifelse(data.id$tmt[data.id$id==i]<3,0,1)*Exp.scorebeta))/Deno
      
      sbetat=sbetat+ssbeta
      
      ibeta=tcrossprod(ssbeta)
      
      Ibetat=Ibetat+ibeta
      
      ### part of score of beta1; part of Fisher of beta1
      ssbeta1=ifelse(data.id$tmt[data.id$id==i]>1,numeric(q),alpha1*desT+alpha2*misqrbetaKbKb)-
        (c(t(d*ifelse(data.id$tmt[data.id$id==i]>1,0,1))%*%Exp.f)*alpha1*desT+colSums(d*ifelse(data.id$tmt[data.id$id==i]>1,0,1)*Exp.scorebeta))/Deno
      
      sbeta1t=sbeta1t+ssbeta1
      
      ibeta1=tcrossprod(ssbeta1)
      
      Ibeta1t=Ibeta1t+ibeta1
      
      ### part of score of beta2; part of Fisher of beta2
      ssbeta2=ifelse(data.id$tmt[data.id$id==i]!=2,numeric(q),alpha1*desT+alpha2*misqrbetaKbKb)-
        (c(t(d*ifelse(data.id$tmt[data.id$id==i]!=2,0,1))%*%Exp.f)*alpha1*desT+colSums(d*ifelse(data.id$tmt[data.id$id==i]!=2,0,1)*Exp.scorebeta))/Deno
      
      sbeta2t=sbeta2t+ssbeta2
      
      ibeta2=tcrossprod(ssbeta2)
      
      Ibeta2t=Ibeta2t+ibeta2
      
      
      ######################################
      igamalp1=colSums(d*((Bbeta+c(X1[data.id$id%in%rs,]%*%eta1))*Exp.f+Exp.fBb)*W[data.id$id%in%rs,])/Deno-
        c(t(d)%*%((Bbeta+c(X1[data.id$id%in%rs,]%*%eta1))*Exp.f+Exp.fBb))*(colSums((d*Exp.f)*W[data.id$id%in%rs,]))/(Deno^2)
      Igamalp1=Igamalp1+igamalp1
      
      igamalp2=colSums(d*Exp.sqrfbetaKb*W[data.id$id%in%rs,])/Deno-c(t(d)%*%Exp.sqrfbetaKb)*(colSums((d*Exp.f)*W[data.id$id%in%rs,]))/(Deno^2)
      Igamalp2=Igamalp2+igamalp2
      
      ialp12=c(t(d)%*%(Exp.fBbbetaKb+(Bbeta+c(X1[data.id$id%in%rs,]%*%eta1))*Exp.sqrfbetaKb)/Deno-(t(d)%*%Exp.sqrfbetaKb)*(t(d)%*%(Exp.fBb+(Bbeta+c(X1[data.id$id%in%rs,]%*%eta1))*Exp.f))/(Deno^2))
      Ialp12=Ialp12+ialp12
      
      
      
    }
    
    ### score of beta; Fisher of beta
    Beta=matrix(0,ncol=q,nrow=m)
    for (i in 1:m){
      if(data.id$tmt[i]==1){
        Beta[i,]=beta1
      }else{
        if(data.id$tmt[i]==2){
          Beta[i,]=beta2
        }else{
          Beta[i,]=beta
        }
      }
    }
    
    res=c(data$Y-apply(des.Y*Beta[rep(c(1:m),l),],1,sum)-X1[rep(c(1:m),l),]%*%eta1-X2*eta2-apply(des.Y[,-1]*b_set[rep(c(1:m),l),],1,sum))# vector of length N
    sigma.2=1/N*(t(data$Y-apply(des.Y*Beta[rep(c(1:m),l),],1,sum)-X1[rep(c(1:m),l),]%*%eta1-X2*eta2)%*%(data$Y-apply(des.Y*Beta[rep(c(1:m),l),],1,sum)-X1[rep(c(1:m),l),]%*%eta1-X2*eta2-2*apply(des.Y[,-1]*b_set[rep(c(1:m),l),],1,sum))+
                   tr+sum(apply(des.Y[,-1]*b_set[rep(c(1:m),l),],1,sum)^2))
    sigma.2=c(sigma.2)
    
    sbeta=c(sbetat)+1/sigma.2*t(des.Y[data$tmt==3,])%*%res[data$tmt==3] 
    Ibeta=Ibetat+(-(sigma.2)^(-2)*2/N)*t(des.Y[data$tmt==3,])%*%res[data$tmt==3]%*%t(res[data$tmt==3])%*%des.Y[data$tmt==3,]+(sigma.2)^(-1)*t(des.Y[data$tmt==3,])%*%des.Y[data$tmt==3,]
    
    sbeta1=c(sbeta1t)+1/sigma.2*t(des.Y[data$tmt==1,])%*%res[data$tmt==1] 
    Ibeta1=Ibeta1t+(-(sigma.2)^(-2)*2/N)*t(des.Y[data$tmt==1,])%*%res[data$tmt==1]%*%t(res[data$tmt==1])%*%des.Y[data$tmt==1,]+(sigma.2)^(-1)*t(des.Y[data$tmt==1,])%*%des.Y[data$tmt==1,]
    
    sbeta2=c(sbeta2t)+1/sigma.2*t(des.Y[data$tmt==2,])%*%res[data$tmt==2] 
    Ibeta2=Ibeta2t+(-(sigma.2)^(-2)*2/N)*t(des.Y[data$tmt==2,])%*%res[data$tmt==2]%*%t(res[data$tmt==2])%*%des.Y[data$tmt==2,]+(sigma.2)^(-1)*t(des.Y[data$tmt==2,])%*%des.Y[data$tmt==2,]
    
    seta1=c(seta1t)+1/sigma.2*t(X1[rep(c(1:m),l),])%*%res
    Ieta1=Ieta1t+(-(sigma.2)^(-2)*2/N)*t(X1[rep(c(1:m),l),])%*%res%*%t(res)%*%X1[rep(c(1:m),l),]+(sigma.2)^(-1)*t(X1[rep(c(1:m),l),])%*%X1[rep(c(1:m),l),]
    
    
    ## Update gamma,alpha1,alpha2,eta1,eta2,beta1,beta2,beta
    Igamalp12=matrix(0,ncol=7,nrow=7)
    Igamalp12[1:5,1:5]=Igamma
    Igamalp12[6,6]=Ialpha1
    Igamalp12[7,7]=Ialpha2
    Igamalp12[1:5,6]=Igamalp1
    Igamalp12[6,1:5]=Igamalp1
    Igamalp12[1:5,7]=Igamalp2
    Igamalp12[7,1:5]=Igamalp2
    Igamalp12[6,7]=Ialp12
    Igamalp12[7,6]=Ialp12
    
    for (v in c(0,2,4,6,8,10)){
      stepsize=2^(-v)
      par.sur=c(gamma,alpha1,alpha2)+stepsize*solve(Igamalp12,c(sgamma,salpha1,salpha2))
      gammanew=par.sur[1:5]
      alpha1new=par.sur[6]
      alpha2new=par.sur[7]
      eta1new=c(eta1+stepsize*solve(Ieta1)%*%seta1)
      betanew=c(beta+stepsize*solve(Ibeta)%*%sbeta)
      beta1new=c(beta1+stepsize*solve(Ibeta1)%*%sbeta1)
      beta2new=c(beta2+stepsize*solve(Ibeta2)%*%sbeta2)
      
      Betanew=matrix(0,ncol=q,nrow=m)
      for (i in 1:m){
        if(data.id$tmt[i]==1){
          Betanew[i,]=beta1new
        }else{
          if(data.id$tmt[i]==2){
            Betanew[i,]=beta2new
          }else{
            Betanew[i,]=betanew
          }
        }
      }
      
      eta2new=c(1/c(t(X2)%*%X2)*t(X2)%*%c(data$Y-apply(des.Y*Betanew[rep(c(1:m),l),],1,sum)-X1[rep(c(1:m),l),]%*%eta1new-apply(des.Y[,-1]*b_set[rep(c(1:m),l),],1,sum)))
      Dnew=1/m*(apply(inv.Fish_set,c(1,2),sum)+t(b_set)%*%b_set)
      Dnew=diag(diag(Dnew))
      sigma2new=1/N*(t(data$Y-apply(des.Y*Betanew[rep(c(1:m),l),],1,sum)-X1[rep(c(1:m),l),]%*%eta1new-X2*eta2new)%*%(data$Y-apply(des.Y*Betanew[rep(c(1:m),l),],1,sum)-X1[rep(c(1:m),l),]%*%eta1new-X2*eta2new-2*apply(des.Y[,-1]*b_set[rep(c(1:m),l),],1,sum))+
                       tr+sum(apply(des.Y[,-1]*b_set[rep(c(1:m),l),],1,sum)^2))
      sigma2new=c(sigma2new)
      
      cc=c()
      for(t in sort(unique(data$Time[data$delta==1]))){
        rs=riskset(t) #id
        Exp.f=numeric(length(rs))
        for(j in rs){
          betanewj=Betanew[data.id$id==j,]
          f=c(exp(alpha1new*des.T[match(t,data.id$Time[data.id$delta==1]),-1]%*%b_set[data.id$id==j,]+
                    alpha2new*sqrt(t(betanewj+c(0,b_set[data.id$id==j,]))%*%K2%*%(betanewj+c(0,b_set[data.id$id==j,])))))
          
          inv.Fishj=inv.Fish_set[,,data.id$id==j]
          fg2=f*(tcrossprod(alpha1new*desT[-1]+alpha2new*c(t(betanewj+c(0,b_set[data.id$id==j,]))%*%K2%*%(betanewj+c(0,b_set[data.id$id==j,])))^(-1/2)*(t(K2[,-1])%*%betanewj+K2[-1,-1]%*%b_set[data.id$id==j,]))-
                   alpha2new*c(t(betanewj+c(0,b_set[data.id$id==j,]))%*%K2%*%(betanewj+c(0,b_set[data.id$id==j,])))^(-3/2)*tcrossprod(t(K2[,-1])%*%betanewj+K2[-1,-1]%*%b_set[data.id$id==j,])+
                   alpha2new*c(t(betanewj+c(0,b_set[data.id$id==j,]))%*%K2%*%(betanewj+c(0,b_set[data.id$id==j,])))^(-1/2)*K2[-1,-1])
          
          Exp.f[which(rs==j)]=f+sum(diag(inv.Fishj%*%fg2))/2
        }
        xbeta=apply(matrix(rep(des.T[match(t,data.id$Time[data.id$delta==1]),],length(rs)),ncol=q,byrow=TRUE)*Betanew[data.id$id%in%rs,],1,sum)#collection of Xj(t)%*%beta1
        cc=c(cc,length(data.id$id[(data.id$Time==t)&(data.id$delta==1)])/c(t(exp(W[data.id$id%in%rs,]%*%gammanew+alpha1new*(xbeta+X1[data.id$id%in%rs,]%*%eta1new)))%*%Exp.f))
      }
      
      cumbase=data.frame(hazard=0,time=c(0,sort(unique(data$Time[data$delta==1]))))
      cumbase[,1]=cumsum(c(0,cc))
      

      
      ### Log-likelihood 
      set.seed(1)
      L=2000
      log.p.Y=numeric(m)
      log.s_set=matrix(0,nrow=m,ncol=L)
      ZTimeb=matrix(0,nrow=sum(delta==1),ncol=L) ##Z_Time%*%b
      bKb=matrix(0,nrow=sum(delta==1),ncol=L)
      rmc=rmvnorm(L,mean=rep(0,q-1))
      for(i in data.id$id){
        Zi=des.Y[data$id==i,-1]
        Xi=cbind(matrix(rep(X1[data.id$id==i,],l[data.id$id==i]),ncol=dim(X1)[2],byrow=T),
                 des.Y[data$id==i,],data$WC[data$id==i])
        Yi=data$Y[data$id==i]
        
        chol_invSig=chol(solve(Dnew)+c(sigma2new^(-1))*t(Zi)%*%Zi)
        Q=solve(chol_invSig)
        if(data.id$tmt[data.id$id==i]==1){
          betai=beta1new
        }else{
          if(data.id$tmt[data.id$id==i]==2){
            betai=beta2new
          }else{
            betai=betanew
          }
        }
        
        mu=c(Q%*%t(Q)%*%(t(Zi)/c(sigma2new))%*%(Yi-Xi%*%c(eta1new,betai,eta2new))) #E(bi|Yi,theta),initial value of bi
        
        Sigma=solve(solve(Dnew)+c(sigma2new^(-1))*t(Zi)%*%Zi)
        
        log.p.Y[data.id$id==i]=dmvnorm(Yi,mean=c(Xi%*%c(eta1new,betai,eta2new)),sigma=sigma2new*diag(length(Yi))+Zi%*%Dnew%*%t(Zi),log=TRUE)
        rmci=t(mu+t(chol(Sigma))%*%t(rmc))## matrix with ncol=q,nrow=L, MC sample from bi|yi
        if(delta[data.id$id==i]){
          ZTimeb[match(i,data.id$id[data.id$delta==1]),]=des.T[match(i,data.id$id[data.id$delta==1]),-1]%*%t(rmci)
          bKb[match(i,data.id$id[data.id$delta==1]),]=colSums((t(cbind(numeric(L),rmci))+betai)*(K.2[,,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%(t(cbind(numeric(L),rmci))+betai)))
        }
        
        a=unique(data.id$Time[(data.id$Time<=Time[data.id$id==i])&(data.id$delta==1)])
        if(length(a)>0){
          if(length(a)==1){
            btKb=colSums((t(cbind(numeric(L),rmci))+betai)*(K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(t(cbind(numeric(L),rmci))+betai)))
          }else{
            btKb=apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) colSums((t(cbind(numeric(L),rmci))+betai)*(y%*%(t(cbind(numeric(L),rmci))+betai))))
          }
          
          log.s_set[data.id$id==i,]=apply(sapply(a,lambda0)*exp(c(gammanew%*%W[data.id$id==i,]))*
                                            exp(alpha1new*(c(X1[data.id$id==i,]%*%eta1new)+des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(t(cbind(numeric(L),rmci))+betai))+
                                                  alpha2new*t(sqrt(btKb))),2,sum)
          
        }
        
      }
      
      log.hazard=matrix(0,ncol=L,nrow=m)
      log.hazard[which(delta==1),]=log(sapply(Time[delta==1],lambda0))+c(W[delta==1,]%*%gammanew)+alpha1new*(c(X1[delta==1,]%*%eta1new+apply(des.T*Beta[delta==1,],1,sum))+ZTimeb)+
        alpha2new*sqrt(bKb)
      log.p.Tb=log.hazard-log.s_set
      p.Tb=exp(log.p.Tb)
      logLikmc=sum(log.p.Y+log(rowMeans(p.Tb)))
      
      if(logLikmc>logLik_set[k-1]){
        print(v)
        break
      }
    }
    
    gamma=gammanew
    alpha1=alpha1new
    alpha2=alpha2new
    eta1=eta1new
    eta2=eta2new
    D=Dnew
    beta=betanew
    beta1=beta1new
    beta2=beta2new
    sigma2=sigma2new
    
    logLik_set[k]=logLikmc
    gamma_set[k,]=gamma
    alpha1_set[k]=alpha1
    alpha2_set[k]=alpha2
    eta1_set[k,]=eta1
    eta2_set[k]=eta2
    beta_set[k,]=beta
    beta1_set[k,]=beta1
    beta2_set[k,]=beta2
    sigma2_set[k]=sigma2
    D_set[,,k]=D
    #Q.fun_set[k]=Q.fun.up
    
    if((abs(logLikmc-logLik_set[k-1])/abs(logLik_set[k-1]))<10^(-5)){
      if(logLikmc>logLik_set[k-1]){
        return(c(k,logLik_set[k],gamma_set[k,],alpha1_set[k],alpha2_set[k],beta_set[k,],beta1_set[k,],beta2_set[k,],
                          eta1_set[k,],eta2_set[k],sigma2_set[k],diag(D_set[,,k])))
      }else{
        return(c(k-1,logLik_set[k-1],gamma_set[k-1,],alpha1_set[k-1],alpha2_set[k-1],beta_set[k-1,],beta1_set[k-1,],beta2_set[k-1,],
                          eta1_set[k-1,],eta2_set[k-1],sigma2_set[k-1],diag(D_set[,,k-1])))
      }
    }
  }
  
  
  
}


