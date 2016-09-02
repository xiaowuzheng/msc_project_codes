#reward function 1 + strategy 1 (static)
#upper and lower bound
ub<-function(x){
  rep(100000*0.5,10)/abs(Predict.price.model3(x)[1,])
}
lb<-function(x){
  -rep(100000*0.5,10)/abs(Predict.price.model3(x)[1,])
}
#----------------correct----------------------------------------------------------#
opt.initial.R1.S1<-function(gam,iw,time,samind){
  #objective function
  eval_f<-function(x){
    fun<-0
    for(i in 1:10){
      fun<-fun-x[i]*Mean.price.model3[time,i]+gam*(x[i]^2)*Var.price.model3[time,i]
    }
    g<-c()
    for(i in 1:10){
      g[i]<- -Mean.price.model3[time,i]+gam*2*x[i]*Var.price.model3[time,i]
    }
    return(list("objective"=fun,"gradient"=g))
  }
  
  #constriant functions
  eval_g_eq<-function(x){
    constr<-c(x%*%Predict.price.model3(samind)[1,]-iw)
    grad<-Predict.price.model3(samind)[1,]
    return(list("constraints"=constr,"jacobian"=grad))
  }
  
  #set.seed(123)
  x0<-runif(10,lb(samind),ub(samind))
  
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 1000,
                "local_opts" = local_opts )
  res<-nloptr(x0=x0,
              eval_f=eval_f,
              lb=lb(samind),
              ub=ub(samind),
              eval_g_eq=eval_g_eq,
              opts=opts)
  res$solution
}

#-------------for the consecutive time intervals-------------#
opt.consecutive.R1.S1<-function(gam,w.prev,L.length,intervalind,samind){
  #end of the interval
  time<-L.length*intervalind
  #objective function
  eval_f<-function(x){
    fun<-0
    for(i in 1:10){
      fun<-fun-x[i]*Mean.price.model3[time,i]+gam*x[i]^2*Var.price.model3[time,i]
    }
    g<-c()
    for(i in 1:10){
      g[i]<- -Mean.price.model3[time,i]+gam*2*x[i]*Var.price.model3[time,i]
    }
    return(list("objective"=fun,"gradient"=g))
  }
  #constriant functions
  eval_g_eq<-function(x){
    constr<-0
    pri<-Predict.price.model3(samind)[L.length*(intervalind-1)+1,]
    for(i in 1:10){
      if(x[i]>w.prev[i]){
        constr<-constr+(1+0.003)*(x[i]-w.prev[i])*pri[i]
      }
      else{constr<-constr+(1-0.003)*(x[i]-w.prev[i])*pri[i]
      }
    }
    g<-c()
    for(i in 1:10){
      if(x[i]>w.prev[i]){
        g[i]<-(1+0.003)*pri[i]
      }
      else{g[i]<-(1-0.003)*pri[i]}
    }
    return(list("constraints"=constr,"jacobian"=g))
  }
  #set.seed(789)
  x0<-runif(10,lb(samind),ub(samind))
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 1000,
                "local_opts" = local_opts )
  res<-nloptr(x0=x0,
              eval_f=eval_f,
              lb=lb(samind),
              ub=ub(samind),
              eval_g_eq=eval_g_eq,
              opts=opts)
  res$solution
}


#-----------------------------------------------------------------------------------------------------------------#
RES.model3.R1.S1<-function(ntime,gammm,intervallength,sam){
  #RES.list<-list()
  res<-matrix(0,nrow=ntime,ncol=10)
  res[1,]<-opt.initial.R1.S1(gammm,100000,intervallength,sam)
  for(i in 2:ntime){
    res[i,]<- opt.consecutive.R1.S1(gammm,res[i-1,],intervallength,i,sam)
  }
  #RES.list[[1]]<-res
  #tc<-0.003*sapply(2:ntime,function(n){sum(abs(res[n,]-res[n-1,])%*%Predict.price.model3(sam)[intervallength*(n-1)+1,])})
  #RES.list[[2]]<-tc
  #RES.list
  res
}
#----------------final answer for the optimised weight matrix
#just input 
AW.model3.R1.S1<-function(ga,totalinte,intelength,whichsam){
  #res.list<-list()
  res<-matrix(0,nrow=100,ncol=10)
  res.inter<-RES.model3.R1.S1(totalinte,ga,intelength,whichsam)
  for(i1 in 1:totalinte){
    for(i2 in ((i1-1)*intelength+1):(i1*intelength)){
      #res[i2,]<-res.inter[[1]][i1,]
      res[i2,]<-res.inter[i1,]
    }
  }
  #res.list[[1]]<-res
  #res.list[[2]]<-res.inter[[2]]
  #res.list
  res
}


#-------------------write a function, only have to input frequency L,interval length,gamma and total number of samples
#-------------------and the return is the sequence of average value of the portfolio, ie sequence of length 100
#-----------answer to the value of the portfolio
AV.model3.R1.S1<-function(L,paragamma,nosample){
  #each row is a series of value of portfolio over 100 trading days for one sample
  res.matrix<-matrix(0,nrow=nosample,ncol=100)
  for (i in 1:nosample){
    AW.sam<-AW.model3.R1.S1(paragamma,100/L,L,i)
    tc<-c(0,0.003*sapply(2:100,function(n){sum(abs(AW.sam[n,]-AW.sam[n-1,])%*%DATA.test[n,2:11])}))
    res.matrix[i,]<-sapply(1:100,function(x){AW.sam[x,]%*%DATA.test[x,2:11]})
    res.matrix[i,]<-res.matrix[i,]-tc
    print(i)
  }
  res.matrix
}

#------------------------------------------------L=20---------------------------------------------------#
#--------------gamma=10-----------------#
model3.L20.R1.S1.gam1<-AV.model3.R1.S1(L=20,paragamma=10,nosample=1000)
#--------------gamma=5-----------------#
model3.L20.R1.S1.gam2<-AV.model3.R1.S1(L=20,paragamma=5,nosample=1000)
#--------------gamma=1-----------------#
model3.L20.R1.S1.gam3<-AV.model3.R1.S1(L=20,paragamma=1,nosample=1000)
#--------------gamma=0.1-----------------#
model3.L20.R1.S1.gam4<-AV.model3.R1.S1(L=20,paragamma=0.1,nosample=1000)
#--------------gamma=0.05-----------------#
model3.L20.R1.S1.gam5<-AV.model3.R1.S1(L=20,paragamma=0.05,nosample=1000)
#--------------gamma=0.01-----------------#
model3.L20.R1.S1.gam6<-AV.model3.R1.S1(L=20,paragamma=0.01,nosample=1000)
#--------------gamma=0.005-----------------#
model3.L20.R1.S1.gam7<-AV.model3.R1.S1(L=20,paragamma=0.005,nosample=1000)
#--------------gamma=0.001-----------------#
model3.L20.R1.S1.gam8<-AV.model3.R1.S1(L=20,paragamma=0.001,nosample=1000)
#--------------gamma=0.0005-----------------#
model3.L20.R1.S1.gam9<-AV.model3.R1.S1(L=20,paragamma=0.0005,nosample=1000)
#--------------gamma=0.0004-----------------#
model3.L20.R1.S1.gam10<-AV.model3.R1.S1(L=20,paragamma=0.0004,nosample=1000)
#--------------gamma=0.0003-----------------#
model3.L20.R1.S1.gam11<-AV.model3.R1.S1(L=20,paragamma=0.0003,nosample=1000)
#--------------gamma=0.0002-----------------#
model3.L20.R1.S1.gam12<-AV.model3.R1.S1(L=20,paragamma=0.0002,nosample=1000)
#--------------gamma=0.0001-----------------#
model3.L20.R1.S1.gam13<-AV.model3.R1.S1(L=20,paragamma=0.0001,nosample=1000)
#--------------gamma=0.00008-----------------#
model3.L20.R1.S1.gam14<-AV.model3.R1.S1(L=20,paragamma=0.00008,nosample=1000)
#--------------gamma=0.00007-----------------#
model3.L20.R1.S1.gam15<-AV.model3.R1.S1(L=20,paragamma=0.00007,nosample=1000)
#--------------gamma=0.00005-----------------#
model3.L20.R1.S1.gam16<-AV.model3.R1.S1(L=20,paragamma=0.00005,nosample=1000)
#--------------gamma=0.00003-----------------#
model3.L20.R1.S1.gam17<-AV.model3.R1.S1(L=20,paragamma=0.00003,nosample=1000)
#--------------gamma=0.00001-----------------#
model3.L20.R1.S1.gam18<-AV.model3.R1.S1(L=20,paragamma=0.00001,nosample=1000)
#----------------------------------------------------------------------------------------#
#---------------L=20 performance--------------------#
per.model3.L20.R1.S1.gam1<-performance(model3.L20.R1.S1.gam1)
per.model3.L20.R1.S1.gam2<-performance(model3.L20.R1.S1.gam2)
per.model3.L20.R1.S1.gam3<-performance(model3.L20.R1.S1.gam3)
per.model3.L20.R1.S1.gam4<-performance(model3.L20.R1.S1.gam4)
per.model3.L20.R1.S1.gam5<-performance(model3.L20.R1.S1.gam5)
per.model3.L20.R1.S1.gam6<-performance(model3.L20.R1.S1.gam6)
per.model3.L20.R1.S1.gam7<-performance(model3.L20.R1.S1.gam7)
per.model3.L20.R1.S1.gam8<-performance(model3.L20.R1.S1.gam8)
per.model3.L20.R1.S1.gam9<-performance(model3.L20.R1.S1.gam9)
per.model3.L20.R1.S1.gam10<-performance(model3.L20.R1.S1.gam10)
per.model3.L20.R1.S1.gam11<-performance(model3.L20.R1.S1.gam11)
per.model3.L20.R1.S1.gam12<-performance(model3.L20.R1.S1.gam12)
per.model3.L20.R1.S1.gam13<-performance(model3.L20.R1.S1.gam13)
per.model3.L20.R1.S1.gam14<-performance(model3.L20.R1.S1.gam14)
per.model3.L20.R1.S1.gam15<-performance(model3.L20.R1.S1.gam15)
per.model3.L20.R1.S1.gam16<-performance(model3.L20.R1.S1.gam16)
per.model3.L20.R1.S1.gam17<-performance(model3.L20.R1.S1.gam17)
per.model3.L20.R1.S1.gam18<-performance(model3.L20.R1.S1.gam18)
##############profit and loss function#####################################
plot(seq(1,100,by=1),colMeans(model3.L20.R1.S1.gam1)/100000-rep(1,100),lty=2,"l",col=1,xlab="time",ylab="profit and loss",ylim=c(-0.03,0.65),main="Reward 1 Static 1 profit and loss plot")
lines(seq(1,100,by=1),colMeans(model3.L20.R1.S1.gam2)/100000-rep(1,100),lty=2,col=2)
lines(seq(1,100,by=1),colMeans(model3.L20.R1.S1.gam3)/100000-rep(1,100),lty=2,col=3)
lines(seq(1,100,by=1),colMeans(model3.L20.R1.S1.gam4)/100000-rep(1,100),lty=2,col=4)
lines(seq(1,100,by=1),colMeans(model3.L20.R1.S1.gam5)/100000-rep(1,100),lty=2,col="forestgreen")
lines(seq(1,100,by=1),colMeans(model3.L20.R1.S1.gam6)/100000-rep(1,100),lty=2,col=6)
lines(seq(1,100,by=1),colMeans(model3.L20.R1.S1.gam7)/100000-rep(1,100),lty=3,col="aquamarine2")
lines(seq(1,100,by=1),colMeans(model3.L20.R1.S1.gam8)/100000-rep(1,100),lty=2,col=8)
lines(seq(1,100,by=1),colMeans(model3.L20.R1.S1.gam9)/100000-rep(1,100),lty=2,col="tomato1")
lines(seq(1,100,by=1),colMeans(model3.L20.R1.S1.gam10)/100000-rep(1,100),lty=3,col="steelblue3")
legend("topleft",c(expression(gamma==10),expression(gamma==1),expression(gamma==0.1),expression(gamma==0.01),expression(gamma==0.005),expression(gamma==0.003),expression(gamma==0.001),expression(gamma==0.0005),expression(gamma==0.0001),expression(gamma==0.00001)),col=c(1:4,"forestgreen",6,"aquamarine2",8,"tomato1","steelblue3"),lty=rep(2,10))

#--------------------------------#final wealth density plot for L=20------------------------------------------#
plot(density(model3.L20.R1.S1.gam1[,100]),lty=1,col=1,xlab=expression(V[T]),ylab="density",main="Reward 1 Static 1 Density Plot",xlim=c(115000,210000),ylim=c(0,0.00018))
abline(v=mean(model3.L20.R1.S1.gam1[,100]),lty=2,col=1)
lines(density(model3.L20.R1.S1.gam2[,100]),lty=1,col=2)
abline(v=mean(model3.L20.R1.S1.gam2[,100]),lty=2,col=2)
lines(density(model3.L20.R1.S1.gam3[,100]),lty=1,col=3)
abline(v=mean(model3.L20.R1.S1.gam3[,100]),lty=2,col=3)
lines(density(model3.L20.R1.S1.gam4[,100]),lty=1,col=4)
abline(v=mean(model3.L20.R1.S1.gam4[,100]),lty=2,col=4)
lines(density(model3.L20.R1.S1.gam5[,100]),lty=1,col="forestgreen")
abline(v=mean(model3.L20.R1.S1.gam5[,100]),lty=2,col="forestgreen")
lines(density(model3.L20.R1.S1.gam6[,100]),lty=1,col=6)
abline(v=mean(model3.L20.R1.S1.gam6[,100]),lty=2,col=6)
lines(density(model3.L20.R1.S1.gam7[,100]),lty=1,col="aquamarine2")
abline(v=mean(model3.L20.R1.S1.gam7[,100]),lty=2,col="aquamarine2")
lines(density(model3.L20.R1.S1.gam8[,100]),lty=1,col=8)
abline(v=mean(model3.L20.R1.S1.gam8[,100]),lty=2,col=8)
lines(density(model3.L20.R1.S1.gam9[,100]),lty=1,col="tomato1")
abline(v=mean(model3.L20.R1.S1.gam9[,100]),lty=2,col="tomato1")
lines(density(model3.L20.R1.S1.gam10[,100]),lty=1,col="steelblue3")
abline(v=mean(model3.L20.R1.S1.gam10[,100]),lty=2,col="steelblue3")
legend("topright",c(expression(gamma==10),expression(gamma==1),expression(gamma==0.1),expression(gamma==0.01),expression(gamma==0.005),expression(gamma==0.003),expression(gamma==0.001),expression(gamma==0.0005),expression(gamma==0.0001),expression(gamma==0.00001)),col=c(1:4,"forestgreen",6,"aquamarine2",8,"tomato1","steelblue3"),lty=rep(2,10))
#----------------------------------------------------------------------------------------------------------#

x.R1.S1.gam1<-c(per.model3.L20.R1.S1.gam9[[2]],per.model3.L20.R1.S1.gam10[[2]],per.model3.L20.R1.S1.gam11[[2]],per.model3.L20.R1.S1.gam12[[2]],per.model3.L20.R1.S1.gam13[[2]],per.model3.L20.R1.S1.gam14[[2]],per.model3.L20.R1.S1.gam15[[2]],per.model3.L20.R1.S1.gam16[[2]],per.model3.L20.R1.S1.gam17[[2]],per.model3.L20.R1.S1.gam18[[2]])
y.R1.S1.gam1<-c(per.model3.L20.R1.S1.gam9[[1]],per.model3.L20.R1.S1.gam10[[1]],per.model3.L20.R1.S1.gam11[[1]],per.model3.L20.R1.S1.gam12[[1]],per.model3.L20.R1.S1.gam13[[1]],per.model3.L20.R1.S1.gam14[[1]],per.model3.L20.R1.S1.gam15[[1]],per.model3.L20.R1.S1.gam16[[1]],per.model3.L20.R1.S1.gam17[[1]],per.model3.L20.R1.S1.gam18[[1]])
smoothingSpline.R1.S1.gam1 <- smooth.spline(x.R1.S1.gam1, y.R1.S1.gam1, spar=0.35)
lines(x.R1.S1.gam1,y.R1.S1.gam1,xlab="Maximum drawdown",ylab="Sharpe Ratio",main="Reward 1 Strategy 2 SR vs MDD")
lines(smoothingSpline.R1.S1.gam1)
plot(seq(1,10,1),x.R1.S1.gam1,"l",xlab="gamma index",ylab="Maximum Drawdown")
plot(seq(1,10,1),y.R1.S1.gam1,"l",xlab="gamma index",ylab="Sharpe Ratio")
#smoothingSpline.R1.S1.gam1 <- smooth.spline(x.R1.S1.gam1, y.R1.S1.gam1, spar=0.35)
#lines(smoothingSpline.R1.S1.gam1)
text(x=0.0665, y =0.0318, labels = expression(gamma==10) , cex = 1,col=2)
text(x=0.0664, y =0.0328, labels = expression(gamma==1) , cex = 1,col=2)
text(x=0.0663, y =0.0331, labels = expression(theta==0.1) , cex = 1,col=2)
text(x=0.0667, y =0.0361, labels = expression(theta==0.01) , cex = 1,col=2)
text(x=0.0662, y =0.0393, labels = expression(theta==0.005) , cex = 1,col=2)
text(x=0.0656, y =0.0433, labels = expression(theta==0.003) , cex = 1,col=2)
text(x=0.0620, y =0.0609, labels = expression(theta==0.001) , cex = 1,col=2)
text(x=0.0592, y =0.0768, labels = expression(theta==0.0005) , cex = 1,col=2)
text(x=0.0572, y =0.0993, labels = expression(theta==0.0001) , cex = 1,col=2)
text(x=0.0586, y =0.0983, labels = expression(theta==0.00001) , cex = 1,col=2)


#--------------------------------ggplot for profit and loss----------------------------#
pandl.R1.S1<-c(colMeans(model3.L20.R1.S1.gam9)/100000-rep(1,100),colMeans(model3.L20.R1.S1.gam10)/100000-rep(1,100),colMeans(model3.L20.R1.S1.gam11)/100000-rep(1,100),colMeans(model3.L20.R1.S1.gam12)/100000-rep(1,100),colMeans(model3.L20.R1.S1.gam13)/100000-rep(1,100),colMeans(model3.L20.R1.S1.gam14)/100000-rep(1,100),colMeans(model3.L20.R1.S1.gam15)/100000-rep(1,100),colMeans(model3.L20.R1.S1.gam16)/100000-rep(1,100),colMeans(model3.L20.R1.S1.gam17)/100000-rep(1,100),colMeans(model3.L20.R1.S1.gam18)/100000-rep(1,100))
para.R1.S1<-c(rep("0.0005",100),rep("0.0004",100),rep("0.0003",100),rep("0.0002",100),rep("0.0001",100),rep("0.00008",100),rep("0.00007",100),rep("0.00005",100),rep("0.00003",100),rep("0.00001",100))
df.R1.S1<-data.frame(time=rep(1:100,10),return=(pandl.R1.S1),gamma=para.R1.S1)
ggplot(data=df.R1.S1,aes(x=time,y=return))+geom_line(aes(colour=gamma))+ scale_colour_hue()

#---------------------------------ggplot for density of final wealth----------------------------#
den.R1.S1<-c(model3.L20.R1.S1.gam9[,100],model3.L20.R1.S1.gam10[,100],model3.L20.R1.S1.gam11[,100],model3.L20.R1.S1.gam12[,100],model3.L20.R1.S1.gam13[,100],model3.L20.R1.S1.gam14[,100],model3.L20.R1.S1.gam15[,100],model3.L20.R1.S1.gam16[,100],model3.L20.R1.S1.gam17[,100],model3.L20.R1.S1.gam18[,100])
para.sam.R1.S1<-c(rep("0.0005",1000),rep("0.0004",1000),rep("0.0003",1000),rep("0.0002",1000),rep("0.0001",1000),rep("0.00008",1000),rep("0.00007",1000),rep("0.00005",1000),rep("0.00003",1000),rep("0.00001",1000))
df.den.R1.S1<-data.frame(wealth=den.R1.S1,gamma=para.sam.R1.S1)
mu<-ddply(df.den.R1.S1,"gamma",summarise,grp.mean=mean(wealth))
ggplot(data=df.den.R1.S1,aes(x=wealth,color=gamma))+stat_density(geom='line',position = 'identity')+coord_cartesian(xlim=c(110000,215000))+geom_vline(data=mu,aes(xintercept=grp.mean,color=gamma),linetype="dashed")
