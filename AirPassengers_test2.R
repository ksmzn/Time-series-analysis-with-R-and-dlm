###################################################
#
# AirPassengersデータはある航空会社の1949年から1960年までの国際旅客数に関する時系列データ
#
###################################################

install.packages("TTR")
install.packages("tseries")
install.packages("lawstat")
library(dlm)
library(TTR)
library(tseries)
library(lawstat)
plot(decompose(data))

#AirPassengers#
data=log(AirPassengers)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE)) 
plot(data,main="Time seriese plot")
lines(SMA(data, 12),lty="dotdash")
lines(SMA(data,12*4),lty="longdash")
legend("bottomright", col = c("black", "black","black"),
       lty = c("solid", "dotdash","longdash"), pch = c(NA, NA,NA), bty = "n",
       legend = c("log(AirPassengers)", "12 simple moving average","12 * 4 simple moving average"))
acf(data,main="acf")
pacf(data,main="pacf")
dev.copy2eps(file="AirPassengersplot.eps")

#model1 ローカル線形レベル＋季節成分
par(mfrow=c(1,1))
freq=12
if(0){
mod1 <- function(u) {
  trend <- dlmModPoly(order=2,
                      dV = exp(u[1]),
                      dW = exp(u[2:3]),
                      m0 = c(level0,slope0),
                      C0 = u[4] * diag(2)
  )
  mod <- dlmModTrig(s=12,q=3,dV = exp(u[4]),dW = exp(u[5]))
  return(trend+mod)
}
init_par=log(var(data))
outMLE1 <- dlmMLE(data,dlmMLE(data,c(0,.1,1,5,.1), mod1)$par  , mod1)
}
mod1 <- function(u) {
  trend <- dlmModPoly(order=2,
                      dV = exp(u[1]),
                      dW = exp(u[2:3]),
                      #m0 = c(level0,slope0),
                      #C0 = exp(u[5]) * diag(2)
  )
  mod <- dlmModSeas(frequency = freq,dV = 0,dW = c(exp(u[4]),rep(0,freq - 2)))
  return(trend+mod)
}
init_par=log(var(data))
outMLE1 <- dlmMLE(data,dlmMLE(data,c(0,1,1,1,0), mod1)$par  , mod1)
outMLE1$par
##  [1]  -8.909802  -7.259852 -26.922506  -9.698092  -5.348046
dlm1 <- mod1(outMLE1$par)
dlm1$FF
res1 <- dlmFilter(data,dlm1)
ress1<- dlmSmooth(data,dlm1)
par(mfrow=c(2,1))
plot(data,xlab="",ylab="value",type="o",col=("darkgrey"))
lines(res1$f,col=2,lty="dotdash")       
legend("bottomright", col = c("black","black"),
       lty = c("solid", "dotdash"), pch = c(NA,NA), bty = "n",
       legend = c("log(AirPassengers)", "filtering value"))
ts.plot(cbind(data,ress1$s[-1,1]+ress1$s[-1,3]+ress1$s[-1,5]),col=c(1,4),lty=c("solid","dotdash"))
legend("bottomright", col = c(1, 4),
       lty = c("solid", "dotdash"), pch = c(NA,NA), bty = "n",
       legend = c("log(AirPassengers)", "smoothing value"))       
par(mfrow=c(4,1))
plot(data)
plot.ts(ts(ress1$s[-1,1]))
plot.ts(ts(ress1$s[-1,3]))
plot.ts(ts(ress1$s[-1,5]))
par(mfrow=c(1,1))

#model2 AR(12)
par(mfrow=c(1,1))
slope0=mean(diff(data))
mod2 <- function(u) {
  trend <- dlmModPoly(order=2,
                      dV = 1e-7,
                      dW = exp(u[6:7]),
                      m0 = c(data[1],slope0),
                      C0 = 2*diag(2)
  )
  ar<- dlmModARMA(
    ar=c(u[1:2],rep(0,8),u[3:4]),
    sigma2 = exp(u[5]),
    #C0 = 2*diag(11))
    #ar=c(u[3:4]),
    #sigma2 = exp(u[5]),
    #C0 = 10^(-10)*diag(11)
    ) 
  return(trend+ar)
}
if(0){
outMLE2 <- dlmMLE(data,dlmMLE(data,
                              c(ar(data, AIC=T, order.max=12)$ar[1:3],
                                ar(data, AIC=T, order.max=12)$ar[10:12],
                                log(var(ar(data, AIC=F, order.max=12)$resid[-c(1:12)]))
                                ,init_par
                              ),
                              mod2)$par,mod2)
}
c(ar(data, AIC=T, order.max=12)$ar[1:2],
  ar(data, AIC=T, order.max=12)$ar[11:12],
  log(var(ar(data, AIC=F, order.max=12)$resid[-c(1:12)]))
  ,init_par,init_par
),
outMLE2 <- dlmMLE(data,dlmMLE(data,
                              c(1,1,1,1,1,1,1),
                              mod2,method="Nelder-Mead")$par,mod2,method="BFGS")

outMLE2
# [1]  1.16911527 -0.39201165 -0.04414805  0.11077001  0.25857009 -0.34678019 -1.53413076 -8.43936547
dlm2 <- mod2(outMLE2$par)
dlm2$F
res2 <- dlmFilter(data,dlm2)
par(mfrow=c(1,1)) 
plot(data,xlab="",ylab="value",type="o",col=("darkgrey"))
lines(res2$f,lty="dotdash")
dev.copy2eps(file="AirPassengersmod2filter.eps",width=10)
ress2<- dlmSmooth(data,dlm2)
par(mfrow=c(2,1), mar=c(1.5,4,0.5,0.5) + 0.1, cex=0.6)
#ts.plot(cbind(data,ress2$s[-1,1]),col=c("darkgrey",4))
plot.ts(dropFirst(ress2$s)[,1],ylab="level+trend",ann=F,yax.flip=TRUE)
legend("bottomright",col=c(1),lty="solid",legend=c("level+trend"))
plot.ts(dropFirst(ress2$s)[,3],main="",ann=F,yax.flip=TRUE, type="o")
legend("bottomright",col=c(1),lty="solid",legend=c("seasonal"))
dev.copy2eps(file="AirPassengersmod2smooth.eps",width=10)

#残差1
par(mfrow=c(2, 1))
zansa1=residuals(res1)$res
qqnorm(zansa1)
qqline(zansa1)
acf(zansa1)
dev.copy2eps(file="AirPassengerszansamod1.eps",width=10)

shapiro.test(zansa1)
ks.test(zansa1, "pnorm", mean=mean(zansa1), sd=sqrt(var(zansa1)))
jarque.bera.test(zansa1)
rjb.test(zansa1)
Box.test(zansa1,type= "Ljung-Box")
Box.test(zansa1,type= "Box-Pierce")
#MAD
sum(abs(zansa1))/length(zansa1)
#MSE
sum(zansa1^2)/length(zansa1)
#MAPE
sum(abs(zansa1)/data)/length(zansa1)
#U
sum((data-res1$f)^2)/sum((data[-1]-data[-length(data)])^2)




#残差2
zansa2=residuals(res2)$res
qqnorm(zansa2, main="")
qqline(zansa2)
acf(zansa2)
dev.copy2eps(file="AirPassengerszansamod2.eps",width=10)
#pacf(zansa2)
shapiro.test(zansa2)
ks.test(zansa2, "pnorm", mean=mean(zansa2), sd=sqrt(var(zansa2)))
jarque.bera.test(zansa2)
rjb.test(zansa2)
Box.test(zansa2,type= "Ljung-Box")
Box.test(zansa2,type= "Box-Pierce")

#MAD
sum(abs(zansa2))/length(zansa2)
#MSE
sum(zansa2^2)/length(zansa2)
#MAPE
sum(abs(zansa2)/data)/length(zansa2)
#U
sum((data-res2$f)^2)/sum((data[-1]-data[-length(data)])^2)


#MAD
mean(abs(res1$f-data))
mean(abs(res2$f-data))
#MSE
mean((res1$f-data)^2)
mean((res2$f-data)^2)
#MAPE
mean(abs(res1$f-data)/data)
mean(abs(res2$f-data)/data)
#U
sqrt(sum((res1$f-data)[-(1:5)]^2)/
       sum(diff(data[-(1:4)])^2))
sqrt(sum((res2$f-data)[-(1:5)]^2)/
       sum(diff(data[-(1:4)])^2))

#予測
par(mfrow=c(2,1),mar=c(2,2,2,0.5))
fore1 <- dlmForecast(res1, nAhead=12)
ciTheory1 <- (outer(sapply(fore1$Q, FUN=function(x) sqrt(diag(x))), qnorm(c(0.25,0.75))) +
                as.vector(t(fore1$f)))
plot(data,xlim=c(1821,1944),type="o",col="darkgrey",main="model1",ylab="")
lines(cbind(ciTheory1,fore1$f[,1])[,1], lty="dotdash")
lines(cbind(ciTheory1,fore1$f[,1])[,2], lty="dotdash")
lines(cbind(ciTheory1,fore1$f[,1])[,3], lty="solid")
legend("bottomright", col = c("black","black"),
       lty = c("solid", "dotdash"), pch = c(NA,NA), bty = "n",
       legend = c("one-step-ahead forecast","95% Confidence Interval"))
fore2 <- dlmForecast(res2, nAhead=12)
ciTheory2 <- (outer(sapply(fore2$Q, FUN=function(x) sqrt(diag(x))), qnorm(c(0.25,0.75))) +
                as.vector(t(fore2$f)))
plot(data,xlim=c(1821,1944),type="o",col="darkgrey",main="model2")
lines(cbind(ciTheory2,fore2$f[,1])[,1], lty="dotdash")
lines(cbind(ciTheory2,fore2$f[,1])[,2], lty="dotdash")
lines(cbind(ciTheory2,fore2$f[,1])[,3], lty="solid")
legend("bottomright", col = c("black","black"),
       lty = c("solid", "dotdash"), pch = c(NA,NA), bty = "n",
       legend = c("one-step-ahead forecast","95% Confidence Interval"))
