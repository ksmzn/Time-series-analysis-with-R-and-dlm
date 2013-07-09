###################################################
#
# USAccDeathsデータは1973~1978年の米国の月別事故死亡者数の1年ずつの時系列データ
#
###################################################

install.packages("TTR")
install.packages("tseries")
install.packages("lawstat")
library(dlm)
library(TTR)
library(tseries)
library(lawstat)
par(mfrow=c(1,1))

#USAccDeaths#
help(USAccDeaths)
USAccDeaths
data=log(USAccDeaths)
ar(data)
par(mfrow=c(3,1))
plot(data,main="Time seriese plot")
plot(diff(data),main="diff(data)")
plot(diff(diff(data)),main="diff(diff(data))")
dev.copy2eps(file="USAccDeathsdiff.eps")

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE)) 
#par(mfrow=c(2,1))
plot(data,main="Time seriese plot")
lines(SMA(data, 12),lty="dotdash")
legend("bottomright", col = c("black","black"),
       lty = c("solid", "dotdash"), pch = c(NA,NA), bty = "n",
       legend = c("log(USAccDeaths)", "12 simple moving average"))
acf(data,main="acf")
pacf(data,main="pacf")

#model1 ローカル線形トレンドモデル＋季節成分
par(mfrow=c(1,1))
slope0 <- mean(diff(data))
freq=12
init_par=log(var(data))
mod1 <- function(u) {
  trend <- dlmModPoly(order=2,
                      #dV = 1e-7,
                      #dV = exp(u[4]),
                      dW = c(exp(u[2]),exp(u[3])),
                      dV = exp(u[1]),
                      #dW = 0,
                      #m0 = c(data[1])
                      m0 = c(data[1],slope0),
                      C0 = (1e-7)*diag(2)
  )
  mod <- dlmModSeas(f=freq,
                    #dV = 0,
                    dV = exp(u[4]),
                    dW = c(exp(u[5]),rep(0,freq-2)),
                    #dW = c(0,0)
                    m0 = rep(0, freq - 1),
                    C0 = 1e+07 * diag(nrow = freq - 1)
  )
  return(trend+mod)
}
init_par=log(var(data))
outMLE1 <- dlmMLE(data,dlmMLE(data,c(init_par,init_par,init_par,init_par,init_par), mod1)$par, mod1)
outMLE1$par
## [1]  -8.431674  -8.331776 -14.144867  -8.431674 -10.801126
dlm1 <- mod1(outMLE1$par)
dlm1$FF
dlm1$V
# Filtering
res1 <- dlmFilter(data,dlm1)
# Smoothing
ress1<- dlmSmooth(data,dlm1)
par(mfrow=c(1,1))
plot(data,xlab="",ylab="value",type="o",col=("darkgrey"))
lines(res1$f,col=2,lty="dotdash")       
legend("bottomleft", col = c("black","black"),
       lty = c("solid", "dotdash"), pch = c(NA,NA), bty = "n",
       legend = c("log(USAccDeaths)", "filtering value"))

par(mfrow=c(3,1))
plot(data,xlab="",ylab="value",type="o",col=("darkgrey"))
ts.plot(ress1$s[-1,1],col=c(1),lty=c("solid"))
ts.plot(ress1$s[-1,3],col=c(1),lty=c("solid"))

#model2 AR(11)
par(mfrow=c(1,1))
ar(data, AIC=T)$aic
slope0=mean(diff(data))
mod2 <- function(u) {
  trend <- dlmModPoly(order=2,
                      dV = 0.000001,
                      dW = exp(u[1:2]),
                      m0 = c(data[1],slope0),
                      C0 = (1e-7)*diag(2)
  )
  ar<- dlmModARMA(
    ar=c(u[4:5],rep(0,9),u[6]),
    #ar=c(u[4:5]),
    #ma=c(u[1]),
    sigma2 = exp(u[3]),
    #C0 = 10^(-10)*diag(11)
    )
  
  return(trend+ar)
}
outMLE2 <- dlmMLE(data,dlmMLE(data,
                              c(init_par,init_par,
                                log(var(ar(data, AIC=F, order.max=13)$resid[-c(1:13)])),
                                ar(data, AIC=T, order.max=13)$ar[1:2],
                                ar(data, AIC=T, order.max=13)$ar[13]
                              ),
                              mod2)$par,mod2)
outMLE2$par
# [1] -5.844469e+00 -6.984875e+00 -6.187222e+00  9.860160e-01 -5.929461e-05
dlm2 <- mod2(outMLE2$par)
dlm2$FF
res2 <- dlmFilter(data,dlm2)
par(mfrow=c(1,1)) 
plot(data,xlab="",ylab="value",type="o",col=("darkgrey"))
lines(res2$f,lty="dotdash")
ress2<- dlmSmooth(data,dlm2)
par(mfrow=c(2,1), mar=c(1.5,4,0.5,0.5) + 0.1, cex=0.6)
plot.ts(dropFirst(ress2$s)[,1],ylab="level+trend",ann=F,yax.flip=TRUE)
legend("bottomright",col=c(1),lty="solid",legend=c("level+trend"))
plot.ts(dropFirst(ress2$s)[,3],main="",ann=F,yax.flip=TRUE, type="o")
legend("bottomright",col=c(1),lty="solid",legend=c("seasonal"))

#残差1
par(mfrow=c(2, 1))
zansa1=residuals(res1)$res
qqnorm(zansa1)
qqline(zansa1)
acf(zansa1)
summary(zansa1)
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
par(mfrow=c(2, 1))
zansa2=residuals(res2)$res
qqnorm(zansa2, main="")
qqline(zansa2)
acf(zansa2)
shapiro.test(zansa2)
ks.test(zansa2, "pnorm", mean=mean(zansa2), sd=sqrt(var(zansa2)))
jarque.bera.test(zansa2)
rjb.test(zansa2)
Box.test(zansa2,type= "Ljung-Box")
Box.test(zansa2,type= "Box-Pierce")

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
plot(data,xlim=c(1973,1980),type="o",col="darkgrey",main="model1",ylab="")
lines(cbind(ciTheory1,fore1$f[,1])[,1], lty="dotdash")
lines(cbind(ciTheory1,fore1$f[,1])[,2], lty="dotdash")
lines(cbind(ciTheory1,fore1$f[,1])[,3], lty="solid")
legend("bottomleft", col = c("black","black"),
       lty = c("solid", "dotdash"), pch = c(NA,NA), bty = "n",
       legend = c("one-step-ahead forecast","theoretical bounds"))
fore2 <- dlmForecast(res2, nAhead=12)
ciTheory2 <- (outer(sapply(fore2$Q, FUN=function(x) sqrt(diag(x))), qnorm(c(0.25,0.75))) +
                as.vector(t(fore2$f)))
plot(data,xlim=c(1973,1980),type="o",col="darkgrey",main="model2")
lines(cbind(ciTheory2,fore2$f[,1])[,1], lty="dotdash")
lines(cbind(ciTheory2,fore2$f[,1])[,2], lty="dotdash")
lines(cbind(ciTheory2,fore2$f[,1])[,3], lty="solid")
legend("bottomleft", col = c("black","black"),
       lty = c("solid", "dotdash"), pch = c(NA,NA), bty = "n",
       legend = c("one-step-ahead forecast","95% Confidence Interval"))
dev.copy2eps(file="USAccDeathsforecast.eps",width=10)
