###################################################
#
# Lynx データは 1821-1934 年に罠にかかったカナダ・オオヤマネコの年別の頭数の, 一年ずつの 時系列データ
#
###################################################

install.packages("TTR")
install.packages("tseries")
install.packages("lawstat")
library(dlm)
library(TTR)
library(tseries)
library(lawstat)

#lynx#
data=log(lynx)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE)) 
plot(data,main="Time seriese plot")
#移動平均のプロット
lines(SMA(data, 12),lty="dotdash")
lines(SMA(data,12*4),lty="longdash")
legend("bottomleft", col = c("black", "black","black"),
       lty = c("solid", "dotdash","longdash"), pch = c(NA, NA,NA), bty = "n",
       legend = c("log(lynx)", "12 simple moving average","12 * 4 simple moving average"))
#自己相関・偏自己相関
acf(data,main="acf")
pacf(data,main="pacf")

#model1 ランダムウォークプラスノイズモデルと季節成分#
par(mfrow=c(1,1))
mod1 <- function(u) {
  trend <- dlmModPoly(order=1,
                      dV = 0,
                      dW = exp(u[1]),
                      m0 = c(data[1])
                      )
  mod <- dlmModTrig(s=10,
                    q=2,
                    dV = exp(u[2]),
                    dW = c(exp(u[3]),0)
                    )
  return(trend+mod)
}
#最尤推定でパラメータを求める
init_par=log(var(data))
outMLE1 <- dlmMLE(data,dlmMLE(data,c(init_par,init_par,init_par), mod1)$par, mod1)
outMLE1
##  [1]  -2.141477 -16.078952  -3.797000
dlm1 <- mod1(outMLE1$par)
#フィルタリング
res1 <- dlmFilter(data,dlm1)
#スムージング
ress1<- dlmSmooth(data,dlm1)
par(mfrow=c(2,1))
plot(data,xlim=c(1818, 1950),xlab="",ylab="value",type="o",col=("darkgrey"))
lines(res1$f,col=2,lty="dotdash")       
legend("bottomleft", col = c("black","black"),
       lty = c("solid", "dotdash"), pch = c(NA,NA), bty = "n",
       legend = c("log(lynx)", "filtering value"))
ts.plot(cbind(data,ress1$s[-1,1]),col=c(1,4),lty=c("solid","dotdash"))
legend("bottomleft", col = c(1, 4),
       lty = c("solid", "dotdash"), pch = c(NA,NA), bty = "n",
       legend = c("log(lynx)", "smoothing value"))       


#model2 AR(11)を考える
par(mfrow=c(1,1))
mod2 <- function(u) {
  trend <- dlmModPoly(order=1,
                      dV = 0.000001,
                      dW = exp(u[8]),
                      m0 = c(data[1]),
                      C0 = 10^5*diag(1)
                      )
  ar<- dlmModARMA(ar=c(u[1:3],rep(0,5),u[4:6]),
                  sigma2 = exp(u[7]),
                  C0 = 10^(-10)*diag(11))
                  
  return(trend+ar)
}
#最尤推定
outMLE2 <- dlmMLE(data,dlmMLE(data,
                             c(ar(data, AIC=T, order.max=11)$ar[1:3],
                               ar(data, AIC=T, order.max=11)$ar[9:11],
                               log(var(ar(data, AIC=F, order.max=11)$resid[-c(1:11)]))
                               ,init_par
                               ),
                             mod2
                             )$par,
                 mod2)

outMLE2
# [1]  1.16911527 -0.39201165 -0.04414805  0.11077001  0.25857009 -0.34678019 -1.53413076 -8.43936547
dlm2 <- mod2(outMLE2$par)
dlm2$F
#フィルタリング
res2 <- dlmFilter(data,dlm2)
par(mfrow=c(1,1)) 
plot(data,xlab="",ylab="value",type="o",col=("darkgrey"))
lines(res2$f,lty="dotdash")
#スムージング
ress2<- dlmSmooth(data,dlm2)
par(mfrow=c(2,1), mar=c(1.5,4,0.5,0.5) + 0.1, cex=0.6)
plot.ts(dropFirst(ress2$s)[,1],ylab="level+trend",ann=F,yax.flip=TRUE)
legend("bottomright",col=c(1),lty="solid",legend=c("level+trend"))
plot.ts(dropFirst(ress2$s)[,3],main="",ann=F,yax.flip=TRUE, type="o")
legend("bottomright",col=c(1),lty="solid",legend=c("seasonal"))
ts.plot(cbind(data,dropFirst(ress2$s)[,1]),col=c("darkgrey",4))

#残差の自己相関や正規性をチェック
#model1
par(mfrow=c(2, 1))
zansa1=residuals(res1)$res
qqnorm(zansa1)
qqline(zansa1)
acf(zansa1)

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




#model2
zansa2=residuals(res2)$res
qqnorm(zansa2, main="")
qqline(zansa2)
acf(zansa2)
dev.copy2eps(file="lynxzansamod2.eps",width=10)
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
