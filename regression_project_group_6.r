rm(list = ls())

#required packages
library(MASS)
library(lattice)
library(olsrr)
library(ggplot2) 
library(lmtest)
library(L1pack)
library(car)
library(ggplot2)
library(leaps)
index=1:31

#data description
fitness <- read.csv(file = "D:/Regression Project/fitness.csv")
attach(fitness)
data <- cbind.data.frame(index,oxy, age, weight, runtime, 
                         restpulse, runpulse, maxpulse)
mydata <- data[,-1]
n=31 #total number of observations
p=7  #number of parameters with the intercept term

#Explorartory data analysis
str(mydata)
summary(mydata)
par(mfrow=c(2,3))
for(i in 2:7){
plot(mydata[,i],mydata[,1],xlab=colnames(mydata)[i],ylab=colnames(mydata)[1],main=paste(colnames(mydata)[1], "vs", colnames(mydata)[i]))
}

#pairwise scatterplot
pairs(mydata[,-1])
par(mfrow=c(1,1))

#changing column names
colnames(mydata) <- c("y", "x1", "x2", "x3", "x4", "x5", "x6")   
#y = oxy, x1 = age, x2 = weight, x3 = runtime, 
#x4 = restpulse, x5 = runpulse, x6 = maxpulse 
head(mydata)
attach(mydata)


#boxplot of covariates and response
boxplot(mydata)


#------------------------------------------------------------------------

#regression analysis
#multiple linear regression
model <- lm(y ~ x1 + x2 + x3 + x4 + x5 + x6)
summary(model) 
plot(1:nrow(mydata),mydata$y,type = "o",pch = 20,
     ylab = "observed & predicted values",xlab = "Index", col = "blue")
lines(1:nrow(mydata),model$fitted.values,type = "o",pch = 20,
      col = "red",lty = 2)
abline(v = 1:nrow(mydata),lty = 2,col = rgb(0,1,1,alpha = 0.4))
legend("topright",legend = c("Observed","Predicted"),fill = c("blue","red"))

#model considering quadratic relation with weight
model1 <- lm(y ~ x1 + poly(x2,2) + x3 + x4 + x5 + x6) 
summary(model1)
plot(1:nrow(mydata),mydata$y,type = "o",pch = 20,
     ylab = "observed & predicted values",xlab = "Index", col = "blue")
lines(1:nrow(mydata),model1$fitted.values,type = "o",pch = 20,
      col = "red",lty = 2)
abline(v = 1:nrow(mydata),lty = 2,col = rgb(0,1,1,alpha = 0.4))
legend("topright",legend = c("Observed","Predicted"),fill = c("blue","red"))

#comparing the models
plot(1:nrow(mydata),mydata$y,type = "o",pch = 20,
     ylab = "observed & predicted values",xlab = "Index", col = "black")
lines(1:nrow(mydata),model$fitted.values,type = "o",pch = 20,
      col = "blue", lwd=2)
lines(1:nrow(mydata),model1$fitted.values,type = "o",pch = 20,
      col = "red",lty = 2, lwd=2)
abline(v = 1:nrow(mydata),lty = 2,col = rgb(0,1,1,alpha = 0.4))
legend("topright",
       legend = c("Observed","Predictedfrom model-1","Predictedfrom model-2"),
       fill = c("black","blue","red"))


#----------------------------------------------------------------------------


#regression diagnostics
#checking for normality
resi<-residuals(model)
with(ToothGrowth, qqPlot(resi,pch=20, col="darkred", 
                         xlab = "Theoretical quantiles"
     ,ylab ="Sample quantiles", main="Quantile-Quantile Plot"))

#shapiro-wilk test
shapiro.test(resi)

#checking homoscedasticity
#residual plot
ggplot(mydata, aes(fitted(model), residuals(model))) +            
  geom_point(col="blue")+geom_abline(intercept=0,slope=0,col="red")+
  labs(x="Fitted Values", y="Residuals", 
       title="Residual vs Fitted values plot")+
  theme_bw()


#residuals vs covariates plot
par(mfrow=c(2,3))
for(i in 3:8)
{
plot(data[,i], residuals(model), xlab=colnames(data)[i], ylab="Residuals",
main=paste("residuals v/s", colnames(data)[i]), pch=19, col= "darkgreen", grid=TRUE)
abline(h=0, col="red")
}
par(mfrow = c(1,1))

#bi plot
A = as.matrix(mydata[,-1])
H = A%*%solve(t(A)%*%A)%*%t(A)
H_i = diag(H)
e_i = residuals(model)
b_i = (e_i^2)/(1-H_i) 
ggplot(mydata, aes(fitted(model), b_i)) +            
  geom_point(col="blue")+
  labs(x="Fitted Values", y="b_i", title="b_i vs Fitted values plot")+
  theme_bw()


#Breusch-pegan test
bptest(model)


#checking for autocorrelation
xyplot(residuals(model)[2:31]~residuals(model)[1:30], 
grid= TRUE, col = "blue", xlab="e_t-1", ylab="e_t", main="Residual vs Residual")

#acf plot
acf(resi,ylab = "",xlab = "",main = "")
title(xlab="Lag", ylab="ACF", line=2)

#pacf plot
pacf(resi,ylab = "",xlab = "",main = "")
title(xlab="Lag", ylab="PACF", line=2)

#durbin-watson test
dwtest(y~x1+x2+x3+x4+x5+x6,data=mydata, alternative="two.sided")

#breusch-godfrey test
bgtest(model,order = 20)


#detecting influencial points
#hat matrix diagonals
hat_d <- hatvalues(model)
head(sort(hat_d,dec=TRUE))
hat_d[hat_d > 2*(p/n)]

#studentized residual plot
stud <- rstudent(model)
ggplot(data, aes(x=index,y=stud),ylim=c(-3,3)) +            
  geom_point(col="blue")+geom_abline(intercept=0,slope=0,col="red")+
  labs(x="Index", y="Studentized Residuals", 
       title="Studentized Residuals Plot")+
  theme_bw()
ols_plot_resid_stud_fit(model)

#dfbetas
par(mfrow = c(3,3))
DFBETAS = dfbetas(model)
for(i in 1:7)
{
  plot(DFBETAS[,i],main=bquote("DFBETAS for" ~ beta[.(i-1)]),
       ylab="",ylim=c(-1.5,1.5),xlab="",pch=20)
  abline(h=c(-2/sqrt(n),2/sqrt(n)))
  ind = which(abs(DFBETAS[,i]) > 2/sqrt(n)) # beta_0
  text = text(ind,DFBETAS[ind,i],pos = 3,labels = ind)
}
par(mfrow=c(1,1))

#dffits
ols_plot_dffits(model) 

#covratio
COVRATIO = covratio(model)
plot(abs(COVRATIO-1),ylab=expression(abs(COVRATIO-1)),pch = 20)
abline(h = 3*p/n)
ind = which(abs(COVRATIO-1) > 3*p/n)
text(ind,(abs(COVRATIO-1))[ind],pos = 1,labels = ind)
influence.measures(model)

#cook's-D
COOKSD = cooks.distance(model)
plot(COOKSD,pch=20)
abline(h = 3*mean(COOKSD))
ind = which(COOKSD > 3*mean(COOKSD)) # beta_0
text = text(ind,COOKSD[ind],pos = 1,labels = ind)


#checking for multicollinearity
#pairwise correlation
round(cor(mydata[,-1]),2)

#eigen-decomposition and condition number
X_mdl = model.matrix(model)[,-1]
e <- eigen(t(X_mdl) %*% X_mdl)
e$val
kappa(scale(X_mdl))

#VIF
vif(lm(y~x1+x2+x3+x4+x5+x6,data = mydata))#high implies collinearity 

#with influencial points
vif(lm(y~x1+x2+x3+x4+x5+x6,data = mydata[]))

#without influencial points
vif(lm(y~x1+x2+x3+x4+x5+x6,data = mydata[-c(10,15),]))

#---------------------------------------------------------------------------


#remedies
#removing influential points
model3=lm(y~x1+x2+x3+x4+x5+x6,data=mydata[-c(10,15),])
summary(model3)
summary(model)

#improvments of the model
hist(rstandard(model3),breaks = 10,main = "Histogram Of Residual Values")
with(ToothGrowth,qqPlot(residuals(model3),pch=20,col="darkred", 
                        xlab = "Theoretical    quantiles",
                        ylab ="Sample quantiles", 
                        main="Quantile-Quantile Plot"))


#remedies for collinearity
#vif
model3=lm(y~x1+x2+x3+x4+x5+x6,data=mydata[-c(10,15),])
X_mdl3 = scale(model.matrix(model3)[ ,-1])
kappa(X_mdl3)
summary(model3)
model4=lm(y~x1+x2+x3+x4+x5,data=mydata[-c(10,15),])
X_mdl4 = scale(model.matrix(model4)[ ,-1])
kappa(X_mdl4)
summary(model4)
model5=lm(y~x1+x2+x3+x4+x6,data=mydata[-c(10,15),])
X_mdl5 = scale(model.matrix(model5)[ ,-1])
kappa(X_mdl5)
summary(model5)


vif(model3)

#added variable plot
avPlots(lm(y~x1+x2+x3+x4+x5+x6,data = mydata[-c(10,15),]))


#remedies for non-normality
#box-cox transformation
boxcox(model3) 
par(mfrow=c(1,1)) 
bc <- boxcox(y~., data=mydata[-c(10,15),], plotit=TRUE) 
lambda_boxcox <- bc$x[which.max(bc$y)] 
ynew_boxcox <- (y^lambda_boxcox-1)/lambda_boxcox
mydata1 = cbind.data.frame(ynew_boxcox, x1, x2, x3, x4, x5, x6) 
model6 = lm(ynew_boxcox~., data = mydata1[-c(10,15),]) 
summary(model6)

resinew_boxcox <-residuals(model6) 
with(ToothGrowth, qqPlot(resinew_boxcox ,pch=20, col="darkred",xlab = "Theoretical quantiles",ylab ="Sample quantiles", main="Quantile-Quantile Plot after applying Box-Cox transformation"
,ylim=c(-0.04,0.04)))
shapiro.test(resinew_boxcox)

#draper-john transformation
l= seq(from=-10,to=10,by=0.3)
like=NULL
for(i in 1:length(l))
{
ynew = sign(y)*((abs(y)+1)^l[i]-1)/l[i]
mydata2 = cbind.data.frame(ynew, x1, x2, x3, x4, x5, x6)
mod=lm(ynew~., data=mydata2[-c(10,15),])
s=sum(residuals(mod)^2)
t=s/23
like[i]=(2*pi*t)^(-29/2)*exp(-1/(2*t)*s)*prod(y^(l[i]-1))
}
xyplot(log(like)~l, grid = TRUE, col = "blue", type = "l")
lambda_john=l[which.max(like)]
lambda_john

ynew_john =sign(y)*((abs(y)+1)^lambda_john-1)/lambda_john
mydata3 = cbind.data.frame(ynew_john, x1, x2, x3, x4, x5, x6)
model7 = lm(ynew_john ~., data = mydata3[-c(10,15),])
summary(model7)

resinew_john <- residuals(model7)
with(ToothGrowth, qqPlot(resinew_john,pch=20, col="darkred", 
                         xlab = "Theoretical quantiles"
     ,ylab ="Sample quantiles", main="Quantile-Quantile Plot after applying Draper-John Transformation", ylim= c(-8*10^10,8*10^10)))
shapiro.test(resinew_john)

#-----------------------------------------------------------------------------

#model selection
#stepwise selection based on aic
mod_wt_out = lm(y~x1+x2+x3+x4+x5+x6, data = mydata[-c(10,15),])
ols_step_both_aic(mod_wt_out,details = TRUE)
AIC(model3)

#stepwise selection based on adjusted r^2
library(leaps)
bestsubset <- regsubsets(y~., data=mydata[-c(10,15),], method="seqrep",nbest=1)
summary(bestsubset)
summary(bestsubset)$adjr2
CP=summary(bestsubset)$cp;
P=1:6
abs(CP-P)
model0=lm(y~1,data=mydata[c(-10,-15), ])
adj0=1/29
CP0=sum(residuals(model0)^2)/((sum(residuals(model3)^2))/23)-29

#plots
P=0:6
aic <- c(176.940, 137.813, 136.140, 130.659, 128.991, 128.894, 130.806)
adj_r <- c(0.034, 0.749, 0.770, 0.815, 0.831, 0.836, 0.829)
mal_p <- c(141.7312, 13.547, 9.842, 2.919, 0.721, 0.067, 1.000)

par(mfrow=c(1,1))

ind = which.min(aic)
plot(P, aic, pch = 20,type="b",ylab="AIC",xlab="No. of covariates",main="AIC vs No. of covariates")
points(ind-1,aic[ind],col="red",pch=20)

ind = which.max(adj_r)
plot(P, adj_r, pch = 20,type="b",ylab="Adjusted R^2",xlab="No. of covariates",main="Adjusted R^2 vs No. of covariates")
points(ind-1,adj_r[ind],col="red",pch=20)

ind = which.min(mal_p)
plot(P, mal_p, pch = 20,type="b",ylab="|Cp-p|",xlab="No. of covariates",main="|Cp-p| vs No. of covariates")
points(ind-1,mal_p[ind],col="red",pch=20)

#best subset selection
X1 <- mydata[-c(10,15),]
names(X1) <- c("R","V1","V2","V3","V4","V5","V6")
models<-list()
models[["V1"]]<-lm(R~V1,X1)
models[["V2"]]<-lm(R~V2,X1)
models[["V3"]]<-lm(R~V3,X1)
models[["V4"]]<-lm(R~V4,X1)
models[["V5"]]<-lm(R~V5,X1)
models[["V6"]]<-lm(R~V6,X1)
models[["V12"]]<-lm(R~V1+V2,X1)
models[["V13"]]<-lm(R~V1+V3,X1)
models[["V14"]]<-lm(R~V1+V4,X1)
models[["V15"]]<-lm(R~V1+V5,X1)
models[["V16"]]<-lm(R~V1+V6,X1)
models[["V23"]]<-lm(R~V2+V3,X1)
models[["V24"]]<-lm(R~V2+V4,X1)
models[["V25"]]<-lm(R~V2+V5,X1)
models[["V26"]]<-lm(R~V2+V6,X1)
models[["V34"]]<-lm(R~V3+V4,X1)
models[["V35"]]<-lm(R~V3+V5,X1)
models[["V36"]]<-lm(R~V3+V6,X1)
models[["V45"]]<-lm(R~V4+V5,X1)
models[["V46"]]<-lm(R~V4+V6,X1)
models[["V56"]]<-lm(R~V5+V6,X1)
models[["V123"]]<-lm(R~V1+V2+V3,X1)
models[["V124"]]<-lm(R~V1+V2+V4,X1)
models[["V125"]]<-lm(R~V1+V2+V5,X1)
models[["V126"]]<-lm(R~V1+V2+V6,X1)
models[["V134"]]<-lm(R~V1+V3+V4,X1)
models[["V135"]]<-lm(R~V1+V3+V5,X1)
models[["V136"]]<-lm(R~V1+V3+V6,X1)
models[["V145"]]<-lm(R~V1+V4+V5,X1)
models[["V146"]]<-lm(R~V1+V4+V6,X1)
models[["V156"]]<-lm(R~V1+V5+V6,X1)
models[["V234"]]<-lm(R~V2+V3+V4,X1)
models[["V235"]]<-lm(R~V2+V3+V5,X1)
models[["V236"]]<-lm(R~V2+V3+V6,X1)
models[["V245"]]<-lm(R~V2+V4+V5,X1)
models[["V246"]]<-lm(R~V2+V4+V6,X1)
models[["V256"]]<-lm(R~V2+V5+V6,X1)
models[["V345"]]<-lm(R~V3+V4+V5,X1)
models[["V346"]]<-lm(R~V3+V4+V6,X1)
models[["V356"]]<-lm(R~V3+V5+V6,X1)
models[["V456"]]<-lm(R~V4+V5+V6,X1)
models[["V1234"]]<-lm(R~V1+V2+V3+V4,X1)
models[["V1235"]]<-lm(R~V1+V2+V3+V5,X1)
models[["V1236"]]<-lm(R~V1+V2+V3+V6,X1)
models[["V1245"]]<-lm(R~V1+V2+V4+V5,X1)
models[["V1246"]]<-lm(R~V1+V2+V4+V6,X1)
models[["V1256"]]<-lm(R~V1+V2+V5+V6,X1)
models[["V1345"]]<-lm(R~V1+V3+V4+V5,X1)
models[["V1346"]]<-lm(R~V1+V3+V4+V6,X1)
models[["V1356"]]<-lm(R~V1+V3+V5+V6,X1)
models[["V1456"]]<-lm(R~V1+V4+V5+V6,X1)
models[["V2345"]]<-lm(R~V2+V3+V4+V5,X1)
models[["V2346"]]<-lm(R~V2+V3+V4+V6,X1)
models[["V2356"]]<-lm(R~V2+V3+V5+V6,X1)
models[["V2456"]]<-lm(R~V2+V4+V5+V6,X1)
models[["V3456"]]<-lm(R~V3+V4+V5+V6,X1)
models[["V12345"]]<-lm(R~V1+V2+V3+V4+V5,X1)
models[["V12346"]]<-lm(R~V1+V2+V3+V4+V6,X1)
models[["V12456"]]<-lm(R~V1+V2+V4+V5+V6,X1)
models[["V13456"]]<-lm(R~V1+V3+V4+V5+V6,X1)
models[["V12356"]]<-lm(R~V1+V2+V3+V5+V6,X1)
models[["V23456"]]<-lm(R~V2+V3+V4+V5+V6,X1)
models[["V123456"]]<-lm(R~V1+V2+V3+V4+V5+V6,X1)

mnames<- factor(names(models),levels = names(models))

#values of different criteria for model selection

adj.R2 <- sapply(models, function(fit) summary(fit)$adj.r.squared)
adj.R2[which.max(adj.R2)]

sigma.sq <- summary(models[["V123456"]])$sigma^2
Cp <- sapply(models, function(fit) extractAIC(fit, scale = sigma.sq)[2])
p <- rep(c(1,2,3,4,5,6), c(6, 15, 20, 15, 6, 1))
mal = Cp-p
mal[which.min(Cp-p)]

AIC <- sapply(models, function(fit) AIC(fit))
AIC[which.min(AIC)]

CV = NULL
for(i in 1:15)
{
  X_mdl_mat = model.matrix(models[[i]])
  head(X_mdl_mat)
  Y_vec = X1$V5
  H = X_mdl_mat%*%solve(t(X_mdl_mat)%*%X_mdl_mat)%*%t(X_mdl_mat)
  h = diag(H)
  n_h = nrow(X1)
  CV[i] = (1/n_h)*sum((Y_vec-H%*%Y_vec)^2/(1-h)^2)
}
mnames[which.max(CV)]
CV[which.max(CV)]

#checking of error assumptions with optimal model

##homoscedasticity
mnew= lm(y~x1+x2+x3+x5+x6, data=mydata[-c(10,15),])
plot(fitted(mnew), residuals(mnew), col = "blue", xlab="Fitted values", ylab="Residuals", main="Residuals vs Fitted values")
abline(h=0)
bptest(mnew)

##autocorrelation
xyplot(as.vector(residuals(mnew))[-29]~as.vector(residuals(mnew))[-1],pch = 20,ylab = "e_t",xlab = "e_t-1", col = "blue", grid=TRUE)
bgtest(mnew)
dwtest(mnew)

##normality
with(ToothGrowth, qqPlot(residuals(mnew),pch=20, col="darkred", xlab = "Theoretical quantiles"
     ,ylab ="Sample quantiles", main="Quantile-Quantile Plot"))
shapiro.test(residuals(mnew))

#checking with randomly generated data
d=rnorm(30)
qqnorm(d)
qqline(d)

#alternative estimation methods
#lasso

library(glmnet)
yn=mydata[-c(10,15),1]
xn=data.matrix(mydata[-c(10,15),-c(1,5)])
model7=glmnet(xn, yn,alpha=1)
summary(model7)
cvmodel=cv.glmnet(xn,yn,alpha=1)
bl=cvmodel$lambda.min
bl
plot(cvmodel)
bmodel=glmnet(xn,yn,alpha=1,lambda=bl)
coef(bmodel)
plot(model7,xvar="lambda")
y_pr=predict(model7,s=bl,newx=xn)
sst=sum((yn-mean(yn))^2)
sse=sum((yn-y_pr)^2)
rsq=1-sse/sst; rsq
adjr2=1-(28)*(1-rsq)/24; adjr2
summary(lm(y~x1+x2+x3+x5+x6, data=mydata[-c(10,15),-5]))

#ridge

yp=mydata[-c(10,15),1]
xp=data.matrix(mydata[-c(10,15),-c(1,5)])
model8=glmnet(xp, yp,alpha=0)
summary(model8)
cvmodel1=cv.glmnet(xp,yp,alpha=0)
bl1=cvmodel1$lambda.min
bl1
plot(cvmodel1)
bmodel1=glmnet(xp,yp,alpha=0,lambda=bl1)
coef(bmodel1)
plot(model8,xvar="lambda")
y_pr1=predict(model8,s=bl1,newx=xp)
sst1=sum((yp-mean(yp))^2)
sse1=sum((yp-y_pr1)^2)
rsq1=1-sse1/sst1; rsq1
adjr21=1-(28)*(1-rsq1)/24; adjr21