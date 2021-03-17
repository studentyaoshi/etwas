args= commandArgs(TRUE)
library(glmnet)
data <- read.csv(args[1], header=T)
data<-na.omit(data) # remove missing
data<-data[-1]
num<-length(names(data))-1
x<-as.matrix(data[-dim(data)[2]])
y<-as.matrix(data[dim(data)[2]])
fit <- cv.glmnet(x,y,family='gaussian',alpha=args[2],nlambda=20,nfolds=10)
lam<-fit$lambda.min
lastm<-glmnet(x,y,alpha=args[2],lambda=lam)
weights<-coef(lastm,s=lam)
write.table(as.data.frame(as.matrix(weights)),file=args[3])
