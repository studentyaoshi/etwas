args = commandArgs(TRUE)
library(glmnet)
train <- read.csv(args[1], header=TRUE)
train <- na.omit(train) # remove missing
train<-train[-1]
num<-length(names(train))-1
x<-as.matrix(train[-dim(train)[2]])
y<-as.matrix(train[dim(train)[2]])
library(doParallel)
registerDoParallel(10)
fit <- cv.glmnet(x,y,family='gaussian',alpha=1,nlambda=20,nfolds=5)
fit2 <- cv.glmnet(x,y,family='gaussian',alpha=0.5,nlambda=20,nfolds=5)
registerDoSEQ();

lassob<-fit$lambda.min
enetb<-fit2$lambda.min

lassom<-glmnet(x,y,alpha=1,lambda=lassob)
enetm<-glmnet(x,y,alpha=0.5,lambda=enetb)

test <- read.csv(args[2], header=TRUE)
test<-na.omit(test)
test<-test[-1]
testx<-as.matrix(test[-dim(test)[2]])
prediction<-predict(lassom,newx=testx,s=lassob)
prediction2<-predict(enetm,newx=testx,s=enetb)
# mark r<-cor(prediction,test$expression)
# mark r2<-cor(prediction2,test$expression)
# mark R<-r*r
# mark R2<-r2*r2
write.table(R,file=paste('../tem/',args[3],'.','lasso.r2',sep=''),row.names=F,col.names=F,quote=F,append=T,sep='\t')
write.table(R2,file=paste('../tem/',args[3],'.','enet.r2',sep=''),row.names=F,col.names=F,quote=F,append=T,sep='\t')
