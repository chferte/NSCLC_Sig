library(caret)

setwd("/home/cferte/FELLOW/cferte/MATRIX_RESP_OBJECTS/")
load("MATRIX_TS.Rdata")
load("MATRIX_VS.Rdata")
load("MATRIX_VS2.Rdata")

load("y_TS.Rdata")
load("y_VS.Rdata")
load("y_VS2.Rdata")

normal <-function(vec){
  return((vec-mean(vec))/sd(vec))
}

# DATA normalization
featureData <-t(scale(t(MATRIX_TS)))
testFeature1 <- t(scale(t(MATRIX_VS)))
testFeature2 <- t(scale(t(MATRIX_VS2)))

# response normalization
NresponseData <- normal(y_TS)
NtestResponse1 <- normal(y_VS)
NtestResponse2 <- normal(y_VS2)
# response need not normalization
responseData <- (y_TS)
testResponse1 <- (y_VS)
testResponse2 <- (y_VS2)


#################### GLMNET /lasso
cvglmnet <-cv.glmnet(t(featureData),responseData,nfolds=3)
predResponse1<-predict(cvglmnet,t(testFeature1),s="lambda.min")
predResponse2<-predict(cvglmnet,t(testFeature2),s="lambda.min")


png("/home/cferte/FELLOW/cferte/IN_SOCK/plot_glmnet.png")
par(mfrow = c(1,2))
plot(testResponse1,predResponse1,main = "VS")
plot(testResponse2,predResponse2,main = "VS2")
dev.off()

#################### Partial Least Square
fit_pls<-train(t(featureData),responseData,method = "pls")
predResponse1<-predict(fit_pls,t(testFeature1))
predResponse2<-predict(fit_pls,t(testFeature2))

png("/home/cferte/FELLOW/cferte/IN_SOCK/plot_PLS.png")
par(mfrow = c(1,2))
plot(testResponse1,predResponse1,main = "VS")
plot(testResponse2,predResponse2,main = "VS2")
dev.off()

#################### Principal Component Regression
fit_pcr<-train(t(featureData),responseData,method = "pcr")
predResponse1<-predict(fit_pcr,t(testFeature1))
predResponse2<-predict(fit_pcr,t(testFeature2))

png("/home/cferte/FELLOW/cferte/IN_SOCK/plot_PCR.png")
par(mfrow = c(1,2))
plot(testResponse1,predResponse1,main = "VS")
plot(testResponse2,predResponse2,main = "VS2")
dev.off()


#################### Linear Regression
fit_lm<-train(t(featureData),responseData,method = "lm")
predResponse1<-predict(fit_lm,t(testFeature1))
predResponse2<-predict(fit_lm,t(testFeature2))

png("/home/cferte/FELLOW/cferte/IN_SOCK/plot_LM.png")
par(mfrow = c(1,2))
plot(testResponse1,predResponse1,main = "VS")
plot(testResponse2,predResponse2,main = "VS2")
dev.off()

#################### SVM linear
fit_svm<-train(t(featureData),responseData,method = "svmLinear")
predResponse1<-predict(fit_svm,t(testFeature1))
predResponse2<-predict(fit_svm,t(testFeature2))

png("/home/cferte/FELLOW/cferte/IN_SOCK/plot_SVMLinear.png")
par(mfrow = c(1,2))
plot(testResponse1,predResponse1,main = "VS")
plot(testResponse2,predResponse2,main = "VS2")
dev.off()


##############
train()