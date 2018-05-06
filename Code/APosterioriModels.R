#A POSTERIORI PREDICTION

library(performanceEstimation)
library(UBL)
library(uba)
library(dplyr)
library(bbmle)

source("./Code/EvalFramework.R")

news <- read.csv("./Data/News_Final.csv") 
news.economy <- read.csv("./Data/Facebook_Economy.csv")
news.microsoft <- read.csv("./Data/Facebook_Microsoft.csv")
news.obama <- read.csv("./Data/Facebook_Obama.csv")
news.palestine <- read.csv("./Data/Facebook_Palestine.csv")

#news.economy <- read.csv("./Data/GooglePlus_Economy.csv")
#news.microsoft <- read.csv("./Data/GooglePlus_Microsoft.csv")
#news.obama <- read.csv("./Data/GooglePlus_Obama.csv")
#news.palestine <- read.csv("./Data/GooglePlus_Palestine.csv")

#news.economy <- read.csv("./Data/LinkedIn_Economy.csv")
#news.microsoft <- read.csv("./Data/LinkedIn_Microsoft.csv")
#news.obama <- read.csv("./Data/LinkedIn_Obama.csv")
#news.palestine <- read.csv("./Data/LinkedIn_Palestine.csv")


##CONSTANT SCALING
const.scale <- function(form,train,test,...) {
  
  #Reduce to those which have a time difference bigger than 2 days with the lowest timestamp in the testset
  train <- train[difftime(as.POSIXct(news[news$IDLink==test[1,]$IDLink,]$PublishDate),news[match(train$IDLink,news$IDLink),]$PublishDate,unit="hours")>48,]
  
  model.results <- data.frame(TimeSlice=numeric(0),
                              mse=numeric(0),mse_phi=numeric(0),
                              prec=numeric(0),rec=numeric(0),
                              F1=numeric(0))
  
  preds <- data.frame(IDLink=numeric(0),TimeSlice=numeric(0),Prediction=numeric(0))
  
  for(row in 1:nrow(test)) {
    
    row.explore <- test[row,]
    
    for(ts.exp in 1:3) { #Prediction for the first 3 timeslices
      
      if(row.explore[,ts.exp+1]==-1) {
        
      } else {
        
        x <- train[,paste0("TS",ts.exp)]
        y <- train$TS144
        
        y <- y[x>=0]
        x <- x[x>=0]
        
        x_test <- row.explore[,paste0("TS",ts.exp)]
        
        alpha <- sum(x/y,na.rm=TRUE)/sum((x/y)^2,na.rm=TRUE)
        
        y_pred <- alpha*x_test
        
        preds.row <- data.frame(IDLink=test[row,]$IDLink,TimeSlice=ts.exp,Prediction=y_pred)
        if(!identical(names(preds.row), names(preds))) { names(preds.row) <- names(preds) }
        preds <- rbind(preds,preds.row)
        
      }
      
    }
    
  }
  
  preds["Target"] <- test[match(preds$IDLink,test$IDLink),]$TS144
  preds$Prediction <- unlist(preds$Prediction)
  
  model.results <- c()
  
  for(a in 1:3) {
    
    aux <- preds[preds$TimeSlice==a,]
    aux <- aux[aux$Target>=0 & !is.na(aux$Target),]
    aux <- aux[aux$Prediction>=0,]
    
    utility.stats <- eval.stats(train[train$TS144>0,]$TS144,aux$Target,aux$Prediction)
    
    model.results <- rbind(model.results,utility.stats)
    
  }
  
  #INFORMATION RETURN
  res <- list(evaluation=model.results)
  res
  
}

##LINEAR LOG
linear.log <- function(form,train,test,...) {
  
  #Reduce to those which have a time difference bigger than 2 days with the lowest timestamp in the testset
  train <- train[difftime(as.POSIXct(news[news$IDLink==test[1,]$IDLink,]$PublishDate),news[match(train$IDLink,news$IDLink),]$PublishDate,unit="hours")>48,]
  
  model.results <- data.frame(TimeSlice=numeric(0),
                              mse=numeric(0),mse_phi=numeric(0),
                              prec=numeric(0),rec=numeric(0),
                              F1=numeric(0))
  
  preds <- data.frame(IDLink=numeric(0),TimeSlice=numeric(0),Prediction=numeric(0))
  
  for(row in 1:nrow(test)) {
    
    row.explore <- test[row,]
    
    for(ts.exp in 1:3) { #Prediction for the first 3 timeslices
      
      if(row.explore[,ts.exp+1]==-1) {
        
      } else {
        
        x <- train[,paste0("TS",ts.exp)]
        y <- train$TS144
        
        y <- y[x>0]
        x <- x[x>0]
        
        x <- log(x)
        y <- log(y)
        
        x_test <- row.explore[,paste0("TS",ts.exp)]
        if(x_test!=0) { x_test <- log(x_test) }
        
        
        LL <- function(beta0, beta1, mu, sigma) {
          R = y - x * beta1 - beta0
          #
          R = suppressWarnings(dnorm(R, mu, sigma, log = TRUE))
          #
          -sum(R)
        }
        
        fit <- mle2(LL, start = list(beta0 = 2, beta1 = 1, mu = 0, sigma=1))
        beta0 <- fit@details$par[1]
        sigma2 <- (fit@details$par[4]^2)
        
        y_pred <- exp(as.numeric(x_test+beta0+sigma2/2))
        
        if(is.na(y_pred)) { y_pred <- 0 }
        
        preds.row <- data.frame(IDLink=test[row,]$IDLink,TimeSlice=ts.exp,Prediction=y_pred)
        if(!identical(names(preds.row), names(preds))) { names(preds.row) <- names(preds) }
        preds <- rbind(preds,preds.row)
        
      }
      
    }
    
  }
  
  preds["Target"] <- test[match(preds$IDLink,test$IDLink),]$TS144
  preds$Prediction <- unlist(preds$Prediction)
  
  model.results <- c()
  
  for(a in 1:3) {
    
    aux <- preds[preds$TimeSlice==a,]
    aux <- aux[aux$Target>=0 & !is.na(aux$Target),]
    aux <- aux[aux$Prediction>=0,]
    
    utility.stats <- eval.stats(train[train$TS144>0,]$TS144,aux$Target,aux$Prediction)
    
    model.results <- rbind(model.results,utility.stats)
    
  }
  
  #INFORMATION RETURN
  res <- list(evaluation=model.results)
  res
  
}

#KERNEL(IQR)
mc.kernel <- function(form,train,test,...) {
  
  #Reduce to those which have a time difference bigger than 2 days with the lowest timestamp in the testset
  train <- train[difftime(as.POSIXct(news[news$IDLink==test[1,]$IDLink,]$PublishDate),news[match(train$IDLink,news$IDLink),]$PublishDate,unit="hours")>48,]
  
  iqr <- c()
  
  for(i in 2:145) {
    values <- train[train[,i]>=0,][,i]
    iqr <- c(iqr,IQR(values[!(values %in% boxplot.stats(values)$out)]))
  }
  
  model.results <- data.frame(TimeSlice=numeric(0),
                              mse=numeric(0),mse_phi=numeric(0),
                              prec=numeric(0),rec=numeric(0),
                              F1=numeric(0))
  
  preds <- data.frame(IDLink=numeric(0),TimeSlice=numeric(0),Prediction=numeric(0))
  
  for(ts.exp in 1:3) { #Prediction for the first 3 timeslices
        
    ys <- test[,ts.exp+1]
    ys <- unique(ys[ys>=0])
    
    for(i in 1:length(ys)) {
      
      upperbound <- as.numeric(ys[i] + iqr[ts.exp])
      lowerbound <- as.numeric(ys[i] - iqr[ts.exp])
      lowerbound[lowerbound<0] <- 0
      
      ids.over <- train[train[,ts.exp+1]<=upperbound,]$IDLink
      ids.under <- train[train[,ts.exp+1]>=lowerbound,]$IDLink
      
      ids.intersect <- intersect(ids.over,ids.under)
      
      if(length(ids.intersect)==0) {
        iqr_i <- 1
        if(iqr[ts.exp]==0) iqr[ts.exp] <- 1
        while(length(ids.intersect)==0) {
          
          upperbound <- as.numeric(ys[i] + (iqr[ts.exp]*iqr_i))
          lowerbound <- as.numeric(ys[i] - (iqr[ts.exp]*iqr_i))
          lowerbound[lowerbound<0] <- 0
          
          ids.over <- train[train[,ts.exp+1]<=upperbound,]$IDLink
          ids.under <- train[train[,ts.exp+1]>=lowerbound,]$IDLink
          
          ids.intersect <- intersect(ids.over,ids.under)
          
          iqr_i <- iqr_i + 1
          
        }
      }
      
      ids.4analysis <- train[train$IDLink %in% ids.intersect & train[,ts.exp+1]>=0,]$IDLink
      
      if(length(ids.4analysis)>0) {
        
        results <- data.frame(TargetNow=numeric(0),DiffWeight=double(0),Target=numeric(0),KMeansPred=double(0))
        
        for(a in 1:length(ids.4analysis)) {
          
          row.analysis <- train[train$IDLink==ids.4analysis[a],]
          
          diff <- as.numeric(abs(row.analysis[ts.exp+1] - ys[i]))
          
          final <- as.numeric(row.analysis[ts.exp+1])
          
          target.now <- as.numeric(row.analysis[ts.exp+1])
          weight <- as.numeric(1-(diff/final))
          if(is.na(weight)) {weight<-1}
          target.final <- as.numeric(row.analysis$TS144)
          
          if(weight>0 & weight<=1) {
            
            kmeanspred <- as.numeric(weight * target.final)
            
            row.aux <- data.frame(TargetNow=target.now,DiffWeight=weight,Target=target.final,KMeansPred=kmeanspred)
            if(!identical(names(row.aux), names(results))) { names(row.aux) <- names(results) }
            results <- rbind(results,row.aux)
            
          }
          
        }
        
        weighted.pred <- sum(results$KMeansPred)/sum(results$DiffWeight)
        if(is.na(weighted.pred)) {weighted.pred<-0}
        now.target <- ys[i]
        if(weighted.pred<now.target) {
          weighted.pred <- now.target
        }
        
        
        test.ids <- test[test[,ts.exp+1]==now.target,]$IDLink
        preds.row <- data.frame(IDLink=test.ids,TimeSlice=ts.exp,Prediction=weighted.pred)
        if(!identical(names(preds.row), names(preds))) { names(preds.row) <- names(preds) }
        preds <- rbind(preds,preds.row)
        
      } else {
        
        #if there is no evidence to support the prediction the prediction will be equal to the present amount
        test.ids <- test[test[,ts.exp+1]==ys[i],]$IDLink
        preds.row <- data.frame(IDLink==test.ids,TimeSlice=ts.exp,Prediction=ys[i])
        print(preds.row)
        if(!identical(names(preds.row), names(preds))) { names(preds.row) <- names(preds) }
        preds <- rbind(preds,preds.row)
        
      }
      
    }
    
  }
  
  preds["Target"] <- test[match(preds$IDLink,test$IDLink),]$TS144
  preds$Prediction <- unlist(preds$Prediction)
  
  model.results <- c()
  
  for(a in 1:3) {
    
    aux <- preds[preds$TimeSlice==a,]
    aux <- aux[aux$Target>=0 & !is.na(aux$Target),]
    aux <- aux[aux$Prediction>=0,]
    
    utility.stats <- eval.stats(train[train$TS144>0,]$TS144,aux$Target,aux$Prediction)
    
    model.results <- rbind(model.results,utility.stats)
    
  }
  
  #INFORMATION RETURN
  res <- list(evaluation=model.results)
  res
  
}

#KNN(IQR)
mc.knn <- function(form,train,test,...) {
  
  #Reduce to those which have a time difference bigger than 2 days with the lowest timestamp in the testset
  train <- train[difftime(as.POSIXct(news[news$IDLink==test[1,]$IDLink,]$PublishDate),news[match(train$IDLink,news$IDLink),]$PublishDate,unit="hours")>48,]
  
  iqr <- c()
  
  for(i in 2:145) {
    values <- train[train[,i]>=0,][,i]
    iqr <- c(iqr,IQR(values[!(values %in% boxplot.stats(values)$out)]))
  }
  
  preds <- data.frame(IDLink=numeric(0),TimeSlice=numeric(0),Prediction=numeric(0))
  
  for(ts.exp in 1:3) { #Prediction for the first 3 timeslices
    
    ys <- test[,ts.exp+1]
    ys <- unique(ys[ys>=0])
    
    for(i in 1:length(ys)) {
      
      upperbound <- as.numeric(ys[i] + iqr[ts.exp])
      lowerbound <- as.numeric(ys[i] - iqr[ts.exp])
      lowerbound[lowerbound<0] <- 0
      
      ids.over <- train[train[,ts.exp+1]<=upperbound,]$IDLink
      ids.under <- train[train[,ts.exp+1]>=lowerbound,]$IDLink
      
      ids.intersect <- intersect(ids.over,ids.under)
      
      if(length(ids.intersect)==0) {
        iqr_i <- 1
        if(iqr[ts.exp]==0) iqr[ts.exp] <- 1
        while(length(ids.intersect)==0) {
          
          upperbound <- as.numeric(ys[i] + (iqr[ts.exp]*iqr_i))
          lowerbound <- as.numeric(ys[i] - (iqr[ts.exp]*iqr_i))
          lowerbound[lowerbound<0] <- 0
          
          ids.over <- train[train[,ts.exp+1]<=upperbound,]$IDLink
          ids.under <- train[train[,ts.exp+1]>=lowerbound,]$IDLink
          
          ids.intersect <- intersect(ids.over,ids.under)
          
          iqr_i <- iqr_i + 1
          
        }
      }
      
      ids.4analysis <- train[train$IDLink %in% ids.intersect & train[,ts.exp+1]>=0,]$IDLink
      
      if(length(ids.4analysis)>0) {
        
        rows.analysis <- train[train$IDLink %in% ids.intersect & train[,ts.exp+1]>=0,]
        
        slope <- numeric(0)
        slope.results <- (rows.analysis$TS144-rows.analysis[,ts.exp+1])/(144-ts.exp)
        #slope <- IQR(slope.results)
        slope <- median(slope.results)
        
        now <- ys[i]
        
        scalar.pred <- now + (144-ts.exp)*slope
        
        test.ids <- test[test[,ts.exp+1]==now,]$IDLink
        preds.row <- data.frame(IDLink=test.ids,TimeSlice=ts.exp,Prediction=scalar.pred)
        if(!identical(names(preds.row), names(preds))) { names(preds.row) <- names(preds) }
        preds <- rbind(preds,preds.row)
        
      } else {
        
        #if there is no evidence to support the prediction the prediction will be equal to the present amount
        test.ids <- test[test[,ts.exp+1]==ys[i],]$IDLink
        preds.row <- data.frame(IDLink==test.ids,TimeSlice=ts.exp,Prediction=ys[i])
        print(preds.row)
        if(!identical(names(preds.row), names(preds))) { names(preds.row) <- names(preds) }
        preds <- rbind(preds,preds.row)
        
      }
      
    }
    
  }
  
  preds["Target"] <- test[match(preds$IDLink,test$IDLink),]$TS144
  preds$Prediction <- unlist(preds$Prediction)
  
  model.results <- c()
  
  for(a in 1:3) {
    
    aux <- preds[preds$TimeSlice==a,]
    aux <- aux[aux$Target>=0 & !is.na(aux$Target),]
    aux <- aux[aux$Prediction>=0,]
    
    utility.stats <- eval.stats(train[train$TS144>0,]$TS144,aux$Target,aux$Prediction)
    
    model.results <- rbind(model.results,utility.stats)
    
  }
  
  #INFORMATION RETURN
  res <- list(evaluation=model.results)
  res
  
}

#ASUR AND HUBERMAN
asur <- function(form,train,test,...) {
  
  #Reduce to those which have a time difference bigger than 2 days with the lowest timestamp in the testset
  train <- train[difftime(as.POSIXct(news[news$IDLink==test[1,]$IDLink,]$PublishDate),news[match(train$IDLink,news$IDLink),]$PublishDate,unit="hours")>48,]
  
  preds <- data.frame(TimeSlice=numeric(0),Prediction=numeric(0),Target=numeric(0))
  
  for(ts.exp in 1:3) { #Prediction for the first 3 timeslices
    
    new_train <- data.frame(Rate=train[,ts.exp+1]/ts.exp)
    new_train["SentimentTitle"] <- news[match(train$IDLink,news$IDLink),]$SentimentTitle
    new_train["SentimentHeadline"] <- news[match(train$IDLink,news$IDLink),]$SentimentHeadline
    new_train["Now"] <- train[,ts.exp+1]
    new_train["Target"] <- train[,144+1]
    new_train <- new_train[new_train$Target>=0,]
    if(any(!complete.cases(new_train))) new_train[is.na(new_train)] <- 0
    
    new_test <- data.frame(Rate=test[,ts.exp+1]/ts.exp)
    new_test["SentimentTitle"] <- news[match(test$IDLink,news$IDLink),]$SentimentTitle
    new_test["SentimentHeadline"] <- news[match(test$IDLink,news$IDLink),]$SentimentHeadline
    new_test["Now"] <- test[,ts.exp+1]
    new_test["Target"] <- test[,144+1]
    new_test <- new_test[new_test$Target>=0,]
    if(any(!complete.cases(new_test))) new_test[is.na(new_test)] <- 0
    
    m <- lm(Target ~ .,new_train)
    p <- predict(m,new_test)
    
    preds.row <- data.frame(TimeSlice=ts.exp,Prediction=as.numeric(p),Target=new_test$Target)
    if(!identical(names(preds.row), names(preds))) { names(preds.row) <- names(preds) }
    preds <- rbind(preds,preds.row)
    
  }
  
  preds$Prediction <- unlist(preds$Prediction)
  if(any(is.na(preds$Prediction))) preds[is.na(preds$Prediction),]$Prediction <- 0
  model.results <- c()
  
  for(a in 1:3) {
    
    aux <- preds[preds$TimeSlice==a,]
    aux <- aux[aux$Target>=0 & !is.na(aux$Target),]
    aux <- aux[aux$Prediction>=0,]
    
    utility.stats <- eval.stats(train[train$TS144>0,]$TS144,aux$Target,aux$Prediction)
    
    model.results <- rbind(model.results,utility.stats)
    
  }
  
  #INFORMATION RETURN
  res <- list(evaluation=model.results)
  res
  
}

# Pinto et al. (2013)
ML <- function(form,train,test,...) {
  
  #Reduce to those which have a time difference bigger than 2 days with the lowest timestamp in the testset
  train <- train[difftime(as.POSIXct(news[news$IDLink==test[1,]$IDLink,]$PublishDate),news[match(train$IDLink,news$IDLink),]$PublishDate,unit="hours")>48,]
  
  preds <- data.frame(TimeSlice=numeric(0),Prediction=numeric(0),Target=numeric(0))
  
  for(ts.exp in 1:3) { #Prediction for the first 3 timeslices
    
    new_train <- as.data.frame(train[,1:(ts.exp+1)])
    new_train[new_train==-1] <- 0
    
    if(ts.exp>1) {
      for(it in 2:ts.exp) {
        new_train[,it+1] <- new_train[,it+1] - new_train[,it]
      }
    }
    
    new_train["Target"] <- train[,144+1]
    new_train <- new_train[new_train$Target>=0,]
    if(any(!complete.cases(new_train))) new_train[is.na(new_train)] <- 0
    
    new_test <- as.data.frame(test[,1:(ts.exp+1)])
    new_test[new_test==-1] <- 0
    
    if(ts.exp>1) {
      for(it in 2:ts.exp) {
        new_test[,it+1] <- new_test[,it+1] - new_test[,it]
      }
    }
    
    new_test["Target"] <- test[,144+1]
    new_test <- new_test[new_test$Target>=0,]
    if(any(!complete.cases(new_test))) new_test[is.na(new_test)] <- 0
    
    new_train$IDLink <- NULL
    new_test$IDLink <- NULL
    
    m <- lm(Target ~ .,new_train)
    p <- predict(m,new_test)
    
    preds.row <- data.frame(TimeSlice=ts.exp,Prediction=as.numeric(p),Target=new_test$Target)
    if(!identical(names(preds.row), names(preds))) { names(preds.row) <- names(preds) }
    preds <- rbind(preds,preds.row)
    
  }
  
  preds$Prediction <- unlist(preds$Prediction)
  if(any(is.na(preds$Prediction))) preds[is.na(preds$Prediction),]$Prediction <- 0
  model.results <- c()
  
  for(a in 1:3) {
    
    aux <- preds[preds$TimeSlice==a,]
    aux <- aux[aux$Target>=0 & !is.na(aux$Target),]
    aux <- aux[aux$Prediction>=0,]
    
    utility.stats <- eval.stats(train[train$TS144>0,]$TS144,aux$Target,aux$Prediction)
    
    model.results <- rbind(model.results,utility.stats)
    
  }
  
  #INFORMATION RETURN
  res <- list(evaluation=model.results)
  res
  
}

# Pinto et al. (2013) (Radial Basis Functions)
MRBF <- function(form,train,test,...) {
  
  require(RSNNS)
  
  #Reduce to those which have a time difference bigger than 2 days with the lowest timestamp in the testset
  train <- train[difftime(as.POSIXct(news[news$IDLink==test[1,]$IDLink,]$PublishDate),news[match(train$IDLink,news$IDLink),]$PublishDate,unit="hours")>48,]
  
  preds <- data.frame(TimeSlice=numeric(0),Prediction=numeric(0),Target=numeric(0))
  
  for(ts.exp in 1:3) { #Prediction for the first 3 timeslices
    
    new_train <- as.data.frame(train[,1:(ts.exp+1)])
    new_train[new_train==-1] <- 0
    
    if(ts.exp>1) {
      for(it in 2:ts.exp) {
        new_train[,it+1] <- new_train[,it+1] - new_train[,it]
      }
    }
    
    new_train["Target"] <- train[,144+1]
    new_train <- new_train[new_train$Target>=0,]
    if(any(!complete.cases(new_train))) new_train[is.na(new_train)] <- 0
    
    new_test <- as.data.frame(test[,1:(ts.exp+1)])
    new_test[new_test==-1] <- 0
    
    if(ts.exp>1) {
      for(it in 2:ts.exp) {
        new_test[,it+1] <- new_test[,it+1] - new_test[,it]
      }
    }
    
    new_test["Target"] <- test[,144+1]
    new_test <- new_test[new_test$Target>=0,]
    if(any(!complete.cases(new_test))) new_test[is.na(new_test)] <- 0
    
    new_train$IDLink <- NULL
    new_test$IDLink <- NULL
    
    m <- rbf(as.data.frame(new_train[,-which(colnames(new_train)=="Target")]),new_train$Target,size=5,maxit=200,linOut = TRUE)
    p <- predict(m,as.data.frame(new_train[,-which(colnames(new_train)=="Target")]))
    new_train["RBF"] <- p
    p <- predict(m,as.data.frame(new_test[,-which(colnames(new_test)=="Target")]))
    new_test["RBF"] <- p
    
    if(any(!complete.cases(new_train))) new_train[is.na(new_train)] <- 0
    if(any(!complete.cases(new_test))) new_test[is.na(new_test)] <- 0
    
    m <- lm(Target ~ .,new_train)
    p <- predict(m,new_test)
    
    preds.row <- data.frame(TimeSlice=ts.exp,Prediction=as.numeric(p),Target=new_test$Target)
    if(!identical(names(preds.row), names(preds))) { names(preds.row) <- names(preds) }
    preds <- rbind(preds,preds.row)
    
  }
  
  preds$Prediction <- unlist(preds$Prediction)
  if(any(is.na(preds$Prediction))) preds[is.na(preds$Prediction),]$Prediction <- 0
  model.results <- c()
  
  for(a in 1:18) {
    
    aux <- preds[preds$TimeSlice==a,]
    aux <- aux[aux$Target>=0 & !is.na(aux$Target),]
    aux <- aux[aux$Prediction>=0,]
    
    utility.stats <- eval.stats(train[train$TS144>0,]$TS144,aux$Target,aux$Prediction)
    
    model.results <- rbind(model.results,utility.stats)
    
  }
  
  #INFORMATION RETURN
  res <- list(evaluation=model.results)
  res
  
}

#############

exp <- performanceEstimation(c(PredTask(TS144 ~ .,news.economy),
                              PredTask(TS144 ~ .,news.microsoft),
                              PredTask(TS144 ~ .,news.obama),
                              PredTask(TS144 ~ .,news.palestine)),
                             c(Workflow("const.scale"),
                               Workflow("linear.log"),
                               Workflow("mc.kernel"),
                               Workflow("mc.knn"),
                               Workflow("asur"),
                               Workflow("ML"),
                               Workflow("MRBF")),
                             EstimationTask("totTime",method=MonteCarlo(nReps=20,szTrain=.5,szTest=.25))
)

#############

#Example to obtain results for the first task, concerning the first iteration and the first workflow.
getIterationsInfo(exp,workflow=1,task=1,it=1)$evaluation

