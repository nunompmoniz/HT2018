#UTILITY-SURFACE
surface.colors <- colorRampPalette(c('#FF6600','#FFCC00','#FFFF66','#CCFFFF','#99FFFF','#33FFFF','#0099FF'))(64)



#UTIL FOR ONE BUMP
utilOneBump <- function(train.y,y,preds=NULL,cf=1.5,thr=0.9,makePlot=FALSE) {
  
  require(UBL)  
  require(akima)
  
  if(is.null(preds)) preds <- y;
  if(any(is.na(preds))) preds[is.na(preds)] <- 0
  
  PHIs <- UBL::phi.control(train.y,method="extremes",coef=cf)
  ph <- UBL::phi(train.y,PHIs)
  ph.test <- UBL::phi(y,PHIs)
  
  phi1 <- min(train.y[ph>=thr])

  phi3 <- 0
  if(length(PHIs$control.pts)>6) {
    phi3 <- PHIs$control.pts[7]
  } else {
    phi3 <- PHIs$control.pts[4]
  }
  
  max.x <- ceiling(max(preds))
  max.y <- max(y)
  
  values <- sort(unique(y))
  if(any(preds>max(values)*3)) {
    preds[preds>max(values)*3] <- max(values)*3
  }
  
  trues <- as.numeric(y)
  if(any(trues>max(preds)*3)) {
    trues[trues>max(preds)*3] <- max(preds)*3
  }
  
  values.df <- NULL
  
  if(makePlot) {
    values.df <- data.frame(trues=min(trues):max(trues),preds=min(preds):max(preds))
  } else {
    values.df <- data.frame(trues=trues,preds=as.numeric(preds))
  }
  
  
  #############
  #FALSE POSITIVE
  
  x <- 0; y <- 0; u <- 0
  
  if(max.x > max.y) {
    
    x <- c(phi1,
           phi1,
           (phi1+phi1),
           max.x,
           max.x,
           phi1,
           (max.y/2),
           seq(max.y/2,max.x,by=5),
           seq(max.y/2,max.x,by=5))
    y <- c(phi1,
           0,
           phi1,
           phi1,
           0,
           (phi1/2),
           0,
           rep(phi1,length(seq(max.y/2,max.x,by=5))),
           rep(0,length(seq(max.y/2,max.x,by=5))))
    u <- c(thr,
           0, 
           0,
           0,
           -1,
           0,
           -1,
           rep(0,length(seq(max.y/2,max.x,by=5))),
           rep(-1,length(seq(max.y/2,max.x,by=5))))    
    
  } else {
    
    x <- c(phi1,
           phi1,
           (phi1+phi1),
           max.y,
           max.y,
           phi1,
           (max.y/2),
           seq(max.y/2,max.y,by=5),
           seq(max.y/2,max.y,by=5))
    y <- c(phi1,
           0,
           phi1,
           phi1,
           0,
           (phi1/2),
           0,
           rep(phi1,length(seq(max.y/2,max.y,by=5))),
           rep(0,length(seq(max.y/2,max.y,by=5))))
    u <- c(thr,
           0, 
           0,
           0,
           -1,
           0,
           -1,
           rep(0,length(seq(max.y/2,max.y,by=5))),
           rep(-1,length(seq(max.y/2,max.y,by=5))))
    
  }
  
  dat <- data.frame(x,y,u)
  dat <- unique(dat)
  
  x <- dat$x; y <- dat$y; u <- dat$u
  
  xo.values <- unique(sort(values.df[values.df$preds>=phi1,]$preds))
  yo.values <- unique(sort(values.df[values.df$trues<=phi1,]$trues))
  if(length(xo.values)==0) { if(max.x>max.y) { xo.values <- c(seq(phi1,max.x,by=5)) } else { xo.values <- c(seq(phi1,max.y,by=5)) } }
  if(length(yo.values)==0) { yo.values <- c(seq(0,phi1,by=5)) }
  
  fp <- interp(x,y,u,duplicate="mean",xo=xo.values,yo=yo.values)
  fp$z[which(fp$z < -1)] <- -1
  fp$z[which(fp$z > 1)] <- 1
  
  #FALSE NEGATIVE
  
  x <- 0; y <- 0; u <- 0
  
  if(max.x > max.y) {
    
    x <- c(0,
           0,
           phi1,
           phi1,
           phi1,
           (phi1/2))
    y <- c(phi1,
           max.x,
           phi1,
           max.x,
           (phi1+phi1),
           phi1)
    u <- c(0, 
           -1,
           thr,
           0,
           0,
           0)
    
  } else {
    
    x <- c(0,
           0,
           phi1,
           phi1,
           phi1,
           (phi1/2))
    y <- c(phi1,
           max.y,
           phi1,
           max.y,
           (phi1+phi1),
           phi1)
    u <- c(0, 
           -1,
           thr,
           0,
           0,
           0)
    
  }
  
  dat <- data.frame(x,y,u)
  dat <- unique(dat)
  
  x <- dat$x; y <- dat$y; u <- dat$u
  
  xo.values <- unique(sort(values.df[values.df$preds<=phi1,]$preds))
  yo.values <- unique(sort(values.df[values.df$trues>=phi1,]$trues))
  if(length(xo.values)==0) { xo.values <- seq(0,phi1,by=5) }
  if(length(yo.values)==0) { if(max.x>max.y) { yo.values <- seq(phi1,max.x,by=5) } else { yo.values <- seq(phi1,max.y,by=5) } }
  
  fn <- interp(x,y,u,duplicate="mean",xo=xo.values,yo=yo.values)
  fn$z[which(fn$z < -1)] <- -1
  fn$z[which(fn$z > 1)] <- 1
  
  #TRUE NEGATIVE
  phi_aux <- data.frame(points=trues[which(ph.test<=thr)], p=ph.test[ph.test<=thr])
  phi_aux <- phi_aux[with(phi_aux,order(phi_aux$points)),]
  phi_aux <- unique(phi_aux)
  
  x <- c(rep(phi1,length(fp$y)),
         0,
         phi_aux$points)
  
  y <- c(fp$y,
         phi1,
         phi_aux$points)
  
  u <- c(fp$z[1,],
         0,
         phi_aux$p)
  
  dat <- data.frame(x,y,u)
  dat <- unique(dat)
  
  x <- dat$x; y <- dat$y; u <- dat$u
  
  xo.values <- unique(sort(values.df[values.df$preds<=phi1,]$preds))
  yo.values <- unique(sort(values.df[values.df$trues<=phi1,]$trues))
  if(length(xo.values)==0) { xo.values <- seq(0,phi1,by=5) }
  if(length(yo.values)==0) { yo.values <- seq(0,phi1,by=5) }
  
  tn <- interp(x,y,u,duplicate="mean",xo=xo.values,yo=yo.values)
  tn$z[which(tn$z < -1)] <- -1
  tn$z[which(tn$z > 1)] <- 1
  
  #TRUE POSITIVE
  phi_aux <- data.frame(points=trues[which(ph.test>=thr)], p=ph.test[ph.test>=thr])
  phi_aux <- phi_aux[with(phi_aux,order(phi_aux$points)),]
  phi_aux <- unique(phi_aux)
  
  x <- 0; y <- 0; u <- 0
  
  if(max.x>max.y) {
    
    x <- c(rep(phi1,length(fn$y)),
           phi_aux$points,
           seq(phi3,max.x,by=5),
           max.x,
           phi1,
           (max.x/2),
           seq(max.x/2,max.x,by=5),
           (phi1+phi1))
    
    y <- c(fn$y,
           phi_aux$points,
           seq(phi3,max.x,by=5),
           phi1,
           max.x,
           phi1,
           rep(phi1,length(seq(max.x/2,max.x,by=5))),
           phi1)
    
    u <- c(fn$z[nrow(fn$z),],
           phi_aux$p,
           rep(1,length(seq(phi3,max.x,by=5))),
           0,
           0,
           0,
           rep(0,length(seq(max.x/2,max.x,by=5))),
           0)
    
    
  } else {
    
    x <- c(rep(phi1,length(fn$y)),
           phi_aux$points,
           seq(phi3,max.y,by=5),
           max.y,
           phi1,
           (max.y/2),
           seq(max.y/2,max.y,by=5),
           (phi1+phi1))
    
    y <- c(fn$y,
           phi_aux$points,
           seq(phi3,max.y,by=5),
           phi1,
           max.y,
           phi1,
           rep(phi1,length(seq(max.y/2,max.y,by=5))),
           phi1)
    
    u <- c(fn$z[nrow(fn$z),],
           phi_aux$p,
           rep(1,length(seq(phi3,max.y,by=5))),
           0,
           0,
           0,
           rep(0,length(seq(max.y/2,max.y,by=5))),
           0)
    
    
  }
  
  dat <- data.frame(x,y,u)
  
  xo.values <- unique(sort(values.df[values.df$preds>=phi1,]$preds))
  yo.values <- unique(sort(values.df[values.df$trues>=phi1,]$trues))
  if(length(xo.values)==0) { if(max.x>max.y) { xo.values <- seq(phi1,max.x,by=5) } else { xo.values <- seq(phi1,max.y,by=5) } }
  if(length(yo.values)==0) { if(max.x>max.y) { yo.values <- seq(phi1,max.x,by=5) } else { yo.values <- seq(phi1,max.y,by=5) } }
  
  dat <- unique(dat)
  x <- dat$x; y <- dat$y; u <- dat$u
  
  tp <- interp(x,y,u,duplicate="mean",xo=xo.values,yo=yo.values)
  if(any(is.na(tp$z[,1]))) { tp$z[,1] <- fp$z[,ncol(fp$z)] }
  tp$z[which(tp$z < -1)] <- -1
  tp$z[which(tp$z > 1)] <- 1
  
  ####
  
  tests <- fp
  tests$x <- c(tn$x,tp$x)
  tests$y <- c(tn$y,tp$y)
  tst.m <- cbind(tn$z,fn$z)
  tst.m2 <- cbind(fp$z,tp$z)
  tst <- rbind(tst.m,tst.m2)
  tests$z <- tst
  
  if(any(table(tests$x)>1)) { 
    tests$z <- tests$z[!duplicated(tests$x),]
    tests$x <- unique(tests$x)
  }
  
  if(any(table(tests$y)>1)) { 
    tests$z <- tests$z[,!duplicated(tests$y)]
    tests$y <- unique(tests$y)
  }
  
  if(makePlot) {
    xl = expression(hat(Y))
    yl = expression(Y)
    image.plot(tests, col=surface.colors,xlab=xl,ylab=yl,legend.shrink = 1,legend.width = 0.8,cex.axis = 1,cex.main = 1.2)
    points(preds,trues,pch=10,cex=.5)
    contour(tests,lwd=0.6,drawlabels=T,labcex=0.7,add=T,cex=0.3,levels = seq(-1,1,by=0.2), nlevels = 1)
    abline(h=phi1,lty=2)
    abline(v=phi1,lty=2)
    box()
  }
  
  #####
  values.df$Utility <- 0
  for(i in 1:nrow(values.df)) {
    values.df[i,]$Utility <- tests$z[which(tests$x==values.df[i,]$preds),which(tests$y==values.df[i,]$trues)]
  }
  
  if(any(is.na(values.df$Utility))) {
    values.df[is.na(values.df$Utility),]$Utility <- 0
  }
  
  detach("package:akima", unload=TRUE)
  
  values.df
  
}

#EVALUATION FRAMEWORK
eval.stats <- function(trues.y,trues,preds,cf=3,thr=0.9,beta=1) {
  
  require(UBL)
  require(uba)
  
  ph <- UBL::phi.control(trues.y, method="extremes",coef=cf)
  ls <- uba::loss.control(trues.y)
  
  preds[preds<0] <- 0 #SOME ALGORITHMS COULD PRODUCE NEGATIVE PREDICTIONS ALTHOUGH THERE ARE NO NEGATIVE VALUES IN THE DATA SETS.
  preds[is.na(preds)] <- 0
  
  u_new <- utilOneBump(train.y=trues.y,y=trues,preds=preds,cf=cf,thr=thr)
  
  phi.trues <- UBL::phi(trues,control.parms = ph)
  phi.preds <- UBL::phi(preds,control.parms = ph)
  
  pr_frm <- u_new
  pr_frm["phiTrues"] <- phi.trues
  pr_frm["phiPreds"] <- phi.preds
  
  rmse= sqrt(mean((trues-preds)^2))
  
  rmse_phi= sqrt(mean(phi.trues[phi.trues>thr]*(trues[phi.trues>thr]-preds[phi.trues>thr])^2))
  
  prec <- utilOB(preds,trues,ph,ls,util.control(umetric="P",event.thr=thr),ut=u_new)
  rec  <- utilOB(preds,trues,ph,ls,util.control(umetric="R",event.thr=thr),ut=u_new)
  F1   <- utilOB(preds,trues,ph,ls,util.control(umetric="Fm",beta=beta,event.thr=thr),ut=u_new)
  
  prec.u <- sum(1+pr_frm[pr_frm$phiTrues>thr & pr_frm$phiPreds>thr,]$Utility)/sum(1+pr_frm[pr_frm$phiPreds>thr,]$phiPreds)
  rec.u <- sum(1+pr_frm[pr_frm$phiTrues>thr & pr_frm$phiPreds>thr,]$Utility)/sum(1+pr_frm[pr_frm$phiTrues>thr,]$phiTrues)
  F1.u <- (1+beta) * prec.u * rec.u / ( beta^2 * prec.u + rec.u)
  
  prec.c <- nrow(pr_frm[pr_frm$phiTrues>thr & pr_frm$phiPreds>thr,])/nrow(pr_frm[pr_frm$phiPreds>thr,])
  rec.c <- nrow(pr_frm[pr_frm$phiTrues>thr & pr_frm$phiPreds>thr,])/nrow(pr_frm[pr_frm$phiTrues>thr,])
  F1.c <- (1+beta) * prec.c * rec.c / ( beta^2 * prec.c + rec.c)
  
  c(
    rmse=rmse, rmse_phi=rmse_phi,prec.u=prec.u,rec.u=rec.u,F1.u=F1.u
  )
  
}
