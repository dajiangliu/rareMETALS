rvmeta.CMH.wss <- function(X.T.times.Y.centered.list,X.T.times.X.list,maf.vec,cov.mat.list,var.Y.list,N.list,alternative=c('two.sided','greater','less'),no.boot=0,alpha=0.05,weight=c('MB','MZ'))
  {
    if(length(alternative)>1) alternative <- "two.sided";
    if(length(weight)>1) weight <- "MB";
    res.list <- list();
    score.stat.numer <- 0;
    score.stat.denom.sq <- 0;
    if(weight=='MZ') w <- rep(1,length(X.T.times.Y.centered.list[[1]]));
    if(weight=='MB') {
      q <- ((maf.vec)+1/(2*N))*(2*N)/(2*N+2);
      w <- 1/sqrt(N*q*(1-q));
    }
    cov.mat.numer <- 0*diag(length(maf.vec));#this is to be used with parametric bootstrap:
    for(ii in 1:length(X.T.times.X.list))
      {
        score.stat.numer <- score.stat.numer+sum(w*X.T.times.Y.centered.list[[ii]]);
        score.stat.denom.sq <- score.stat.denom.sq+as.numeric(t(w)%*%(N.list[[ii]]*cov.mat.list[[ii]])%*%w)*(var.Y.list[[ii]]);
        cov.mat.numer <- cov.mat.numer+(cov.mat.list[[ii]])*(N.list[[ii]])*(var.Y.list[[ii]]);        
      }
    
    wss.stat <- as.numeric(score.stat.numer/sqrt(score.stat.denom.sq));
    if(no.boot==0) {
      if(alternative=='two.sided')
        {
          statistic <- wss.stat^2;
          p.value <- pnorm(-sqrt(statistic))*2;
        }
      if(alternative=='greater')
        {
          statistic <- wss.stat;
          p.value <- pnorm(-statistic);
        }
      if(alternative=='less')
        {
          statistic <- wss.stat;
          p.value <- pnorm(statistic);
        }
    }
    wss.stat.perm <- 0;
    if(no.boot>0)
      {
        X.T.times.Y.perm <- rmvnorm(no.boot,mean=rep(0,length(maf.vec)),sigma=cov.mat.numer,method='svd');
        wss.stat.perm <- as.vector(X.T.times.Y.perm%*%w)/as.numeric(score.stat.denom.sq);
        if(alternative=='two.sided')
          {p.value <- mean((wss.stat.perm^2)>(wss.stat^2))/2+mean((wss.stat.perm^2)>=(wss.stat^2))/2;
           statistic <- wss.stat^2;
         }
        if(alternative=='greater')
          {
            p.value <- mean((wss.stat.perm)>(wss.stat))/2+mean((wss.stat.perm)>=(wss.stat))/2;
            statistic <- wss.stat;
          }
        if(alternative=='less')
          {
            p.value <- mean(wss.stat.perm<wss.stat)/2+mean(wss.stat.perm<=wss.stat)/2;
            statistic <- wss.stat;
          }
      }
    no.site <- length(maf.vec);
    return(list(p.value=p.value,
                statistic=statistic,
                no.site=no.site));
  }
