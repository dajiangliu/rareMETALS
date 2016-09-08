library(getopt);
opt.spec <- matrix(c('file','a',2,'character',
                     'col','b',2,'integer',
                     'maf','d',2,'integer',
                     'cutoff','e',2,'double',
                     'out','c',2,'character'),
                   ncol=4,byrow=TRUE);
A <- getopt(opt.spec)
fname.out <- A$out;
if(length(A$cutoff)==0) A$cutoff <- 1;
res <- read.table(A$file,as.is=TRUE,fill=TRUE,header=TRUE);
##pvals <- as.numeric(res[,A$col]);
pvals <- as.numeric(res$PVALUE);
ix.na <- which(is.na(pvals));
print(length(ix.na));
if(length(ix.na)>0) {pvals <- pvals[-ix.na];res <- res[-ix.na,];}
summary(pvals);
if(length(A$maf)>0) {
  maf <- as.numeric(unlist(res[,A$maf]));
  ix.maf <- which(maf<A$cutoff | (1-maf)<A$cutoff);
  pvals <- pvals[ix.maf];
}
if(length(fname.out)==0) fname.out <- paste(A$file,'maf',A$cutoff,'tiff',sep='.',collapse='.');
main.txt <- paste(A$file,'maf',A$cutoff,sep="_",collapse="_");
qq.conf.beta <- function(pvals, num.disp = length(pvals), sz.conf = 0.95, xlim = NULL, ylim = NULL, drawaxis = TRUE, ptcolor = "black" ,main.txt='') {
 N <- length(pvals)
 Z <- qchisq(1-pvals,df=1);
 Z.median <- median(Z,na.rm=TRUE);
 lambda <- format(Z.median/qchisq(0.5,df=1),digits=3);
 M <- num.disp
 if ( M > N ) {
   stop("num.disp must be equal or smaller than length(pvals)")
 }
 if ( ( sz.conf <= 0 ) || ( sz.conf >= 1 ) ) {
   stop("sz.conf must be between 0 and 1")
 }
 pvals.sorted <- sort(pvals)
 lmean<-(1:N)/(N+1)
 lc <- qbeta(rep((1-sz.conf)/2,M),1:M,N+1-(1:M))
 lmax<-qbeta(rep((1+sz.conf)/2,M),1:M,N+1-(1:M))
 sub.txt <- paste('lambda=',lambda,sep="",collapse="");
 if ( drawaxis ) {
   plot(-log10(lmean)[M:1],-log10(pvals.sorted[M:1]),type="p",xlab="Expected ordered -log10(pvalue)",ylab="Observed ordered -log10(pvalue)", xlim = xlim, ylim = ylim,cex.lab = 1.5,pch = 10,cex.axis = 1.5,main=main.txt,sub=sub.txt);
   polygon(c(-log10(lmean)[M:1],-log10(lmean[1:M])),c(-log10(lc[M:1]),-log10(lmax[1:M])),col="gray",border="gray",lwd = 2)
   #polygon(c(-log10(lmean)[M:1],-log10(lmean[1:M])),c(-log10(lc[M:1]),-log10(lmax[1:M])),col="red",border=NA)
   abline(0,1)
 }
 else {
   par(new=TRUE)
 }
 points(-log10(lmean)[1:M],-log10(pvals.sorted[1:M]),pch=19,cex=.8,col=ptcolor)
}
tiff(file=fname.out,height=6,width=6,units='in',res=300,compression='lzw');
qq.conf.beta(pvals,main.txt=main.txt);
dev.off();


