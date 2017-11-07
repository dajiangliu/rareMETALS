
QC <- function(raw.data,QC.par,cov=1)
  {
    
      hwe.cutoff <- QC.par$hwe.cutoff;
      callrate.cutoff <- QC.par$callrate.cutoff;
      if(length(hwe.cutoff)==0) hwe.cutoff <- 0;
      if(length(callrate.cutoff)==0) callrate.cutoff <- 0;
      log.mat <- matrix("Untyped",nrow=length(raw.data$ref[[1]]),ncol=length(raw.data$ref));
      for(ii in 1:length(raw.data$hwe))
          {
              hwe.ii <- raw.data$hwe[[ii]];
              callrate.ii <- raw.data$callrate[[ii]];
              af.ii <- raw.data$af[[ii]];
              ix.rm <- which((callrate.ii<callrate.cutoff | hwe.ii<hwe.cutoff ) & af.ii!=0 & af.ii!=1);
              ix.hwe <- which(hwe.ii<hwe.cutoff & af.ii!=0 & af.ii!=1);
              ix.callrate <- which(callrate.ii<callrate.cutoff & af.ii!=0 & af.ii!=1);
              if(length(ix.hwe)>0) {
                  log.mat[ix.hwe,ii] <- "HWE";

              }
              if(length(ix.callrate)>0) {
                  log.mat[ix.callrate,ii] <- "CallRate";
              }
              ix.bug <- integer(0);
              if(cov==1)
                  {

                      res.diag <- rm.na(diag(raw.data$cov[[ii]])*raw.data$nSample[[ii]]);
                      res.vstat <- rm.na((raw.data$vstat[[ii]])^2)
                      diag.diff <- abs(res.diag-res.vstat);
                      ix.bug <- which(diag.diff>0.5);
                      msg <- paste(c('study',ii,'different missingness among variants'),sep=" ",collapse=' ');
                      
                      if(length(ix.bug)>0)
                          {
                              log.mat[ix.bug,ii] <- "bug";
                              warning(msg);
                          }
                      
                  }
              ix.rm <- unique(ix.rm);              
              if(length(ix.rm)>0)
                  {
                      raw.data$ustat[[ii]][ix.rm] <- NA;
                      raw.data$vstat[[ii]][ix.rm] <- NA;
                      if(cov==1)
                          {
                              raw.data$cov[[ii]][ix.rm,] <- NA;
                              raw.data$cov[[ii]][,ix.rm] <- NA;
                          }
                      raw.data$ref[[ii]][ix.rm] <- NA;
                      raw.data$alt[[ii]][ix.rm] <- NA;
                      raw.data$nSample[[ii]][ix.rm] <- NA;
                      raw.data$af[[ii]][ix.rm] <- NA;
                      raw.data$ac[[ii]][ix.rm] <- NA;
                      raw.data$nref[[ii]][ix.rm] <- NA;
                      raw.data$nhet[[ii]][ix.rm] <- NA;
                      raw.data$nalt[[ii]][ix.rm] <- NA;
                      raw.data$effect[[ii]][ix.rm] <- NA;
                      raw.data$pVal[[ii]][ix.rm] <- NA; 
                  }
          }
      raw.data$log.mat <- log.mat;
      return(raw.data);
  }

#' This is the function for flipping alleles
#'
#' @param raw.data The input datasets to be considered flipped
#' @param raw.data.ori The input datasets to be considered flipped
#' @param refaltList The list consists of ref, alt, pos, af and af.diff.max, as well as the option of whether throw away sites with large af.differences checkAF;
#' @param ix.pop index of the population
#' @param ix.var index of the variant;
#' @param log.mat.var The log for QC procedure;
#' @param correctFlip Correct for score and covariance matrices for flipped alleles;
#' @return A list consist of modified raw.data, ix.include and log.mat.var
#' @export
flipAllele <- function(raw.data,raw.data.ori,refaltList,ix.pop,ix.var,log.mat.var,correctFlip=TRUE,analyzeRefAltListOnly=TRUE)
    {
        ii <- ix.pop;
        if(is.na(raw.data$af[[ii]][ix.var])) {
            raw.data$af[[ii]][ix.var] <- (raw.data$nalt[[ii]][ix.var]+1/2*raw.data$nhet[[ii]][ix.var])/raw.data$nSample[[ii]][ix.var];
        }
        if((is.na(refaltList$ref[ix.var]) | is.na(refaltList$alt[ix.var])) & analyzeRefAltListOnly )
            {
                ii <- ix.pop;
                
                ix.include <- rep(0,length(raw.data$ustat));
                
                if(length(raw.data$cov)>0)
                    {
                        raw.data$cov[[ii]][ix.var,] <- NA;
                        raw.data$cov[[ii]][,ix.var] <- NA;
                    }
                log.mat.var[ii] <- "NotInRefAltList";
                raw.data$nSample[[ii]][ix.var] <- NA;
                raw.data$af[[ii]][ix.var] <- NA;
                raw.data$ac[[ii]][ix.var] <- NA;
                raw.data$ustat[[ii]][ix.var] <- NA;
                raw.data$vstat[[ii]][ix.var] <- NA;
                raw.data$nref[[ii]][ix.var] <- NA;
                raw.data$nhet[[ii]][ix.var] <- NA;
                raw.data$nalt[[ii]][ix.var] <- NA;                
                return(list(raw.data=raw.data,
                            log.mat.var=log.mat.var,
                            ix.include=ix.include));
            }

        if((is.na(refaltList$ref[ix.var]) | is.na(refaltList$alt[ix.var])) & !analyzeRefAltListOnly )
            {
                ii <- ix.pop;                
                ix.include <- rep(0,length(raw.data$ustat));
                log.mat.var[ii] <- "NotInRefAltList";
                return(list(raw.data=raw.data,
                            log.mat.var=log.mat.var,
                            ix.include=ix.include));
            }


        
        ii <- ix.pop;
        ref.gold <- refaltList$ref;alt.gold <- refaltList$alt;af.gold <- refaltList$af;af.diff.max <- refaltList$af.diff.max;checkAF <- refaltList$checkAF;
        if(length(checkAF)==0) checkAF <- FALSE;
        
        af.diff <- abs(raw.data$af[[ii]][ix.var]-af.gold[ix.var]);
        ix.include <- rep(0,length(raw.data$ustat));
        if(is.na(af.diff)) af.diff <- 0;
        if(!is.na(raw.data.ori$ref[[ii]][ix.var]) | !is.na(raw.data.ori$alt[[ii]][ix.var]))
            {
                if(rm.na(raw.data$af[[ii]][ix.var])==0 | rm.na(raw.data$af[[ii]][ix.var])==1)
                    {

                        if(rm.na(af.gold[ix.var])<.5)
                            {
                                ix.include[ii] <- 1;
                                raw.data$ustat[[ii]][ix.var] <- 0;
                                raw.data$af[[ii]][ix.var] <- 0;
                                raw.data$vstat[[ii]][ix.var] <- 0;
                                log.mat.var[ii] <- "Monomorphic";
                            }
                        if(rm.na(af.gold[ix.var])>=.5)
                            {
                                raw.data$af[[ii]][ix.var] <- 1;
                                if(length(raw.data$cov)>0)
                                    {
                                        raw.data$cov[[ii]][ix.var,] <- (-1)*raw.data$cov[[ii]][ix.var,];
                                        raw.data$cov[[ii]][,ix.var] <- (-1)*raw.data$cov[[ii]][,ix.var]
                                    }
                                raw.data$ustat[[ii]][ix.var] <- 0;
                                raw.data$vstat[[ii]][ix.var] <- 0;
                                nref.tmp <- (raw.data$nalt[[ii]][ix.var]);
                                nalt.tmp <- (raw.data$nref[[ii]][ix.var]);
                                nhet.tmp <- (raw.data$nhet[[ii]][ix.var]);
                                raw.data$nref[[ii]][ix.var] <- nref.tmp;
                                raw.data$nalt[[ii]][ix.var] <- nalt.tmp;
                                raw.data$nhet[[ii]][ix.var] <- nhet.tmp;
                                ix.include[ii] <- 1;
                                log.mat.var[ii] <- "Monomorphic";

                            }
                    }           
            }
        
        if(!is.na(raw.data$ref[[ii]][ix.var]) & !is.na(raw.data$alt[[ii]][ix.var]) & !(rm.na(raw.data$af[[ii]][ix.var])==0 | rm.na(raw.data$af[[ii]][ix.var])==1))
            {
                strandAmbiguous <- (((ref.gold[ix.var]=="A") & (alt.gold[ix.var]=="T")) | ((ref.gold[ix.var]=="T") & (alt.gold[ix.var]=="A")) | ((ref.gold[ix.var]=="C") & (alt.gold[ix.var]=="G")) | ((ref.gold[ix.var]=="G") & (alt.gold[ix.var]=="C")));
                flip.ref.alt <- (raw.data$ref[[ii]][ix.var]==refaltList$alt[ix.var] & raw.data$alt[[ii]][ix.var]==refaltList$ref[ix.var]);
                match.ref.alt <- (refaltList$ref[ix.var]==(raw.data$ref[[ii]][ix.var]) & refaltList$alt[ix.var]==raw.data$alt[[ii]][ix.var]);
                mono <- (((raw.data$ref[[ii]][ix.var]==".") | (raw.data$ref[[ii]][ix.var]==0) | (raw.data$alt[[ii]][ix.var]==".") | (raw.data$alt[[ii]][ix.var]==0) ) & (rm.na(raw.data$af[[ii]][ix.var])==0 | rm.na(raw.data$af[[ii]][ix.var])==1) & !is.na(raw.data.ori$nSample[[ii]][ix.var]))
                if(!match.ref.alt & !correctFlip &!mono)
                    {
                        
                        ix.include[ii] <- 1;
                        if(length(raw.data$cov)>0)
                            {
                                raw.data$cov[[ii]][ix.var,] <- NA;
                                raw.data$cov[[ii]][,ix.var] <- NA;
                            }
                        
                        log.mat.var[ii] <- "MismatchRemove";
                        
                        raw.data$nSample[[ii]][ix.var] <- NA;
                        raw.data$af[[ii]][ix.var] <- NA;
                        raw.data$ac[[ii]][ix.var] <- NA;
                        raw.data$ustat[[ii]][ix.var] <- NA;
                        raw.data$vstat[[ii]][ix.var] <- NA;
                        raw.data$nref[[ii]][ix.var] <- NA;
                        raw.data$nhet[[ii]][ix.var] <- NA;
                        raw.data$nalt[[ii]][ix.var] <- NA;
                    }

                if(match.ref.alt & (!strandAmbiguous))
                    {
                        ix.include[ii] <- 1;
                        log.mat.var[ii] <- "Match";
                    }

                if(flip.ref.alt & (!strandAmbiguous))
                {
                    raw.data$af[[ii]][ix.var] <- 1-raw.data$af[[ii]][ix.var];
                    raw.data$ustat[[ii]][ix.var] <- (-1)*(raw.data$ustat[[ii]][ix.var]);
                    if(length(raw.data$cov)>0)
                    {
                        raw.data$cov[[ii]][ix.var,] <- (-1)*raw.data$cov[[ii]][ix.var,];
                        raw.data$cov[[ii]][,ix.var] <- (-1)*raw.data$cov[[ii]][,ix.var]
                    }
                    
                    nref.tmp <- (raw.data$nalt[[ii]][ix.var]);
                    nalt.tmp <- (raw.data$nref[[ii]][ix.var]);
                    nhet.tmp <- (raw.data$nhet[[ii]][ix.var]);
                    raw.data$nref[[ii]][ix.var] <- nref.tmp;
                    raw.data$nalt[[ii]][ix.var] <- nalt.tmp;
                    raw.data$nhet[[ii]][ix.var] <- nhet.tmp;
                    ix.include[ii] <- 1;
                    
                    log.mat.var[ii] <- "FlipRefAlt";
                    
                }
                af.diff.min <- 0.05;
                if(flip.ref.alt & strandAmbiguous)
                    {                        
                        if(af.diff<=af.diff.min)
                            {
                                ix.include[ii] <- 1;
                                ## if(length(raw.data$cov)>0)
                                ##     {
                                ##         raw.data$cov[[ii]][ix.var,] <- NA;
                                ##         raw.data$cov[[ii]][,ix.var] <- NA;
                                ##     }
                                
                                log.mat.var[ii] <- "FlipStrand";
                                
                                ## raw.data$nSample[[ii]][ix.var] <- NA;
                                ## raw.data$af[[ii]][ix.var] <- NA;
                                ## raw.data$ac[[ii]][ix.var] <- NA;
                                ## raw.data$ustat[[ii]][ix.var] <- NA;
                                ## raw.data$vstat[[ii]][ix.var] <- NA;
                                ## raw.data$nref[[ii]][ix.var] <- NA;
                                ## raw.data$nhet[[ii]][ix.var] <- NA;
                                ## raw.data$nalt[[ii]][ix.var] <- NA;                                
                                
                            }
                        if(af.diff>af.diff.max)
                            {
                                raw.data$af[[ii]][ix.var] <- 1-raw.data$af[[ii]][ix.var];
                                raw.data$ustat[[ii]][ix.var] <- (-1)*(raw.data$ustat[[ii]][ix.var]);
                                if(length(raw.data$cov)>0)
                                    {
                                        
                                        raw.data$cov[[ii]][ix.var,] <- (-1)*raw.data$cov[[ii]][ix.var,];
                                        raw.data$cov[[ii]][,ix.var] <- (-1)*raw.data$cov[[ii]][,ix.var]

                                    }
                                nref.tmp <- (raw.data$nalt[[ii]][ix.var]);
                                nalt.tmp <- (raw.data$nref[[ii]][ix.var]);
                                nhet.tmp <- (raw.data$nhet[[ii]][ix.var]);
                                raw.data$nref[[ii]][ix.var] <- nref.tmp;
                                raw.data$nalt[[ii]][ix.var] <- nalt.tmp;
                                raw.data$nhet[[ii]][ix.var] <- nhet.tmp;
                                ix.include[ii] <- 1;
                                
                                log.mat.var[ii] <- "FlipRefAlt";
                            }
                    }

                if(match.ref.alt & strandAmbiguous)
                    {
                        if(af.diff<=af.diff.min)
                            {
                                ix.include[ii] <- 1;
                                
                                log.mat.var[ii] <- "Match";
                            }
                        if(af.diff>af.diff.max)
                            {
                                log.mat.var[ii] <- "FlipStrand";
                                if(length(raw.data$cov)>0)
                                    {
                                        raw.data$cov[[ii]][ix.var,] <- NA;
                                        raw.data$cov[[ii]][,ix.var] <- NA;
                                    }
                                
                                raw.data$af[[ii]][ix.var] <- 1-raw.data$af[[ii]][ix.var]
                                raw.data$ac[[ii]][ix.var] <- 2*raw.data$nSample[[ii]][ix.var]-raw.data$ac[[ii]][ix.var];
                                raw.data$ustat[[ii]][ix.var] <- raw.data$ustat[[ii]][ix.var]*(-1);
                                tmp <- raw.data$nref[[ii]][ix.var];
                                raw.data$nref[[ii]][ix.var] <- 2*raw.data$nSample[[ii]][ix.var]-raw.data$nref[[ii]][ix.var];
                                raw.data$nalt[[ii]][ix.var] <- tmp;
                            }
                    }                
                if(mono)
                    {
                        if(rm.na(af.gold[ix.var])<.5)
                            {
                                ix.include[ii] <- 1;
                                log.mat.var[ii] <- "Monomorphic";
                                raw.data$ustat[[ii]][ix.var] <- 0;
                                raw.data$af[[ii]][ix.var] <- 0;
                            }
                        if(rm.na(af.gold[ix.var])>=.5)
                            {
                                raw.data$af[[ii]][ix.var] <- 1-raw.data$af[[ii]][ix.var];
                                raw.data$ustat[[ii]][ix.var] <- 0;
                                if(length(raw.data$cov)>0)
                                    {
                                        raw.data$cov[[ii]][ix.var,] <- (-1)*raw.data$cov[[ii]][ix.var,];
                                        raw.data$cov[[ii]][,ix.var] <- (-1)*raw.data$cov[[ii]][,ix.var]
                                    }
                                
                                nref.tmp <- (raw.data$nalt[[ii]][ix.var]);
                                nalt.tmp <- (raw.data$nref[[ii]][ix.var]);
                                nhet.tmp <- (raw.data$nhet[[ii]][ix.var]);
                                raw.data$nalt[[ii]][ix.var] <- 2*raw.data$nSample[[ii]][ix.var];
                                raw.data$nref[[ii]][ix.var] <- 0;
                                raw.data$nhet[[ii]][ix.var] <- 0;
                                ix.include[ii] <- 1;
                                log.mat.var[ii] <- "Monomorphic";

                            }                        
                    }
                if(!match.ref.alt & !flip.ref.alt & !mono)
                    {
                        ix.include[ii] <- 1;
                        if(length(raw.data$cov)>0)
                            {
                                raw.data$cov[[ii]][ix.var,] <- NA;
                                raw.data$cov[[ii]][,ix.var] <- NA;
                            }
                        
                        log.mat.var[ii] <- "Unmatched";
                        raw.data$nSample[[ii]][ix.var] <- NA;
                        raw.data$af[[ii]][ix.var] <- NA;
                        raw.data$ac[[ii]][ix.var] <- NA;
                        raw.data$ustat[[ii]][ix.var] <- NA;
                        raw.data$vstat[[ii]][ix.var] <- NA;
                        raw.data$nref[[ii]][ix.var] <- NA;
                        raw.data$nhet[[ii]][ix.var] <- NA;
                        raw.data$nalt[[ii]][ix.var] <- NA;
                    }
                if(checkAF==TRUE)
                    {
                        af.diff.new <- abs(raw.data$af[[ii]][ix.var]-af.gold[ix.var]);
                        if(is.na(af.diff.new)) af.diff.new <- 0;
                        if(af.diff.new>af.diff.max)
                            {
                                ix.include[ii] <- 1;
                                log.mat.var[ii] <- "DiffAF";
                                raw.data$nSample[[ii]][ix.var] <- NA;
                                raw.data$af[[ii]][ix.var] <- NA;
                                raw.data$ac[[ii]][ix.var] <- NA;
                                raw.data$ustat[[ii]][ix.var] <- NA;
                                raw.data$vstat[[ii]][ix.var] <- NA;
                                
                                raw.data$nref[[ii]][ix.var] <- NA;
                                raw.data$nhet[[ii]][ix.var] <- NA;
                                raw.data$nalt[[ii]][ix.var] <- NA;                                
                            }
                    }
     
            }
        return(list(raw.data=raw.data,
                    log.mat.var=log.mat.var,
                    ix.include=ix.include));
    }


#' Impute missing summary association statistics assuming
#'
#' @param ustat.list the score statistics;
#' @param vstat.list the vstat list;
#' @param cov.mat.list the list of the covariance matrix
#' @param N.mat the matrix of sample sizes;each row for a study and each column for a variant site;
#' @export
imputeConditional <- function(ustat.list,vstat.list,cov.mat.list,N.mat,beta.vec=NULL,ix.candidate,ix.known) {
    U.imp <- 0;nSample.U <- 0;
    covG <- matrix(0,nrow=nrow(cov.mat.list[[1]]),ncol=ncol(cov.mat.list[[1]]));
    nSample.covG <- matrix(1,nrow=nrow(cov.mat.list[[1]]),ncol=ncol(cov.mat.list[[1]]));
    N.mat.imp <- N.mat;
    U.meta <- 0;
    
    for(ii in 1:length(ustat.list))
    {
        N.mat.imp[ii,] <- max(rm.na(N.mat[ii,]));
        U.meta <- U.meta+rm.na(ustat.list[[ii]]);
        nSample.U <- nSample.U+rm.na(N.mat[ii,]);
        for(jj in 1:length(ustat.list[[1]]))
        {
            for(kk in 1:jj)
            {
                covG[jj,kk] <- covG[jj,kk]+rm.na(sqrt(as.numeric(N.mat[ii,jj])*as.numeric(N.mat[ii,kk]))*cov.mat.list[[ii]][jj,kk]);
                covG[kk,jj] <- covG[jj,kk];
                nSample.covG[jj,kk] <- nSample.covG[jj,kk]+sqrt(as.numeric(rm.na(N.mat[ii,jj]))*as.numeric(rm.na(N.mat[ii,kk])));
                nSample.covG[kk,jj] <- nSample.covG[jj,kk];
            }
        }
    }
    U.meta <- U.meta/nSample.U;
    U.XY <- U.meta[ix.candidate];
    U.ZY <- U.meta[ix.known];
    covG.ori <- covG;
    covG <- rm.na(covG/nSample.covG);
    V.XZ <- matrix(covG[ix.candidate,ix.known],nrow=length(ix.candidate),ncol=length(ix.known));
    V.ZZ <- matrix(covG[ix.known,ix.known],nrow=length(ix.known),ncol=length(ix.known));
    V.XX <- matrix(covG[ix.candidate,ix.candidate],nrow=length(ix.candidate),ncol=length(ix.candidate));
    conditional.ustat <- U.XY-V.XZ%*%ginv(V.ZZ)%*%U.ZY;
    beta.ZY <- ginv(V.ZZ)%*%U.ZY;
    scaleMat <- as.matrix(diag(as.matrix(nSample.covG))%*%t(diag(as.matrix(nSample.covG))));

    var.U.XY <- covG.ori[ix.candidate,ix.candidate]/(scaleMat[ix.candidate,ix.candidate]);
    var.U.ZY <- covG.ori[ix.known,ix.known]/(scaleMat[ix.known,ix.known]);
    cov.U.XY.U.ZY <- covG.ori[ix.candidate,ix.known]/matrix(scaleMat[ix.candidate,ix.known],nrow=length(ix.candidate),ncol=length(ix.known));   
    conditional.V <- var.U.XY+V.XZ%*%ginv(V.ZZ)%*%var.U.ZY%*%ginv(V.ZZ)%*%t(V.XZ)-cov.U.XY.U.ZY%*%t(V.XZ%*%ginv(V.ZZ))-(V.XZ%*%ginv(V.ZZ))%*%t(cov.U.XY.U.ZY);
    
    sigma.sq.est <- abs(1-(t(U.ZY)%*%ginv(V.ZZ)%*%U.ZY));
    conditional.V <- conditional.V*as.numeric(sigma.sq.est);
    lambda <- .1; 
    conditional.V <- regMat(conditional.V,lambda);
    N.out <- sum(apply(N.mat.imp,1,max,na.rm=T));
    conditional.ustat <- conditional.ustat*N.out;
    conditional.V <- conditional.V*N.out^2;
    return(list(conditional.ustat=conditional.ustat,
                conditional.V=conditional.V,
                U.ZY=U.ZY,
                var.U.ZY=var.U.ZY,
                N=N.out,
                sigma.sq.est=sigma.sq.est));
}
regMat <- function(M,lambda) {
    cor.tmp <- rm.na(cov2cor(M));
    sd.mat <- matrix(0,nrow=nrow(M),ncol=ncol(M));
    id.mat <- matrix(0,nrow=nrow(M),ncol=ncol(M));
    diag(id.mat) <- 1;
    diag(sd.mat) <- sqrt(abs(diag(M)));
    cor.tmp <- cor.tmp+lambda*id.mat;
    M.reg <- sd.mat%*%(cor.tmp)%*%sd.mat;
    return(M.reg);
}
#' read rvtests and raremetal format;
#'
#' @param score.stat.file The score statistics file;
#' @return the formatted stat;
#' @export
readData <- function(score.stat.file) {
    res.tmp <- read.table(file=gzfile(score.stat.file),header=TRUE,as.is=TRUE,comment.char="",fill=TRUE,nrows=1000,sep='\n');
    ix.header <- grep("CHROM",res.tmp[,1]);
    res.header <- unlist(strsplit(res.tmp[ix.header,],split="\t"));
    if(res.header[1]=="#CHROM") {
            res <- read.table(file=gzfile(score.stat.file),header=FALSE,as.is=TRUE,fill=TRUE);
            colnames(res) <- gsub("#","",res.header);
        }
    if(res.header[1]=="CHROM") {
        res <- read.table(file=gzfile(score.stat.file),header=TRUE,as.is=TRUE,fill=TRUE);
        
    }
    ix.af <- which(res.header=="ALL_AF");
    if(length(ix.af)>0) {
        res.header[ix.af] <- "ALT_FREQ";
        colnames(res) <- gsub("#","",res.header);
    }
    return(res);
}

#' check sample overlaps between cohorts;
#'
#' @param score.stat.file
#' @return a list with estimated lambda values;
#' @export
checkOverlap <- function(score.stat.file1,score.stat.file2,...) {
    extraPar <- list(...);
    maf.cutoff <- extraPar$maf.cutoff;
    nVar <- extraPar$nVar;
    if(is.null(maf.cutoff)) maf.cutoff <- .1;
    if(is.null(nVar)) nVar <- 10000;
    res1 <- readData(score.stat.file1);
    res2 <- readData(score.stat.file2);
    pos1 <- paste(res1$CHROM,res1$POS,sep=":");
    pos2 <- paste(res2$CHROM,res2$POS,sep=":");
    pos.both <- intersect(pos1,pos2);
    res1 <- res1[match(pos.both,pos1),];
    res2 <- res2[match(pos.both,pos2),];
    maf1 <- res1$ALT_FREQ;
    maf1[maf1>.5] <- 1-maf1[maf1>.5];
    ix.common <- which(maf1>maf.cutoff);
    res1 <- res1[ix.common,];
    res2 <- res2[ix.common,];
    ixVar <- as.integer(seq(1,nrow(res1),length=nVar));
    stat1 <- res1$U_STAT[ixVar]/res1$SQRT_V_STAT[ixVar];
    stat2 <- res2$U_STAT[ixVar]/res2$SQRT_V_STAT[ixVar];
    stat.out <- median((stat1-stat2)^2,na.rm=TRUE)/qchisq(.5,df=1,lower.tail=FALSE)/2
    return(stat.out);
    
    
}


        
    
    


    
