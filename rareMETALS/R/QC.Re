
QC <- function(raw.data,QC.par,cov=1)
  {
    
      ######################################print(c("QC.par",QC.par));
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
                              ## ######################################print("bug");
                              ## ix.rm <- c(ix.rm,ix.bug);
                              ## ######################################print(c("study",ii));
                              ## ######################################print(res.diag);
                              ## ######################################print(res.vstat)
                              log.mat[ix.bug,ii] <- "bug";
                              warning(msg);
                          }
                      
                  }
              ix.rm <- unique(ix.rm);              
              if(length(ix.rm)>0)
                  {
                      raw.data$ustat[[ii]][ix.rm] <- NA;
                      raw.data$vstat[[ii]][ix.rm] <- NA;
                      ######################################print("out");
                      if(cov==1)
                          {
                              ######################################print("this is exec");
                              raw.data$cov[[ii]][ix.rm,ix.rm] <- NA;
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
        ############################print(c("REF","ALT",refaltList$ref[ix.var],refaltList$alt[ix.var]))
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
                ############################print("NotInRefAltList");
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

                        if(af.diff<=af.diff.max)
                            {
                                ix.include[ii] <- 1;
                                raw.data$ustat[[ii]][ix.var] <- 0;
                                raw.data$vstat[[ii]][ix.var] <- 0;
                                log.mat.var[ii] <- "Monomorphic";
                            }
                        if(af.diff>=af.diff.max)
                            {
                                raw.data$af[[ii]][ix.var] <- 1-raw.data$af[[ii]][ix.var];
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
                ##mono <- (((raw.data$ref[[ii]][ix.var]==".") | (raw.data$ref[[ii]][ix.var]==0) | (raw.data$alt[[ii]][ix.var]==".") | (raw.data$alt[[ii]][ix.var]==0) | (rm.na(raw.data$af[[ii]][ix.var])==0)| rm.na(raw.data$af[[ii]][ix.var])==1) & !is.na(raw.data.ori$nSample[[ii]][ix.var]))
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
                if(match.ref.alt & (!strandAmbiguous))
                    {
                        ix.include[ii] <- 1;
                        log.mat.var[ii] <- "Match";
                    }
                af.diff.min <- 0.05;
                if(flip.ref.alt & strandAmbiguous)
                    {                        
                        if(af.diff<=af.diff.min)
                            {
                                ix.include[ii] <- 1;
                                if(length(raw.data$cov)>0)
                                    {
                                        raw.data$cov[[ii]][ix.var,] <- NA;
                                        raw.data$cov[[ii]][,ix.var] <- NA;
                                    }
                                
                                log.mat.var[ii] <- "FlipStrand";
                                
                                raw.data$nSample[[ii]][ix.var] <- NA;
                                raw.data$af[[ii]][ix.var] <- NA;
                                raw.data$ac[[ii]][ix.var] <- NA;
                                raw.data$ustat[[ii]][ix.var] <- NA;
                                raw.data$vstat[[ii]][ix.var] <- NA;
                                raw.data$nref[[ii]][ix.var] <- NA;
                                raw.data$nhet[[ii]][ix.var] <- NA;
                                raw.data$nalt[[ii]][ix.var] <- NA;
                                
                                
                            }
                        if(af.diff>af.diff.min)
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
                        if(af.diff<=af.diff.max)
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
                
                if(mono)
                    {
                        if(af.diff<=af.diff.max)
                            {
                                ix.include[ii] <- 1;
                                log.mat.var[ii] <- "Monomorphic";
                                raw.data$ustat[[ii]][ix.var] <- 0;
                            }
                        if(af.diff>=af.diff.max)
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
                                raw.data$nref[[ii]][ix.var] <- nref.tmp;
                                raw.data$nalt[[ii]][ix.var] <- nalt.tmp;
                                raw.data$nhet[[ii]][ix.var] <- nhet.tmp;
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
                        ##########################print(c("raw.data$af[[ii]][ix.var]",raw.data$af[[ii]][ix.var]));
                        ##########################print(c("af.gold[ix.var]",af.gold[ix.var]));
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
imputeMeta <- function(ustat.list,vstat.list,cov.mat.list,N.mat,beta.vec=NULL) {
    U.imp <- 0;nSample.U <- 0;
    covG <- matrix(0,nrow=nrow(cov.mat.list[[1]]),ncol=ncol(cov.mat.list[[1]]));
    nSample.covG <- covG;
    
    for(ii in 1:length(ustat.list))
    {
        
        U.imp <- U.imp+rm.na(ustat.list[[ii]]);
        nSample.U <- nSample.U+rm.na(N.mat[ii,]);
        for(jj in 1:length(ustat.list[[1]]))
        {
            for(kk in 1:jj)
            {
                covG[jj,kk] <- covG[jj,kk]+rm.na(sqrt(N.mat[ii,jj]*N.mat[ii,kk])*cov.mat.list[[ii]][jj,kk]);
                covG[kk,jj] <- covG[jj,kk];
                nSample.covG[jj,kk] <- nSample.covG[jj,kk]+sqrt(rm.na(N.mat[ii,jj])*rm.na(N.mat[ii,kk]));
                nSample.covG[kk,jj] <- nSample.covG[jj,kk];
                ##if(is.na(covG[jj,kk])) stop();
            }
            
        }
    }
    covG.ori <- covG;    
    covG <- (covG/nSample.covG);
    ix.missing <- which(is.na(covG),arr.ind=TRUE);
    covG[ix.missing] <- median(covG,na.rm=TRUE);
    N.meta <- apply(nSample.covG,1,max,na.rm=T);
    N.meta <- rep(max(N.meta),length(U.imp));
    V.imp <- covG*max(N.meta);
    ustat.list.imp <- ustat.list;
    vstat.list.imp <- vstat.list;
    cov.mat.list.imp <- cov.mat.list;
    V.list <- list();
    impState <- matrix(0,nrow=length(ustat.list),ncol=length(ustat.list[[1]]));
    for(ii in 1:length(ustat.list)) {
        ix.miss <- which(is.na(ustat.list[[ii]]));
        if(length(ix.miss)>0) {
            impState[ii,ix.miss] <- 1;
            cov.mat.list.imp[[ii]][ix.miss,] <- covG[ix.miss,];
            cov.mat.list.imp[[ii]][,ix.miss] <- covG[,ix.miss];
            
            vstat.list.imp[[ii]][ix.miss] <- sqrt(diag(cov.mat.list.imp[[ii]])[ix.miss]*max(N.mat[ii,],na.rm=T));
            
            if(is.null(beta.vec)) {
                ##ustat.list.imp[[ii]][ix.miss] <- (ginv(covG)%*%(U.imp/nSample.U))[ix.miss]*(vstat.list.imp[[ii]][ix.miss])^2;
                         ustat.list.imp[[ii]][ix.miss] <- (U.imp[ix.miss])/(diag(covG.ori)[ix.miss])*(vstat.list.imp[[ii]][ix.miss])^2;
            }
            if(!is.null(beta.vec)) {
                ustat.list.imp[[ii]][ix.miss] <- beta.vec[ix.miss]*vstat.list.imp[[ii]][ix.miss]^2;
            }
        }
        V.list[[ii]] <- cov.mat.list.imp[[ii]]*max(N.mat[ii,],na.rm=T);
    }
    
    return(list(covG=covG,
                nSample.covG=nSample.covG,
                V.imp=V.imp,
                U.imp=U.imp,
                V.list=V.list,
                impState=impState,
                ustat.list.imp=ustat.list.imp,
                vstat.list.imp=vstat.list.imp,
                cov.mat.list.imp=cov.mat.list.imp,
                N.vec=apply(N.mat,1,max,na.rm=T),
                N.meta=N.meta));
}
 
