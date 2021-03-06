#' Format data into analyzable format;
#'
#' @param raw.data The data that is imported using seqminer;
#' @return a list of formatted data;
#' @export
rareMETALS.formatData <- function(raw.data) {
    pos <- raw.data$pos;
    pos <- gsub("_.*","",pos);
    pos.tab <- table(pos);
    posMulti <- names(pos.tab)[which(pos.tab>1)];
    ref.mat <- matrix(unlist(raw.data$ref),ncol=length(raw.data$ref),nrow=length(raw.data$ref[[1]]));
    alt.mat <- matrix(unlist(raw.data$alt),ncol=length(raw.data$ref),nrow=length(raw.data$ref[[1]]));
    uniq.allele <- function(x) {x.tab <- table(x);return(paste(names(x.tab),sep=',',collapse=','))}
    ref.tab <- apply(ref.mat,1,uniq.allele);
    alt.tab <- apply(alt.mat,1,uniq.allele);
    nSample.mat <- matrix(unlist(raw.data$nSample),ncol=length(raw.data$ref),nrow=length(raw.data$ref[[1]]));
    nref.mat <- matrix(unlist(raw.data$nref),ncol=length(raw.data$ref),nrow=length(raw.data$ref[[1]]));
    nalt.mat <- matrix(unlist(raw.data$nalt),ncol=length(raw.data$ref),nrow=length(raw.data$ref[[1]]));
    nhet.mat <- matrix(unlist(raw.data$nhet),ncol=length(raw.data$ref),nrow=length(raw.data$ref[[1]]));
    af.mat <- matrix(unlist(raw.data$af),ncol=length(raw.data$ref),nrow=length(raw.data$ref[[1]]));
    ustat.mat <- matrix(unlist(raw.data$ustat),ncol=length(raw.data$ref),nrow=length(raw.data$ref[[1]]));
    vstat.mat <- matrix(unlist(raw.data$vstat),ncol=length(raw.data$ref),nrow=length(raw.data$ref[[1]]));
    ref.mat <- matrix(unlist(raw.data$ref),ncol=length(raw.data$ref),nrow=length(raw.data$ref[[1]]));
    alt.mat <- matrix(unlist(raw.data$alt),ncol=length(raw.data$ref),nrow=length(raw.data$ref[[1]]));
    
    direction.mat <- ustat.mat;
    direction.mat[which(ustat.mat>0,arr.ind=T)] <- "+";
    direction.mat[which(ustat.mat<0,arr.ind=T)] <- "-";
    direction.mat[which(ustat.mat==0,arr.ind=T)] <- "=";
    direction.mat[which(is.na(ustat.mat),arr.ind=T)] <- "X";
    ref.alt.mat <- matrix(paste(ref.mat,alt.mat,sep=','),nrow=nrow(ref.mat),ncol=ncol(alt.mat));
    pos.ref.alt.mat <- apply(ref.alt.mat,2,function(vec1,pos) return(paste(pos,vec1,sep=',')),pos=pos);
    
    V.list <- list();
    for(ii in 1:length(raw.data$cov)) {
        V.list[[ii]] <- (raw.data$cov[[ii]])*mean(raw.data$nSample[[ii]]);
    }
    
    return(list(nref.mat=nref.mat,
                nalt.mat=nalt.mat,
                ref.tab=ref.tab,
                alt.tab=alt.tab,
                af.mat=af.mat,
                ref.mat=ref.mat,
                alt.mat=alt.mat,
                nhet.mat=nhet.mat,
                ustat.mat=ustat.mat,
                pos=pos,
                pos.ref.alt=paste(pos,paste(ref.tab,alt.tab,sep="/"),sep="_"),
                posMulti=posMulti,
                V.list=V.list,             
                nSample.mat=nSample.mat,
                ref.alt.mat=ref.alt.mat,
                pos.ref.alt.mat=pos.ref.alt.mat,
                direction.mat=direction.mat,
                vstat.mat=vstat.mat));
}
#' Function to perform multi-allelic association test
#'
#' @param score.stat.file The file names for score statistics files;
#' @param cov.file The file names for covariance matrix files;
#' @param tabixRange The vectors of tabix ranges;
#' @param rangeName The name of the genes that are represented by the corresponding tabixRange;
#' @param maf.cutoff The MAF cutoffs;
#' @param rvTest The rare variant association tests to be performed. Options include GRANVIL, SKAT or VT;
#' @return A matrix of results; each row being the result for a gene;
#' @export
rareMETALS.range.multiAllele <- function(score.stat.file,cov.file,tabixRange,rangeName,maf.cutoff,rvTest) {
    capture.output(raw.data.all <- rvmeta.readDataByRange( score.stat.file, cov.file, tabixRange));  
    
    res.list <- lapply(raw.data.all,assocGene,maf.cutoff=maf.cutoff,rvTest=rvTest);
    res.tmp <- do.call(rbind,res.list);
    res.out <- cbind(rangeName,res.tmp);
    colnames(res.out) <- c("GENE","STAT","PVALUE","NO_SITE");
    return(res.out);
}
#' Association analysis of gene-level test with multi-allele sites;
#'
#' @param raw.data raw.data from seqminer output;
#' @param maf.cutoff The minor allele freq cutoff
#' @param rvTest The rare variant association test
#' @return A list of results;
#' @export
assocGene <- function(raw.data,maf.cutoff,rvTest) {
    formattedData <- rareMETALS.formatData(raw.data);
    ustat.meta <- rowSums(formattedData$ustat.mat,na.rm=TRUE);
    vstat.sq.meta <- rowSums(formattedData$vstat.mat^2,na.rm=TRUE)
    nSample.meta <- rowSums((formattedData$nSample.mat),na.rm=TRUE);
    af.meta <- rowSums((formattedData$af.mat)*(formattedData$nSample.mat),na.rm=TRUE)/nSample.meta;
    V.meta <- 0;
    for(ii in 1:length(formattedData$V.list)) {
        V.meta <- V.meta+rm.na(formattedData$V.list[[ii]]);
    }
    maf <- af.meta;
    maf[which(maf>.5)] <- 1-maf[which(maf>.5)];
    ix.rare <- which(maf>0 & maf<maf.cutoff)
    pos.rare <- formattedData$pos[ix.rare];
    pos.ref.alt.rare <- formattedData$pos.ref.alt[ix.rare];
    pos.ref.alt.x2 <- (formattedData$pos.ref.alt)[(formattedData$pos %in% pos.rare) & !(formattedData$pos.ref.alt %in% pos.ref.alt.rare)];
    if(length(pos.ref.alt.rare)>0) {
        if(length(pos.ref.alt.x2)==0) {
            
            ix.x1 <- (1:length(ustat.meta))[formattedData$pos.ref.alt %in% pos.ref.alt.rare];
            ustat.meta.rare <- ustat.meta[ix.x1];
            V.meta.rare <- matrix(V.meta[ix.x1,ix.x1],nrow=length(ix.x1),ncol=length(ix.x1));
            maf.rare <- maf[ix.x1];
            ##calcGeneLevelTest <- function(ustat,V,maf.vec,N,test="GRANVIL") 
            res.gene <- calcGeneLevelTest(ustat.meta.rare,V.meta.rare,maf.rare,N=mean(nSample.meta),test=rvTest);
            
        }
        if(length(pos.ref.alt.x2)>0) {
            ##get.conditional.score.stat.RE <- function(ustat,X.T.times.X,N,ix.candidate,ix.known);
            ix.x1.tmp <- (1:length(ustat.meta))[formattedData$pos.ref.alt %in% pos.ref.alt.rare];
            ix.x2.tmp <- (1:length(ustat.meta))[formattedData$pos.ref.alt %in% pos.ref.alt.x2];
            ix.tmp <- c(ix.x1.tmp,ix.x2.tmp);
            ustat.meta.tmp <- ustat.meta[ix.tmp];
            V.meta.tmp <- matrix(V.meta[ix.tmp,ix.tmp],nrow=length(ix.tmp),ncol=length(ix.tmp));
            res.cond <- get.conditional.score.stat(ustat.meta.tmp,V.meta.tmp,mean(nSample.meta),1:length(ix.x1.tmp),length(ix.x1.tmp)+(1:length(ix.x2.tmp)));
            ustat.meta <- res.cond$conditional.ustat;
            V.meta <- res.cond$conditional.V;
            
            maf.meta <- maf[ix.x1.tmp];
            res.gene <- calcGeneLevelTest(ustat.meta,V.meta,maf.meta,N=mean(nSample.meta[ix.x1.tmp]),test=rvTest);
        }
        return(matrix(c(res.gene$statistic,res.gene$p.value,res.gene$no.site),nrow=1));
    }
    if(length(pos.ref.alt.rare)==0) return(matrix(rep(NA,3),nrow=1));
}
