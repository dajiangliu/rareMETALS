cov.file <- c("study1.MetaCov.assoc.gz","study2.MetaCov.assoc.gz");
score.stat.file <- c("study1.MetaScore.assoc.gz","study2.MetaScore.assoc.gz")
res <- rareMETALS.single(score.stat.file,
                         cov.file=NULL,
                         range="19:11200093-11201275",
                         alternative="two.sided",
                         ix.gold=1,
                         callrate.cutoff=0,
                         hwe.cutoff=0)
res <- rareMETALS.range(score.stat.file,
                        cov.file,
                        range="19:11200093-11201275",
                        range.name="LDLR",
                        test = "GRANVIL",
                        maf.cutoff = 0.05,
                        alternative = c("two.sided"),
                        ix.gold = 1,
                        out.digits = 4,
                        callrate.cutoff = 0,
                        hwe.cutoff = 0,
                        max.VT = NULL)
 
