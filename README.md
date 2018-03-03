# rareMETALS
**rareMETALS** is an R-package for performing single or gene-level tests for detecting rare variant associations. For questions regarding the use of this package, please contact Dajiang Liu([dajiang.liu at outlook dot com](mailto:dajiang.liu@outlook.com)) or Gonçalo Abecasis([goncalo at umich dot edu](mailto:goncalo@umich.edu)), or follow [Feedback/Contact](#Feedback-Contact).   
The same methodology is also implemented in command line tools. Please see [here](https://genome.sph.umich.edu/wiki/RAREMETAL_Documentation)

## **Table of Contents**

1. [Change Log](#change-log)
2. [Download and Installation](#download-and-installation)
3. [Documentation](#documatation)
4. [Forum](#forum) 
5. [Supported Functionalities](#Supported-Functionalities)  
6. [Exemplar Datase](#exemplar)
7. [How to Generate Summary Association Statistics and Prepare Them for Meta-analysis](#HTSASPTMA) 
8. [A Simple Tutorial for Using functions](#simple-tutorial-functions)  
  - [Using the rareMETALS.single function]  
  - [Using the rareMETALS.single.group function]  
  - [Using the rareMETALS.range function]  
  - [Using the rareMETALS.range.group function]  
  - [Using the conditional.rareMETALS.single]  
  - [Using the conditional.rareMETALS.range]  
9. [Feedback/Contact](#Feedback-Contact)

### Change Log <a name="change-log"></a>
- 05/18/2017 Version 6.8 incorporates a number of new features and bug fixes. We included support for multi-allelic variants, the support for a new conditional analysis method, the support for cohort level genomic controls, and the bug fixes for calculating heterogeneity statistics such as Q and I2.
- 04/09/2016 Version 6.3 is released. Minor bug fix: Due to different level of missingness of variants in the gene, the single variant association statistics calculated using the covariance matrices of score statistics can be different than single variant association statistics calculated using vstat. This has lead to confusions. It has been fixed in version 6.3. The primary results from version 6.2 should be correct.
- 09/25/2015 Version 6.2 is released. Minor bug fix: Removed the incorrect warning information in version 6.1 when quantitative traits are meta-analyzed. The software incorrectly consider it as binary trait and suggested the use of rareMETALS2.
- 07/23/2015 Version 6.1 is released. Minor feature changes include output for VT the sites where the statistics are maximized; fixed a bug for determining monomorphic sites. Issue warnings when rareMETALS is used to analyze binary trait for meta-analysis.
- 05/19/2015 Version 6.0 is released. Minor feature addition: rareMETALS can now output of the set of variants that are analyzed in VT (i.e. the set of variants with MAF < the threshold where the VT statistic is maximized).
- 04/01/2015 Version 5.9 is released (not a April's fool joke)! A bug in calculating Cochran-Q statistic is fixed. A bug in conditional.rareMETALS.range.group is also fixed. No other analyses are affected.
- 01/24/2015 Version 5.8 is released, which fixed a serious bug for single variant unconditional association tests with group file. If you happen to run the analyses using rareMETALS.single.group() in version 5.7, the results are likely to be incorrect. Please rerun using version 5.8. Please note only rareMETALS.single.group function is affected. All other functions should not be affected by this error.
- 01/04/2015 Version 5.7 is released, which added metrics for heterogeneity of genetic effects, including I2 and Q for single variant association statistics
- 12/09/2014 Version 5.6 is released, which added function conditional.rareMETALS.range.group, and fixed a minor issue for estimating sample sizes.
- 11/19/2014 Version 5.5 is released, which fixes a few bugs on the version 5.4.
- 11/09/2014 Version 5.4 is posted with the following change 1.) Allowing for performing conditional analysis for multiple candidate variants 2.) add option correctFlip to rareMETALS.single.group, rareMETALS.range.group allowing for options to discard sites with non-matching ref or alt alleles. Default is TRUE
- 09/08/2014 Version 5.2 is posted. One change in version 5.0 and 5.1 is reverted, which could lead to undesirable effect. It improves on some border line cases as compared to Versions 4.7 - 4.9. But in general, version 5.2 and 4.7-4.9 should give very comparable results. Please update to the latest version. I would expect that version 5.2 should run stably for all models under all circumstances.
- 08/21/2014 Version 4.9 is posted. A bug is fixed for VT test. While the p-values and statistics were correct, the number of sites and the beta estimate could sometimes be incorrect in version 4.8. Now it is fixed. Please download the newest version. Thanks!
- 08/18/2014 Version 4.8 is posted. A bug for recessive model analysis is fixed. Additive and dominant models should remain unaffected. Thanks!
- 08/06/2014 Version 4.7 is posted, where a few minor bugs were fixed. Thanks to Heather Highland and Xueling Sim for careful testing!! Please update. Thanks!
- 07/15/2014 Fixed a bug in conditional.rareMETALS.single and conditional.rareMETALS.range; Please update. Thanks!
- 06/27/2014 Updated to version 4.0: Many updates are implemented, including support for group files in both single variant and gene-level association test; checks for allele flips based upon variant frequency, the detection of possible allele flips using a novel statistic based upon variations of allele frequency between studies;

### Downlaod and installation <a name="download-and-installation"></a>
- [Option one] The R package can be downloaded from [rareMETALS_6.8.tar.gz](https://genome.sph.umich.edu/w/images/b/b2/RareMETALS_6.8.tar.gz). It will be eventually released on the Comprehensive R-archive Network. If you want to perform gene-level association test using automatically generated annotations, you will also need [refFlat_hg19.txt.gz](https://genome.sph.umich.edu/w/images/0/04/RefFlat_hg19.txt.gz), which is the gene definition modified from refFlat.
   - To install the package, please use `R CMD INSTALL rareMETALS_XXX.tar.gz` command, where **XXX** is the version number for rareMETALS

- [Option two] Right now, the package is hosted on github, which allows installation and update to be very easy. 
  - First, make sure you have the `mvtnorm` and `data.table` packages installed:
  ```
  install.packages("devtools")
  ```
  - Then you could use:
   ```
   library(devtools)   
   install_github("dajiangliu/rareMETALS")
   ```
   - With `library(rareMETALS)`, your are ready to go!


### Documentation <a name=documatation></a>
An R automatically generated documentation is available here: [rareMETALS-manual.pdf](https://genome.sph.umich.edu/w/images/4/44/RareMETALS-manual.pdf). Please note that it is still rough in places. Please let us know if you see any problems. Thanks!


### Forum <a name=forum></a>
I have created a google group for discussion on the usage and for bug reports etc. As you can see, there are numerous updates to the package since its release, thanks to the valuable suggestions from many users. We are committed to continue to update the package and improve its functionalities. If you find any issues to the package and think that the discussions may also benefit others, please post them to the user group. Here is the link to the [discussion group](https://groups.google.com/forum/#!forum/raremetals)


### Supported Functionalities <a name=Supported-Functionalities></a>  
- Marginal meta-analysis of single variant or gene-level association test
- Conditional analysis of single variant or gene-level association, for variants (gene) where there are covariance information available between candidate variants and known variants.
- Estimates of genetic effects and locus genetic variance
- Estimate measures of genetic effect heterogeneities between studies


### Exemplar Dataset <a name=exemplar></a>
Four datasets are useful to get you started on how to use rareMETALS R package for meta-analyses of gene-level association test
- Media:[study1.MetaScore.assoc.gz](https://genome.sph.umich.edu/w/images/0/0e/Study1.MetaScore.assoc.gz)
- Media:[study2.MetaScore.assoc.gz](https://genome.sph.umich.edu/w/images/e/e5/Study2.MetaScore.assoc.gz)
- Media:[study1.MetaCov.assoc.gz](https://genome.sph.umich.edu/w/images/d/da/Study1.MetaCov.assoc.gz) 
- Media:[study2.MetaCov.assoc.gz](https://genome.sph.umich.edu/w/images/a/a5/Study2.MetaCov.assoc.gz)


### How to Generate Summary Association Statistics and Prepare Them for Meta-analysis <a name=HTSASPTMA></a>
Meta-analysis summary association statistics can be generated by both RVTESTS and RAREMETALWORKER. Please refer to their documentations for generating summary association statistics.  
Once you have generated summary association statistics, you need to compress them with bgzip, and index them with tabix. If you use RAREMETALWORKER, the command should be like:   
If you use RVTESTS, your command should be  
```
bgzip study1.MetaScore.assoc
tabix -s 1 -b 2 -e 2 -S 1 study1.MetaScore.assoc.gz
tabix -s 1 -b 2 -e 2 -S 1 study1.MetaCov.assoc.gz
```
:point_right: **NOTE: Tabix 1.X does not seem to support the indexing for generic tab-delimited files. To index the file, please use tabix 0.2.5 or earlier versions**. 


### A Simple Tutorial for Using functions <a name=simple-tutorial-functions></a>
  - Using the rareMETALS.single function 
rareMETALS.single function allow you to perform meta-analyses for single variant association tests. The summary association statistics are combined using Mantel Haenszel test statistic. The details are described in our method paper: **Meta-analysis of gene-level tests of rare variant association, Nature Genetics, 46, 200–204 (2014)** [doi: 10.1038/ng.2852.](https://www.nature.com/articles/ng.2852)
Assume that you have a set of single variant score statistics and their covariance matrices.
Example:
```
cov.file <- c("study1.MetaCov.assoc.gz","study2.MetaCov.assoc.gz")
score.stat.file <- c("study1.MetaScore.assoc.gz","study2.MetaScore.assoc.gz")
```

```
library(rareMETALS)
res <- rareMETALS.single(score.stat.file,cov.file=NULL,range="19:11200093-11201275",alternative="two.sided",ix.gold=1,callrate.cutoff=0,hwe.cutoff=0)

###result can be explored as below###
 > names(res)
 [1] "p.value"            "ref"                "alt"                "integratedData"     "raw.data"          
 [6] "clean.data"         "statistic"          "direction.by.study" "anno"               "maf"               
 [11] "QC.by.study"        "no.sample"          "beta1.est"          "beta1.sd"           "hsq.est"           
 [16] "nearby"             "pos"    
 > print(res$pos)
 [1] "19:11200093" "19:11200213" "19:11200235" "19:11200272" "19:11200282" "19:11200309" "19:11200412" "19:11200419"
 [9] "19:11200431" "19:11200442" "19:11200475" "19:11200508" "19:11200514" "19:11200557" "19:11200579" "19:11200728"
 [17] "19:11200753" "19:11200754" "19:11200806" "19:11200839" "19:11200840" "19:11200896" "19:11201124" "19:11201259"
 [25] "19:11201274" "19:11201275"
 > print(res$p.value)
 [1] 0.551263675 0.056308558 0.172481571 0.734935815 0.922326732 0.053804524 0.886985353 0.903835162 0.005280228 0.266575301
 [11] 0.196122312 0.157114376 0.951477852 0.840523624 0.759482777 0.112743041 0.414147263 0.825877149 0.006090142 0.096474975
[21] 0.096474975 0.956407850 0.038234190 0.253512486 0.550935361 0.482315038```
```


  
  - [Using the rareMETALS.single.group function]  
  - [Using the rareMETALS.range function]  
  - [Using the rareMETALS.range.group function]  
  - [Using the conditional.rareMETALS.single]  
  - [Using the conditional.rareMETALS.range]  
  

### Feedback/Contact <a name="Feedback-Contact"></a>
Questions and requests can be sent to
Github issue page ([link](https://github.com/dajiangliu/rareMETALS/issues))
or
Dajiang Liu ([dajiang.liu@outlook.com](mailto:dajiang.liu@outlook.com "mailto:dajiang.liu@outlook.com")) and Fang Chen([fchen1@hmc.psu.edu](mailto:fchen1@hmc.psu.edu))
