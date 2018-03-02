# rareMETALS

## **Table of Contents**

1. [Change Log](#change-log)
2. [Download and Installation](#download-and-installation)
3. [Documentation](#documatation)
4. [Forum]() 
5. [Supported Functionalities]()  
6. [Exemplar Datase]()
7. [How to Generate Summary Association Statistics and Prepare Them for Meta-analysis]  
8. [A Simple Tutorial for Using functions]  
  - [A Simple Tutorial for Using the rareMETALS.single function]  
  - [A Simple Tutorial for Using the rareMETALS.single.group function]  
  - [A Simple Tutorial for Using the rareMETALS.range function]  
  - [A Simple Tutorial for Using the rareMETALS.range.group function]  
  - [A Simple Tutorial for Using the conditional.rareMETALS.single]  
  - [A Simple Tutorial for Using the conditional.rareMETALS.range]  


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
