Rscript roxygenise.R
R CMD build rareMETALS;
R CMD INSTALL rareMETALS_6.7.tar.gz
git commit -am "message"
git push 
