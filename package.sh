Rscript roxygenise.R
R CMD build rareMETALS;
R CMD INSTALL rareMETALS_7.1.tar.gz
git commit -am "message"
git push 
