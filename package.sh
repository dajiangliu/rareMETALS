Rscript roxygenise.R
R CMD build rareMETALS;
R CMD INSTALL rareMETALS_6.9.tar.gz
git commit -am "message"
git push 
