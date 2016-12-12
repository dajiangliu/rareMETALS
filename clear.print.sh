##sed -ie s/print/##print/g rareMETALS/R/*.R;
sed -ie '/[^\t]\+print/d' rareMETALS/R/*.R

./package.sh
