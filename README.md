# MetaGLIMPSE2
Meta Imputation of Low Coverage Sequencing


This method takes the results from two or more single panel GLIMPSE2 imputations and combines the output using weights estimated via HMM tailored for each individual and marker. 

The output of this method is a vcf file at the union set of markers in each input file with the estimated genotype dosage for each sample and marker.

We have shown this method is better than the best single panel imputation in the cases where imputation with the mega panel is better than the other single panel imputations (which is most cases). See pre-print for more information: 
