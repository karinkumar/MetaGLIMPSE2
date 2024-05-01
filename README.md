# MetaGLIMPSE2
Meta Imputation of Low Coverage Sequencing

***Overview***

This method takes the results from two or more single panel GLIMPSE2 imputations and combines the output using weights estimated via HMM tailored for each individual and marker. 

The output of this method is a vcf file at the union set of markers in each input file with the estimated genotype dosage for each sample and marker.

We have shown this method is better than the best single panel imputation in the cases where imputation with the mega panel is better than the other single panel imputations (which is most cases). See pre-print for more information: TBA

***1. Installation***

git clone https://github.com/karinkumar/MetaGLIMPSE2.git

Once installed, the main file is RunMetaGLIMPSE.py. The following options are required:


-- dosages:  paths of imputed genotypes vcf files (output from GLIMPSE2), 

-- gl:  vcf file with genotype likelihoods for each position in the union set of markers of the dosage files, 

-- out:  prefix of outfiles. 

***2. Run Example***

See the example files and run the following code once you have installed the program.

To run GLIMPSE2 please check out the GLIMPSE tutorial https://odelaneau.github.io/GLIMPSE/ 

***3 Ligate*** 

MetaGLIMPSE2 produces meta-imputed chunks. In order to be turned into one vcf file for an entire chromosome, they need to be ligated used bcftools. The following code ligates the chunks in the example. 

ls -v out*.vcf.gz > list.txt

bcftools concat -f list.txt -Oz -o fullchrom.vcf.gz

bcftools index fullchrom.vcf.gz



