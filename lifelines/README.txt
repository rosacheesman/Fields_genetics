README

Place 10 gzipped files for 10 eduFields GWAS in ../GWASSummary
Then run the following bash scripts in order:

CleanEducFields #a bash script that cleans each Education Field GWAS
EducFieldSbayesR.sh #a bash script which uses SbayesR to create PGIs for each of the ten fields. Here, the genotyped data is stored in plink binary format 
					#+ ../../INPUT/PLINKFILES/UGLICyto, which is cleaned according to the protocol outlined in the paper
EduFieldsMIToPGI.sh #a bash script which creates PGIs for each of the ten fields. Here, the genotyped data is the imputed genotyped data of one's parents, which 
					was created using SNIPar, and is stored chromosome-by-chromosome in plink binary format ../../TEMP/MI1, ../../TEMP/MI2, etc.
					#+ The corresponding pedigrees are stred in ../../TEMP/MIPedigree1.csv, ../../TEMP/MIPedigree2.csv, etc.


The scores are next linked to administrative data linkage into the secured microdata enclave of Statistics Netherlands. Within this enclave, the association analyses between each score and educational fields are conducted by running the STATA dofiles LifelinesEducationFields.do and LifelinesEducationFieldsSpouses.do
