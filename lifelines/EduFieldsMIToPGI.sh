#!/bin/bash
#SBATCH -t 2-00:00:00 ## WALL CLOCK TIME
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --output=../logs/EduFieldsMIToPGI.log 
#SBATCH --export=NONE
#SBATCH --get-user-env=60L
#SBATCH --begin=now+0hour


module load R/4.2.2-foss-2022a-bare  

P=1
echo $P
for var in {edu_field_EA,arts_field_EA,social_field_EA,natural_sci_field_EA,ict_field_EA,engineering_field_EA,agri_field_EA,health_field_EA,services_field_EA}
	do 
echo $var
rm -f ../../INPUT/SCORES/MI/${var}_${P}_MI.csv
rm -fr ../../TEMP/ImputedParentalScoresChr*
Rscript ../Rscripts/MIToPGI.R ../../TEMP/SbayesScoreWeights.${var}.txt ../../TEMP/SbayesScorePValues.${var}.txt ../../TEMP/MI ../../TEMP/MIPedigree ${var} $P
wait 
	done 
	wait
	

for var in {edu_field_EA,arts_field_EA,social_field_EA,natural_sci_field_EA,ict_field_EA,engineering_field_EA,agri_field_EA,health_field_EA,services_field_EA}
do
Rscript ../Rscripts/MI_PGI_QC.R ../../TEMP/MIPedigree ${var} 1
wait
done 
