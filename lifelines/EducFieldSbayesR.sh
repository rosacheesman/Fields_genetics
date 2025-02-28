#!/bin/bash
#SBATCH -t 2-00:00:00 ## WALL CLOCK TIME
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=40G
#SBATCH --output=../logs/EducFieldSbayesR.log 
#SBATCH --export=NONE
#SBATCH --get-user-env=60L
#SBATCH --begin=now+0hour


#cd ../SOFTWARE/gctb
#wget https://cnsgenomics.com/software/gctb/download/gctb_2.05beta_Linux.zip
#unzip gctb_2.05beta_Linux.zip
#unzip band_ukb_10k_hm3.zip
#wget https://cnsgenomics.com/data/GCTB/band_ukb_10k_hm3.zip #use banded LD matrix 

cd ../../CODE


gctbpath=../SOFTWARE/gctb

module load R/4.2.2-foss-2022a-bare  

Rscript ../Rscripts/gwas_to_ma.R ../../INPUT/GWASClean/Cheesman_meta_edu.txt.clean edu_field_EA SNP Allele1 Allele2 Freq1 Beta SE P-value NULL 460000
Rscript ../Rscripts/gwas_to_ma.R ../../INPUT/GWASClean/Cheesman_meta_arts.txt.clean arts_field_EA SNP Allele1 Allele2 Freq1 Beta SE P-value NULL 460000
Rscript ../Rscripts/gwas_to_ma.R ../../INPUT/GWASClean/Cheesman_meta_social.txt.clean social_field_EA SNP Allele1 Allele2 Freq1 Beta SE P-value NULL 460000
Rscript ../Rscripts/gwas_to_ma.R ../../INPUT/GWASClean/Cheesman_meta_business.txt.clean business_field_EA SNP Allele1 Allele2 Freq1 Beta SE P-value NULL 460000
Rscript ../Rscripts/gwas_to_ma.R ../../INPUT/GWASClean/Cheesman_meta_natural_sci.txt.clean natural_sci_field_EA SNP Allele1 Allele2 Freq1 Beta SE P-value NULL 460000
Rscript ../Rscripts/gwas_to_ma.R ../../INPUT/GWASClean/Cheesman_meta_ict.txt.clean ict_field_EA SNP Allele1 Allele2 Freq1 Beta SE P-value NULL 460000
Rscript ../Rscripts/gwas_to_ma.R ../../INPUT/GWASClean/Cheesman_meta_engineering.txt.clean engineering_field_EA SNP Allele1 Allele2 Freq1 Beta SE P-value NULL 460000
Rscript ../Rscripts/gwas_to_ma.R ../../INPUT/GWASClean/Cheesman_meta_agri.txt.clean agri_field_EA SNP Allele1 Allele2 Freq1 Beta SE P-value NULL 460000
Rscript ../Rscripts/gwas_to_ma.R ../../INPUT/GWASClean/Cheesman_meta_health.txt.clean health_field_EA SNP Allele1 Allele2 Freq1 Beta SE P-value NULL 460000
Rscript ../Rscripts/gwas_to_ma.R ../../INPUT/GWASClean/Cheesman_meta_services.txt.clean services_field_EA SNP Allele1 Allele2 Freq1 Beta SE P-value NULL 460000

for trait in {edu_field_EA,arts_field_EA,social_field_EA,natural_sci_field_EA,ict_field_EA,engineering_field_EA,agri_field_EA,health_field_EA,services_field_EA}
do
for c in `seq 1 22`
do
./${gctbpath}/gctb_2.05beta_Linux/gctb --sbayes R --ldm ${gctbpath}/band_ukb_10k_hm3/band_chr${c}.ldm.sparse --pi 0.95,0.02,0.02,0.01 --gamma 0.0,0.01,0.1,1 --gwas-summary ../../TEMP/${trait}.ma --chain-length 10000  --burn-in 2000 --out-freq 10 --out ../../TEMP/sbayeschr/${trait}_${c} 
done 
wait
done
wait 

for trait in {edu_field_EA,arts_field_EA,social_field_EA,natural_sci_field_EA,ict_field_EA,engineering_field_EA,agri_field_EA,health_field_EA,services_field_EA}
do
Rscript ../Rscripts/sbayesr_concatenate_weights.R ${trait} ${trait} &
done
wait 

module load PLINK/1.9-beta6-20190617 

#Create score, p-value threshold of 1 (lower not necessary because SbayesR has been used):
echo "1 0 1" >> range_list 

for trait in {edu_field_EA,arts_field_EA,social_field_EA,natural_sci_field_EA,ict_field_EA,engineering_field_EA,agri_field_EA,health_field_EA,services_field_EA}

do

filter duplicates
wc -l  ../../OUTPUT/${trait}sbayes.txt

awk -v c1=2 '{print $c1}'  ../../OUTPUT/${trait}sbayes.txt |\
sort |\
uniq -d |grep -vF .  > ../../TEMP/duplicate.snp 

echo "found"
wc -l ../../TEMP/duplicate.snp 
echo "duplicates"

if test -f "../../TEMP/duplicate.snp"; then
LC_ALL=C sort -u -f ../../OUTPUT/${trait}sbayes.txt | uniq > ../../TEMP/${trait}sbayes_nodup.txt
else 
fi

awk -v c1=2 -v c2=6 'FNR>1{print $c1,$c2}' ../../OUTPUT/${trait}sbayes.txt > ../../TEMP/SbayesScorePValues.${trait}.txt
awk -v c1=2 -v c2=3 -v c3=5 'FNR>1{print $c1,$c2,$c3}' ../../OUTPUT/${trait}sbayes.txt  > ../../TEMP/SbayesScoreWeights.${trait}.txt

wc -l  ../../TEMP/SbayesScoreWeights.${trait}.txt

awk -v c1=2 '{print $c1}'  ../../TEMP/SbayesScoreWeights.${trait}.txt |\
sort |\
uniq -d |grep -vF .  > ../../TEMP/duplicate.snp 

echo "filtering"
wc -l ../../TEMP/duplicate.snp 
echo "SNPs"

if test -f "../TEMP/duplicate.snp"; then
grep -vFf ../../TEMP/duplicate.snp ../../TEMP/SbayesScoreWeights.${trait}.txt > ../../TEMP/SbayesScoreWeights.${trait}.txt.nodup
rm ../../TEMP/SbayesScoreWeights.${trait}.txt
mv ../../TEMP/SbayesScoreWeights.${trait}.txt.nodup ../../TEMP/SbayesScoreWeights.${trait}.txt

else 
    echo "none"
fi


plink --bfile ../../INPUT/PLINKFILES/UGLICyto --score ../../TEMP/SbayesScoreWeights.${trait}.txt --q-score-range range_list ../../TEMP/SbayesScorePValues.${trait}.txt --out ../../INPUT/SCORES/Sbayes${trait}

done 
wait 


