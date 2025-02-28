#!/bin/bash
#SBATCH -t 24:00:00 ## WALL CLOCK TIME
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --output=../logs/CleanEducFields.log 
#SBATCH --export=NONE
#SBATCH --get-user-env=60L
#SBATCH --begin=now+0hour


module load R/4.2.2-foss-2022a-bare  

bash GWASSumCleanProtocolEducFields.sh PC1.N_weighted_GWAMA.results.txt 2 5 6 12 15 eduFactor1 no
bash GWASSumCleanProtocolEducFields.sh PC2.N_weighted_GWAMA.results.txt 2 5 6 12 15 eduFactor2 no
bash GWASSumCleanProtocolEducFields.sh Cheesman_meta_edu.txt 12 2 3 9 10 edu_field_EA yes
bash GWASSumCleanProtocolEducFields.sh Cheesman_meta_arts.txt 12 2 3 9 10 arts_field_EA yes
bash GWASSumCleanProtocolEducFields.sh Cheesman_meta_social.txt 12 2 3 9 10 social_field_EA yes
bash GWASSumCleanProtocolEducFields.sh Cheesman_meta_business.txt 12 2 3 9 10 business_field_EA yes
bash GWASSumCleanProtocolEducFields.sh Cheesman_meta_natural_sci.txt 12 2 3 9 10 natural_sci_field_EA yes
bash GWASSumCleanProtocolEducFields.sh Cheesman_meta_ict.txt 12 2 3 9 10 ict_field_EA yes
bash GWASSumCleanProtocolEducFields.sh Cheesman_meta_engineering.txt 12 2 3 9 10 engineering_field_EA yes
bash GWASSumCleanProtocolEducFields.sh Cheesman_meta_agri.txt 12 2 3 9 10 agri_field_EA yes
bash GWASSumCleanProtocolEducFields.sh Cheesman_meta_health.txt 12 2 3 9 10 health_field_EA yes
bash GWASSumCleanProtocolEducFields.sh Cheesman_meta_services.txt 12 2 3 9 10 services_field_EA yes


