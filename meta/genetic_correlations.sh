
cd /gpfs/work5/0/gusr0487/abdel/rosa/workingdir/revision/

# load LDSC regression
module load 2021
module load Anaconda3/2021.05
conda init
source ~/.bashrc
conda activate /gpfs/work5/0/gusr0487/anaconda2/envs/ldsc

zcat PC1.N_weighted_GWAMA.results.txt.gz | sed '1d' | \
awk 'BEGIN{print "SNPID CHR BP A1 A2 EAF N N_obs Z PVAL"} \
{print $2, $3, $4, $5, $6, $7, $9, $10, $14, $15}' > PC1.rosa_edu.sumstats
head PC1.rosa_edu.sumstats

zcat PC2.N_weighted_GWAMA.results.txt.gz | sed '1d' | \
awk 'NR==1{print "SNPID CHR BP A1 A2 EAF N N_obs Z PVAL"} \
{print $2, $3, $4, $5, $6, $7, $9, $10, $14, $15}' > PC2.rosa_edu.sumstats
head PC2.rosa_edu.sumstats

# munge
/projects/0/gusr0487/UKBAUMC/software/ldsc/munge_sumstats.py \
--sumstats PC1.rosa_edu.sumstats \
--chunksize 1000000 \
--N 10413 \
--out PC1.rosa_edu.LDSC \
--merge-alleles /projects/0/gusr0487/UKBAUMC/software/ldsc/w_hm3.snplist

/projects/0/gusr0487/UKBAUMC/software/ldsc/munge_sumstats.py \
--sumstats PC2.rosa_edu.sumstats \
--chunksize 1000000 \
--ignore beta \
--N 7353 \
--out PC2.rosa_edu.LDSC \
--merge-alleles /projects/0/gusr0487/UKBAUMC/software/ldsc/w_hm3.snplist

#########################################################

# rgs

cd /gpfs/work5/0/gusr0487/lisa_ukbaumc_backup/munged_sumstats_rgs
# cp /gpfs/work5/0/gusr0487/lisa_abdel_backup/OCD/workingdir/munged_sumstats/*.* .
cp /gpfs/work5/0/gusr0487/abdel/rosa/workingdir/revision/PC1.rosa_edu.LDSC.sumstats.gz .
cp /gpfs/work5/0/gusr0487/abdel/rosa/workingdir/revision/PC2.rosa_edu.LDSC.sumstats.gz .

/projects/0/gusr0487/UKBAUMC/software/ldsc/ldsc.py \
--rg PC1.rosa_edu.LDSC.sumstats.gz,PC2.rosa_edu.LDSC.sumstats.gz,ukb.GWAS.memory.399.HRC.25PCs.LDSC.sumstats.gz,raymond.F5.LDSC.sumstats.gz,raymond.F10.LDSC.sumstats.gz,raymond.F15.LDSC.sumstats.gz,F1.LDSC.sumstats.gz,F2.LDSC.sumstats.gz,occupational_status_siops.LDSC.sumstats.gz,occupational_status_isei.LDSC.sumstats.gz,occupational_status_camsis.LDSC.sumstats.gz,creativity_kim_et_al_2024.LDSC.sumstats.gz,OPENN.rosa.LDSC.sumstats.gz,CONSC.rosa.LDSC.sumstats.gz,EXTRA.rosa.LDSC.sumstats.gz,AGREE.rosa.LDSC.sumstats.gz,NEURO.rosa.LDSC.sumstats.gz,bip.2021.LDSC.sumstats.gz,SCZ3.EUR.LDSC.sumstats.gz,cocdep.sumstats.gz,intelligence.savage.2018.LDSC.sumstats.gz,ENIGMA.surface_area.LDSC.sumstats.gz,ENIGMA.thickness.LDSC.sumstats.gz,ENIGMA2.Accumbens.LDSC.sumstats.gz,ENIGMA2.Amygdala.LDSC.sumstats.gz,ENIGMA2.Caudate.LDSC.sumstats.gz,ENIGMA2.Hippocampus.LDSC.sumstats.gz,ENIGMA2.Pallidum.LDSC.sumstats.gz,ENIGMA2.Thalamus.LDSC.sumstats.gz,COVID19_HGI_A2_R6.EUR.LDSC.sumstats.gz,COVID19_HGI_B1_R6.EUR.LDSC.sumstats.gz,COVID19_HGI_B2_R6.EUR.LDSC.sumstats.gz,COVID19_HGI_C2_R6.EUR.LDSC.sumstats.gz,vit_D.BMI_cov.revezetal_2020.LDSC.sumstats.gz,vit_D.revezetal_2020.LDSC.sumstats.gz,EA.non_cog.LDSC.sumstats.gz,EA.cog.LDSC.sumstats.gz,schizophreniaPGC3.sumstats.gz,UKB.Memory.LDSC.sumstats.gz,SmokingCessation.NG_2019.LDSC.sumstats.gz,UKB.Income.LDSC.sumstats.gz,number_sexual_partners.LDSC.sumstats.gz,DrinksPerWeek.NG_2019.LDSC.sumstats.gz,Alzheimer.2019.LDSC.sumstats.gz,risk.linner.2019.LDSC.sumstats.gz,NumberChildrenEverBorn_Pooled.LDSC.sumstats.gz,PTSD.2019.LDSC.sumstats.gz,BIP.2018.LDSC.sumstats.gz,ukb.GWAS.bmi.25PCs.LDSC.sumstats.gz,ukb.GWAS.bodyfat.25PCs.LDSC.sumstats.gz,bmi_combined.LDSC.sumstats.gz,UKB.tiredness.LDSC.sumstats.gz,parkinsons.2019.LDSC.sumstats.gz,cad.add.160614.LDSC.sumstats.gz,MIS.sumstats.gz,ASD.LDSC.sumstats.gz,OCD.14042020.LDSC.sumstats.gz,inflammatory_bowel.LDSC.sumstats.gz,crohns_disease.LDSC.sumstats.gz,ulcerative_colitis.LDSC.sumstats.gz,suicidality.2019.LDSC.sumstats.gz,ADHD.EUR.LDSC.sumstats.gz,childhood_maltreatment.2020.LDSC.sumstats.gz,atopic_dermatitis.2015.LDSC.sumstats.gz,ALS.2020.LDSC.sumstats.gz,AN2.2019.LDSC.sumstats.gz,AgeFirstBirth_Pooled.LDSC.sumstats.gz,AgeOfInitiation.NG_2019.LDSC.sumstats.gz,CD.LDSC.sumstats.gz,CigarettesPerDay.NG_2019.LDSC.sumstats.gz,EA3_excl_23andMe.LDSC.sumstats.gz,HDL.LDSC.sumstats.gz,IQ.LDSC.sumstats.gz,LDL.LDSC.sumstats.gz,MDD.2018.LDSC.sumstats.gz,MDD.2019.LDSC.sumstats.gz,ND.LDSC.sumstats.gz,RA.2013.LDSC.sumstats.gz,SCZ.pardinas.2018.LDSC.sumstats.gz,SWB_Full.LDSC.sumstats.gz,SmokingInitiation.NG_2019.LDSC.sumstats.gz,TS.2018.LDSC.sumstats.gz,UKB.Reaction_time.LDSC.sumstats.gz,UKB.Townsend.LDSC.sumstats.gz,UKB.VNR.LDSC.sumstats.gz,UKB.self_rated_health.LDSC.sumstats.gz,age_at_menarche.LDSC.sumstats.gz,age_at_menopauze.LDSC.sumstats.gz,agreeableness.GPC.23andme.LDSC.sumstats.gz,alcdep.eur_discovery.aug20.LDSC.sumstats.gz,alcohol_clarke.LDSC.sumstats.gz,anxiety.UKB_iPSYCH.LDSC.sumstats.gz,asthma.adult_onset.LDSC.sumstats.gz,asthma.child_onset.LDSC.sumstats.gz,birth_weight.LDSC.sumstats.gz,body_fat.LDSC.sumstats.gz,cad.add.160614.LDSC.sumstats.gz,caffeine.LDSC.sumstats.gz,cannabis_ever_2018.revision.LDSC.sumstats.gz,childhoodIQ.CHIC_Benyamin_2014.LDSC.sumstats.gz,conscientiousness.GPC.23andme.LDSC.sumstats.gz,diagram_T2D.LDSC.sumstats.gz,extraversion.GPC.23andme.LDSC.sumstats.gz,family_satisfaction.LDSC.sumstats.gz,father_death.LDSC.sumstats.gz,focal_epilepsy.LDSC.sumstats.gz,freq_friend_visit.LDSC.sumstats.gz,friend_satisfaction.LDSC.sumstats.gz,generalised_epilepsy.LDSC.sumstats.gz,harm_avoidance.2012.LDSC.sumstats.gz,height_combined.LDSC.sumstats.gz,hip_combined.LDSC.sumstats.gz,intelligence.savage.2018.LDSC.sumstats.gz,job_satisfaction.rapid_UKB.LDSC.sumstats.gz,loneliness_HMG.LDSC.sumstats.gz,lupus.2018.LDSC.sumstats.gz,mi.add.030315.LDSC.sumstats.gz,morningness.2019.LDSC.sumstats.gz,mother_death.LDSC.sumstats.gz,neuroticism.nagel.2018.LDSC.sumstats.gz,neuroticism.GPC.23andme.LDSC.sumstats.gz,openness.GPC.23andme.LDSC.sumstats.gz,parents_death.LDSC.sumstats.gz,sleep_duration.janssen.2019.LDSC.sumstats.gz,total_cholesterol.LDSC.sumstats.gz,triglycerides.LDSC.sumstats.gz,ukb.GWAS.age_at_first_sex.25PCs.LDSC.sumstats.gz,ukb.GWAS.glasses.25PCs.LDSC.sumstats.gz,ukb.GWAS.meaningful_life.25PCs.LDSC.sumstats.gz,morningness.2019.LDSC.sumstats.gz,insomnia.2019.LDSC.sumstats.gz,job_satisfaction.rapid_UKB.LDSC.sumstats.gz,wc_combined.LDSC.sumstats.gz,whr_combined.LDSC.sumstats.gz \
--ref-ld-chr /projects/0/gusr0487/UKBAUMC/software/ldsc/eur_w_ld_chr/ \
--w-ld-chr /projects/0/gusr0487/UKBAUMC/software/ldsc/eur_w_ld_chr/ \
--out /gpfs/work5/0/gusr0487/abdel/rosa/workingdir/revision/PC1.24102024.LDSC.rgs &

/projects/0/gusr0487/UKBAUMC/software/ldsc/ldsc.py \
--rg PC2.rosa_edu.LDSC.sumstats.gz,PC1.rosa_edu.LDSC.sumstats.gz,ukb.GWAS.memory.399.HRC.25PCs.LDSC.sumstats.gz,raymond.F5.LDSC.sumstats.gz,raymond.F10.LDSC.sumstats.gz,raymond.F15.LDSC.sumstats.gz,F1.LDSC.sumstats.gz,F2.LDSC.sumstats.gz,occupational_status_siops.LDSC.sumstats.gz,occupational_status_isei.LDSC.sumstats.gz,occupational_status_camsis.LDSC.sumstats.gz,creativity_kim_et_al_2024.LDSC.sumstats.gz,OPENN.rosa.LDSC.sumstats.gz,CONSC.rosa.LDSC.sumstats.gz,EXTRA.rosa.LDSC.sumstats.gz,AGREE.rosa.LDSC.sumstats.gz,NEURO.rosa.LDSC.sumstats.gz,bip.2021.LDSC.sumstats.gz,SCZ3.EUR.LDSC.sumstats.gz,cocdep.sumstats.gz,intelligence.savage.2018.LDSC.sumstats.gz,ENIGMA.surface_area.LDSC.sumstats.gz,ENIGMA.thickness.LDSC.sumstats.gz,ENIGMA2.Accumbens.LDSC.sumstats.gz,ENIGMA2.Amygdala.LDSC.sumstats.gz,ENIGMA2.Caudate.LDSC.sumstats.gz,ENIGMA2.Hippocampus.LDSC.sumstats.gz,ENIGMA2.Pallidum.LDSC.sumstats.gz,ENIGMA2.Thalamus.LDSC.sumstats.gz,COVID19_HGI_A2_R6.EUR.LDSC.sumstats.gz,COVID19_HGI_B1_R6.EUR.LDSC.sumstats.gz,COVID19_HGI_B2_R6.EUR.LDSC.sumstats.gz,COVID19_HGI_C2_R6.EUR.LDSC.sumstats.gz,vit_D.BMI_cov.revezetal_2020.LDSC.sumstats.gz,vit_D.revezetal_2020.LDSC.sumstats.gz,EA.non_cog.LDSC.sumstats.gz,EA.cog.LDSC.sumstats.gz,schizophreniaPGC3.sumstats.gz,UKB.Memory.LDSC.sumstats.gz,SmokingCessation.NG_2019.LDSC.sumstats.gz,UKB.Income.LDSC.sumstats.gz,number_sexual_partners.LDSC.sumstats.gz,DrinksPerWeek.NG_2019.LDSC.sumstats.gz,Alzheimer.2019.LDSC.sumstats.gz,risk.linner.2019.LDSC.sumstats.gz,NumberChildrenEverBorn_Pooled.LDSC.sumstats.gz,PTSD.2019.LDSC.sumstats.gz,BIP.2018.LDSC.sumstats.gz,ukb.GWAS.bmi.25PCs.LDSC.sumstats.gz,ukb.GWAS.bodyfat.25PCs.LDSC.sumstats.gz,bmi_combined.LDSC.sumstats.gz,UKB.tiredness.LDSC.sumstats.gz,parkinsons.2019.LDSC.sumstats.gz,cad.add.160614.LDSC.sumstats.gz,MIS.sumstats.gz,ASD.LDSC.sumstats.gz,OCD.14042020.LDSC.sumstats.gz,inflammatory_bowel.LDSC.sumstats.gz,crohns_disease.LDSC.sumstats.gz,ulcerative_colitis.LDSC.sumstats.gz,suicidality.2019.LDSC.sumstats.gz,ADHD.EUR.LDSC.sumstats.gz,childhood_maltreatment.2020.LDSC.sumstats.gz,atopic_dermatitis.2015.LDSC.sumstats.gz,ALS.2020.LDSC.sumstats.gz,AN2.2019.LDSC.sumstats.gz,AgeFirstBirth_Pooled.LDSC.sumstats.gz,AgeOfInitiation.NG_2019.LDSC.sumstats.gz,CD.LDSC.sumstats.gz,CigarettesPerDay.NG_2019.LDSC.sumstats.gz,EA3_excl_23andMe.LDSC.sumstats.gz,HDL.LDSC.sumstats.gz,IQ.LDSC.sumstats.gz,LDL.LDSC.sumstats.gz,MDD.2018.LDSC.sumstats.gz,MDD.2019.LDSC.sumstats.gz,ND.LDSC.sumstats.gz,RA.2013.LDSC.sumstats.gz,SCZ.pardinas.2018.LDSC.sumstats.gz,SWB_Full.LDSC.sumstats.gz,SmokingInitiation.NG_2019.LDSC.sumstats.gz,TS.2018.LDSC.sumstats.gz,UKB.Reaction_time.LDSC.sumstats.gz,UKB.Townsend.LDSC.sumstats.gz,UKB.VNR.LDSC.sumstats.gz,UKB.self_rated_health.LDSC.sumstats.gz,age_at_menarche.LDSC.sumstats.gz,age_at_menopauze.LDSC.sumstats.gz,agreeableness.GPC.23andme.LDSC.sumstats.gz,alcdep.eur_discovery.aug20.LDSC.sumstats.gz,alcohol_clarke.LDSC.sumstats.gz,anxiety.UKB_iPSYCH.LDSC.sumstats.gz,asthma.adult_onset.LDSC.sumstats.gz,asthma.child_onset.LDSC.sumstats.gz,birth_weight.LDSC.sumstats.gz,body_fat.LDSC.sumstats.gz,cad.add.160614.LDSC.sumstats.gz,caffeine.LDSC.sumstats.gz,cannabis_ever_2018.revision.LDSC.sumstats.gz,childhoodIQ.CHIC_Benyamin_2014.LDSC.sumstats.gz,conscientiousness.GPC.23andme.LDSC.sumstats.gz,diagram_T2D.LDSC.sumstats.gz,extraversion.GPC.23andme.LDSC.sumstats.gz,family_satisfaction.LDSC.sumstats.gz,father_death.LDSC.sumstats.gz,focal_epilepsy.LDSC.sumstats.gz,freq_friend_visit.LDSC.sumstats.gz,friend_satisfaction.LDSC.sumstats.gz,generalised_epilepsy.LDSC.sumstats.gz,harm_avoidance.2012.LDSC.sumstats.gz,height_combined.LDSC.sumstats.gz,hip_combined.LDSC.sumstats.gz,intelligence.savage.2018.LDSC.sumstats.gz,job_satisfaction.rapid_UKB.LDSC.sumstats.gz,loneliness_HMG.LDSC.sumstats.gz,lupus.2018.LDSC.sumstats.gz,mi.add.030315.LDSC.sumstats.gz,morningness.2019.LDSC.sumstats.gz,mother_death.LDSC.sumstats.gz,neuroticism.nagel.2018.LDSC.sumstats.gz,neuroticism.GPC.23andme.LDSC.sumstats.gz,openness.GPC.23andme.LDSC.sumstats.gz,parents_death.LDSC.sumstats.gz,sleep_duration.janssen.2019.LDSC.sumstats.gz,total_cholesterol.LDSC.sumstats.gz,triglycerides.LDSC.sumstats.gz,ukb.GWAS.age_at_first_sex.25PCs.LDSC.sumstats.gz,ukb.GWAS.glasses.25PCs.LDSC.sumstats.gz,ukb.GWAS.meaningful_life.25PCs.LDSC.sumstats.gz,morningness.2019.LDSC.sumstats.gz,insomnia.2019.LDSC.sumstats.gz,job_satisfaction.rapid_UKB.LDSC.sumstats.gz,wc_combined.LDSC.sumstats.gz,whr_combined.LDSC.sumstats.gz \
--ref-ld-chr /projects/0/gusr0487/UKBAUMC/software/ldsc/eur_w_ld_chr/ \
--w-ld-chr /projects/0/gusr0487/UKBAUMC/software/ldsc/eur_w_ld_chr/ \
--out /gpfs/work5/0/gusr0487/abdel/rosa/workingdir/revision/PC2.24102024.LDSC.rgs &


