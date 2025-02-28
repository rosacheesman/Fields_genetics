use "EduFieldData.dta", clear

*merge to spousal ID:
drop _merge
merge 1:1 rinpersoon using "../GenesSES/INPUT/FamilyRelations.dta", nogenerate keep(match master)

drop rinpersoon
rename rinpersoon_spouse1 rinpersoon 

drop e_field* 
drop education_field 

merge m:1 rinpersoon using "G:/Onderwijs/HOOGSTEOPLTAB/2022/Geconverteerde data/HOOGSTEOPL2022TABV1", keepusing(RICHTSOI2021SCEDF2013HGNIRWO) keep(master match)


rename RICHTSOI2021SCEDF2013HGNIRWO education_field_sp

replace education_field_sp = substr(education_field_sp,1,2)

tab education_field_sp, generate(e_field)

label variable e_field2 "Education"
label variable e_field3 "Arts"
label variable e_field4 "Social"
label variable e_field5 "Business"
label variable e_field6 "Natural Sciences"
label variable e_field7 "ICT"
label variable e_field8 "Engineering"
label variable e_field9 "Agriculture"
label variable e_field10 "Health"
label variable e_field11 "Services"

drop if edu_field_EA_PGI == . 

global demcontrols PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10  male age age_sq age_cube c.male#c.age c.male#c.age_sq c.male#c.age_cube

local i = 2
foreach score in edu_field_EA_PGI arts_field_EA_PGI social_field_EA_PGI business_field_EA_PGI natural_sci_field_EA_PGI ict_field_EA_PGI engineering_field_EA_PGI agri_field_EA_PGI health_field_EA_PGI services_field_EA_PGI{
  logit e_field`i' $demcontrols
	local r2null = e(r2_p)
	eststo reg`i': logit e_field`i' `score' $demcontrols
	local r2add = e(r2_p)- `r2null' 
	
	if `i' ==2{
	outreg2 using  "${output}\EducFieldsSP", drop(e_field* $demcontrols) addnote("All specifications control for the first 20 PCs, gender, a cube in age, and interactions between gender and a cube in age") label excel replace noaster stats(coef se pval) paren(se) bracket(pval) addstat("Pseudo-R_squared", e(r2_p),"Incremental Pseudo-R_squared", `r2add')
	}
	else{
	 	outreg2 using  "${output}\EducFieldsSP", drop(e_field* $demcontrols) label excel append noaster stats(coef se pval) paren(se) bracket(pval) addstat("Pseudo-R_squared", e(r2_p),"Incremental Pseudo-R_squared", `r2add')	   
	}
	
	local i = `i'+1
}


local i = 2
foreach score in edu_field_EA_PGI arts_field_EA_PGI social_field_EA_PGI business_field_EA_PGI natural_sci_field_EA_PGI ict_field_EA_PGI engineering_field_EA_PGI agri_field_EA_PGI health_field_EA_PGI services_field_EA_PGI{
  
  if "`score'" == "natural_sci_field_EA_PGI"{
  local parscore natural_field_EA_PGI_Parental    
  }
  else{
  if "`score'" == "engineering_field_EA_PGI" {
  local parscore engineering_EA_Parental   
  }
  else{
      local parscore `score'_Parental
  }
  }
  
  logit e_field`i' $demcontrols `parscore'  
	local r2null = e(r2_p)
	eststo reg`i': logit e_field`i' `score' `parscore' $demcontrols 
	local r2add = e(r2_p)- `r2null' 
	
	if `i' ==2{
	outreg2 using  "${output}\EducFieldsMISP", drop(e_field* $demcontrols) addnote("All specifications control for the first 20 PCs, gender, a cube in age, and interactions between gender and a cube in age") label excel replace noaster stats(coef se pval) paren(se) bracket(pval)  addstat("Incremental R_squared", `r2add')	
	}
	else{
	 	outreg2 using  "${output}\EducFieldsMISP", drop(e_field* $demcontrols) label excel append noaster stats(coef se pval) paren(se) bracket(pval) addstat("Pseudo-R_squared", e(r2_p),"Incremental Pseudo-R_squared", `r2add')	  
	}
	
	local i = `i'+1
}
