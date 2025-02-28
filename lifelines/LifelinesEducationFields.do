*use "H:\CBSLifelines\INPUT\Sample0_sibs", clear 
cd "H:/EducationFields"
global output "H:/EducationFields"
use "H:\GenesSES\INPUT\LifelinesMergedWealth.dta", clear 

*some variable names were cut off because they were too lengthy:
rename _v49 natural_field_EA_PGI_Parental

rename _v50 engineering_EA_Parental
*destring edu_field_EA_PGI arts_field_EA_PGI social_field_EA_PGI business_field_EA_PGI natural_sci_field_EA_PGI ict_field_EA_PGI engineering_field_EA_PGI agri_field_EA_PGI health_field_EA_PGI services_field_EA_PGI, replace

sum age 

keep if age >= 25

merge 1:1 rinpersoon using "G:/Onderwijs/HOOGSTEOPLTAB/2022/Geconverteerde data/HOOGSTEOPL2022TABV1", keepusing(RICHTSOI2021SCEDF2013HGNIRWO) keep(master match)

rename RICHTSOI2021SCEDF2013HGNIRWO education_field 

replace education_field = substr(education_field,1,2)

tab education_field, generate(e_field)

forvalues i = 1/12 {
    logit e_field`i' EA4_PGI
	logit e_field`i' EA4_PGI EA4_PGI_Parental
}

global destringvar edu_field_EA_PGI arts_field_EA_PGI arts_field_EA_PGI_Parental social_field_EA_PGI business_field_EA_PGI natural_sci_field_EA_PGI agri_field_EA_PGI engineering_field_EA_PGI health_field_EA_PGI ict_field_EA_PGI services_field_EA_PGI ///
 edu_field_EA_PGI_Parental social_field_EA_PGI_Parental business_field_EA_PGI_Parental natural_field_EA_PGI_Parental agri_field_EA_PGI_Parental engineering_EA_Parental health_field_EA_PGI_Parental ict_field_EA_PGI_Parental services_field_EA_PGI_Parental 
foreach v of varlist $destringvar { 
replace `v' = "" if `v'=="NA" 
display "`v'"
destring `v', force replace 
}

gen age_sq = age^2
gen age_cube = age^3	

global demcontrols PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10  male age age_sq age_cube c.male#c.age c.male#c.age_sq c.male#c.age_cube
global famcontrols male age age_sq age_cube c.male#c.age c.male#c.age_sq c.male#c.age_cube

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
*eststo reg1: reg e_field1 gen_prog_field_EA_PGI $demcontrols

save "EduFieldData.dta", replace 

drop if edu_field_EA_PGI == . 
destring eduFactor1_PGI eduFactor2_PGI eduFactor1_PGI_Parental eduFactor2_PGI_Parental, replace force 
forval i = 2/11 {
	logit e_field`i' $demcontrols
	local r2null = e(r2_p)
	eststo reg`i': logit e_field`i' eduFactor1_PGI eduFactor2_PGI $demcontrols
	local r2add = e(r2_p)- `r2null'
	if `i' ==2{
	outreg2 using  "${output}\EducFieldsFactor", drop(e_field* $demcontrols) addnote("All specifications control for the first 20 PCs, gender, a cube in age, and interactions between gender and a cube in age") label excel replace noaster stats(coef se pval) paren(se) bracket(pval) addstat("Pseudo-R_squared", e(r2_p),"Incremental Pseudo-R_squared", `r2add')
	}
	else{
	 	outreg2 using  "${output}\EducFieldsFactor", drop(e_field* $demcontrols) label excel append noaster stats(coef se pval) paren(se) bracket(pval) addstat("Pseudo-R_squared", e(r2_p),"Incremental Pseudo-R_squared", `r2add')	   
	}
	
	
	eststo reg`i': logit e_field`i' eduFactor1_PGI eduFactor2_PGI eduFactor1_PGI_Parental eduFactor2_PGI_Parental $demcontrols
	local r2add = e(r2_p)- `r2null'
	if `i' ==2{
	outreg2 using  "${output}\EducFieldsFactor_par", drop(e_field* $demcontrols) addnote("All specifications control for the first 20 PCs, gender, a cube in age, and interactions between gender and a cube in age") label excel replace noaster stats(coef se pval) paren(se) bracket(pval) addstat("Pseudo-R_squared", e(r2_p),"Incremental Pseudo-R_squared", `r2add')
	}
	else{
	 	outreg2 using  "${output}\EducFieldsFactor_par", drop(e_field* $demcontrols) label excel append noaster stats(coef se pval) paren(se) bracket(pval) addstat("Pseudo-R_squared", e(r2_p),"Incremental Pseudo-R_squared", `r2add')	   
	}
	
	local i = `i'+1
}



local i = 2
foreach score in edu_field_EA_PGI arts_field_EA_PGI social_field_EA_PGI business_field_EA_PGI natural_sci_field_EA_PGI ict_field_EA_PGI engineering_field_EA_PGI agri_field_EA_PGI health_field_EA_PGI services_field_EA_PGI{
  logit e_field`i' $demcontrols
	local r2null = e(r2_p)
	eststo reg`i': logit e_field`i' `score' $demcontrols
	local r2add = e(r2_p)- `r2null' 
	
	if `i' ==2{
	outreg2 using  "${output}\EducFields", drop(e_field* $demcontrols) addnote("All specifications control for the first 20 PCs, gender, a cube in age, and interactions between gender and a cube in age") label excel replace noaster stats(coef se pval) paren(se) bracket(pval) addstat("Pseudo-R_squared", e(r2_p),"Incremental Pseudo-R_squared", `r2add')
	}
	else{
	 	outreg2 using  "${output}\EducFields", drop(e_field* $demcontrols) label excel append noaster stats(coef se pval) paren(se) bracket(pval) addstat("Pseudo-R_squared", e(r2_p),"Incremental Pseudo-R_squared", `r2add')	   
	}
	
	local i = `i'+1
}


local i = 2
foreach score in edu_field_EA_PGI arts_field_EA_PGI social_field_EA_PGI business_field_EA_PGI natural_sci_field_EA_PGI ict_field_EA_PGI engineering_field_EA_PGI agri_field_EA_PGI health_field_EA_PGI services_field_EA_PGI{
  logit e_field`i' $demcontrols
	local r2null = e(r2_p)
	eststo reg`i': logit e_field`i' `score' c.`score'#c.male $demcontrols
	local r2add = e(r2_p)- `r2null' 
	
	if `i' ==2{
	outreg2 using  "${output}\EducFieldsBySex", drop(e_field* $demcontrols) addnote("All specifications control for the first 20 PCs, gender, a cube in age, and interactions between gender and a cube in age") label excel replace noaster stats(coef se pval) paren(se) bracket(pval)  addstat("Incremental Pseudo-R_squared", `r2add')	
	}
	else{
	 	outreg2 using  "${output}\EducFieldsBySex", drop(e_field* $demcontrols) label excel append noaster stats(coef se pval) paren(se) bracket(pval) addstat("Pseudo-R_squared", e(r2_p),"Incremental Pseudo-R_squared", `r2add')		   
	}
	
	local i = `i'+1
}


local i = 2
foreach score in edu_field_EA_PGI arts_field_EA_PGI social_field_EA_PGI business_field_EA_PGI natural_sci_field_EA_PGI ict_field_EA_PGI engineering_field_EA_PGI agri_field_EA_PGI health_field_EA_PGI services_field_EA_PGI{
  logit e_field`i' $demcontrols YearsEducation
	local r2null = e(r2_p)
	eststo reg`i': logit e_field`i' `score' $demcontrols YearsEducation
	local r2add = e(r2_p)- `r2null' 
	
	if `i' ==2{
	outreg2 using  "${output}\EducFieldsYearsEducation", drop(e_field* $demcontrols) addnote("All specifications control for the first 20 PCs, gender, a cube in age, and interactions between gender and a cube in age") label excel replace noaster stats(coef se pval) paren(se) bracket(pval)  addstat("Incremental Pseudo-R_squared", `r2add')	
	}
	else{
	 	outreg2 using  "${output}\EducFieldsYearsEducation", drop(e_field* $demcontrols) label excel append noaster stats(coef se pval) paren(se) bracket(pval) addstat("Pseudo-R_squared", e(r2_p),"Incremental Pseudo-R_squared", `r2add')		   
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
	outreg2 using  "${output}\EducFieldsMI", drop(e_field* $demcontrols) addnote("All specifications control for the first 20 PCs, gender, a cube in age, and interactions between gender and a cube in age") label excel replace noaster stats(coef se pval) paren(se) bracket(pval)  addstat("Incremental Pseudo-R_squared", `r2add')	
	}
	else{
	 	outreg2 using  "${output}\EducFieldsMI", drop(e_field* $demcontrols) label excel append noaster stats(coef se pval) paren(se) bracket(pval) addstat("Pseudo-R_squared", e(r2_p),"Incremental Pseudo-R_squared", `r2add')	  
	}
	
	logit e_field`i' $demcontrols `parscore'  
	local r2null = e(r2_p)
	eststo reg`i': logit e_field`i' `score' `parscore' c.`score'#c.male c.`parscore'#c.male $demcontrols 
	local r2add = e(r2_p)- `r2null' 
	if `i' ==2{
	outreg2 using  "${output}\EducFieldsMIBySex", drop(e_field* $demcontrols) addnote("All specifications control for the first 20 PCs, gender, a cube in age, and interactions between gender and a cube in age") label excel replace noaster stats(coef se pval) paren(se) bracket(pval)  addstat("Incremental R_squared", `r2add')	
	}
	else{
	 	outreg2 using  "${output}\EducFieldsMIBySex", drop(e_field* $demcontrols) label excel append noaster stats(coef se pval) paren(se) bracket(pval) addstat("Pseudo-R_squared", e(r2_p),"Incremental Pseudo-R_squared", `r2add')		   
	}	
	
	local i = `i'+1
}


*bootstrap the difference 
capture program drop diff_coef_logit
program define diff_coef_logit, rclass
	quietly logit e_field2 edu_field_EA_PGI PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 male age age_sq age_cube c.male#c.age c.male#c.age_sq c.male#c.age_cube
	scalar b1 = _b[edu_field_EA_PGI]
	quietly logit e_field2 edu_field_EA_PGI edu_field_EA_PGI_Parental male age age_sq age_cube c.male#c.age c.male#c.age_sq c.male#c.age_cube
	scalar b2 = _b[edu_field_EA_PGI]
	
	scalar differ = b1 - b2
	
	return scalar diff = differ
end

bootstrap r(diff), reps(10): diff_coef_logit 

capture program drop diff_coef_logit
program define diff_coef_logit, rclass
	syntax, yvar(name) xvar(name) controls(name)
	quietly logit `yvar' `xvar' PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 male age age_sq age_cube c.male#c.age c.male#c.age_sq c.male#c.age_cube
	scalar b1 = _b[`xvar']
	quietly logit `yvar' `xvar' `controls' male age age_sq age_cube c.male#c.age c.male#c.age_sq c.male#c.age_cube
	scalar b2 = _b[`xvar']
	
	scalar differ = b1 - b2
	
	return scalar diff = differ
end


bootstrap r(diff), reps(10): diff_coef_logit, yvar(e_field2) xvar(edu_field_EA_PGI) controls(edu_field_EA_PGI_Parental) 

mat results = J(9, 1, .)  

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
bootstrap r(diff), reps(1000): diff_coef_logit, yvar(e_field`i') xvar(`score') controls(`parscore') 

mat temp = r(table)
mat results = results, temp

mat rownames results = rnames(r(diff))

local i = `i' + 1
  
}

putexcel set "bootstrap_results.xlsx", replace
putexcel A1=matrix(results)

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
  
  
  
  logit e_field`i' $demcontrols `parscore' YearsEducation 
	local r2null = e(r2_p)
	eststo reg`i': logit e_field`i' `score' `parscore' $demcontrols YearsEducation
	local r2add = e(r2_p)- `r2null' 
	
	if `i' ==2{
	outreg2 using  "${output}\EducFieldsMIYearsEducation", drop(e_field* $demcontrols) addnote("All specifications control for the first 20 PCs, gender, a cube in age, and interactions between gender and a cube in age") label excel replace noaster stats(coef se pval) paren(se) bracket(pval) addstat("Pseudo-R_squared", e(r2_p),"Incremental Pseudo-R_squared", `r2add')
	}
	else{
	 	outreg2 using  "${output}\EducFieldsMIYearsEducation", drop(e_field* $demcontrols) label excel append noaster stats(coef se pval) paren(se) bracket(pval) addstat("Pseudo-R_squared", e(r2_p),"Incremental Pseudo-R_squared", `r2add')		   
	}
	
	local i = `i'+1
}
