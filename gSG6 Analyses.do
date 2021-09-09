// Author: Ellen Kearney
//Project: Anopheles Salivary Biomarkers - Systematic Review with multilevel modelling
//Contact: ellen.kearney@burnet.edu.au
*note: all analyses are performed with Stata 15.1

cd "C:\Users\ellen.kearney\Desktop\eLife\Do and Dta files" 
use "Anopheles salivary biomarkers review.dta", clear
************

gen log10_hbr_estimate = log10(hbr_estimate)

gen ln_hbr_estimate = ln(hbr_estimate)
lab var ln_hbr_estimate "log (HBR)"

gen ln_eir_estimate = ln(eir_estimate)
lab var ln_eir_estimate "log (EIR)"

gen ln_hbr1_estimate = ln(hbr_estimate +1)

********************************************************
*generate a variable that is both gSG6 and gSG6-P1 IgG*
gen all_gsg6_igg_seroprevalence=gsg6_igg_seroprevalence

replace all_gsg6_igg_seroprevalence = gsg6p1_igg_seroprevalence if gsg6p1_igg_seroprevalence!=.
lab var all_gsg6_igg_seroprevalence "IgG response to both recombinant and peptide gSG6"
gen all_gsg6_antigen=.
recode all_gsg6_antigen .=1 if gsg6p1_igg_seroprevalence!=.
recode all_gsg6_antigen .=0 if gsg6_igg_seroprevalence!=.
lab var all_gsg6_antigen "IgG gSG6 Antigen"
lab define gSG6 0"gSG6" 1"gSG6-P1"
lab val all_gsg6_antigen gSG6
	
gen all_gsg6_igg_nsubpop=gsg6_igg_nsubpop
replace all_gsg6_igg_nsubpop = gsg6p1_igg_nsubpop if gsg6p1_igg_nsubpop!=.
lab var all_gsg6_igg_nsubpop "Sample size of subpop used to measure gSG6"

*gen integer for successes
gen  successes_all_gsg6 = all_gsg6_igg_seroprevalence*all_gsg6_igg_nsubpop
replace  successes_all_gsg6 = round(successes_all_gsg6,1)
lab var successes_all_gsg6 "# seropositive individuals"

*gen method for gSG6
gen  method_all_gsg6 = gsg6_igg_method
replace method_all_gsg6 = gsg6p1_igg_method if method_all_gsg6==""
encode(method_all_gsg6), gen(method_all_gsg6_n)


encode study_type, gen(study_type_id)
encode country, gen(country_id)
tab endemicity_map, gen(endemicity_map_d)
tab country_id, gen(country_id_d)
tab study_type_id, gen(study_type_id_d)
gen hbr_ln_sq = ln_hbr_estimate^2
encode participants_checked,  gen(participants_checked_n)
encode hbr_casedetect,  gen(hbr_casedetect_n)
encode eir_casedetect,  gen(eir_casedetect_n)


gen c1_hbr = country_id_d1*ln_hbr_estimate
gen c2_hbr = country_id_d2*ln_hbr_estimate
gen c3_hbr = country_id_d3*ln_hbr_estimate
gen c4_hbr = country_id_d4*ln_hbr_estimate
gen c5_hbr = country_id_d5*ln_hbr_estimate
gen c6_hbr = country_id_d6*ln_hbr_estimate
gen c7_hbr = country_id_d7*ln_hbr_estimate
gen c8_hbr = country_id_d8*ln_hbr_estimate
gen c9_hbr = country_id_d9*ln_hbr_estimate
gen c10_hbr = country_id_d10*ln_hbr_estimate
gen c11_hbr = country_id_d11*ln_hbr_estimate
gen c12_hbr = country_id_d12*ln_hbr_estimate
gen c13_hbr = country_id_d13*ln_hbr_estimate


* deriving logits again
gen n1 = successes_all_gsg6
gen d1 = all_gsg6_igg_nsubpop
gen sg6_prob_logit = logit((n1)/(d1))
gen prob = n1/d1

gen sample_se = sqrt((prob * (1- prob)) / d1)
gen wt1 = 1 / sample_se

gen pfpos_estimate = pfpositive_pcr_estimate
replace pfpos_estimate = pfpositive_micro_estimate if pfpos_estimate==.
gen pfpos_method =.
recode pfpos_method .=0 if pfpositive_micro_estimate!=.
recode pfpos_method .=1 if pfpositive_pcr_estimate!=.
recode pfpos_method 0=1 if pfpositive_pcr_estimate!=.
lab define malaria_detection 0"Microscopy" 1"PCR"
lab val pfpos_method malaria_detection

gen parapos_estimate = parapos_pcr_estimate
replace parapos_estimate = parapos_micro_estimate if parapos_estimate==.
gen parapos_method =.
recode parapos_method .=0 if parapos_micro_estimate!=.
recode parapos_method .=1 if parapos_pcr_estimate!=.
recode parapos_method 0=1 if parapos_pcr_estimate!=.
lab val parapos_method malaria_detection

gen pvpos_estimate = pvpos_pcr_estimate
replace pvpos_estimate = pvpos_micro_estimate if pvpos_estimate==.
gen pvpos_method =.
recode pvpos_method .=0 if pvpos_micro_estimate!=.
recode pvpos_method .=1 if pvpos_pcr_estimate!=.
recode pvpos_method 0=1 if pvpos_pcr_estimate!=.
lab val pvpos_method malaria_detection

gen anypspppos_estimate = parapos_estimate
replace anypspppos_estimate = pfpos_estimate if anypspppos_estimate==.
replace anypspppos_estimate = pvpos_estimate if anypspppos_estimate==.
gen anypspppos_method =.
recode anypspppos_method .=0 if parapos_method==0 | pfpos_method==0 | pvpos_method==0
recode anypspppos_method .=1 if parapos_method==1 | pfpos_method==1 | pvpos_method==1
lab val anypspppos_method malaria_detection
gen anypspppos_spp = .
recode anypspppos_spp .=0 if parapos_estimate!=.
recode anypspppos_spp .=1 if pfpos_estimate!=.
recode anypspppos_spp .=2 if pvpos_estimate!=.
lab define malaria_spp 0"Plasmodium spp." 1"Pf only" 2"Pv only"
lab val anypspppos_spp malaria_spp

gen anypspppos_micro_estimate = parapos_micro_estimate
replace anypspppos_micro_estimate = pfpositive_micro_estimate if anypspppos_micro_estimate==.
replace anypspppos_micro_estimate = pvpos_micro_estimate if anypspppos_micro_estimate==.
gen anypspppos_pcr_estimate = parapos_pcr_estimate
replace anypspppos_pcr_estimate = pfpositive_pcr_estimate if anypspppos_pcr_estimate==.
replace anypspppos_pcr_estimate = pvpos_pcr_estimate if anypspppos_pcr_estimate==.
gen anypspppos_micro_spp = .
recode anypspppos_micro_spp .=0 if parapos_micro_estimate!=.
recode anypspppos_micro_spp .=1 if pfpositive_micro_estimate!=.
recode anypspppos_micro_spp .=2 if pvpos_micro_estimate!=.
lab val anypspppos_micro_spp malaria_spp

gen anypspppos_pcr_spp = .
recode anypspppos_pcr_spp .=0 if parapos_pcr_estimate!=.
recode anypspppos_pcr_spp .=1 if pfpositive_pcr_estimate!=.
recode anypspppos_pcr_spp .=2 if pvpos_pcr_estimate!=.
lab val anypspppos_pcr_spp malaria_spp

foreach var of varlist *pos_estimate *_pcr_estimate *_micro_estimate {
                gen p_`var' = `var'*100
        }
		
lab var p_pfpositive_pcr_estimate "P. falciparum Prevalence by PCR (%)"
lab var  p_parapos_pcr_estimate "Plasmodium spp. Prevalence by PCR (%)"
lab var  p_pvpos_pcr_estimate "P. vivax Prevalence by PCR (%)"
lab var  p_pfpositive_micro_estimate "P. falciparum Prevalence by Microscopy (%)"
lab var  p_parapos_micro_estimate "Plasmodium spp. Prevalence by Microscopy (%)"
lab var  p_pvpos_micro_estimate "P. vivax Prevalence by Microscopy (%)"
lab var  p_anypspppos_micro_estimate "Any Plasmodium spp. Prevalence by Microscopy (%)"
lab var  p_anypspppos_pcr_estimate "Any Plasmodium spp. Prevalence by PCR (%)"

rename p_pfpositive_pcr_estimate p_pfpcr
rename  p_parapos_pcr_estimate p_pspppcr
rename  p_pvpos_pcr_estimate p_pvpcr
rename  p_pfpositive_micro_estimate p_pfmicro
rename  p_parapos_micro_estimate p_psppmicro
rename p_pvpos_micro_estimate p_pvmicro
rename  p_anypspppos_micro_estimate p_apsppmicro
rename  p_anypspppos_pcr_estimate p_apspppcr

foreach var of varlist *_igg_estimate *_igg1_estimate *_igg3_estimate {
                gen p_`var' = `var'*100
        }
		
lab var p_pfama1_igg_estimate "PfAMA1 IgG Seroprevalence (%)"
lab var p_pfmsp119_igg_estimate "PfMSP1(19) IgG Seroprevalence (%)"
lab var p_pfmsp2_igg_estimate "PfMSP2 IgG Seroprevalence (%)"
lab var p_pfcsp_igg_estimate "PfCSP IgG Seroprevalence (%)"
lab var p_pfse_igg_estimate "PfSchizont Extract IgG Seroprevalence (%)"
lab var p_pfglurp_igg_estimate "PfGLURP IgG Seroprevalence (%)"
lab var p_pfmsp3_igg_estimate "PfMSP3 IgG Seroprevalence (%)"
lab var p_pvama1_igg_estimate "PvAMA1 IgG Seroprevalence (%)"
lab var p_pvmsp119_igg_estimate "PvMSP1(19) IgG Seroprevalence (%)"
lab var p_pfcsp_igg1_estimate "PfCSP IgG1 Seroprevalence (%)"
lab var p_pfse_igg1_estimate "PfSchizont Extract IgG1 Seroprevalence (%)"
lab var p_pfcsp_igg3_estimate "PfCSP IgG3 Seroprevalence (%)"
lab var p_pfse_igg3_estimate "PfSchizont Extract IgG3 Seroprevalence (%)"

rename p_pfama1_igg_estimate p_pfama1_igg
rename p_pfmsp119_igg_estimate p_pfmsp119_igg
rename p_pfmsp2_igg_estimate p_pfmsp2_igg
rename p_pfcsp_igg_estimate p_pfcsp_igg
rename p_pfse_igg_estimate p_pfse_igg
rename p_pfglurp_igg_estimate p_pfglurp_igg
rename p_pfmsp3_igg_estimate p_pfmsp3_igg
rename p_pvama1_igg_estimate p_pvama1_igg
rename p_pvmsp119_igg_estimate p_pvmsp119_igg
rename p_pfcsp_igg1_estimate p_pfcsp_igg1
rename p_pfse_igg1_estimate p_pfse_igg1
rename p_pfcsp_igg3_estimate p_pfcsp_igg3
rename p_pfse_igg3_estimate p_pfse_igg3


gen dvs_albi = 1 if dvs_full==1
gen dvs_arab = 1 if dvs_full==2 | dvs_full==3
gen dvs_diru = 1 if dvs_full==4 | dvs_full==11
gen dvs_fara = 1 if dvs_full==5 
gen dvs_fune = 1 if dvs_full==3 | dvs_full==6 | dvs_full==8
gen dvs_gamb = 1 if dvs_full==3 |  dvs_full==7 | dvs_full==8
gen dvs_mini = 1 if dvs_full==9 |  dvs_full==10 | dvs_full==11
gen dvs_macu = 1 if dvs_full==10 |  dvs_full==11 
gen dvs_phar = 1 if dvs_full==12 
gen dvs_gambsl = 1  if dvs_full==2 |dvs_full==3 |  dvs_full==7 | dvs_full==8

recode dvs_albi .=0 if dvs_albi==.
recode dvs_arab .=0 if dvs_arab ==. 
recode dvs_diru .=0 if dvs_diru ==. 
recode dvs_fara .=0 if dvs_fara ==. 
recode dvs_fune .=0 if dvs_fune ==.
recode dvs_gamb .=0 if dvs_gamb ==. 
recode dvs_mini .=0 if dvs_mini==.
recode dvs_macu .=0 if dvs_macu==. 
recode dvs_phar .=0 if dvs_phar==.
recode dvs_gambsl .=0 if dvs_gambsl ==. 

gen dvs_2 =.
recode dvs_2 .=0 if dvs_gambsl ==0
recode dvs_2 .=1 if dvs_gambsl ==1

lab def DVS_2 0"Non-An. gambiae s.l." 1"An. gambiae s.l." 
lab val dvs_2 DVS_2

gen dvs_3 =.
recode dvs_3 .=0 if dvs_full ==11 | dvs_full==10 | dvs_full ==9 | dvs_full ==5 | dvs_full ==1 | dvs_full ==4
recode dvs_3 .=1 if dvs_full ==3 | dvs_full ==6  | dvs_full ==8 | dvs_full ==12
recode dvs_3 .=2 if dvs_full ==2 | dvs_full ==7

lab def DVS_3 0"Non-An. gambiae s.l." 1"Non-dominant An. gambiae s.l." 2"Dominant An. gambiae s.l."
lab val dvs_3 DVS_3

gen dvs_3_2 =.
recode dvs_3_2 .=3 if dvs_full ==11 | dvs_full==10 | dvs_full ==9 | dvs_full ==5 | dvs_full ==1 | dvs_full ==4
recode dvs_3_2 .=2 if dvs_full ==3 | dvs_full ==6  | dvs_full ==8 | dvs_full ==12
recode dvs_3_2 .=1 if dvs_full ==2 | dvs_full ==7

lab def DVS_3_2 3"Non-An. gambiae s.l." 2"Non-dominant An. gambiae s.l." 1"Dominant An. gambiae s.l." , replace
lab val dvs_3_2 DVS_3_2

/*
* plots of observed data
scatter sg6_prob_logit hbr_estimate
scatter sg6_prob_logit ln_hbr_estimate


* study as level 2
xtmelogit successes_all_gsg6 hbr_estimate, binomial(all_gsg6_igg_nsubpop) || id: , 	 var 
estimates store id
estat ic

* study  as level 3 and sub population as level 2
xtmelogit successes_all_gsg6 hbr_estimate, binomial(all_gsg6_igg_nsubpop) || id: || subpopulation: , 	 var 
estimates store id_subpop
estat ic
lrtest id id_subpop

* country as level 3 and study as level 2 
xtmelogit successes_all_gsg6 hbr_estimate, binomial(all_gsg6_igg_nsubpop) || country: || id: , 	 var iter(100)
estimates store country_id
lrtest id country_id
estat ic

* country as level 3 and study as level 2 log hbr - linear
xtmelogit successes_all_gsg6 ln_hbr_estimate, binomial(all_gsg6_igg_nsubpop) || country: || id: , 	 var iter(100)
estimates store a
*lrtest id country_id
estat ic
foreach i of numlist 1.01 1.05 1.10 1.20 1.50 1.75 2.0  {
 nlcom `i'^_b[ln_hbr_estimate]
}


* country as level 3 and study as level 2 log hbr - linear. Random slope for study. Emperical support shown
xtmelogit successes_all_gsg6 ln_hbr_estimate, binomial(all_gsg6_igg_nsubpop) || country: || id: ln_hbr_estimate, 	 var 
estimates store b
*lrtest id country_id
estat ic
foreach i of numlist 1.01 1.05 1.10 1.20 1.50 1.75 2.0  {
 nlcom `i'^_b[ln_hbr_estimate]
}
lrtest a b 

* country as level 3 and study as level 2 log hbr - linear. Random slope for study, covariance. No emperical support shown
xtmelogit successes_all_gsg6 ln_hbr_estimate, binomial(all_gsg6_igg_nsubpop) || country: || id: ln_hbr_estimate, cov(unstructured)	 var 
estimates store c
*lrtest id country_id
estat ic
foreach i of numlist 1.01 1.05 1.10 1.20 1.50 1.75 2.0  {
 nlcom `i'^_b[ln_hbr_estimate]
}
lrtest b c 

* country as level 3 and study as level 2 log hbr - linear. Random slope for country. 
xtmelogit successes_all_gsg6 ln_hbr_estimate, binomial(all_gsg6_igg_nsubpop) || country: ln_hbr_estimate || id: ,  var 
estimates store d
*lrtest id country_id
estat ic
foreach i of numlist 1.01 1.05 1.10 1.20 1.50 1.75 2.0  {
 nlcom `i'^_b[ln_hbr_estimate]
}
lrtest a d 

* country as level 3 and study as level 2 log hbr - linear. Random slope for study and country, covariance. No emperical support shown
xtmelogit successes_all_gsg6 ln_hbr_estimate, binomial(all_gsg6_igg_nsubpop) || country: ln_hbr_estimate || id: ln_hbr_estimate,  var 
estimates store e
*lrtest id country_id
estat ic
foreach i of numlist 1.01 1.05 1.10 1.20 1.50 1.75 2.0  {
 nlcom `i'^_b[ln_hbr_estimate]
}
lrtest b e 

* country as level 3 and study as level 2 log hbr - quadratic
xtmelogit successes_all_gsg6 c.ln_hbr_estimate##c.ln_hbr_estimate , binomial(all_gsg6_igg_nsubpop) || country: || id: , 	 var 
estimates store f
estat ic

* country as level 3 and study as level 2 log hbr - quadratic. Random slope for study (linear and square term)
*xtmelogit successes_all_gsg6 c.ln_hbr_estimate##c.ln_hbr_estimate , binomial(all_gsg6_igg_nsubpop) || country: || id: ln_hbr_estimate hbr_ln_sq, 	 var 
*estimates store g
*estat ic


/*
xtmelogit successes_all_gsg6 hbr_estimate, binomial(all_gsg6_igg_nsubpop) diff || country: || id: || subpopulation:, 	 var  
estimates store country_id_sub
lrtest id country_id_sub
*/ 



* country as level 3 and study as level 2 log hbr - cubic
xtmelogit successes_all_gsg6 c.ln_hbr_estimate#c.ln_hbr_estimate#c.ln_hbr_estimate c.ln_hbr_estimate##c.ln_hbr_estimate   , binomial(all_gsg6_igg_nsubpop) || country: || id: , 	 var 
estimates store g
lrtest f g
estat ic


* using a linear model with single linear term for log hbr
xtmixed sg6_prob_logit c.ln_hbr_estimate , || country: || id: , 	 var 
est store a
estat ic
lincom _b[ln_hbr_estimate]
foreach i of numlist 1.01 1.05 1.10 1.20 1.50 1.75 2.0  {
 nlcom `i'^_b[ln_hbr_estimate]
}

* constraining residual term for study specific standard error
constraint 1 _b[/var(sample_se[country>id])] = 1
meglm sg6_prob_logit c.ln_hbr_estimate , || country: || id: sample_se, 	 nocons constraints(1) 
meglm sg6_prob_logit c.ln_hbr_estimate##c.ln_hbr_estimate , || country:   || id: sample_se, 	 nocons constraints(1)

* using a linear model with single linear term for log hbr. random slope for study
xtmixed sg6_prob_logit c.ln_hbr_estimate , || country: || id: ln_hbr_estimate, 	 var  cov(unstructured) 
est store b
estat ic
lincom _b[ln_hbr_estimate]
foreach i of numlist 1.01 1.05 1.10 1.20 1.50 1.75 2.0  {
 nlcom `i'^_b[ln_hbr_estimate]
}
lrtest a b

* using a linear model with single linear term for log hbr. random slope for study
xtmixed sg6_prob_logit c.ln_hbr_estimate , || country: ln_hbr_estimate,  var cov(unstructured)  || id: , 	
est store c
estat ic
lincom _b[ln_hbr_estimate]
foreach i of numlist 1.01 1.05 1.10 1.20 1.50 1.75 2.0  {
 nlcom `i'^_b[ln_hbr_estimate]
}
lrtest a c



* using a linear model with single linear term for log hbr. random slope for study
xtmixed sg6_prob_logit c.ln_hbr_estimate , || country: ln_hbr_estimate,  || id: ln_hbr_estimate, 	 var  
est store d
estat ic
lincom _b[ln_hbr_estimate]
foreach i of numlist 1.01 1.05 1.10 1.20 1.50 1.75 2.0  {
 nlcom `i'^_b[ln_hbr_estimate]
}
lrtest a d


* using a linear model with a quadratic function for log hbr 
xtmixed sg6_prob_logit c.ln_hbr_estimate##c.ln_hbr_estimate , || country:  || id:  , 	 var 
estat ic 
lincom _b[ln_hbr_estimate]
lincom _b[ln_hbr_estimate#ln_hbr_estimate]


*estimates splines for ln_hbr
*estimate splines*
mkspline ln_hbr2_1 2.3025851 ln_hbr2_2  = ln_hbr_estimate,


*include it in the model
xtmelogit successes_all_gsg6 ln_hbr2_1-ln_hbr2_2 , binomial(all_gsg6_igg_nsubpop) ///
|| id: || subpopulation: , or var 
estimates store spline2_ln
estat ic

*predicted fitted values incl. random effects*
predict all_splhbr_mu, mu
predict all_splhbr_fixed, xb
predict all_splhbr_fixed_se, stdp 
predict all_splhbr_reeffects*, reffects

replace all_splhbr_fixed=. if all_gsg6_igg_seroprevalence==.
replace all_splhbr_fixed_se=. if all_gsg6_igg_seroprevalence==.

gen all_splhbr_fitted = (all_splhbr_fixed + all_splhbr_reeffects1 + all_splhbr_reeffects2)

gen all_splhbr_fixed_lci = all_splhbr_fixed - (1.96*all_splhbr_fixed_se)
gen all_splhbr_fixed_uci = all_splhbr_fixed + (1.96*all_splhbr_fixed_se)

gen invlogit_all_splhbr_fitted = invlogit(all_splhbr_fitted)
gen invlogit_all_splhbr_fixed = invlogit(all_splhbr_fixed)
gen invlogit_all_splhbr_fixed_lci = invlogit(all_splhbr_fixed_lci)
gen invlogit_all_splhbr_fixed_uci = invlogit(all_splhbr_fixed_uci)

*observed vs predicted 
twoway (rarea invlogit_all_splhbr_fixed_lci invlogit_all_splhbr_fixed_uci ///
ln_hbr_estimate, sort col(cranberry) fintensity(30)) ///
(scatter all_gsg6_igg_seroprevalence ln_hbr_estimate [aweight=all_gsg6_igg_nsubpop] , col(navy)) ///
(line invlogit_all_splhbr_fixed ln_hbr_estimate, col(red) sort) , legend(order(2 3 1)) 



*********************estimate splines*
mkspline ln_hbr3_1 1.6094379 ln_hbr3_2 2.3025851 ln_hbr3_3 = ln_hbr_estimate,
*include it in the model
xtmelogit successes_all_gsg6 ln_hbr3_1-ln_hbr3_3 , binomial(all_gsg6_igg_nsubpop) ///
|| id: || subpopulation: , or var 
estimates store spline3_ln
estat ic

*predicted fitted values incl. random effects*
predict all_spl3hbr_mu, mu
predict all_spl3hbr_fixed, xb
predict all_spl3hbr_fixed_se, stdp 
predict all_spl3hbr_reeffects*, reffects

replace all_spl3hbr_fixed=. if all_gsg6_igg_seroprevalence==.
replace all_spl3hbr_fixed_se=. if all_gsg6_igg_seroprevalence==.

gen all_spl3hbr_fitted = (all_spl3hbr_fixed + all_spl3hbr_reeffects1 + all_spl3hbr_reeffects2)

gen all_spl3hbr_fixed_lci = all_spl3hbr_fixed - (1.96*all_spl3hbr_fixed_se)
gen all_spl3hbr_fixed_uci = all_spl3hbr_fixed + (1.96*all_spl3hbr_fixed_se)

gen invlogit_all_spl3hbr_fitted = invlogit(all_spl3hbr_fitted)
gen invlogit_all_spl3hbr_fixed = invlogit(all_spl3hbr_fixed)
gen invlogit_all_spl3hbr_fixed_lci = invlogit(all_spl3hbr_fixed_lci)
gen invlogit_all_spl3hbr_fixed_uci = invlogit(all_spl3hbr_fixed_uci)

*observed vs predicted 
twoway (rarea invlogit_all_spl3hbr_fixed_lci invlogit_all_spl3hbr_fixed_uci ///
ln_hbr_estimate, sort col(cranberry) fintensity(30)) ///
(scatter all_gsg6_igg_seroprevalence ln_hbr_estimate [aweight=all_gsg6_igg_nsubpop] , col(navy)) ///
(line invlogit_all_spl3hbr_fixed ln_hbr_estimate, col(red) sort) , legend(order(2 3 1)) 

*observed vs predicted including random effects
twoway (scatter all_gsg6_igg_seroprevalence ln_hbr_estimate, col(navy)) ///
(scatter invlogit_all_spl3hbr_fitted ln_hbr_estimate, col(cranberry)) ///
(line invlogit_all_spl3hbr_fixed ln_hbr_estimate, col(red) sort)  if ///
hbr_estimate!=. & all_gsg6_igg_seroprevalence!=. ,  by (id) 

*move into GSEM as some models no longer conver in xtmelogit
gsem (ln_hbr_estimate -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) (M1[id] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) , covstruct(_lexogenous, diagonal) latent(M1   ) nocapslatent noest
matrix e = e(b)
gsem (ln_hbr_estimate -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) (M1[id] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) , covstruct(_lexogenous, diagonal) latent(M1   ) nocapslatent from(e)
est store twolvl
estat ic

gsem (ln_hbr_estimate -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) (M1[id] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit))  (M2[id]#c.ln_hbr_estimate -> successes_all_gsg6), covstruct(_lexogenous, diagonal) latent(M1 M2  ) nocapslatent
est store twolvl_rslope
estat ic

gsem (ln_hbr_estimate -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) (M1[country] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) (M3[country>id] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) , covstruct(_lexogenous, diagonal) latent(M1 M3 ) nocapslatent noest
matrix c = e(b)
gsem (ln_hbr_estimate -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) (M1[country] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) (M3[country>id] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) , covstruct(_lexogenous, diagonal) latent(M1 M3 ) nocapslatent from(c)
est store threelvl
estat ic

gsem (ln_hbr_estimate -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) (M1[country] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) (M3[country>id] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) (M2[country>id]#c.ln_hbr_estimate -> successes_all_gsg6), covstruct(_lexogenous, diagonal) latent(M1 M2 M3 ) nocapslatent
est store threelvl_rslope


lrtest twolvl twolvl_rslope // support for random slope
lrtest twolvl threelvl // no support for inclusion of country
lrtest threelvl_rslope twolvl_rslope // no support for inclusion of country
*/

// run models for given k-fold changes in x use this to estimate the correct 95% CI for log transformed exposures
gsem (ln_hbr_estimate -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) (M1[country] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) (M3[country>id] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) (M2[country>id]#c.ln_hbr_estimate -> successes_all_gsg6), covstruct(_lexogenous, diagonal) latent(M1 M2 M3 ) nocapslatent
estimates store hbr
estimates save hbr, replace

//Generate data for Figure 3
foreach i of numlist 1.01 1.05 1.10 1.20 1.50 1.75 2.0  {

estimates use hbr.ster

  nlcom `i'^_b[ln_hbr_estimate] // CIs are not correct with this syntax

  nlcom log(`i'^_b[ln_hbr_estimate]), post // need to work with the log

  matrix b = e(b)

  matrix x = e(V)

  local se = sqrt(x[1,1])

  local b = b[1,1]

* Standard error 

  disp "SE: "  "`se'"

*** 95% CIs

display "lower limit: " exp(`b' -invnormal(0.975)* `se')

display "upper limit: " exp(`b' +invnormal(0.975)* `se')

}
est restore hbr
* incorporate HBR random slope variance to get 95% reference range
global orhbr = 1.5^(_b[ln_hbr_estimate])
global rrhbrll95 = 1.5^(_b[ln_hbr_estimate] - (sqrt(_b[/var(M2[country>id])]) * 1.96))
global rrhbrul95 = 1.5^(_b[ln_hbr_estimate] + (sqrt(_b[/var(M2[country>id])]) * 1.96))

disp "odds ratio " $orhbr 
disp "lower 95% reference range " $rrhbrll95  
disp "upper 95% reference range " $rrhbrul95 


//now lets do the same for the hbrXdvs interaction estimates
gsem (ln_hbr_estimate -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(M1[country] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(M3[country>id] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(2.dvs_3_2 -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(3.dvs_3_2 -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(2.dvs_3_2#c.ln_hbr_estimate -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(3.dvs_3_2#c.ln_hbr_estimate -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(M2[country>id]#c.ln_hbr_estimate -> successes_all_gsg6), ///
covstruct(_lexogenous, diagonal) latent(M1 M2 M3 ) nocapslatent 
est store hbrXdvs
est save hbrXdvs , replace

*corrected 95%CIs
foreach i of numlist 2 3  {

estimates use hbrXdvs.ster
gsem

  nlcom 1.5^(_b[ln_hbr_estimate] + _b[`i'.dvs_3_2#c.ln_hbr_estimate]) // CIs are not correct with this syntax

  nlcom log(1.5^(_b[ln_hbr_estimate] + _b[`i'.dvs_3_2#c.ln_hbr_estimate])), post // need to work with the log
  matrix b = e(b)
  matrix x = e(V)
  local se = sqrt(x[1,1])
  local b = b[1,1]

* Standard error 
  disp "SE: "  "`se'"

*** 95% CIs
display "lower limit: " exp(`b' -invnormal(0.975)* `se')
display "upper limit: " exp(`b' +invnormal(0.975)* `se')

}

estimates use hbrXdvs.ster 
gsem
//use this to calculate the conditional icc for each model
disp "*** level-3 - country (same country, different study)"
disp _b[/var(M1[country])]  / (_b[/var(M1[country])] + _b[/var(M3[country>id])] + _b[/var(M2[country>id])]+ ((3.14^2)/3))
disp "*** level-2 - study  (same country, same study)"
disp (_b[/var(M1[country])] + _b[/var(M3[country>id])]) / (_b[/var(M1[country])] + _b[/var(M3[country>id])] + _b[/var(M2[country>id])]+ ((3.14^2)/3))
disp "*** level-2 - study + HBR (same country, same study, same effect of HBR)"
disp (_b[/var(M1[country])] + _b[/var(M3[country>id])] + _b[/var(M2[country>id])]) / (_b[/var(M1[country])] + _b[/var(M3[country>id])] + _b[/var(M2[country>id])]+ ((3.14^2)/3))

//interact with entomological detection method
gen hbr_entomethod = 1 if hbr_casedetect_n==1 | hbr_casedetect_n==2 | hbr_casedetect_n==4
recode hbr_entomethod .= 0 if hbr_casedetect_n==3

gsem (ln_hbr_estimate -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(M1[country] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(M3[country>id] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(1.hbr_entomethod -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(1.hbr_entomethod#c.ln_hbr_estimate -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(M2[country>id]#c.ln_hbr_estimate -> successes_all_gsg6), ///
covstruct(_lexogenous, diagonal) latent(M1 M2 M3 ) nocapslatent 
est store hbrXento
est save hbrXento

//Figure 2
/*
preserve
set obs 540
replace id = 9999 if _n > 387
egen obs2 = seq() if _n > 387
replace obs = obs2 if id==9999

 

egen hbr_estimate_2 = seq() if _n > 389
replace hbr_estimate = hbr_estimate_2 if id==9999
replace  hbr_estimate = 0.01 if obs==1 & id==9999
replace  hbr_estimate = 0.1 if obs==2 & id==9999


est restore hbr
predictnl sg6_hbr_rc_prob_predictnl = invlogit(_b[_cons] + _b[ln_hbr_estimate]*log(hbr_estimate)) if id == 9999, se(sg6_hbr_rc_prob_predictnl_se) ci(sg6_hbr_rc_prob_predictnl_lci sg6_hbr_rc_prob_predictnl_uci)

twoway rarea sg6_hbr_rc_prob_predictnl_uci sg6_hbr_rc_prob_predictnl_lci hbr_estimate if id == 9999 & hbr_estimate < 14, color(red%6) lw(none) || ///
scatter prob hbr_estimate if hbr_estimate < 12.2  [aweight=all_gsg6_igg_nsubpop] , mcolor(white) mlcolor(black)  msize(small) || ///
function y=invlogit(_b[_cons] + _b[ln_hbr_estimate]*(log(x))), range(0 13) n(300) lcolor(red)  graphregion(fcolor(white))  graphregion(lcolor(white))  graphregion(ilc(white))   plotregion(ilc(white))  ///
legend(order( 2 "Observed study anti-gSG6 IgG seroprevalence" - " " 3 "Average anti-gSG6 IgG probability (95%CI) by HBR"  " " 1 "" ) size(small)) ytitle(anti-gSG6 IgG seroprevalence (%), size(small)) xtitle(HBR (bites per person per night), size(small)) ///
xlabel(,labs(small)) ylabel(,labs(small)) play(HBR_yaxis)

*data for figure 3a + 4a
br hbr_estimate  sg6_hbr_rc_prob_predictnl  sg6_hbr_rc_prob_predictnl_lci  sg6_hbr_rc_prob_predictnl_uci   if obs==1|obs==2 | hbr_estimate==0.01 | hbr_estimate==0.1  | hbr_estimate==1 | hbr_estimate==5 | hbr_estimate==10 | hbr_estimate==50 | hbr_estimate==100 | hbr_estimate==150 


restore
*/

//Figure 4 
/*
//recode for dummy indicators
gen dvs_angam = 1 if dvs_3_2==1
recode dvs_angam .=0 if dvs_3_2==2 | dvs_3_2==3
gen dvs_angam_others = 1 if dvs_3_2==2
recode dvs_angam_others .=0 if dvs_3_2==3 | dvs_3_2==1
gen dvs_nonangam = 1 if dvs_3_2==3
recode dvs_nonangam .=0 if dvs_3_2==2 | dvs_3_2==1


//rerun model with dummy indicators
gsem (ln_hbr_estimate -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(M1[country] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(M3[country>id] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(1.dvs_angam_others -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(1.dvs_nonangam -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(1.dvs_angam_others#c.ln_hbr_estimate -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(1.dvs_nonangam#c.ln_hbr_estimate -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(M2[country>id]#c.ln_hbr_estimate -> successes_all_gsg6), ///
covstruct(_lexogenous, diagonal) latent(M1 M2 M3 ) nocapslatent 
est store hbrXdvs_cat

*******************************************************************************************
//produce dummy dataset
preserve
set obs 540
egen obs2 = seq() if _n > 387
egen hbr_estimate_2 = seq() if _n > 389
recode dvs_angam .= 1 if _n > 387

set obs 693
egen obs3 = seq() if _n > 540
egen hbr_estimate_3 = seq() if _n > 542
recode dvs_angam_others .= 1 if _n > 540

set obs 846
egen obs4 = seq() if _n > 693
egen hbr_estimate_4 = seq() if _n > 695
recode dvs_nonangam .= 1 if _n > 693

replace obs = obs2
replace obs = obs3 if obs==.
replace obs = obs4 if obs==.
 replace id = 9999 if _n > 387

replace hbr_estimate = hbr_estimate_2 if id == 9999
replace hbr_estimate = hbr_estimate_3 if id == 9999 & hbr_estimate==.
replace hbr_estimate = hbr_estimate_4 if id == 9999 & hbr_estimate==.

replace  hbr_estimate = 0.01 if obs==1 & id==9999
replace  hbr_estimate = 0.1 if obs==2 & id==9999



recode dvs_angam .=0 if id==9999
recode dvs_angam_others .=0 if id==9999
recode dvs_nonangam .=0 if id==9999

*******************************************************************************************
//predict out using dummy dataset 
est restore hbrXdvs_cat
predictnl sg6_hbrxdvs_prob_predictnl = invlogit(_b[_cons] + _b[ln_hbr_estimate]*log(hbr_estimate) + _b[1.dvs_angam_others#c.ln_hbr_estimate]*log(hbr_estimate)*dvs_angam_others + _b[1.dvs_angam_others]*dvs_angam_others + _b[1.dvs_nonangam#c.ln_hbr_estimate]*log(hbr_estimate)*dvs_nonangam + _b[1.dvs_nonangam]*dvs_nonangam) if id == 9999, se(sg6_hbrxdvs_prob_predictnl_se) ci(sg6_hbrxdvs_prob_predictnl_lci sg6_hbrxdvs_prob_predictnl_uci)

***********************************************************************************
//data for inclusion in Figure 4

br hbr_estimate  sg6_hbrxdvs_prob_predictnl  sg6_hbrxdvs_prob_predictnl_lci  sg6_hbrxdvs_prob_predictnl_uci  sg6_hbrxdvs_prob_predictnl_uci if obs==1|obs==2 | hbr_estimate==0.01 | hbr_estimate==0.1  | hbr_estimate==1 | hbr_estimate==5 | hbr_estimate==10 | hbr_estimate==50 | hbr_estimate==100 | hbr_estimate==150 


restore
*/

***********************************************************************************
**************************************EIR*****************************************
***********************************************************************************

gsem (ln_eir_estimate -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) (M1[country] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) (M3[country>id] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) (M2[country>id]#c.ln_eir_estimate -> successes_all_gsg6), covstruct(_lexogenous, diagonal) latent(M1 M2 M3 ) nocapslatent
estimates store eir
estimates save eir, replace


foreach i of numlist 1.01 1.05 1.10 1.20 1.50 1.75 2.0  {

estimates use eir.ster

  nlcom `i'^_b[ln_eir_estimate] // CIs are not correct with this syntax

  nlcom log(`i'^_b[ln_eir_estimate]), post // need to work with the log

  matrix b = e(b)

  matrix x = e(V)

  local se = sqrt(x[1,1])

  local b = b[1,1]

* Standard error 

  disp "SE: "  "`se'"

*** 95% CIs

display "lower limit: " exp(`b' -invnormal(0.975)* `se')

display "upper limit: " exp(`b' +invnormal(0.975)* `se')

}

est restore eir

* incorporate EIR random slope variance to get 95% reference range
global oreir = 1.5^(_b[ln_eir_estimate])
global rreirll95 = 1.5^(_b[ln_eir_estimate] - (sqrt(_b[/var(M2[country>id])]) * 1.96))
global rreirul95 = 1.5^(_b[ln_eir_estimate] + (sqrt(_b[/var(M2[country>id])]) * 1.96))

disp "odds ratio " $oreir 
disp "lower 95% reference range " $rreirll95  
disp "upper 95% reference range " $rreirul95


//EIRxDVS
gsem (ln_eir_estimate -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(M1[country] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(M3[country>id] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(2.dvs_3_2 -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(3.dvs_3_2 -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(2.dvs_3_2#c.ln_eir_estimate -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(3.dvs_3_2#c.ln_eir_estimate -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) ///
(M2[country>id]#c.ln_eir_estimate -> successes_all_gsg6), ///
covstruct(_lexogenous, diagonal) latent(M1 M2 M3 ) nocapslatent 
est store eirXdvs
est save eirXdvs , replace

*corrected 95%CIs
foreach i of numlist 2 3  {

estimates use eirXdvs.ster
gsem

  nlcom 1.5^(_b[ln_eir_estimate] + _b[`i'.dvs_3_2#c.ln_eir_estimate]) // CIs are not correct with this syntax

  nlcom log(1.5^(_b[ln_eir_estimate] + _b[`i'.dvs_3_2#c.ln_eir_estimate])), post // need to work with the log
  matrix b = e(b)
  matrix x = e(V)
  local se = sqrt(x[1,1])
  local b = b[1,1]

* Standard error 
  disp "SE: "  "`se'"

*** 95% CIs
display "lower limit: " exp(`b' -invnormal(0.975)* `se')
display "upper limit: " exp(`b' +invnormal(0.975)* `se')

}


//interact with entomological detection method
gen eir_entomethod = 1 if eir_casedetect_n==1 | eir_casedetect_n==4
recode eir_entomethod .= 0 if eir_casedetect_n==2 | eir_casedetect_n==3

xtmelogit successes_all_gsg6 c.ln_eir_estimate##i.eir_entomethod, binomial(all_gsg6_igg_nsubpop) ///
|| country: || id: ln_eir_estimate, 	 var iter(30)


est store eirXento
est save eirXento
estat ic



//Figure EIR (not included in final manuscript)
/*
preserve
set obs 538
replace id = 9999 if _n > 387
egen obs2 = seq() if _n > 387
replace obs = obs2

egen eir_estimate_2 = seq() if _n > 387
replace eir_estimate = eir_estimate_2 - 1 if id == 9999
replace eir_estimate = .01 if eir_estimate_2 == 1

est restore eir 
predictnl sg6_eir_rc_prob_predictnl = invlogit(_b[_cons] + _b[ln_eir_estimate]*log(eir_estimate)) if id == 9999, se(sg6_eir_rc_prob_predictnl_se) ci(sg6_eir_rc_prob_predictnl_lci sg6_eir_rc_prob_predictnl_uci)

twoway rarea sg6_eir_rc_prob_predictnl_uci sg6_eir_rc_prob_predictnl_lci eir_estimate if id == 9999 & eir_estimate<37, color(red%6) lw(none) || ///
scatter prob eir_estimate if eir_estimate < 36.5  [aweight=all_gsg6_igg_nsubpop] , mcolor(white) mlcolor(black)  msize(small) || ///
function y=invlogit(_b[_cons] + _b[ln_eir_estimate]*(log(x))), range(0 36) n(300) lcolor(red)  graphregion(fcolor(white))  graphregion(lcolor(white))  graphregion(ilc(white))   plotregion(ilc(white))  ///
legend(order( 2 "Observed study anti-gSG6 IgG prevalence" - " " 3 "Average anti-gSG6 IgG probability (95%CI) by EIR " 1 "" ) size(small)) ytitle(anti-gSG6 IgG seroprevalence (%), size(small)) xtitle(EIR (infective bites per person per year), size(small)) ///
xlabel(,labs(small)) ylabel(,labs(small)) play(HBR_yaxis)

restore 
*/
*****************************************************************************************************************************
************************MALARIA PREVALENCE DATA***********************************
**********************************************************************************
foreach var of varlist parapos_estimate pfpos_estimate pvpos_estimate anypspppos_estimate {
gen p10_`var'=`var'*10
}

***********************************************************************************
*LETS EXPLORE ANYPSPP

/*
*exploratory analyses
xtmelogit successes_all_gsg6 p10_anypspppos_estimate, binomial(all_gsg6_igg_nsubpop) ///
|| id: , 	 var 
estat ic
estimates store anypspp_a

xtmelogit successes_all_gsg6 p10_anypspppos_estimate, binomial(all_gsg6_igg_nsubpop) ///
|| id: p10_anypspppos_estimate, 	 var iter(30)
estimates store anypspp_b
estat ic
lrtest anypspp_a anypspp_b // support for random slope

xtmelogit successes_all_gsg6 p10_anypspppos_estimate, binomial(all_gsg6_igg_nsubpop) ||country: || id: p10_anypspppos_estimate, 	 var iter(30)
estimates store anypspp_c
estat ic
lrtest anypspp_b anypspp_c //no support for country

gen ln_prev_10 = log(p10_anypspppos_estimate) 
gen sq_prev = p10_anypspppos_estimate*p10_anypspppos_estimate

gen ln_prev = log(p_anypspppos_estimate) 

xtmelogit successes_all_gsg6 ln_prev, binomial(all_gsg6_igg_nsubpop) ///
|| id: ln_prev, 	 var iter(30)
estat ic 

xtmelogit successes_all_gsg6 p10_anypspppos_estimate sq_prev, binomial(all_gsg6_igg_nsubpop) ///
|| id: p10_anypspppos_estimate, 	 var iter(30)
estat ic

*** THIS LN-PREV MODEL SHOWS BEST FIT (AIC AND BIC). 
*/



gen ln_prev = log(p_anypspppos_estimate) 

xtmelogit successes_all_gsg6 ln_prev, binomial(all_gsg6_igg_nsubpop) ///
|| id: ln_prev, 	 var iter(30)
estat ic
est store ln_prev
est save ln_prev , replace

foreach i of numlist 1.10   {

estimates use ln_prev.ster

  nlcom `i'^_b[ln_prev] // CIs are not correct with this syntax

  nlcom log(`i'^_b[ln_prev]), post // need to work with the log

  matrix b = e(b)

  matrix x = e(V)

  local se = sqrt(x[1,1])

  local b = b[1,1]
  local n = e(N)

* Standard error 

  disp "SE: "  "`se'"

*** 95% CIs

display "lower limit: " exp(`b' -invnormal(0.975)* `se')

display "upper limit: " exp(`b' +invnormal(0.975)* `se')

di "Meta-N:" `n'
}
est restore ln_prev
xtmelogit
*** level-2 - study (same study)
disp exp(2*[lns1_1_1]_cons)  / ((exp(2*[lns1_1_1]_cons))+(exp(2*[lns1_1_2]_cons))+ (3.14^2)/3)
*** level-2 - study + effect (same study, same effect)
disp (exp(2*[lns1_1_1]_cons) +exp(2*[lns1_1_2]_cons)) / ((exp(2*[lns1_1_1]_cons)) +(exp(2*[lns1_1_2]_cons)) + (3.14^2 / 3))




xtmelogit successes_all_gsg6 c.ln_prev##i.anypspppos_method, binomial(all_gsg6_igg_nsubpop) ///
|| id: ln_prev, 	 var iter(30)
*** level-2 - study (same study)
disp exp(2*[lns1_1_1]_cons)  / ((exp(2*[lns1_1_1]_cons))+(exp(2*[lns1_1_2]_cons))+ (3.14^2)/3)
*** level-2 - study + effect (same study, same effect)
disp (exp(2*[lns1_1_1]_cons) +exp(2*[lns1_1_2]_cons)) / ((exp(2*[lns1_1_1]_cons)) +(exp(2*[lns1_1_2]_cons)) + (3.14^2 / 3))

xtmelogit successes_all_gsg6 c.ln_prev##i.anypspppos_spp, binomial(all_gsg6_igg_nsubpop) ///
|| id: ln_prev, 	 var iter(30)

*** level-2 - study (same study)
disp exp(2*[lns1_1_1]_cons)  / ((exp(2*[lns1_1_1]_cons))+(exp(2*[lns1_1_2]_cons))+ (3.14^2)/3)
*** level-2 - study + effect (same study, same effect)
disp (exp(2*[lns1_1_1]_cons) +exp(2*[lns1_1_2]_cons)) / ((exp(2*[lns1_1_1]_cons)) +(exp(2*[lns1_1_2]_cons)) + (3.14^2 / 3))


/* 
//Figure 5
preserve
set obs 488

replace id = 9999 if _n > 387

egen obs2 = seq() if _n > 387

replace obs = obs2

 

egen ln_prev_2 = seq() if _n > 387

replace ln_prev = ln_prev_2 - 1 if id == 9999
replace p_anypspppos_estimate = ln_prev_2 - 1 if id == 9999


predictnl sg6_hbr_rc_prob_predictnl = invlogit(_b[_cons] + _b[ln_prev]*log(ln_prev)) if id == 9999, se(sg6_ln_prev_prob_predictnl_se) ci(sg6_ln_prev_prob_predictnl_lci sg6_ln_prev_prob_predictnl_uci)

replace sg6_ln_prev_prob_predictnl_lci =0 if sg6_ln_prev_prob_predictnl_lci <0
replace sg6_ln_prev_prob_predictnl_uci =1 if sg6_ln_prev_prob_predictnl_uci >=1

replace sg6_ln_prev_prob_predictnl_lci =0 if id!=9999
replace sg6_ln_prev_prob_predictnl_uci =0 if id!=9999

twoway rarea sg6_ln_prev_prob_predictnl_uci sg6_ln_prev_prob_predictnl_lci p_anypspppos_estimate if id == 9999, color(red%6) lw(none) || ///
scatter prob p_anypspppos_estimate if all_gsg6_igg_seroprevalence < 100 & all_gsg6_igg_seroprevalence >0 [aweight=all_gsg6_igg_nsubpop] , mcolor(white) mlcolor(black)  msize(small) || ///
function y=invlogit(_b[_cons] + _b[ln_prev]*(log(x))), range(0 100) n(300) lcolor(red)  graphregion(fcolor(white))  graphregion(lcolor(white))  graphregion(ilc(white))   plotregion(ilc(white))  ///
legend(order( 2 "Observed study anti-gSG6 IgG prevalence" - " " 3 "Average anti-gSG6 IgG probability (95%CI) by Plasmodium spp. prevalence (%)"  " " 1 "" ) size(small)) ytitle(anti-gSG6 IgG seroprevalence (%), size(small)) xtitle(Plasmodium spp. prevalence (%), size(small)) ///
xlabel(,labs(small)) yscale(r(0(.2)1)) ylabel(,labs(small)) play(HBR_yaxis)
restore
*/

*************************************************************************************
****************************MALARIA ENDEMICITY CLASS*********************************
*************************************************************************************


xtmelogit successes_all_gsg6 i.endem_can5 , binomial(all_gsg6_igg_nsubpop) || id: ,  var iter(100) or
*xtmelogit successes_all_gsg6 i.endem_can5 , binomial(all_gsg6_igg_nsubpop) ||country: || id: ,  var  or



margins endem_can5 , expression(invlogit(predict(xb)))

*******************************************************************************************
//malarial seroprevalence
// 10% increase in exposure
foreach var of varlist *igg_estimate{
gen p_`var'= `var'*100

}
lab var p_pfama1_igg_estimate "PfAMA1 IgG Seroprevalence (%)"
lab var p_pfmsp119_igg_estimate "PfMSP1(19) IgG Seroprevalence (%)"
lab var p_pfmsp2_igg_estimate "PfMSP2 IgG Seroprevalence (%)"
lab var p_pfcsp_igg_estimate "PfCSP IgG Seroprevalence (%)"
lab var p_pfse_igg_estimate "PfSchizont Extract IgG Seroprevalence (%)"
lab var p_pfglurp_igg_estimate "PfGLURP IgG Seroprevalence (%)"
lab var p_pfmsp3_igg_estimate "PfMSP3 IgG Seroprevalence (%)"
lab var p_pvama1_igg_estimate "PvAMA1 IgG Seroprevalence (%)"
lab var p_pvmsp119_igg_estimate "PvMSP1(19) IgG Seroprevalence (%)"

foreach var of varlist p_pfama1_igg p_pfmsp119_igg p_pfmsp2_igg p_pfcsp_igg p_pfse_igg p_pfglurp_igg p_pfmsp3_igg p_pvama1_igg p_pvmsp119_igg{
gen ln_`var'= ln(`var')

}

foreach var of varlist p_pfama1_igg p_pfmsp119_igg p_pfmsp2_igg p_pfmsp3_igg p_pfcsp_igg p_pfglurp_igg p_pfse_igg  p_pvama1_igg p_pvmsp119_igg{

xtmelogit successes_all_gsg6 `var' , binomial(all_gsg6_igg_nsubpop) ///
|| id: , 	 var iter(60)
estat ic
est store `var'_a

xtmelogit successes_all_gsg6 `var' , binomial(all_gsg6_igg_nsubpop) ///
|| id: `var', 	 var iter(60)
estat ic
est store `var'_b

lrtest `var'_a `var'_b
}

foreach var of varlist p_pfama1_igg p_pfmsp3_igg { ///pfama1 and msp3 do not run converge with inclusion of a random slope

xtmelogit successes_all_gsg6 ln_`var' , binomial(all_gsg6_igg_nsubpop) ///
|| id: , 	 var iter(60)
estat ic
est store ln_`var'_a

xtmelogit successes_all_gsg6 ln_`var' , binomial(all_gsg6_igg_nsubpop) ///
|| id: ln_`var', 	 var iter(60)
estat ic
est store ln_`var'_b

lrtest ln_`var'_a ln_`var'_b
}

foreach var of varlist  p_pfmsp119_igg p_pfmsp2_igg  p_pfcsp_igg p_pfglurp_igg p_pfse_igg  p_pvama1_igg p_pvmsp119_igg{

xtmelogit successes_all_gsg6 ln_`var' , binomial(all_gsg6_igg_nsubpop) ///
|| id: , 	 var iter(60)
estat ic
est store ln_`var'_a

xtmelogit successes_all_gsg6 ln_`var' , binomial(all_gsg6_igg_nsubpop) ///
|| id: ln_`var', 	 var iter(60)
estat ic
est store ln_`var'_b

lrtest ln_`var'_a ln_`var'_b
}

//Table S8
/*

//10% increase in antimalarial seroprevalence
putexcel set gSG6_malariaantigens20210701.xlsx, sheet(logodds) modify

//give the columns names
putexcel A1 = "Variable" , border(bottom top) vcenter
putexcel B1 = "b" ,  border(bottom top) vcenter hcenter
putexcel C1= "95%CI" , border(bottom top) vcenter hcenter
putexcel E1 = "",  border(bottom top) vcenter hcenter
putexcel F1 = "P value",  border(bottom top) vcenter hcenter
putexcel G1 = "" , border (bottom top) vcenter hcenter
putexcel H1 = "Obs" , border (bottom top) vcenter hcenter
putexcel I1 = "Studies" , border (bottom top) vcenter hcenter
putexcel A2 = "Fixed part" , bold italic vcenter left

local i =3
foreach var of varlist ln_p_pfcsp_igg  ln_p_pfmsp119_igg ln_p_pfmsp2_igg ln_p_pfglurp_igg ln_p_pfse_igg  {
local var_names: variable label `var'
xtmelogit successes_all_gsg6 `var' , binomial(all_gsg6_igg_nsubpop) || id: `var' , 	 var iter(30) 
est save correct_`var' , replace
putexcel A`i' = "`var_names'"

matrix coef_mat1 = r(table)
matlist coef_mat1
matrix est_mat = coef_mat1[1, 1..1]'
matrix se_mat = coef_mat1[2, 1..1]'
matrix ll_mat = coef_mat1[5, 1..1]'
matrix ul_mat = coef_mat1[6, 1..1]'
matrix pvalue_mat = coef_mat1[4, 1..1]'

//export the results to an Excel file.
putexcel B`i' = matrix(est_mat), nformat("0.00") right vcenter 
putexcel C`i' = matrix(se_mat), nformat("(0.00)") left vcenter
putexcel D`i' = matrix(ll_mat), nformat("(0.00; (-0.00") right vcenter
putexcel E`i' = matrix(ul_mat), nformat("0.00)") left vcenter
putexcel F`i' = matrix(pvalue_mat), nformat("0.000") vcenter hcenter	

local n = e(N)
matrix n_g = e(N_g)

putexcel H`i' = `n', nformat("0") vcenter hcenter	
putexcel I`i' = matrix(n_g), nformat("0") vcenter hcenter			
				
				local ++i


} 

local i = 10
foreach var of varlist ln_p_pfama1_igg ln_p_pfmsp3_igg ln_p_pvama1_igg ln_p_pvmsp119_igg {
local var_names: variable label `var'
xtmelogit successes_all_gsg6 `var' , binomial(all_gsg6_igg_nsubpop) || id: , 	 var iter(30) 
est save correct_`var' , replace
putexcel A`i' = "`var_names'"

matrix coef_mat1 = r(table)
matlist coef_mat1
matrix est_mat = coef_mat1[1, 1..1]'
matrix se_mat = coef_mat1[2, 1..1]'
matrix ll_mat = coef_mat1[5, 1..1]'
matrix ul_mat = coef_mat1[6, 1..1]'
matrix pvalue_mat = coef_mat1[4, 1..1]'

//export the results to an Excel file.
putexcel B`i' = matrix(est_mat), nformat("0.00") right vcenter 
putexcel C`i' = matrix(se_mat), nformat("(0.00)") left vcenter
putexcel D`i' = matrix(ll_mat), nformat("(0.00; (-0.00") right vcenter
putexcel E`i' = matrix(ul_mat), nformat("0.00)") left vcenter
putexcel F`i' = matrix(pvalue_mat), nformat("0.000") vcenter hcenter				
			
local n = e(N)
matrix n_g = e(N_g)

putexcel H`i' = `n', nformat("0") vcenter hcenter	
putexcel I`i' = matrix(n_g), nformat("0") vcenter hcenter			
				local ++i


}

foreach i of numlist 1.10   {

foreach var of varlist ln_p_pfcsp_igg ln_p_pfama1_igg ln_p_pfmsp119_igg ln_p_pfmsp2_igg ln_p_pfmsp3_igg ln_p_pfglurp_igg ln_p_pfse_igg ln_p_pvama1_igg ln_p_pvmsp119_igg {
estimates use correct_`var'.ster

  nlcom `i'^_b[`var'] // CIs are not correct with this syntax

  nlcom log(`i'^_b[`var']), post // need to work with the log

  matrix b = e(b)

  matrix x = e(V)

  local se = sqrt(x[1,1])

  local b = b[1,1]

  local n = e(N)
  
* Standard error 

  disp "SE: "  "`se'"

*** 95% CIs

display "lower limit: " exp(`b' -invnormal(0.975)* `se')

display "upper limit: " exp(`b' +invnormal(0.975)* `se')

disp "Meta-N: " `n'

}
}
*/

**********************************************************************************************
*************************EBI of intercept only gsg6 IgG model*********************************
**********************************************************************************************

bysort country: egen country_seq = seq()   
bysort id: egen id_seq = seq()   



gsem (M2[country>id] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) (M1[country] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)), covstruct(_lexogenous, diagonal) latent(M2 M1 ) nocapslatent noestimate
matrix base = e(b)
gsem (M2[country>id] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)) (M1[country] -> successes_all_gsg6, family(binomial all_gsg6_igg_nsubpop) link(logit)), covstruct(_lexogenous, diagonal) latent(M2 M1 ) nocapslatent from(base)



predict gsg6_re*, latent intp(11) //  linear predictor for the random effects        
predict gsg6_xb, eta cond(fixedonly) intp(11) // linear predictor for the fixed effects
gen mu = 1 /(1+exp(-1*(gsg6_xb))) // predict prob in the average country and study
gen mu2 = 1 /(1+exp(-1*(gsg6_xb + gsg6_re1))) // predicted gsg6 prev for country
gen mu3 = 1 /(1+exp(-1*(gsg6_xb + gsg6_re2))) // predicted gsg6 prev for study

predict  ebi*, latent se(ebi_se*) intp(11)
gen pr_se = _se[_cons]
gen pred_se_c = sqrt(pr_se^2 + ebi_se1^2) //standard error for country
gen pred_se_s = sqrt(pr_se^2 + ebi_se2^2) //standard error for study

gen pred_se_prob_prob_c=mu2*(1-mu2)*pred_se_c
gen pred_se_prob_prob_s=mu3*(1-mu3)*pred_se_s


gsort + country_id - ebi1
replace ebi1 = ebi1[_n-1] if country_id == country_id[_n-1]  
replace ebi_se1 = ebi_se1[_n-1] if country_id == country_id[_n-1]  

gsort + id - ebi2
replace ebi2 = ebi2[_n-1] if id == id[_n-1]  
replace ebi_se2 = ebi_se2[_n-1] if id == id[_n-1]  

gsort country_id  - mu2
replace mu2 = mu2[_n-1] if country_id == country_id[_n-1]  
replace pred_se_prob_prob_c = pred_se_prob_prob_c[_n-1] if country_id == country_id[_n-1]  


* Country estimates
gsort + ebi1 - country_seq
generate rank_ebi1 = sum(country_seq) if country_seq ==1 & id ~= . & country_id!=3 & country_id!=6 & country_id!=9 & country_id!=15
gen labpos_ebi1 = ebi1 + 1.96*ebi_se1 + .1


serrbar ebi1 ebi_se1 rank_ebi1 if country_seq ==1 & id ~= ., mvopt(mcolor(red) ms(D) msize(tiny)) lw(thin) lcolor(black) addplot(scatter labpos_ebi1 rank_ebi1, mlabel(country) msymbol(none) mlabpos(o) mlabcol(black)) ///
scale(1.96) xtitle(Rank, size(small)) ytitle(Sg6 - random intercept, size(small)) xscale(range(1 13)) xlabel(1(1)13, labs(small)) ysize(2) ylabel(,nogrid labs(small)) ///
legend(off) graphregion(fcolor(white))  graphregion(lcolor(white))  graphregion(ilc(white))   plotregion(ilc(white)) ///
 yline(0, lp(dot) lcolor(black) ) 

 gsort + mu2 - country_seq
generate rank_pr_c = sum(country_seq) if country_seq ==1 & id ~= . & country_id!=3 & country_id!=6 & country_id!=9 & country_id!=15
gen labpos_prob_c = mu2 + 1.96*(pred_se_prob_prob_c) + .025

serrbar mu2 pred_se_prob_prob_c rank_pr_c if country_seq ==1 & id ~= ., mvopt(mcolor(red) ms(D) msiz(small)) lw(thin) lcolor(black) addplot(scatter labpos_prob_c rank_pr_c, mlabel(country) msymbol(none) mlabpos(o) mlabcol(black)) ///
scale(1.96) xtitle(Rank, size(small)) ytitle(Predicted gSG6 IgG Seroprevalence (%) - Country, size(small)) yscale(range(0 (.2) 1)) ylabel(0 (.2)1) xscale(range(1 13)) xlabel(1(1)13, labs(small)) ysize(2) ylabel(,nogrid labs(small)) ///
legend(off) graphregion(fcolor(white))  graphregion(lcolor(white))  graphregion(ilc(white))   plotregion(ilc(white)) ///
 yline(0, lp(dot) lcolor(black) ) 


  

gsort id  - mu3
replace mu3 = mu3[_n-1] if id == id[_n-1]  
replace pred_se_prob_prob_s = pred_se_prob_prob_s[_n-1] if id == id[_n-1]  

gsort + ebi2 - id_seq
generate rank_ebi2 = sum(id_seq) if id_seq ==1 & id ~= . & all_gsg6_igg_seroprevalence!=.
gen labpos_ebi2 = ebi2 + 1.96*ebi_se2 + .1


* study estimates 
serrbar ebi2 ebi_se2 rank_ebi2 if id_seq ==1 & id ~= ., mvopt(mcolor(red) ms(D) msize(tiny)) lw(thin) lcolor(black) addplot(scatter labpos_ebi2 rank_ebi2, mlabel(id) msymbol(none) mlabpos(o) mlabcol(black)) ///
scale(1.96) xtitle(Rank, size(small)) ytitle(Sg6 - random intercept, size(small)) xscale(range(1 22)) xlabel(1(1)22, labs(small)) ysize(2) ylabel(,nogrid labs(small)) ///
legend(off) graphregion(fcolor(white))  graphregion(lcolor(white))  graphregion(ilc(white))   plotregion(ilc(white)) ///
 yline(0, lp(dot) lcolor(black) )  ///
 title("Figure S1: Caterpillar plot showing study random intercept predictions for Sg6 and 95% confidence intervals versus Study ranking", size(small) position(7) color(black*255)) ///
 saving("random_intercept_predXstudy_clabel", replace) 
graph export random_intercept_predXstudy_clabel.tif ,   width(6000) height(1093) replace
  



 gsort + mu3 - id_seq
generate rank_pr_s = sum(id_seq) if id_seq ==1 & id ~= .  & all_gsg6_igg_seroprevalence!=.
gen labpos_prob_s = mu3 + 1.96*(pred_se_prob_prob_s) + .025
  
 serrbar  mu3 pred_se_prob_prob_s rank_pr_s if id_seq ==1 & id ~= ., mvopt(mcolor(red) ms(D) msize(tiny)) lw(thin) lcolor(black) ///
scale(1.96) xtitle(Rank, size(small)) ytitle(Predicted gSG6 IgG Seroprevalence (%) - Study, size(small)) yscale(range(0 (.2) 1)) ylabel(0 (.2)1) xscale(range(1 22)) xlabel(1(1)22, labs(small)) yscale(range(0 (0.2)1)) ysize(2) ylabel(,nogrid labs(small)) ///
legend(off) graphregion(fcolor(white))  graphregion(lcolor(white))  graphregion(ilc(white))   plotregion(ilc(white)) ///
 yline(0, lp(dot) lcolor(black) )     



  
  
**********************************************************************************************
*************************Other Anopheles salivary antigens************************************
**********************************************************************************************

gen  successes_fsg6_igg = fsg6_igg_seroprevalence*fsg6_igg_nsubpop
replace  successes_fsg6_igg = round(successes_fsg6_igg,1)
lab var successes_fsg6_igg "# seropositive individuals"

gen  successes_gsg6p2_igg = gsg6p2_igg_seroprevalence*gsg6p2_igg_nsubpop
replace  successes_gsg6p2_igg = round(successes_gsg6p2_igg,1)
lab var successes_gsg6p2_igg "# seropositive individuals"

*fsg6 and HBR (data from more than one study)
xtmelogit successes_fsg6_igg ln_hbr_estimate, binomial(fsg6_igg_nsubpop) ///
 || id: , 	 var 
 
 est save fsg6_hbr

 foreach i of numlist  1.50   {

estimates use fsg6_hbr.ster

  nlcom `i'^_b[ln_hbr_estimate] // CIs are not correct with this syntax

  nlcom log(`i'^_b[ln_hbr_estimate]), post // need to work with the log

  matrix b = e(b)

  matrix x = e(V)

  local se = sqrt(x[1,1])

  local b = b[1,1]

* Standard error 

  disp "SE: "  "`se'"

*** 95% CIs

display "lower limit: " exp(`b' -invnormal(0.975)* `se')

display "upper limit: " exp(`b' +invnormal(0.975)* `se')
}

*gSG6-P2 and pfcsp igg (data from more than one study)
xtmelogit successes_gsg6p2_igg ln_p_pfcsp_igg, binomial(gsg6p2_igg_nsubpop) ///
|| id:  , 	 var 
est save gsg6p2_pfcsp

foreach i of numlist  1.10   {

estimates use gsg6p2_pfcsp.ster

  nlcom `i'^_b[ln_p_pfcsp_igg] // CIs are not correct with this syntax

  nlcom log(`i'^_b[ln_p_pfcsp_igg]), post // need to work with the log

  matrix b = e(b)

  matrix x = e(V)

  local se = sqrt(x[1,1])

  local b = b[1,1]

* Standard error 

  disp "SE: "  "`se'"

*** 95% CIs

display "lower limit: " exp(`b' -invnormal(0.975)* `se')

display "upper limit: " exp(`b' +invnormal(0.975)* `se')
}

*gSG6-P2 and pfglurp igg (data from more than one study)
xtmelogit successes_gsg6p2_igg ln_p_pfglurp_igg, binomial(gsg6p2_igg_nsubpop) ///
|| id:  , 	 var 
est save gsg6p2_pfglurp

 foreach i of numlist  1.10   {

estimates use gsg6p2_pfglurp.ster

  nlcom `i'^_b[ln_p_pfglurp_igg] // CIs are not correct with this syntax

  nlcom log(`i'^_b[ln_p_pfglurp_igg]), post // need to work with the log

  matrix b = e(b)

  matrix x = e(V)

  local se = sqrt(x[1,1])

  local b = b[1,1]

* Standard error 

  disp "SE: "  "`se'"

*** 95% CIs

display "lower limit: " exp(`b' -invnormal(0.975)* `se')

display "upper limit: " exp(`b' +invnormal(0.975)* `se')
}
