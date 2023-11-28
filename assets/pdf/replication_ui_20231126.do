/***************************************************************************** 
STATA code to reproduce part of Figure 3 in the paper "Two-Way Fixed 
Effects and Differences-in-Differences with Heterogeneous Treatment Effects: A 
Survey " by C. de Chaisemartin and X. D'Haultfoeuille (dCDH, 2023). 

by Much Rosidi, PhD student in Public Policy and Administration 
@ Martin School of Public Policy and Administration, University of Kentucky 
(https://muchrosidi.com)
******************************************************************************/

/* The original dataset originates from Wolfers (2006) titled "Did Unilateral 
Divorce Laws Raise Divorce Rates? A Reconciliation and New Results." This modified 
dataset version is adapted from the dataset provided in dCDH (2023). For instructional clarity, the dataset has been simplified to include only relevant variables. */


/* Replication Details:
   - Utilizing the provided do file from dCDH (2023) with certain modifications outlined below:
     1. Excluding always-treated states by dropping observations where lfdivlaw<1956, specifically AK and OK states.
     2. Conducting estimation for the entire dataset to ensure balanced groups for each lead post-treatment. Contrary to the paper's restriction of estimating only for year<1989.
     3. Assigning Dt_10=0 for never-treated states (identified by lfdivlaw=2000 for 20 states). This differs from the paper's coding where Dt_10=1 for these states.
     4. Although irrelevant due to exclusion, coding Dt15=0 for always-treated states.
     5. Reducing bootstrap replications to 50 instead of 200 for the dCDH method to expedite the estimation process.
*/


*Install these before running the code, if not yet. 
*ssc install csdid //net install csdid
*ssc install drdid
*ssc install did_multiplegt
*ssc install event_plot
*ssc install reghdfe 
*ssc install ftools 
*ssc install avar 
*ssc install moremata 



*Change the path after "cd" with your path where file *.dta located  

clear 
cd "D:\Bahan Kuliah\Paper to read\Two-way fixed effects and differences-in-differences with heterogeneous treatment effects a survey\Bahan_UI_final"
use "Divorce-Wolfers-AER_ui.dta", clear 



**Preparation

set matsize 800

*Transform string variable of state into numeric (alternative: egen stnum=group(st))
encode st, generate(state)

sum

*drop states with law before the time window (before 1956)
tab st lfdivlaw
drop if lfdivlaw<1956 //2 states dropped: AK and OK


gen Dur=year-lfdivlaw if lfdivlaw<2000 
replace Dur=15 if Dur>15
replace Dur=-10 if Dur<-10
replace Dur=-10 if lfdivlaw==2000

****1. TWFE 

*Generate dummy for lags (pre-treatment) 
forvalues x = 10(-1)1 {
	gen Dt_`x'=(Dur==-`x')
}

*Generate dummy for leads (post-treatment)
forvalues x = 0(1)15 {
	gen Dt`x'=(Dur==`x')
}

replace Dt_10=0 if lfdivlaw==2000

*Replace lag t-1 to be all zero (will be ommitted in the regression)
replace Dt_1=Dt_1*0




*Static Model
reg div_rate unilateral i.state i.year [w=stpop], vce(cluster state)

*Event Study (Dynamic Model)
reg div_rate Dt_10 Dt_9 Dt_8 Dt_7 Dt_6 Dt_5 Dt_4 Dt_3 Dt_2 Dt_1 Dt0 Dt1 Dt2 Dt3 Dt4 Dt5 Dt6 Dt7 Dt8 Dt9 Dt10 Dt11 Dt12 Dt13 Dt14 Dt15 i.state i.year [w=stpop], vce(cluster st)
estimates store twfe //store the regression result

event_plot twfe, graph_opt(xtitle("Relative time to change in law") ytitle("Effect") ///
	title("TWFE estimates") xlabel(-10(5)15) ylab(-1(.5)0.5) name(g_twfe, replace)) stub_lag(Dt#) stub_lead(Dt_#) together ciplottype(rcap) 
	/*plottype(scatter) lag_opt(msymbol(o) color(navy%45)) lag_ci_opt(color(maroon45))*/

	
* Joint test of null for pre-trend 
testparm Dt_*

****2. Callaway-Sant'Anna estimator 

gen cohort=lfdivlaw
replace cohort=0 if cohort==2000

csdid div_rate [weight=stpop], ivar(state) time(year)  ///
	gvar(cohort) notyet agg(event)	
	
*Note: to look like dCDH graph add long2 as an option
estat event, window(-9 15) estore(cs)

event_plot cs, /*default_look*/ graph_opt(xtitle("Relative time to change in law") ytitle("Effect") ///
	title("Callaway & Sant'Anna") xlabel(-10(5)15) ylab(-1(.5)0.5) name(g_cs, replace)) stub_lag(Tp#) stub_lead(Tm#) together	ciplottype(rcap)
	/*plottype(scatter) lag_opt(msymbol(o) color(navy%45)) lag_ci_opt(color(maroon45))*/

*Group-time ATT
csdid div_rate , ivar(state) time(year) ///
	gvar(cohort) notyet 	
	

	


* Calculate group-time ATT manually --> ATT(1970,1970)
preserve

keep if year == 1969 | year == 1970
sort state year
by state: egen missing_treatment = total(missing(div_rate))
drop if missing_treatment >= 1
gen treated = (cohort == 1970)
gen control = (cohort > 1970 | cohort == 0)

/* Try calculate it manually (mynote: do all these lines at one, if not, it will be failed) */
sum div_rate if year == 1969 & control == 1
loc div_rate_n1 = r(mean)
sum div_rate if year == 1970 & control == 1
loc div_rate_n2 = r(mean)
sum div_rate if year == 1969 & treated == 1
loc div_rate_s1 = r(mean)
sum div_rate if year == 1970 & treated == 1
loc div_rate_s2 = r(mean)

loc att_1970_1970 = (`div_rate_s2' - `div_rate_s1') - (`div_rate_n2' - `div_rate_n1')

disp "ATT(1970, 1970): `att_1970_1970'"

restore



* Calculate group-time ATT manually --> ATT(1971,1974)
preserve

keep if year == 1970 | year == 1974
sort state year
by state: egen missing_treatment = total(missing(div_rate))
drop if missing_treatment >= 1
gen treated = (cohort == 1971)
gen control = (cohort > 1974 | cohort == 0)

/* Try calculate it manually (mynote: do all these lines at one, if not, it will be failed) */
sum div_rate if year == 1970 & control == 1
loc div_rate_n1 = r(mean)
sum div_rate if year == 1974 & control == 1
loc div_rate_n2 = r(mean)
sum div_rate if year == 1970 & treated == 1
loc div_rate_s1 = r(mean)
sum div_rate if year == 1974 & treated == 1
loc div_rate_s2 = r(mean)

loc att_1971_1974 = (`div_rate_s2' - `div_rate_s1') - (`div_rate_n2' - `div_rate_n1')

disp "ATT(1971, 1974): `att_1971_1974'"

restore
	
	
****3. dC-DH estimator

did_multiplegt div_rate state year unilateral, average_effect ///
	robust_dynamic dynamic(15) placebo(9) ///
	cluster(state) breps(50) weight(stpop) seed(1) ///
	graphoptions(ylabel(-1(.5)0.5) yscale(range(-1.1 0.8)) legend(off) ///
	xtitle(Relative time to change in law) title(dC&DH w/o lin. trends) /// 
	ytitle(Effect) name(g_dCDH, replace))
	
		
* Joint test of null for pre-trend 
display e(p_jointplacebo)

* Combine all graphs (Figure 3)
graph combine g_twfe g_cs g_dCDH
graph export graphs_combine.pdf, replace

graph drop g*


