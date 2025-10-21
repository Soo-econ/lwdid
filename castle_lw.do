

***************************************************************
* Lee and Wooldridge (2025)

* Title: Simple Approaches to Inference with Difference-in-Differences Estimators with Small Cross-Sectional Sample Sizes
* Replication code for Section 
* Using CASTLE.DTA data
***************************************************************

global date = string(date(c(current_date),"DMY"), "%tdCCYYNNDD")

set more off
eststo clear


use castle, clear
tab year
rename state state_str
gen state = sid


xtset state year


* FIRST_YEAR: the time of first treatment
gen first_temp = year*castle
replace first_temp = . if first_temp==0
bys state: egen first_year = min(first_temp)


drop first_temp // drop temp-var


* RYEAR: relative time
gen ryear = year-first_year
replace first_year = 0 if first_year==.

* Cohort 
egen tag = tag(state)
tab first_year if tag
* 2005 (1), 2006(13)...!

order state year ryear first_year castle lhomicide 



* Find all unique values of "state"
qui levelsof state, local(sids)
qui global sids `sids'
* display "$sids" 

* State-specific, pre-treatment demeaning each cohort:

qui forval year =2005/2009 {
    gen y`year'd = .
	* loop over each state: sids
	foreach i of global sids {
        qui reg lhomicide if state == `i' & year < `year'
        predict double y_bar, xb
        qui replace y`year'd = lhomicide - y_bar if state == `i' 
		
		* clean up temporary variable
        drop y_bar   
		}
}

* State-specific, pre-treatment detrending each cohort:

qui forval year = 2005/2009 {
		        gen y`year'dd = .
			* loop over each state: sids
	foreach i of global sids {
		// Run regression for each cid and predict 
        qui reg lhomicide year if state == `i' & year < `year' 
        predict double yhat, xb
        qui replace y`year'dd = lhomicide - yhat if state==`i' 
        drop yhat
		}
} 



**# Single treatment effect estimated following equations (7.18) and (7.19) in our manuscript.
*! ever-treated indicator
gen d = (first_year !=0)
label var d "ever-treated indicator, = 1 if eventually treated"

**###[1] post-average for treated units
forval year = 2005/2009 {
bys state:	egen double ybar_temp_`year' = mean(cond(year>=`year'& d==1,y`year'd,.))
}

forval year = 2005/2009 {
bys state:	egen double ybart_temp_`year' = mean(cond(year>=`year'& d==1,y`year'dd,.))
}


**###[2] weighted post-average for control units

*! weights in post-period w = N_g/N_treat
replace first_year=. if first_year==0
tab first_year, matcell(freqs) 
mat list freqs
mat w = freqs / r(N)

mat li w

forvalues year = 2005/2009 {
    local i = `year' - 2004
    scalar w_g = w[`i',1]
    bysort state: egen double ybar_cont_`year' = mean(cond(year >= `year' & d == 0, w_g * y`year'd, .))
}

forvalues year = 2005/2009 {
    local i = `year' - 2004
    scalar w_g = w[`i',1]
    bysort state: egen double ybart_cont_`year' = ///
        mean(cond(year >= `year' & d == 0, w_g * y`year'dd, .))
}

**### [3] ydot_bar in eq (7.18)
*! ydot_bar  (demeaning)
gen ydot_bar = .
forval year = 2005/2009 {
replace ydot_bar = ybar_temp_`year' if first_year==`year' 
}

replace ydot_bar =  (ybar_cont_2005 +ybar_cont_2006+ybar_cont_2007 + ybar_cont_2008+ybar_cont_2009) if d ==0

*! ydott_bar (detrending)
gen ydott_bar = .
forval year = 2005/2009 {
replace ydott_bar = ybart_temp_`year' if first_year==`year' 
}

replace ydott_bar =  (ybart_cont_2005 +ybart_cont_2006+ybart_cont_2007 + ybart_cont_2008+ybart_cont_2009) if d ==0

 

**# Cross-sectional regression  (for any year)
*! demeaning
reg ydot_bar d if year==2007  
reg ydot_bar d if year==2007  , vce(hc3) 

*! detrending

reg ydott_bar d if year==2007  

*! SDiD
sdid lhomicide state year castle, vce(placebo)
















