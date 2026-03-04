*! lwdid: Lee & Wooldridge (2025, working, https://dx.doi.org/10.2139/ssrn.5325686) — rolling DID estimator
*! version 1.0 13oct2025
*! authors: Soo Jeong Lee, Jeffrey M. Wooldridge
*! contact: soojeong.lee@siu.edu, wooldri1@msu.edu

*! This command is currently designed for the ****common timing case with small N *****. 
*! A forthcoming version will extend the method to accommodate staggered treatment adoption and {it:large-N} panels (see, e.g., Lee and Wooldridge (2023))
**********************************************************************************************************************
*! SYNTAX
*!     lwdid y d, ivar(id) tvar(year [quarter]) post(post)
*!           rolling({demean|detrend|demeanq|detrendq})
*!           [vce(vartype) table("filename") gid(string) graph gopts(string)
*!            controls(varlist) save(string)
*!            ri rireps(integer) riseed(integer) ]

*! ARGUMENTS
*!   y            outcome variable
*!   d            treatment indicator (non-zero treated, zero otherwise)

*! REQUIRED OPTIONS
*!   ivar(id)     panel identifier (numeric or string; strings are internally grouped)
*!   tvar()       time variable(s):
*!                  • yearly data: tvar(year)
*!                  • year–quarter: tvar(year quarter)  (quarter must be 1–4 or labeled)
*!   post(var)    post-period indicator (1 = post, 0 = pre)
*!   rolling()    type of transformation applied to residualize y:
*!                  - demean    : unit-specific mean removal over pre periods
*!                  - detrend   : unit-specific linear time trend removal
*!                  - demeanq   : as demean, plus quarter fixed effects (requires tvar(year quarter))
*!                  - detrendq  : as detrend, plus quarter fixed effects (requires tvar(year quarter))

*! OPTIONAL [... ]
*!   vce(vartype)   HC3 variance estimator for regressions (e.g., hc3.)
*!   TABLE(string)   Save results to files in the current working directory:
*!                   • Single-effect (overall post average) → Excel file via outreg2
*!                       - Always saved as "name.xls" (same stem as TABLE())
*!                   • Period-by-period ATT estimates → CSV file
*!                       - If TABLE("name") saved as "`name'_byperiod.csv"

*!   graph          Plot the residualized outcome: control-group mean vs. treated units over time
*!   gopts(string)  Additional twoway graph options passed directly to the plot command
*!                  (effective only when the graph option is specified).
*!                  For example:
*!                      gopts(ytitle("Residualized outcome") xtitle("Year") legend(pos(1) ring(0)) title("Control vs Treated"))

*! controls(varlist)  List of additional covariates to include in the regressions
*!                   (e.g., controls(x1 x2 x3)). These variables are used as adjustment
*!                   covariates when estimating ATT_t. As discussed in the paper,
*!                   controls are included in the cross-sectional regressions only when
*!                   N_1 > K+1 and N_0 > K+1, where K is the number of control variables,
*!                   and N_1 (N_0) is the number of treated (control) units.


*!   gid(id)   Choose what to plot on the treated side in the graph.
*!                    - If omitted: the command plots the treated-group average
*!                      of the residualized outcome over time.
*!                    - If provided: plots the series for that specific treated unit.
*!                      Works with numeric or string IDs, e.g.:
*!                          gid(101)   or   gid("CA")
*!                    - The value must correspond to a unit with d==1 in the data;
*!                      otherwise an error is thrown.

*!   ri               Performs a manual randomization inference (RI) test (instead of `ritest`)
*!                    for the estimated ATT and reports the simulated p-value.
*!   rireps(integer)  Number of randomization repetitions (default = 1000).
*!   riseed(integer)  Specifies the random seed used for reproducibility in the RI procedure.
*!                    Although technically parsed as a string, the user should enter a
*!                    numeric value (e.g., riseed(12345)). When omitted,  a random seed is used.

*! NOTES
*!   • For demeanq/detrendq you must supply tvar(year quarter).
*!   • Non-numeric ids are automatically converted to a numeric group.
*!   • The treatment indicator and post variable are coerced to {0,1}.


set more off
set trace off
capture program drop lwdid
program define lwdid, eclass sortpreserve

    *------------------------------------------------------------------------------*
    * (0) Parse syntax & basic checks
    *------------------------------------------------------------------------------*
    /* varlist must be: outcome treatment_indicator (2 variables) */
    syntax varlist(min=2 max=2 numeric) [if] [in] , ///
        IVAR(varname) ///                
        TVAR(varlist min=1 max=2) ///     
        POST(varname) ///
        ROLLING(name) ///
       [ GID(string) ///
  CONTROLS(varlist) ///
  GRAPH ///
  GOPTS(string asis) ///
  TABLE(string) ///
  VCE(string) ///
  RI ///
  RISEED(string) ///
  RIREPS(integer 1000) ]
  
* handle RI options (seed & reps)
if missing(`rireps') local rireps 1000
* convert riseed() if provided
capture confirm number `riseed'
if _rc local riseed = .
if missing(`riseed') local riseed = ceil(runiform() * 1e6 + 1000*runiform())




    marksample touse, novarlist
    gettoken y d : varlist

	* normalize option names into locals used below
	local id `ivar'
	local t  `tvar'

    *-- rolling() check
    local rolling = lower("`rolling'")
    local ok = inlist("`rolling'","demean","detrend","demeanq","detrendq")
    if !`ok' {
        di as err "rolling() must be one of: demean, detrend, demeanq, detrendq"
        exit 198
    }

    *-- basic confirms
    confirm variable `y' `d' `id' `post'
	local id_orig `ivar' 
    tokenize "`t'"

	*---- ensure ID is numeric (convert string IDs to grouped numeric) ----*
	capture confirm numeric variable `ivar'
	if _rc {   // string id
		tempvar id_num
		egen `id_num' = group(`ivar')
		label var `id_num' "`ivar' (grouped from string)"
		local id `id_num'
	}
	else {
		local id `ivar'
	}

    local t1 `1'
    local t2 `2'
    confirm variable `t1'
    if "`t2'"!="" confirm variable `t2'

    *------------------------------------------------------------------------------*
    * (1) Keep analysis sample & create core flags
    *------------------------------------------------------------------------------*
    tempvar post_ tindex tq yhat ydot cs ctrlSum ydot_tr
	
    tempname b V
	

    preserve
        qui keep if `touse'
		
		qui capture drop d_
		qui capture drop firstyear 

        qui gen byte d_    = (`d'   != 0) if !missing(`d')
        qui gen byte `post_' = (`post'!= 0) if !missing(`post')

        qui drop if missing(`y', `id', `post_', d_)

    *------------------------------------------------------------------------------*
    * (2) Build time index: yearly or year-quarter, plus nice tq label
    *------------------------------------------------------------------------------*
        capture drop `tq'
        if "`t2'"=="" {
            su `t1', meanonly
            qui gen long `tindex' = `t1' - r(min) + 1
            label var `tindex' "time index (year-based; min->1)"
        }
        else {
            qui gen double `tq' = yq(`t1', `t2')
            format `tq' %tq
            su `tq', meanonly
            qui gen long `tindex' = `tq' - r(min) + 1
            label var `tindex' "time index (year-quarter; min->1)"
        }

    *------------------------------------------------------------------------------*
    * (3) Identify K (max pre tindex) and first post period (tpost1)
    *------------------------------------------------------------------------------*
        qui su `tindex' if `post_'==0, meanonly
        if r(N)==0 {
            di as err "No pre-treatment observations (post==0)."
            exit 2000
        }
        local K = r(max)

        qui su `tindex' if `post_'==1, meanonly
        if r(N)==0 {
            di as err "No post-treatment observations (post==1)."
            exit 2000
        }
        local tpost1 = r(min)

    *------------------------------------------------------------------------------*
    * (4) Fit pre-period model by ID; predict yhat for ALL periods
    *------------------------------------------------------------------------------*
        qui gen double `yhat' = .
        qui levelsof `id', local(IDlist)

        quietly foreach ii of local IDlist {
            if "`rolling'"=="demean" {
                regress `y' if `id'==`ii' & `post_'==0
                predict double __fit if `id'==`ii', xb
            }
            else if "`rolling'"=="detrend" {
                regress `y' c.`tindex' if `id'==`ii' & `post_'==0
                predict double __fit if `id'==`ii', xb
            }
            else if "`rolling'"=="demeanq" {
                if "`t2'"=="" {
                    di as err "rolling(demeanq) requires tvar(year quarter)."
                    exit 198
                }
                regress `y' i.`t2' if `id'==`ii' & `post_'==0
                predict double __fit if `id'==`ii', xb
            }
            else if "`rolling'"=="detrendq" {
                if "`t2'"=="" {
                    di as err "rolling(detrendq) requires tvar(year quarter)."
                    exit 198
                }
                regress `y' c.`tindex' i.`t2' if `id'==`ii' & `post_'==0
                predict double __fit if `id'==`ii', xb
            }
            replace `yhat' = __fit if `id'==`ii'
            cap drop __fit
        }


	*------------------------------------------------------------------------------*
	* (5) Residualized outcome for ALL periods: ydot = y - yhat
	*------------------------------------------------------------------------------*

* Always also create internal ydot for the subsequent analysis
    quietly gen double `ydot' = `y' - `yhat'
    qui label var `ydot' "y - yhat (residualized via pre-period fit)"

*------------------------------------------------------------------------------*
* (6) Collapse to first post cross-section & overall ATT
*------------------------------------------------------------------------------*

quietly bysort `id': egen double ydot_postavg = mean(cond(`post_'==1, `ydot', .))
qui label var ydot_postavg "mean(ydot) over post periods, by id"

quietly gen byte firstpost = (`tindex'==`tpost1') & !missing(ydot_postavg, d_)

* ---- build Excel filename for single-effect table: <table>.xls ----
local __do_export = ("`table'" != "")
if `__do_export' {
    * If user passed a name with any extension, strip it and add .xls
    mata: st_local("__suf", pathsuffix(st_local("table")))
    if inlist(strlower("`__suf'"), ".xls", ".xlsx", ".csv") {
        mata: st_local("__stem", pathrmsuffix(st_local("table")))
        local __xls "`__stem'.xls"
    }
    else {
        local __xls "`table'.xls"
    }
}


* Check outreg2 availability (Excel export depends on it)
local __has_or2 = 1
capture which outreg2
if _rc {
    di as txt ">>> outreg2 not found; attempting SSC install..."
    capture noisily ssc install outreg2, replace
    capture which outreg2
    if _rc {
        di as err ">>> Auto-install failed. Excel export skipped. Manually run:  ssc install outreg2"
        local __has_or2 = 0
    }
}



* --- build centered X and interactions ONCE (time-invariant X) ---
local RHS
local INTERACTS


local K1 0
if "`controls'" != "" {
    local K1 : word count `controls'
}
local nk = `K1' + 1

quietly count if `tindex'==1 & d_==1
local nt = r(N)
quietly count if `tindex'==1 & d_==0
local nc = r(N)

 
* only if controls specified AND hold the condition.
if "`controls'" != "" {
    if `nt' > `nk' & `nc' > `nk' {
        foreach x of varlist `controls' {
            tempvar xc Dx
            * time-invariant X: center once using treated mean (any consistent rule)
            quietly summarize `x' if d_==1, meanonly
            quietly gen double `xc' = `x' - r(mean)
            quietly gen double `Dx' = d_ * `xc'
            local INTERACTS `INTERACTS' `Dx'
        }
        local RHS "`controls' `INTERACTS'"
    }
    else {
        local RHS ""
        di as txt "------------------------------------------------------------"
        di as txt "Controls not applied: sample does not satisfy N_1 > K+1 and N_0 > K+1"
        di as txt "------------------------------------------------------------"
    }
}
else {
    local RHS ""
}

*------------------------------------------------------------------------------*
* (Main regression and optional randomization inference)
*------------------------------------------------------------------------------*

if "`ri'" != "" {
di as txt "--------------------------------------------------------------"
di as txt " Running manually randomization inference (instead of RITEST) "
di as txt " to ensure fully reproducible RI p-values"
di as txt "--------------------------------------------------------------"


    * --- reproducibility setup ---
    quietly set rng mt64
    quietly set seed `riseed'
    local reps = `rireps'

    * --- save a clean copy of the current dataset (avoid preserve conflict) ---
    tempfile base
    quietly save `base', replace

    * --- baseline regression (observed effect) ---
    quietly regress ydot_postavg d_ `RHS' if firstpost
    scalar __b0 = _b[d_]

    * --- permutation loop ---
    tempname res
    mat `res' = J(`reps', 1, .)

    forvalues r = 1/`reps' {
        quietly use `base', clear
        gen double d_perm = d_[runiformint(1, _N)]   // permute treatment indicator
        quietly regress ydot_postavg d_perm `RHS' if firstpost
        mat `res'[`r', 1] = _b[d_perm]
    }

    * --- compute p-value ---
    mata: st_matrix("abs_res", abs(st_matrix("`res'")))
    mata: st_numscalar("__p_ri", mean(st_matrix("abs_res") :>= abs(st_numscalar("__b0"))))
    ereturn scalar p_ri = __p_ri
}

* --- FYI: OLS (unadjusted) ---
di as txt "------------------------------------------------------------"
di as txt "  (FYI: OLS regression for the single average effect, without variance adjustment) "
di as txt "------------------------------------------------------------"
regress ydot_postavg d_ `RHS' if firstpost

if `__do_export' & `__has_or2' {
    outreg2 using "`__xls'", replace label nose ///
        ctitle("`rolling' (OLS)") addstat(N, e(N))
}

* --- Main regression with/without vce() ---
if "`vce'" != "" {
    di as txt "============================================================"
    di as txt "  Main Results: regression for the single average effect, with `vce' adjustment"
    di as txt "============================================================"
    regress ydot_postavg d_ `RHS' if firstpost, vce(`vce')

    if `__do_export' & `__has_or2' {
        outreg2 using "`__xls'", append label nose ///
            ctitle("`rolling' (OLS) with `vce'") addstat(N, e(N))
    }
}
else {
    di as txt "----------------------------------------------------------------------------"
    di as txt "  OLS regression for the single average effect, without variance adjustment "
    di as txt "----------------------------------------------------------------------------"
    regress ydot_postavg d_ `RHS' if firstpost

    if `__do_export' & `__has_or2' {
        outreg2 using "`__xls'", replace label nose ///
            ctitle("`rolling' (OLS)") addstat(N, e(N))
    }
}

* keep these for e()
matrix `b' = e(b)
matrix `V' = e(V)
scalar __att_overall = _b[d_]
scalar __se_overall  = _se[d_]

*------------------------------------------------------------------------------*
* (6b) Period-by-period post effects table
*------------------------------------------------------------------------------*
tempname __pph
tempfile __ppfile
postfile `__pph' int tindex str20 period double beta se tstat pval N using `__ppfile', replace

* first row = overall 'average' from the single-effect above
local b_avg = _b[d_]
local s_avg = _se[d_]
local t_avg = `b_avg'/`s_avg'
local p_avg = 2*ttail(e(df_r), abs(`t_avg'))
local N_avg = e(N)
post `__pph' (.) ("average") (`b_avg') (`s_avg') (`t_avg') (`p_avg') (`N_avg')

quietly su `tindex', meanonly
local Tmax = r(max)

forvalues tt = `tpost1'/`Tmax' {
    if ("`t2'"=="") {
        quietly su `t1' if `tindex'==`tt', meanonly
        local lab : display %9.0f r(mean)
    }
    else {
        quietly su `tq' if `tindex'==`tt', meanonly
        local lab : display %tq r(mean)
    }

    capture {
        if "`vce'" != "" {
          qui  regress `ydot' d_ `RHS' if `tindex'==`tt', vce(`vce')
        }
        else {
          qui  regress `ydot' d_ `RHS' if `tindex'==`tt'
        }
    }
    if _rc==0 {
        local coef = _b[d_]
        local se1  = _se[d_]
        local t    = `coef'/`se1'
        local p    = 2*ttail(e(df_r), abs(`t'))
        post `__pph' (`tt') ("`lab'") (`coef') (`se1') (`t') (`p') (e(N))
    }
    else {
        post `__pph' (`tt') ("`lab'") (.) (.) (.) (.) (.)
    }
}
postclose `__pph'

* frame-safe load + pretty print + export
local __curframe = c(frame)
if "`__curframe'" == "__lwdid_pp" frame change default
capture confirm frame __lwdid_pp
if !_rc {
    capture frame drop __lwdid_pp
}
frame create __lwdid_pp
frame change  __lwdid_pp

capture confirm file "`__ppfile'"
if _rc {
    di as err "Internal: period-by-period temp file not found."
    frame change default
}
else {
    use "`__ppfile'", clear

    qui gen byte __isavg = period=="average"
    gsort -__isavg tindex
    qui gen str10 tindex_s = cond(__isavg, "-", string(tindex))
    qui drop tindex
    rename tindex_s tindex
    qui drop __isavg

    label var tindex "time index"
    label var period "period"
    label var beta   "ATT_t"
    label var se     "SE(ATT_t)"
    label var tstat  "t"
    label var pval   "p-value"
    qui gen double ci_lw = beta - 1.96*se
    qui gen double ci_up = beta + 1.96*se
    order period tindex beta se ci_lw ci_up tstat pval N

    di as txt "=== Period-by-period post-treatment effects ==="
    list period tindex beta se ci_lw ci_up tstat pval N, noobs abbrev(20)

    * robust export:
    if `__do_export' {
        if inlist("`__ext'", ".xls", ".xlsx") {
            * keep single-effect in Excel (via outreg2 above),
            * and save the by-period table alongside as CSV
            local __csvbp = "`table'_byperiod.csv"
            capture noisily export delimited using "`__csvbp'", replace
            if _rc {
                di as err ">>> export (by-period CSV) failed: rc = " _rc
            }
            else {
                di as txt ">>> Period-by-period results exported to " as res "`__csvbp'"
            }
        }
        else if "`__ext'"==".csv" {
            * user explicitly wants CSV → save the by-period table to that file
            capture noisily export delimited using "`__fname'", replace
            if _rc {
                di as err ">>> export (CSV) failed: rc = " _rc
            }
            else {
                di as txt ">>> Period-by-period results exported to " as res "`__fname'"
            }
        }
        else {
            * unknown extension → default to CSV with suffix
            local __csvbp = "`table'_byperiod.csv"
            capture noisily export delimited using "`__csvbp'", replace
            if _rc {
                di as err ">>> export (by-period CSV) failed: rc = " _rc
            }
            else {
                di as txt ">>> Period-by-period results exported to " as res "`__csvbp'"
            }
        }
    }

    quietly count
    local __Nrow = r(N)
    tempname __ATTt __SEt
    if (`__Nrow' > 0) {
        capture noisily mkmat beta, matrix(`__ATTt')
        capture noisily mkmat se,   matrix(`__SEt')
    }

    frame change default
    capture confirm matrix `__ATTt'
    if !_rc {
        matrix ATTt = `__ATTt'
        capture ereturn matrix ATTt = ATTt
    }
    capture confirm matrix `__SEt'
    if !_rc {
        matrix SEt  = `__SEt'
        capture ereturn matrix SEt  = SEt
    }

    capture frame drop __lwdid_pp
}

    *------------------------------------------------------------------------------*
    * (7) Control mean trajectory & optional graph
    *------------------------------------------------------------------------------*
* Always compute control-group mean over time (d==0)
        qui gen double `ctrlSum' = `ydot' if d_==0
        bysort `tindex': egen double y_dot_cont = mean(`ctrlSum')
        label var y_dot_cont "Control mean of ydot (d==0)"

* Always compute treated-group mean over time (d==1)
tempvar trSum
qui gen double `trSum' = `ydot' if d_==1
bysort `tindex': egen double y_dot_tr = mean(`trSum')
label var y_dot_tr "Treated mean of ydot (d==1)"
		
*---------------- Ensure there is at least one treated obs ----------------*
quietly count if d_==1
if r(N)==0 {
    di as err "No treated units found (d==1)."
    exit 2000
}

* In your syntax: GID(string)
local gid = trim("`gid'")
if "`gid'" != "" {
    * Compare using the ORIGINAL ivar(), not the grouped numeric `id'
    * because the user passed the id value in the ivar's domain.
    capture confirm numeric variable `id_orig'
    if !_rc {
        * ivar() is numeric → gid() must be numeric (or numeric-like string)
        local __gidnum = real("`gid'")
        if missing(`__gidnum') {
            di as err "gid(`gid') must be numeric because ivar() is numeric."
            exit 198
        }
        quietly count if `id_orig'==`__gidnum' & d_==1
        if r(N)==0 {
            di as err "gid(`gid') not found among treated units (d==1)."
            exit 2000
        }
        capture drop `ydot_tr'
        gen double `ydot_tr' = `ydot' if `id_orig'==`__gidnum'
    }
    else {
        * ivar() is string → use gid() literally (must be quoted if it has letters)
        quietly count if `id_orig'=="`gid'" & d_==1
        if r(N)==0 {
            di as err "gid(`gid') not found among treated units (d==1)."
            exit 2000
        }
        capture drop `ydot_tr'
        gen double `ydot_tr' = `ydot' if `id_orig'=="`gid'"
    }
    label var `ydot_tr' "Treated unit (id=`gid')"
}


    *------------------------------------------------------------------------------*
    * (7-1) Graph — Specific treated unit specified by gid()
	*------------------------------------------------------------------------------*
	
	if "`gid'" != "" {
	*====================== YEAR-ONLY PLOT======================*
        if "`t2'"=="" {
		tempvar xyear
		qui gen double `xyear' = `t1'
		label var `xyear' "year"

		quietly su `t1' if `post_'==1, meanonly
		local xpost = r(min)

		if "`graph'" != "" {
	      *---- y-axis ticks: nice step (>=0.1), aligned bounds, hard-limit tick count ----*
			local y_max_ticks = 6      // target max number of labeled ticks (~5–6)
			local y_min_step  = 0.1    // never show spacing smaller than 0.1

			qui su y_dot_cont, meanonly
			local ymin = r(min)
			local ymax = r(max)
			qui su `ydot_tr', meanonly
			local ymin = min(`ymin', r(min))
			local ymax = max(`ymax', r(max))

			local yrange = `ymax' - `ymin'

			if (`yrange' <= 0) {
				local ystep = `y_min_step'
				local ymin2 = `ymin' - 5*`ystep'
				local ymax2 = `ymax' + 5*`ystep'
			}
			else {
				* 1) raw step from range / target count, enforce minimum
				local raw = max(`yrange'/`y_max_ticks', `y_min_step')

				* 2) round UP to nice sequence {1,2,5} * 10^k   (true "ceil" to nice)
				local k     = floor(log10(`raw'))
				local base  = 10^`k'
				local frac  = `raw'/`base'
				local nice  = cond(`frac'<=1, 1, cond(`frac'<=2, 2, cond(`frac'<=5, 5, 10)))
				local ystep = max(`nice'*`base', `y_min_step')

				* 3) align bounds to multiples of step
				local ymin2 = `ystep'*floor(`ymin'/`ystep')
				local ymax2 = `ystep'*ceil(`ymax'/`ystep')

				* 4) if we still have too many ticks, bump step up until <= y_max_ticks
				local nticks = floor((`ymax2' - `ymin2')/`ystep') + 1
				local guard = 0
				while (`nticks' > `y_max_ticks' & `guard' < 10) {
					* multiply step by 2 (you could use 2.5 or 5 if you prefer larger jumps)
					local ystep = `ystep' * 2
					local ymin2 = `ystep'*floor(`ymin'/`ystep')
					local ymax2 = `ystep'*ceil(`ymax'/`ystep')
					local nticks = floor((`ymax2' - `ymin2')/`ystep') + 1
					local guard = `guard' + 1
				}
			}

			* label format: integers if step>=1, else one decimal
			local yfmt = cond(`ystep' >= 1, "%9.0f", "%9.1f")

			* build ylabel option
			local YL "ylabel(`ymin2'(`ystep')`ymax2', format(`yfmt'))"


    *------------------------------- Plot ----------------------------------*
    twoway ///
        (line y_dot_cont `xyear', lpattern(dash)  lcolor(stc1)) ///
        (line `ydot_tr'  `xyear', lpattern(solid) lcolor(stc2)) ///
        , xline(`xpost', lcolor(gray) lpattern(dash)) ///
          `YL' ///
          legend(order(1 "Control" 2 "Treated (ID=`gid')") pos(12) col(2)) ///
          xtitle("year") ytitle("") title("") `gopts'
            }
        }
		
*====================== YEAR-QUARTER PLOT ======================*
else {
    if "`graph'" != "" {

        *---------------- Build quarterly x ----------------*
        tempvar xdate tq_
        qui gen double `tq_' = yq(`t1', `t2')
        format `tq_' %tq
        qui gen double `xdate' = `tq_'
        quietly su `tq_' if `post_'==1, meanonly
        local xpostline = r(min)
        local xfmt %tq

        *---------------- Y-axis: nice ticks (>=0.1), <= ~6 labels ----------------*
        local y_max_ticks = 6
        local y_min_step  = 0.1

        quietly su y_dot_cont, meanonly
        local ymin = r(min)
        local ymax = r(max)
        quietly su `ydot_tr', meanonly
        local ymin = min(`ymin', r(min))
        local ymax = max(`ymax', r(max))
        local yrange = `ymax' - `ymin'

        if (`yrange' <= 0) {
            local ystep = `y_min_step'
            local ymin2 = `ymin' - 5*`ystep'
            local ymax2 = `ymax' + 5*`ystep'
        }
        else {
            local raw   = max(`yrange'/`y_max_ticks', `y_min_step')
            local k     = floor(log10(`raw'))
            local base  = 10^`k'
            local frac  = `raw'/`base'
            local nice  = cond(`frac'<=1,1, cond(`frac'<=2,2, cond(`frac'<=5,5,10)))
            local ystep = max(`nice'*`base', `y_min_step')

            local ymin2 = `ystep'*floor(`ymin'/`ystep')
            local ymax2 = `ystep'*ceil(`ymax'/`ystep')

            local nticks = floor((`ymax2' - `ymin2')/`ystep') + 1
            while (`nticks' > `y_max_ticks') {
                local ystep  = `ystep' * 2
                local ymin2  = `ystep'*floor(`ymin'/`ystep')
                local ymax2  = `ystep'*ceil(`ymax'/`ystep')
                local nticks = floor((`ymax2' - `ymin2')/`ystep') + 1
            }
        }
        local yfmt = cond(`ystep'>=1, "%9.0f", "%9.1f")
        local YL "ylabel(`ymin2'(`ystep')`ymax2', format(`yfmt'))"

        *---------------- X-axis: quarterly ticks, end at last obs ----------------*
        quietly su `xdate', meanonly
        local xmin = floor(r(min))
        local xmax = floor(r(max))

        local rawx     = (`xmax' - `xmin')/10
        local tickstep = cond(`rawx' < 4, 4, 4*ceil(`rawx'/4))
        local firsttick = `tickstep' * ceil(`xmin'/`tickstep')

        local labs
        local lasttick = .
        if (`firsttick' <= `xmax') {
            forvalues t = `firsttick'(`tickstep')`xmax' {
                local labs "`labs' `t'"
                local lasttick = `t'
            }
        }
        if missing(`lasttick') | (`lasttick' < `xmax') local labs "`labs' `xmax'"
        local XL "xlabel(`labs', format(`xfmt') angle(vertical) labsize(small)) xscale(range(`xmin' `xmax'))"

        *---------------- Plot ----------------*
        twoway ///
            (line y_dot_cont `xdate', lpattern(dash)  lcolor(stc1)) ///
            (line `ydot_tr'  `xdate', lpattern(solid) lcolor(stc2)) ///
            , xline(`xpostline', lcolor(gray) lpattern(dash)) ///
              `XL' `YL' ///
              legend(order(1 "Control" 2 "Treated (ID=`gid')") pos(12) col(2)) ///
              ytitle("") xtitle("") title("") `gopts' ///
              plotregion(margin(b=14))
    }
}
}

    *------------------------------------------------------------------------------*
    * (7-2) Graph — Average treated effect across all treated units
	*------------------------------------------------------------------------------*
	

else{
*====================== YEAR-ONLY PLOT======================*
        if "`t2'"=="" {
		tempvar xyear
		qui gen double `xyear' = `t1'
		label var `xyear' "year"

		quietly su `t1' if `post_'==1, meanonly
		local xpost = r(min)

		if "`graph'" != "" {
	      *---- y-axis ticks: nice step (>=0.1), aligned bounds, hard-limit tick count ----*
			local y_max_ticks = 6      // target max number of labeled ticks (~5–6)
			local y_min_step  = 0.1    // never show spacing smaller than 0.1

			qui su y_dot_cont, meanonly
			local ymin = r(min)
			local ymax = r(max)
			qui su y_dot_tr, meanonly
			local ymin = min(`ymin', r(min))
			local ymax = max(`ymax', r(max))

			local yrange = `ymax' - `ymin'

			if (`yrange' <= 0) {
				local ystep = `y_min_step'
				local ymin2 = `ymin' - 5*`ystep'
				local ymax2 = `ymax' + 5*`ystep'
			}
			else {
				* 1) raw step from range / target count, enforce minimum
				local raw = max(`yrange'/`y_max_ticks', `y_min_step')

				* 2) round UP to nice sequence {1,2,5} * 10^k   (true "ceil" to nice)
				local k     = floor(log10(`raw'))
				local base  = 10^`k'
				local frac  = `raw'/`base'
				local nice  = cond(`frac'<=1, 1, cond(`frac'<=2, 2, cond(`frac'<=5, 5, 10)))
				local ystep = max(`nice'*`base', `y_min_step')

				* 3) align bounds to multiples of step
				local ymin2 = `ystep'*floor(`ymin'/`ystep')
				local ymax2 = `ystep'*ceil(`ymax'/`ystep')

				* 4) if we still have too many ticks, bump step up until <= y_max_ticks
				local nticks = floor((`ymax2' - `ymin2')/`ystep') + 1
				local guard = 0
				while (`nticks' > `y_max_ticks' & `guard' < 10) {
					* multiply step by 2 (you could use 2.5 or 5 if you prefer larger jumps)
					local ystep = `ystep' * 2
					local ymin2 = `ystep'*floor(`ymin'/`ystep')
					local ymax2 = `ystep'*ceil(`ymax'/`ystep')
					local nticks = floor((`ymax2' - `ymin2')/`ystep') + 1
					local guard = `guard' + 1
				}
			}

			* label format: integers if step>=1, else one decimal
			local yfmt = cond(`ystep' >= 1, "%9.0f", "%9.1f")

			* build ylabel option
			local YL "ylabel(`ymin2'(`ystep')`ymax2', format(`yfmt'))"


    *------------------------------- Plot ----------------------------------*
    twoway ///
        (line y_dot_cont `xyear', lpattern(dash)  lcolor(stc1)) ///
        (line y_dot_tr  `xyear', lpattern(solid) lcolor(stc2)) ///
        , xline(`xpost', lcolor(gray) lpattern(dash)) ///
          `YL' ///
          legend(order(1 "Control" 2 "Treated") pos(12) col(2)) ///
          xtitle("year") ytitle("") title("") `gopts'
            }
        }
		
*====================== YEAR-QUARTER PLOT ======================*
else {
    if "`graph'" != "" {

        *---------------- Build quarterly x ----------------*
        tempvar xdate tq_
        qui gen double `tq_' = yq(`t1', `t2')
        format `tq_' %tq
        qui gen double `xdate' = `tq_'
        quietly su `tq_' if `post_'==1, meanonly
        local xpostline = r(min)
        local xfmt %tq

        *---------------- Y-axis: nice ticks (>=0.1), <= ~6 labels ----------------*
        local y_max_ticks = 6
        local y_min_step  = 0.1

        quietly su y_dot_cont, meanonly
        local ymin = r(min)
        local ymax = r(max)
        quietly su y_dot_tr, meanonly
        local ymin = min(`ymin', r(min))
        local ymax = max(`ymax', r(max))
        local yrange = `ymax' - `ymin'

        if (`yrange' <= 0) {
            local ystep = `y_min_step'
            local ymin2 = `ymin' - 5*`ystep'
            local ymax2 = `ymax' + 5*`ystep'
        }
        else {
            local raw   = max(`yrange'/`y_max_ticks', `y_min_step')
            local k     = floor(log10(`raw'))
            local base  = 10^`k'
            local frac  = `raw'/`base'
            local nice  = cond(`frac'<=1,1, cond(`frac'<=2,2, cond(`frac'<=5,5,10)))
            local ystep = max(`nice'*`base', `y_min_step')

            local ymin2 = `ystep'*floor(`ymin'/`ystep')
            local ymax2 = `ystep'*ceil(`ymax'/`ystep')

            local nticks = floor((`ymax2' - `ymin2')/`ystep') + 1
            while (`nticks' > `y_max_ticks') {
                local ystep  = `ystep' * 2
                local ymin2  = `ystep'*floor(`ymin'/`ystep')
                local ymax2  = `ystep'*ceil(`ymax'/`ystep')
                local nticks = floor((`ymax2' - `ymin2')/`ystep') + 1
            }
        }
        local yfmt = cond(`ystep'>=1, "%9.0f", "%9.1f")
        local YL "ylabel(`ymin2'(`ystep')`ymax2', format(`yfmt'))"

        *---------------- X-axis: quarterly ticks, end at last obs ----------------*
        quietly su `xdate', meanonly
        local xmin = floor(r(min))
        local xmax = floor(r(max))

        local rawx     = (`xmax' - `xmin')/10
        local tickstep = cond(`rawx' < 4, 4, 4*ceil(`rawx'/4))
        local firsttick = `tickstep' * ceil(`xmin'/`tickstep')

        local labs
        local lasttick = .
        if (`firsttick' <= `xmax') {
            forvalues t = `firsttick'(`tickstep')`xmax' {
                local labs "`labs' `t'"
                local lasttick = `t'
            }
        }
        if missing(`lasttick') | (`lasttick' < `xmax') local labs "`labs' `xmax'"
        local XL "xlabel(`labs', format(`xfmt') angle(vertical) labsize(small)) xscale(range(`xmin' `xmax'))"

        *---------------- Plot ----------------*
        twoway ///
            (line y_dot_cont `xdate', lpattern(dash)  lcolor(stc1)) ///
            (line y_dot_tr  `xdate', lpattern(solid) lcolor(stc2)) ///
            , xline(`xpostline', lcolor(gray) lpattern(dash)) ///
              `XL' `YL' ///
              legend(order(1 "Control" 2 "Treated") pos(12) col(2)) ///
              ytitle("") xtitle("") title("") `gopts' ///
              plotregion(margin(b=14))
    }
}	
	
}
    *------------------------------------------------------------------------------*
    * (8) Store results in e()
    *------------------------------------------------------------------------------*
        ereturn post `b' `V', esample(firstpost)
        ereturn local cmd        "lwdid"
        ereturn local depvar     "`y'"
        ereturn local rolling    "`rolling'"
        ereturn local controls   "`controls'"
        ereturn local gid  "`gid'"

        ereturn scalar att       = __att_overall
        ereturn scalar se_att    = __se_overall
        ereturn scalar K         = `K'
        ereturn scalar tpost1    = `tpost1'

        di as txt "lwdid (rolling: `rolling')  | first post-treatment tindex = `tpost1' "
        di as res "single ATT = " %9.4f e(att) "   SE = " %9.4f e(se_att)
if "`ri'" != "" {
		di as txt "============================================================"
			di as txt ">>> Randomization inference (RI) p-value = " %9.4f __p_ri 
		di as txt "============================================================"
}

        capture confirm matrix ATTt
        if !_rc ereturn matrix ATTt = ATTt
        capture confirm matrix SEt
        if !_rc ereturn matrix SEt  = SEt

    restore
end


