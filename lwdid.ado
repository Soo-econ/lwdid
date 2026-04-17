
*! lwdid - Lee & Wooldridge rolling DID estimator (unified: small-N + large-N)
*! version 2.1 April 16 2026
*! authors: Soo Jeong Lee, Jeffrey M. Wooldridge
*! contact: soojeong.lee@siu.edu, wooldri1@msu.edu
*! https://github.com/Soo-econ/lwdid.git  [Readme]

set more off
*set trace on
set trace off

capture program drop lwdid
capture program drop lwdid_small_single
capture program drop lwdid_small_staggered
capture program drop lwdid_large



**# [1] MAIN PROGRAM: lwdid 
***>> Dispatches to the appropriate subroutine: lwdid_small_single/lwdid_small_staggered/lwdid_large
program define lwdid, eclass sortpreserve
		version 16.0

		syntax varlist(min=1 numeric) [if] [in], ///
			IVAR(name) TVAR(name) GVAR(name)     ///
			ROLLING(name)                        ///
			[METHOD(name)                        ///   
			 Small                               ///    
			 REPS(integer 999)                   ///
			 LEVEL(integer 95)                   ///
			 CLUSTER(name)                       ///
			 SEED(integer -1)                    ///
			 VCE(string)                         ///
			 TABLE(string)                       ///
			 GRAPH                               ///
			 SCHEME(string)						 ///
			 GOPTS(string asis)                  ///
			 SAVE(string)                      ///
			 GID(string)                         ///
			 RI                                  ///
			 RIREPS(integer 999)                ///
			 RISEED(string)                      ///
			]

		marksample touse, novarlist

		* allow missing gvar → convert to never-treated with 0
		quietly replace `gvar' = 0 if missing(`gvar') & `touse'
		
		capture assert `gvar' == floor(`gvar') if `touse' & `gvar' > 0
		if _rc {
			di as err "gvar() must be integer-valued."
			exit 198
		}
		
		*-- parse varlist
		local y    : word 1 of `varlist'
		local xlist: list varlist - y       

			
		*-- rolling() check  (small-N adds demeanq/detrendq)
		local rolling = lower("`rolling'")
		if "`small'" != "" {
			local ok = inlist("`rolling'","demean","detrend","demeanq","detrendq","demeanm","detrendm")
		}
		else {
			local ok = inlist("`rolling'","demean","detrend")
		}
		if !`ok' {
			if "`small'" != "" {
				di as err "rolling() must be: demean | detrend | demeanq | detrendq | demeanm | detrendm"
			}
			else {
				di as err "rolling() must be: demean | detrend  (large-N mode)"
			}
			exit 198
		}

	
	*-- DISPATCH (small option)
		if "`small'" != "" {
			
				tempvar cohort
				quietly egen `cohort' = tag(`gvar') if `gvar' > 0 & `touse'
				quietly count if `cohort'
				local n_cohort = r(N)


				if `n_cohort' == 0 {
					di as err "gvar(): no treated units found."
					exit 2000
				}
				* with single treatment
				else if `n_cohort' == 1 {
					lwdid_small_single `0'            
					exit
				}
				*with staggered adoption
				else {
					lwdid_small_staggered `0'
					exit
				}
			}
				*with Large-N
		else {
				lwdid_large `0'
			}
end

**# [2] lwdid_small_single
**>> Subroutine for the small-N common timing case (single treated unit)
		
program define lwdid_small_single, eclass
		version 16.0

		syntax varlist(min=1 numeric) [if] [in], ///
			IVAR(name) TVAR(name)  GVAR(name)     ///
			ROLLING(name)                        ///
			[METHOD(name)                        ///   
			 Small                               ///    
			 REPS(integer 999)                   ///
			 LEVEL(integer 95)                   ///
			 CLUSTER(name)                       ///
			 SEED(integer -1)                    ///
			 VCE(string)                         ///
			 TABLE(string)                       ///
			 GRAPH                               ///
			 SCHEME(string)						 ///
			 GOPTS(string asis)                  ///
			 SAVE(string)                      ///
			 TITLE(string)                       ///
			 GID(string)                         ///
			 RI                                  ///
			 RIREPS(integer 999)                ///
			 RISEED(string)                      ///
			]

	

		marksample touse, novarlist
		local y    : word 1 of `varlist'
		local xlist: list varlist - y
		if "`controls'" != "" local xlist `controls'
		local rolling = lower("`rolling'")
		local method  = lower("`method'")
		if "`method'" == "" local method "ra"
		
		*-- RI seed
		qui capture confirm number `riseed'
		qui if _rc local riseed = .
		qui if missing(`riseed') local riseed = ceil(runiform()*1e6 + 1000*runiform())

		*-- validate variables
		confirm variable `gvar'
		confirm variable `ivar'
		confirm variable `y'
		confirm variable `tvar'
		capture confirm numeric variable `tvar'
		if _rc {
			di as err "tvar() must be a single numeric time variable."
			exit 198
		}

		*-- build internal time variable (single numeric index supplied in tvar())
		tempvar __timevar __qvar __mvar
		quietly gen double `__timevar' = `tvar' if `touse'
		local timevar `__timevar'
		local __tfmt : format `tvar'
		local __gfmt : format `gvar'

		if inlist("`rolling'","demeanq","detrendq") {
			quietly gen byte `__qvar' = quarter(dofq(`timevar')) if `touse'
			quietly count if `touse' & missing(`__qvar')
			if r(N) > 0 {
				di as err "rolling(`rolling') requires tvar() to be a Stata quarterly date variable created by yq()."
				exit 198
			}
			if strpos("`__gfmt'","%tq")==0 {
				di as err "rolling(`rolling') requires gvar() to be on the same quarterly scale as tvar() (format %tq, created by yq())."
				exit 198
			}
		}
		if inlist("`rolling'","demeanm","detrendm") {
			quietly gen byte `__mvar' = month(dofm(`timevar')) if `touse'
			quietly count if `touse' & missing(`__mvar')
			if r(N) > 0 {
				di as err "rolling(`rolling') requires tvar() to be a Stata monthly date variable created by ym()."
				exit 198
			}
			if strpos("`__gfmt'","%tm")==0 {
				di as err "rolling(`rolling') requires gvar() to be on the same monthly scale as tvar() (format %tm, created by ym())."
				exit 198
			}
		}

		*-- ensure ivar is numeric (internal id) + keep original id for gid()
		local id_orig `ivar'
		capture confirm numeric variable `ivar'
		if _rc {
			tempvar __id
			quietly egen long `__id' = group(`ivar') if `touse'
			local id `__id'
		}
		else {
			local id `ivar'
		}
			

		*-- count cohorts, build d and post
		tempvar cohort post_
		qui egen `cohort' = tag(`gvar') if `gvar' > 0 & `touse'
		qui count if `cohort'
		local n_cohort = r(N)

		* build d_ and post_ from gvar
		
		* treatment year
		quietly summarize `gvar' if `gvar' > 0 & `touse', meanonly
		local gyear = r(min)

		* treatment indicator
		capture drop d_
		gen byte d_ = (`gvar' > 0)

		* post indicator
		gen byte `post_' = (`timevar' >= `gyear') if `touse'
		replace `post_' = . if !`touse'
		di as txt "------------------------------------------------------------"
		di as txt "gvar(): single cohort detected -> common timing"
		di as txt "lwdid [small-N mode]  rolling=`rolling'  method=ols  ci=`ci'"
		di as txt "------------------------------------------------------------"



	* --- Keep analysis sample & create core flags
		tempvar tindex yhat ydot cs ctrlSum ydot_tr
		tempname b V
		
		preserve
		qui keep if `touse'
		qui drop if missing(`y', `id', `post_', d_)

	* ---  (2) Build time index: yearly or year-quarter, plus nice tq label
		quietly su `timevar', meanonly
		qui gen long `tindex' = `timevar' - r(min) + 1 if !missing(`timevar')

		label var `tindex' "time index (min->1)"

    * ---  Identify K (max pre tindex) and first post period (tpost1)
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


    * ---  For each unit, predict yhat for ALL periods (based on selected rolling)
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
                regress `y' i.`__qvar' if `id'==`ii' & `post_'==0
                predict double __fit if `id'==`ii', xb
            }
            else if "`rolling'"=="detrendq" {
                regress `y' c.`tindex' i.`__qvar' if `id'==`ii' & `post_'==0
                predict double __fit if `id'==`ii', xb
            }
            else if "`rolling'"=="demeanm" {
                regress `y' i.`__mvar' if `id'==`ii' & `post_'==0
                predict double __fit if `id'==`ii', xb
            }
            else if "`rolling'"=="detrendm" {
                regress `y' c.`tindex' i.`__mvar' if `id'==`ii' & `post_'==0
                predict double __fit if `id'==`ii', xb
            }
            qui replace `yhat' = __fit if `id'==`ii'
            cap qui drop __fit
        }



	* ---  Residualized outcome for ALL periods: ydot = y - yhat
		* Always also create internal ydot for the subsequent analysis
			quietly gen double `ydot' = `y' - `yhat'
			qui label var `ydot' "y - yhat (residualized via pre-period fit)"


	* ---  Collapse to first post cross-section & overall ATT
		quietly bysort `id': egen double ydot_postavg = mean(cond(`post_'==1, `ydot', .))
		qui label var ydot_postavg "mean(ydot) over post periods, by id"

		quietly gen byte firstpost = (`tindex'==`tpost1') & !missing(ydot_postavg, d_)


		* ---- cf. build Excel filename for single-effect table: <table>.xls ----
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

		* cf. Check outreg2 availability (Excel export depends on it)
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
		if "`xlist'" != "" {
			local K1 : word count `xlist'
		}
		local nk = `K1' + 1

		quietly count if `tindex'==1 & d_==1
		local nt = r(N)
		quietly count if `tindex'==1 & d_==0
		local nc = r(N)

		* only if controls specified AND hold the condition.
		if "`xlist'" != "" {

			if `nt' > `nk' & `nc' > `nk' {

				foreach x of varlist `xlist' {

					tempvar xc 

					* center X using treated mean
				   quietly summarize `x' if d_==1, meanonly
					quietly gen double `xc' = `x' - r(mean)

					gen double D_`x' = d_ * `xc'
					label var D_`x' "d_ × centered(`x')"

					local INTERACTS `INTERACTS' D_`x'
				}

				local RHS "`xlist' `INTERACTS'"
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

	
* --- Main regression & optional randomization inference

		* --- Main regression with/without vce() ---
		if "`vce'" != "" {
			di as txt " -- [1] Single ATT, with `vce' adjustment"
			regress ydot_postavg d_ `RHS' if firstpost, vce(`vce')

			if `__do_export' & `__has_or2' {
				outreg2 using "`__xls'", append label nose ///
					ctitle("`rolling' (OLS) with `vce'") addstat(N, e(N))
			}
		}
		else {
			di as txt " -- [1] Single ATT, without variance adjustment "
			regress ydot_postavg d_ `RHS' if firstpost

			if `__do_export' & `__has_or2' {
				outreg2 using "`__xls'", replace label nose ///
					ctitle("`rolling' (OLS)") addstat(N, e(N))
			}
		}

**# RI
			
		if "`ri'" != "" {
    quietly set rng mt64
    if "`riseed'" != "" quietly set seed `riseed'
    local reps = `rireps'

    scalar __b0 = _b[d_]
    mata: st_numscalar("__p_ri", ///
        lwdid_ri_inline(`reps', st_numscalar("__b0"), "`RHS'"))
}

		* keep these for e()
		matrix `b' = e(b)
		matrix `V' = e(V)
		scalar __att_overall = _b[d_]
		scalar __se_overall  = _se[d_]


		* ---[2] Period-by-period post effects table

		tempname __pph
		tempfile __ppfile
		qui postfile `__pph' int tindex str20 period double beta se tstat pval N using `__ppfile', replace

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
				quietly su `timevar' if `tindex'==`tt', meanonly
				if "`__tfmt'" == "%tq" {
					local lab : display %tq r(mean)
				}
				else if "`__tfmt'" == "%tm" {
					local lab : display %tm r(mean)
				}
				else {
					local lab : display %9.0f r(mean)
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

		capture quietly confirm file "`__ppfile'"
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
			qui rename tindex_s tindex
			qui capture drop __isavg

			label var tindex "time index"
			label var period "period"
			label var beta   "ATT_t"
			label var se     "SE(ATT_t)"
			label var tstat  "t"
			label var pval   "p-value"
			qui scalar cr=invttail(e(df_r), 0.025)
			qui gen double ci_lw = beta - cr*se
			qui gen double ci_up = beta + cr*se
			order period tindex beta se ci_lw ci_up tstat pval N
			
			di as txt " ----[2] Period-by-period ATTs"

			format beta se ci_lw ci_up %9.3f
			format tstat pval %9.3f

			list period beta se ci_lw ci_up tstat pval N, noobs 

			if "`save'" != "" {

			qui keep period beta se ci_lw ci_up tstat pval N
			qui rename beta att
			order period att se ci_lw ci_up tstat pval N

			* reset display format before saving
			format att se ci_lw ci_up %9.0g
			format tstat pval %9.0g

			save "`save'", replace
			qui rename att beta
		}
			
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

	* --- Control mean trajectory & optional graph
		* Always compute control-group mean over time (d==0)
			qui gen double `ctrlSum' = `ydot' if d_==0
			bysort `tindex': egen double y_dot_cont = mean(`ctrlSum')
			label var y_dot_cont "Control mean of ydot (d==0)"

		* Always compute treated-group mean over time (d==1)
			tempvar trSum
			qui gen double `trSum' = `ydot' if d_==1
			bysort `tindex': egen double y_dot_tr = mean(`trSum')
			label var y_dot_tr "Treated mean of ydot (d==1)"
	
	
	* --- Ensure there is at least one treated obs
		quietly count if d_==1
		if r(N)==0 {
			di as err "No treated units found (d==1)."
			exit 2000
		}

		
				
	* --- Graph — Specific treated unit specified by gid()
	* ---(1) normalize gid() ---
		local gid = trim("`gid'")

		* --- choose treated series for plotting ---
		*     - if gid() omitted: use average treated series y_dot_tr
		*     - if gid() specified: create series `ydot_tr' for that unit only
		tempvar ydot_tr
		if "`gid'" == "" {

			* avg treated (already computed elsewhere as y_dot_tr)
			local TRSER  "y_dot_tr"
			local TRLAB  "Treated"
		}
		else {
			*--- Case 1: original ivar is numeric ---
			capture confirm numeric variable `id_orig'
			if !_rc {

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
				qui gen double `ydot_tr' = `ydot' if `id_orig'==`__gidnum'
			}
			else {
			*--- Case 2: original ivar is string ---
				quietly count if `id_orig'=="`gid'" & d_==1
				if r(N)==0 {
					di as err "gid(`gid') not found among treated units (d==1)."
					exit 2000
				}

				capture drop `ydot_tr'
				qui gen double `ydot_tr' = `ydot' if `id_orig'=="`gid'"
			}

			label var `ydot_tr' "Treated unit (id=`gid')"
			local TRSER  "`ydot_tr'"
			local TRLAB  "Treated (ID=`gid')"
		}

	* --- If graph requested, build x-axis + y-axis once, then plot     *

		if "`graph'" != "" {
			local ISQ = strpos("`__tfmt'","%tq")>0
			local ISM = strpos("`__tfmt'","%tm")>0

			tempvar xdate __yearlab
			qui gen double `xdate' = `timevar'
			local XVAR "`xdate'"
			quietly su `timevar' if `post_'==1, meanonly
			local XLINE = r(min)
			local XTITLE "xtitle("Time")"
			quietly su `xdate', meanonly
			local xmin = r(min)
			local xmax = r(max)

			if `ISM' {
				qui gen int `__yearlab' = yofd(dofm(`xdate')) if !missing(`xdate')
			}
			else if `ISQ' {
				qui gen int `__yearlab' = yofd(dofq(`xdate')) if !missing(`xdate')
			}
			else {
				qui gen int `__yearlab' = `xdate' if !missing(`xdate')
			}

			quietly su `__yearlab', meanonly
			local yminx = floor(r(min))
			local ymaxx = floor(r(max))
			local yspan = `ymaxx' - `yminx'

			local XLABS
			local xang = 0
			if `ISM' {
				* keep monthly labels horizontal, matching quarterly style
				local xang = 0
				local mspan = `xmax' - `xmin'
				if `mspan' <= 24       local xstep = 3
				else if `mspan' <= 60  local xstep = 6
				else if `mspan' <= 120 local xstep = 12
				else if `mspan' <= 240 local xstep = 24
				else                   local xstep = 60
				forvalues xx = `=floor(`xmin')'(`xstep')`=floor(`xmax')' {
					local lab : display %tm `xx'
					local XLABS `XLABS' `xx' `"`lab'"'
				}
			}
			else if `ISQ' {
				local qspan = `xmax' - `xmin'
				if `qspan' <= 16       local xstep = 2
				else if `qspan' <= 40  local xstep = 4
				else if `qspan' <= 80  local xstep = 8
				else if `qspan' <= 160 local xstep = 20
				else                   local xstep = 40
				forvalues xx = `=floor(`xmin')'(`xstep')`=floor(`xmax')' {
					local lab : display %tq `xx'
					local XLABS `XLABS' `xx' `"`lab'"'
				}
			}
			else {
				* for yearly data, keep the original axis behavior (no sparse relabeling)
			}

			if `"`XLABS'"' != "" {
				local XL "xlabel(`XLABS', angle(`xang') labsize(vsmall) nogrid) xscale(range(`xmin' `xmax'))"
			}
			else {
				local XL "xscale(range(`xmin' `xmax'))"
			}

			*---------------- Y-axis ticks: nice step (>=0.1), aligned bounds, max tick count ----------------*
			local y_max_ticks = 6
			local y_min_step  = 0.1

			quietly su y_dot_cont, meanonly
			local ymin = r(min)
			local ymax = r(max)

			quietly su `TRSER', meanonly
			local ymin = min(`ymin', r(min))
			local ymax = max(`ymax', r(max))

			local yrange = `ymax' - `ymin'

			if (`yrange' <= 0) {
				local ystep = `y_min_step'
				local ymin2 = `ymin' - 5*`ystep'
				local ymax2 = `ymax' + 5*`ystep'
			}
			else {
				local raw = max(`yrange'/`y_max_ticks', `y_min_step')
				local k     = floor(log10(`raw'))
				local base  = 10^`k'
				local frac  = `raw'/`base'
				local nice  = cond(`frac'<=1, 1, cond(`frac'<=2, 2, cond(`frac'<=5, 5, 10)))
				local ystep = max(`nice'*`base', `y_min_step')
				local ymin2 = `ystep'*floor(`ymin'/`ystep')
				local ymax2 = `ystep'*ceil(`ymax'/`ystep')
				local nticks = floor((`ymax2' - `ymin2')/`ystep') + 1
				local guard  = 0
				while (`nticks' > `y_max_ticks' & `guard' < 10) {
					local ystep  = `ystep' * 2
					local ymin2  = `ystep'*floor(`ymin'/`ystep')
					local ymax2  = `ystep'*ceil(`ymax'/`ystep')
					local nticks = floor((`ymax2' - `ymin2')/`ystep') + 1
					local guard  = `guard' + 1
				}
			}

			local yfmt = cond(`ystep' >= 1, "%9.0f", "%9.1f")
			local YL   "ylabel(`ymin2'(`ystep')`ymax2', format(`yfmt') nogrid)"

			local mono = 0
			if strpos("`scheme'","mono") local mono = 1

			if `mono' {
				local col_tr black%80
				local col_ct black%60
				local pat_tr solid
				local pat_ct dash
				local gsch  "scheme(s1mono)"
				local aheight "aspectratio(0.45)"
			}
			else {
				local col_tr cranberry%90
				local col_ct navy%70
				local pat_tr solid
				local pat_ct dash
				local gsch  ""
				local aheight ""
			}

			twoway ///
				(line `TRSER' `XVAR', lpattern(`pat_tr') lcolor(`col_tr') lwidth(medthick)) ///
				(line y_dot_cont `XVAR', lpattern(`pat_ct') lcolor(`col_ct') lwidth(medthick)) ///
				, ///
				xline(`XLINE', lcolor(gs8) lpattern(dash)) ///
				`XL' ///
				`YL' ///
				`XTITLE' ///
				legend(order(1 "Treated" 2 "Control") ///
					   pos(2) ring(0) col(2) ///
					   region(lcolor(none))) ///
				graphregion(color(white)) ///
				plotregion(color(white)) ///
				`gsch' ///
				`gopts' ///
				`aheight'
		}
			
			
   * Store results in e()
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
		if "`ri'" != "" ereturn scalar p_ri = __p_ri
		
		di as res "Single ATT = " %9.3f e(att) "   SE = " %9.3f e(se_att)
		
		if "`ri'" != "" {
				di as txt "-------------------------------------------"
					di as res "Randomization inference (RI) p-value = " %9.3f e(p_ri) 
				di as txt "-------------------------------------------"
		}

				capture confirm matrix ATTt
				if !_rc ereturn matrix ATTt = ATTt
				capture confirm matrix SEt
				if !_rc ereturn matrix SEt  = SEt
	restore
end


**# [3] lwdid_small_staggered
**>> Subroutine for the small-N staggered adoption (multiple treated units)
program define lwdid_small_staggered, eclass
    version 16.0

    syntax varlist(min=1 numeric) [if] [in], ///
        IVAR(name) TVAR(name) GVAR(name)     ///
        ROLLING(name)                        ///
        [METHOD(name)                        ///
         Small                               ///
         REPS(integer 999)                   ///
         LEVEL(integer 95)                   ///
         CLUSTER(name)                       ///
         SEED(integer -1)                    ///
         VCE(string)                         ///
         TABLE(string)                       ///
         GRAPH                               ///
         SCHEME(string)                      ///
         GOPTS(string asis)                  ///
         SAVE(string)                        ///
         GID(string)                         ///
         RI                                  ///
         RIREPS(integer 999)                 ///
         RISEED(string)                      ///
        ]

    marksample touse, novarlist

	
		* --- ensure ivar is numeric (internal id)
		local id_orig `ivar'

		capture confirm numeric variable `ivar'
		if _rc {
			tempvar __id
			quietly egen long `__id' = group(`ivar') if `touse'
			local id `__id'
		}
		else {
			local id `ivar'
		}
			
			
		local y : word 1 of `varlist'
		local xlist : list varlist - y
		local rolling = lower("`rolling'")

		*-- controls not yet supported
		if "`xlist'" != "" {
			di as err "Controls are not yet supported for staggered treatment timing."
			di as err "Please run lwdid without additional regressors."
			exit 198
		}


		*-- validate variables
		confirm variable `gvar'
		confirm variable `ivar'
		confirm variable `y'
		confirm variable `tvar'
		capture confirm numeric variable `tvar'
		if _rc {
			di as err "tvar() must be a single numeric time variable."
			exit 198
		}

		tempvar timevar __qvar __mvar
		quietly gen double `timevar' = `tvar' if `touse'
		local __tfmt : format `tvar'
		if inlist("`rolling'","demeanq","detrendq") {
			quietly gen byte `__qvar' = quarter(dofq(`timevar')) if `touse'
			quietly count if `touse' & missing(`__qvar')
			if r(N) > 0 {
				di as err "rolling(`rolling') requires tvar() to be a Stata quarterly date variable created by yq()."
				exit 198
			}
		}
		if inlist("`rolling'","demeanm","detrendm") {
			quietly gen byte `__mvar' = month(dofm(`timevar')) if `touse'
			quietly count if `touse' & missing(`__mvar')
			if r(N) > 0 {
				di as err "rolling(`rolling') requires tvar() to be a Stata monthly date variable created by ym()."
				exit 198
			}
		}

		*-- count cohorts for message
		tempvar cohorttag
		quietly egen `cohorttag' = tag(`gvar') if `gvar' > 0 & `touse'
		quietly count if `cohorttag'
		local n_cohort = r(N)

		di as txt "gvar(): `n_cohort' cohorts detected -> Staggered adoption"
		di as txt "lwdid [small-N mode] rolling=`rolling'"

		*-- graph not yet implemented
		if "`graph'" != "" {
			di as txt "------------------------------------------------------------"
			di as err "graph option not yet supported for small-N staggered designs."
			di as txt "This feature will be included in a future update of lwdid."
			di as txt "Please run the command **without the graph option.**"
			di as txt "------------------------------------------------------------"
			exit 198
		}

  quietly {
        capture drop d_
        gen byte d_ = (`gvar' != 0) if `touse'
        replace d_ = 0 if missing(d_)
        label var d_ "ever-treated indicator"

        preserve
            keep if `touse'
            drop if missing(`y', `ivar', `timevar')

			
			/* later
			* --- reproducibility setup (only if RI requested) ---
			if "`ri'" != "" {
				quietly set rng mt64

				capture confirm number `riseed'
				if _rc | missing(real("`riseed'")) {
					local riseed = ceil(runiform()*1e6 + 1000*runiform())
				}

				quietly set seed `riseed'
				local reps = `rireps'
			}
			*/

            *-- cohort list
            levelsof `gvar' if `gvar' > 0, local(Glist)

            *-- create cohort-specific residualized outcomes
            foreach g of local Glist {
                tempvar yhat
                gen double `yhat' = .

                levelsof `id', local(IDlist)

                foreach idval of local IDlist {
                    capture quietly {
                        if "`rolling'" == "demean" {
                            qui regress `y' if `id'==`idval' & `timevar' < `g'
                        }
                        else if "`rolling'" == "detrend" {
                            qui regress `y' c.`timevar' if `id'==`idval' & `timevar' < `g'
                        }
                        else if "`rolling'" == "demeanq" {
                            regress `y' i.`__qvar' if `id'==`idval' & `timevar' < `g'
                        }
                        else if "`rolling'" == "detrendq" {
                            regress `y' c.`timevar' i.`__qvar' if `id'==`idval' & `timevar' < `g'
                        }
                        else if "`rolling'" == "demeanm" {
                            regress `y' i.`__mvar' if `id'==`idval' & `timevar' < `g'
                        }
                        else if "`rolling'" == "detrendm" {
                            regress `y' c.`timevar' i.`__mvar' if `id'==`idval' & `timevar' < `g'
                        }
                        else {
                            di as err "rolling() must be demean, detrend, demeanq, detrendq, demeanm, or detrendm"
                            exit 198
                        }
                    }

                    if _rc == 0 {
                        predict double __fit if `id'==`idval', xb
                        replace `yhat' = __fit if `id'==`idval'
                        cap drop __fit
                    }
                    else {
                        cap drop __fit
                    }
                }

                gen double y`g'd = `y' - `yhat'
                label var y`g'd "(rolling=`rolling') Residualized outcome cohort g=`g' "
                drop `yhat'
            }

            *-- treated post averages
			foreach g of local Glist {
				bysort `id': egen double ybar_temp_`g' = ///
					mean(cond(`timevar' >= `g' & d_==1, y`g'd, .))
			}

            *-- cohort weights
            tempvar gtemp
            gen `gtemp' = `gvar'
            replace `gtemp' = . if `gtemp' == 0
            tab `gtemp', matcell(freqs)
            matrix w = freqs / r(N)

            *-- control weighted averages
            local i = 0
            foreach g of local Glist {
                local ++i
                scalar w_g = w[`i',1]
                bysort `id': egen double ybar_cont_`g' = ///
                    mean(cond(`timevar' >= `g' & d_==0, w_g*y`g'd, .))
            }

            *-- build ydot_bar (equation 7.18)
            gen double ydot_bar = .

            foreach g of local Glist {
                replace ydot_bar = ybar_temp_`g' if `gvar' == `g'
            }

            gen double __contsum = 0
            foreach g of local Glist {
                replace __contsum = __contsum + ybar_cont_`g'
            }

            qui replace ydot_bar = __contsum if d_==0
            drop __contsum

            *-- first treated cohort
            summarize `gvar' if `gvar'>0, meanonly
			local gmin = r(min)

			* first observed post-treatment period
			summarize `timevar' if `timevar'>=`gmin', meanonly
			local tpost1 = r(min)

			qui regress ydot_bar d_ if `timevar'==`tpost1'
    }

    di as txt "*--- Aggregated Single treatment effect (Lee & Wooldridge: equation 7.18)"

    regress ydot_bar d_ if `timevar' == `gmin'
    matrix b = e(b)
    matrix V = e(V)

    ereturn post b V
    ereturn scalar att = _b[d_]
    ereturn scalar se_att = _se[d_]
    ereturn local cmd     "lwdid"
    ereturn local depvar  "`y'"
    ereturn local rolling "`rolling'"

    di as res "lwdid (rolling: `rolling') ATT = " %9.3f e(att) ///
              "   SE = " %9.3f e(se_att)

    restore
end

**# [4] lwdid_large
**>> Subroutine for the LARGE-N (Common&staggered)
program define lwdid_large, eclass
    version 16.0

    syntax varlist(min=1 numeric) [if] [in], ///
        IVAR(name) TVAR(name)  GVAR(name)     ///
        ROLLING(name)                        ///
        [METHOD(name)                        ///   
         Small                               ///    
         REPS(integer 999)                   ///
         LEVEL(integer 95)                   ///
         CLUSTER(name)                       ///
         SEED(integer -1)                    ///
         VCE(string)                         ///
         TABLE(string)                       ///
         GRAPH                               ///
		 SCHEME(string)						 ///
         GOPTS(string asis)                  ///
         SAVE(string)                      ///
         TITLE(string)                       ///
         GID(string)                         ///
         RI                                  ///
         RIREPS(integer 999)                ///
         RISEED(string)                      ///
        ]

		marksample touse, novarlist
		local y    : word 1 of `varlist'
		local xlist: list varlist - y
		local rolling = lower("`rolling'")

		local method  = lower("`method'") 
			local method = lower("`method'")
			*-- large-N: method required
			if "`method'" == "" {
				di as err "method() required for large-N mode: ra, ipw, or ipwra"
				exit 198
			}
			if !inlist("`method'","ra","ipwra","ipw") {
				di as err "method() must be ra, ipw, or ipwra  (large-N mode)"
				exit 198
			}
			if inlist("`method'","ipwra","ipw") & "`xlist'"=="" {
				di as err "method(`method') requires covariates."
				exit 198
			}
			
			*--- Confidence bands
			local ci_type "simultaneous"
			
			
		* --- Internal numeric id for large-N computations (xtset-safe) / keep original id as well.
			local id_orig `ivar'
			tempvar __id
			capture confirm numeric variable `ivar'
			if _rc {
				quietly egen long `__id' = group(`ivar') if `touse'
				local id `__id'
			}
			else {
				local id `ivar'
			}	
				
		* --- Graph-scheme option
			if "`scheme'" == "" {
					local scheme s2color
				}

		* --- LARGE-N PATH  (wild bootstrap, staggered, influence-function)
			di as txt "------------------------------------------------------------"
			di as txt " lwdid [large-N mode]  rolling=`rolling'  method=`method'"

			* --- cluster (default = internal id; cluster() can still be user-specified)
			if "`cluster'" == "" {
				local cluster_var `id'
			}
			else {
				local cluster_var `cluster'
			}
			
			* --- xtset
			quietly xtset `id' `tvar'

			* --- time range
			qui su `tvar' if `touse', meanonly
			local tmin = r(min)
			local tmax = r(max)

			* --- total obs
			qui count if `touse'
			local Nobs = r(N)

		* --- Cohort dummies + centered X from gvar
			tempvar gvar_clean dinf
			qui gen `gvar_clean' = `gvar' if `touse'
			
			qui gen byte `dinf' = (`gvar_clean'==0) if `touse'
			qui replace `dinf' = 0 if missing(`dinf')
			label var `dinf' "Never treated"
			
			qui levelsof `gvar_clean' if `gvar_clean' > 0 & `touse', local(cohorts)
			if "`cohorts'" == "" {
				di as err "No treated cohorts found in gvar()."; exit 2000
				}

			foreach g of local cohorts {
				qui capture drop d`g'
				qui gen byte d`g' = (`gvar_clean' == `g') if `touse'
				qui replace d`g' = 0 if missing(d`g')
				label var d`g' "First treated in `g'"
				}

			if "`xlist'" != "" {
				qui {
						foreach v in `xlist' {
							foreach g of local cohorts {
								su `v' if d`g' == 1 & `touse', meanonly
								local m_`v'_`g' = r(mean)
								capture drop `v'_g`g'
								gen double `v'_g`g' = `v' - `m_`v'_`g'' if `touse'
							}
						}
					}
				}

		* ---  time dummies
			foreach t of numlist `tmin'/`tmax' {
				qui capture drop f`t'
				qui gen byte f`t' = (`tvar' == `t') if `touse'
				qui replace f`t' = 0 if missing(f`t')
			}


		* --- Stage 1: Residualized outcome y{g}d for each cohort
				* >> using the reduced-form representation (computationally efficient)
		
				quietly {
				tempvar yobs cy ct ctt cty cn
				qui gen byte `yobs' = (`touse' & !missing(`y'))
				bys `id' (`tvar'): gen double `cy'  = sum(cond(`yobs', `y',        0))
				bys `id' (`tvar'): gen double `ct'  = sum(cond(`yobs', `tvar',     0))
				bys `id' (`tvar'): gen double `ctt' = sum(cond(`yobs', `tvar'^2,   0))
				bys `id' (`tvar'): gen double `cty' = sum(cond(`yobs', `tvar'*`y', 0))
				bys `id' (`tvar'): gen double `cn'  = sum(cond(`yobs', 1,          0))

				foreach g of local cohorts {
					capture drop y`g'd
					gen double y`g'd = . if `touse'
				if ("`rolling'" == "demean") {
					tempvar Sy_pre n_pre

					* Total sum and count of pre-treatment outcomes for each unit
					bys `id': egen double `Sy_pre' = total(cond(`tvar' < `g' & `yobs', `y', 0)) if `touse'
					bys `id': egen double `n_pre'  = total(cond(`tvar' < `g' & `yobs', 1, 0)) if `touse'

					* Pre- and post-treatment periods:
					* Subtract the mean over all pre-treatment periods
					replace y`g'd = `y' - (`Sy_pre'/`n_pre') ///
						if `touse' & `yobs' & `n_pre' > 0 ///
						& !missing(`Sy_pre', `n_pre')

					drop `Sy_pre' `n_pre'
				}
					else {  // detrend: use one fixed pre-treatment trend for both pre and post
						tempvar SyP StP SttP StyP nP denomP bP aP fitP

						* Pre-treatment totals for each unit (using all periods t < g)
						bys `id': egen double `SyP'  = max(cond(`tvar' < `g', `cy',  .)) if `touse'
						bys `id': egen double `StP'  = max(cond(`tvar' < `g', `ct',  .)) if `touse'
						bys `id': egen double `SttP' = max(cond(`tvar' < `g', `ctt', .)) if `touse'
						bys `id': egen double `StyP' = max(cond(`tvar' < `g', `cty', .)) if `touse'
						bys `id': egen double `nP'   = max(cond(`tvar' < `g', `cn',  .)) if `touse'

						* Slope and intercept from the full pre-treatment sample
						gen double `denomP' = `nP' * `SttP' - (`StP')^2 if `touse'
						gen double `bP' = .
						replace `bP' = (`nP' * `StyP' - `StP' * `SyP') / `denomP' ///
							if `touse' & `nP' > 1 & `denomP' != 0

						gen double `aP' = .
						replace `aP' = (`SyP' - `bP' * `StP') / `nP' ///
							if `touse' & `nP' > 1 & !missing(`bP')

						* One fixed fitted pre-trend for all periods
						gen double `fitP' = .
						replace `fitP' = `aP' + `bP' * `tvar' ///
							if `touse' & !missing(`aP', `bP')

						* Apply the same detrending transformation to both pre and post periods
						replace y`g'd = `y' - `fitP' ///
							if `touse' & `yobs' & !missing(`fitP')

						drop `SyP' `StP' `SttP' `StyP' `nP' `denomP' `bP' `aP' `fitP'
					}
					label var y`g'd "(rolling=`rolling') Residualized outcome cohort g=`g' "
				}
			 } 
			 
	 
	* --- Stage 2: ATT(g,t) point estimates + Influence Functions
			quietly {
					preserve
					keep if `touse' & `gvar' > 0
					keep if `touse' & `gvar' > 0
					bys `gvar' `ivar': gen byte one = (_n==1)
					bys `gvar': egen Ng = total(one)
					keep `gvar' Ng
					duplicates drop
					rename `gvar' cohort
					tempfile COHORTSIZE
					save "`COHORTSIZE'", replace
					restore
				}

				qui tempfile ATTfile
				qui tempname pf
				qui postfile `pf' int cohort int time int ryear double att using "`ATTfile'", replace

				quietly tempvar cont
				qui gen byte `cont' = 0

				mata: IF_mat = J(`Nobs', 0, .)
				mata: cell_g = J(1, 0, .)
				mata: cell_t = J(1, 0, .)
				local cell_count = 0

				quietly foreach g of local cohorts {
					local yvar y`g'd
					local dvar_g d`g'

					local current_x_list ""
					foreach v in `xlist' {
						local current_x_list `current_x_list' `v'_g`g'
					}

			forval t = `tmin'/`tmax' {
				local r = `t' - `g'

				quietly replace `cont' = 0 if `touse'

				if `t' < `g' {
					* pre-period: never-treated + cohort g + cohorts treated after g
					quietly replace `cont' = (`gvar' == 0 | `gvar' == `g' | `gvar' > `g') if `touse'
				}
				else {
					* post-period: never-treated + cohort g + cohorts treated after g
					* and still untreated at time t
					quietly replace `cont' = (`gvar' == 0 | `gvar' == `g' | (`gvar' > `g' & `gvar' > `t')) if `touse'
				}

				quietly count if `touse' & f`t' & `cont'
				if r(N) == 0 continue

				tempvar esamp_gt psamp_gt
				qui gen byte `esamp_gt' = 0 if `touse'
				qui gen byte `psamp_gt' = 0 if `touse'

				* ------------------------------------------------------------
				* Propensity score estimation and ATT weights
				* ------------------------------------------------------------
				if inlist("`method'", "ipw", "ipwra") {
					tempvar p_hat_gt ipw_gt

					sort `id' `tvar'
					cap qui logit `dvar_g' `xlist' if `touse' & f`t' & `cont', nolog
					if _rc != 0 continue

					qui replace `psamp_gt' = e(sample)

					qui predict double `p_hat_gt' if `psamp_gt', pr
					qui gen double `ipw_gt' = cond(`dvar_g' == 1, 1, `p_hat_gt' / (1 - `p_hat_gt')) ///
						if `psamp_gt'
				}

				* ------------------------------------------------------------
				* Point estimation: ATT(g,t)
				* ------------------------------------------------------------
				if "`method'" == "ra" {
					if "`xlist'" == "" {
						qui regress `yvar' `dvar_g' if `touse' & f`t' & `cont'
					}
					else {
						qui reg `yvar' `dvar_g' `xlist' c.`dvar_g'#c.(`current_x_list') ///
							if `touse' & f`t' & `cont'
					}
				}
				else if "`method'" == "ipw" {
					qui reg `yvar' `dvar_g' [aw=`ipw_gt'] if `psamp_gt'
				}
				else if "`method'" == "ipwra" {
					if "`xlist'" == "" {
						qui reg `yvar' `dvar_g' [aw=`ipw_gt'] if `psamp_gt'
					}
					else {
						qui reg `yvar' `dvar_g' `current_x_list' c.`dvar_g'#c.(`current_x_list') ///
							[aw=`ipw_gt'] if `psamp_gt'
					}
				}
				else if "`method'" == "psmatch" {
					teffects psmatch (`yvar') (`dvar_g' `xlist') ///
						if `touse' & f`t' & `cont', atet
				}

				if _rc != 0 {
					if inlist("`method'", "ipw", "ipwra") {
						cap drop `p_hat_gt' `ipw_gt'
					}
					drop `esamp_gt' `psamp_gt'
					continue
				}

				qui replace `esamp_gt' = e(sample)

		* ------------------------------------------------------------
		* Store point estimate
		* ------------------------------------------------------------
		qui tempname b_att
		if inlist("`method'", "ra", "ipw", "ipwra") scalar `b_att' = _b[`dvar_g']
		else                                         scalar `b_att' = _b[ATET:r1vs0.`dvar_g']
		post `pf' (`g') (`t') (`r') (`b_att')

		* ------------------------------------------------------------
		* Influence function
		* ------------------------------------------------------------
		if inlist("`method'", "ra", "ipw", "ipwra") {

			tempvar IF_gt
			qui gen double `IF_gt' = 0 if `touse'

			qui su `dvar_g' if `esamp_gt', meanonly
			local mean_D = r(mean)

			if "`method'" == "ra" {

            tempvar uhat
            qui predict double `uhat' if `esamp_gt', resid

            local zvars `dvar_g' `xlist'
            local intvars
            if "`xlist'" != "" {
                foreach v in `current_x_list' {
                    tempvar intv
                    qui gen double `intv' = `dvar_g' * `v' if `esamp_gt'
                    local intvars `intvars' `intv'
                }
                local zvars `zvars' `intvars'
            }

            local IFname `IF_gt'
            local ESname `esamp_gt'
            local Uname  `uhat'

            mata: es  = st_data(., "`ESname'")
            mata: idx = selectindex(es :== 1)
            mata: Z   = st_data(idx, tokens("`zvars'"))
            mata: Z   = J(rows(Z),1,1), Z
            mata: u   = st_data(idx, "`Uname'")
            mata: ZZinv = invsym(quadcross(Z,Z))
            mata: IFb   = (Z * ZZinv) :* u
            mata: ifvec = J(rows(es),1,0)
            mata: ifvec[idx] = IFb[,2]
            mata: st_store(., "`IFname'", ifvec)

            qui replace `IF_gt' = 0 if missing(`IF_gt')
        }
        else if "`method'" == "ipw" {

            qui count if `esamp_gt'
            local n_eff_gt = r(N)

            tempvar __plug
            qui gen double `__plug' = ///
                ((`dvar_g' / `mean_D') * (`yvar' - `b_att') ///
                - ((1 - `dvar_g') * `p_hat_gt' / ((1 - `p_hat_gt') * `mean_D')) * (`yvar')) ///
                if `esamp_gt'

            local IFname `IF_gt'
            local ESname `esamp_gt'
            local Yname  `yvar'
            local Dname  `dvar_g'
            local Pname  `p_hat_gt'
            local PLname `__plug'

            mata: es   = st_data(., "`ESname'")
            mata: idx  = selectindex(es :== 1)
            mata: y    = st_data(idx, "`Yname'")
            mata: d    = st_data(idx, "`Dname'")
            mata: p    = st_data(idx, "`Pname'")
            mata: plug = st_data(idx, "`PLname'")
            mata: X    = st_data(idx, tokens("`xlist'"))
            mata: X    = J(rows(X),1,1), X
            mata: W    = p :* (1 :- p)
            mata: A    = quadcross(X, X :* W)
            mata: Ainv = invsym(A)
            mata: s    = X :* (d :- p)
            mata: IFg  = s * Ainv
            mata: gvec = ((1 :- d) :* y :* p :/ (1 :- p))
            mata: Gamma = mean(X :* gvec)'
            mata: corr = IFg * Gamma
            mata: ifvec = J(rows(es),1,0)
            mata: ifvec[idx] = (plug :- corr) :/ `n_eff_gt'
            mata: st_store(., "`IFname'", ifvec)

            drop `__plug'
        }
        else if "`method'" == "ipwra" {

            tempvar uhat
            qui predict double `uhat' if `esamp_gt', resid

            if "`xlist'" == "" {

                local IFname `IF_gt'
                local ESname `esamp_gt'
                local Uname  `uhat'
                local Wname  `ipw_gt'

                mata: es   = st_data(., "`ESname'")
                mata: idx  = selectindex(es :== 1)
                mata: Z    = st_data(idx, tokens("`dvar_g'"))
                mata: Z    = J(rows(Z),1,1), Z
                mata: u    = st_data(idx, "`Uname'")
                mata: w    = st_data(idx, "`Wname'")
                mata: Qw   = quadcross(Z, Z :* w)
                mata: Qwinv = invsym(Qw)
                mata: M    = Z :* (w :* u)
                mata: IFb  = M * Qwinv
                mata: ifvec = J(rows(es),1,0)
                mata: ifvec[idx] = IFb[,2]
                mata: st_store(., "`IFname'", ifvec)
            }
            else {

                local zvars `dvar_g' `current_x_list'
                local intvars
                foreach v in `current_x_list' {
                    tempvar intv
                    qui gen double `intv' = `dvar_g' * `v' if `esamp_gt'
                    local intvars `intvars' `intv'
                }
                local zvars `zvars' `intvars'

                local IFname `IF_gt'
                local ESname `esamp_gt'
                local Uname  `uhat'
                local Wname  `ipw_gt'
                local Pname  `p_hat_gt'
                local Dname  `dvar_g'
                local Zname  `zvars'
                local Xps    `xlist'

                mata: es    = st_data(., "`ESname'")
                mata: idx   = selectindex(es :== 1)
                mata: Z     = st_data(idx, tokens("`Zname'"))
                mata: Z     = J(rows(Z),1,1), Z
                mata: X     = st_data(idx, tokens("`Xps'"))
                mata: X     = J(rows(X),1,1), X
                mata: u     = st_data(idx, "`Uname'")
                mata: w     = st_data(idx, "`Wname'")
                mata: p     = st_data(idx, "`Pname'")
                mata: d     = st_data(idx, "`Dname'")
                mata: Qw    = quadcross(Z, Z :* w)
                mata: Qwinv = invsym(Qw)
                mata: M     = Z :* (w :* u)
                mata: A     = quadcross(X, X :* (p :* (1 :- p)))
                mata: Ainv  = invsym(A)
                mata: S     = X :* (d :- p)
                mata: IFg   = S * Ainv
                mata: H     = quadcross(Z :* ((1 :- d) :* w :* u), X)
                mata: IFb   = (M + IFg * H') * Qwinv
                mata: ifvec = J(rows(es),1,0)
                mata: ifvec[idx] = IFb[,2]
                mata: st_store(., "`IFname'", ifvec)
            }

            qui replace `IF_gt' = 0 if missing(`IF_gt')
        }

        qui replace `IF_gt' = 0 if missing(`IF_gt')

        tempvar if_full
        qui gen double `if_full' = 0 if `touse'
        qui replace `if_full' = `IF_gt' if `esamp_gt' == 1 & `touse'

        mata: IF_col = st_data(., "`if_full'", "`touse'")
        mata: IF_mat = IF_mat, IF_col
        mata: cell_g = cell_g, `g'
        mata: cell_t = cell_t, `t'
        local cell_count = `cell_count' + 1

				drop `if_full' `IF_gt'
				if inlist("`method'", "ipw", "ipwra") cap drop `p_hat_gt' `ipw_gt'
			}
			else {
				if inlist("`method'", "ipw", "ipwra") cap drop `p_hat_gt' `ipw_gt'
			}

			drop `esamp_gt' `psamp_gt'
		}
				}
					postclose `pf'


			
	* --- Stage 2b: WATT(r) aggregation
			quietly {
			tempfile WATT_point WATT_weights
			preserve
			use "`ATTfile'", clear
			merge m:1 cohort using "`COHORTSIZE'", nogen
			if ("`rolling'" == "demean") {
				drop if missing(att)
			}
	else if ("`rolling'" == "detrend") {
				replace att = 0 if inlist(ryear,-1,-2)
				drop if missing(att) & !inlist(ryear,-1,-2)
			}
			tempvar tag
			egen byte `tag' = tag(ryear cohort)
			bys ryear: egen N_cohort = total(`tag')
			bys ryear: egen N_units  = total(Ng)
			gen double weight = Ng / N_units
			duplicates drop ryear cohort, force
			keep ryear cohort Ng N_units weight
			save "`WATT_weights'", replace
			restore

			preserve
			use "`ATTfile'", clear
			merge m:1 cohort using "`COHORTSIZE'", nogen
			if ("`rolling'" == "demean") {
				drop if missing(att)
			}
	else if ("`rolling'" == "detrend") {
				replace att = 0 if inlist(ryear,-1,-2)
				drop if missing(att) & !inlist(ryear,-1,-2)
			}
			tempvar tag
			qui egen byte `tag' = tag(ryear cohort)
			bys ryear: egen N_cohort = total(`tag')
			bys ryear: egen N_units  = total(Ng)
			gen double weight = Ng / N_units
			gen double w_att  = weight * att
			duplicates drop ryear cohort, force
			collapse (sum) watt = w_att (sum) weight (first) N_cohort N_units, by(ryear)
			sort ryear
			qui save "`WATT_point'", replace
			restore
			}

			
			
	* --- Stage 3: Wild Bootstrap for WATT(r)

        if inlist("`method'","ra","ipw","ipwra") & `cell_count' > 0 {
            
            qui if `seed' >= 0 set seed `seed'
            quietly {
                preserve
                use "`WATT_point'", clear
                qui drop if missing(watt)
                sort ryear
                mkmat ryear watt, matrix(WATT_pmat)
                restore
                mata: WATT_pmat = st_matrix("WATT_pmat")

                preserve
                use "`WATT_weights'", clear
                mkmat ryear cohort weight, matrix(Wmat)
                restore
                mata: Wmat = st_matrix("Wmat")

                tempvar cl_num unit_num
                qui egen long `cl_num'   = group(`cluster_var') if `touse'
                qui egen long `unit_num' = group(`ivar')        if `touse'

                mata: cl_vec   = st_data(., "`cl_num'",   "`touse'")
                mata: unit_vec = st_data(., "`unit_num'", "`touse'")

                qui levelsof `cluster_var' if `touse', local(cl_vals)
                local n_clusters : word count `cl_vals'

                mata: st_numscalar("n_units_sc", max(unit_vec))
                local n_units = n_units_sc
            }

            mata {
                n_vr      = rows(WATT_pmat)
                Nobs_m    = rows(IF_mat)
                n_cells_m = cols(IF_mat)

             * --- Build aggregated IF for each event time r
                IF_r = J(Nobs_m, n_vr, 0)

                for (rv=1; rv<=n_vr; rv++) {
                    r_val = WATT_pmat[rv, 1]

                    for (k=1; k<=n_cells_m; k++) {
                        g_k = cell_g[1, k]
                        t_k = cell_t[1, k]

                        if ((t_k - g_k) != r_val) continue

                        w_k = 0
                        for (m=1; m<=rows(Wmat); m++) {
                            if (Wmat[m,1] == r_val & Wmat[m,2] == g_k) {
                                w_k = Wmat[m,3]
                                break
                            }
                        }

                        IF_r[., rv] = IF_r[., rv] :+ w_k * IF_mat[., k]
                    }
                }

             * --- Center the aggregated IF column by column
                IF_r = IF_r :- J(rows(IF_r), 1, 1) * (colsum(IF_r) / rows(IF_r))

                reps_m = `reps'
                n_cl   = `n_clusters'

             * --- Star bootstrap draws: only the fluctuation part
                BS_star = J(reps_m, n_vr, .)

                for (rep=1; rep<=reps_m; rep++) {
                    xi_cl = ((runiform(n_cl,1):>0.5) :- 0.5) * 2
                    xi_i  = xi_cl[cl_vec]

                    BS_star[rep,.] = colsum(IF_r :* xi_i)
                }

            * --- What gets plotted / reported depends on rolling()
                WATT_plot = WATT_pmat[,2]
                BS_plot   = BS_star

                if ("`rolling'" == "demean") {
                    base_idx = .
                    for (rv=1; rv<=n_vr; rv++) {
                        if (WATT_pmat[rv,1] == -1) {
                            base_idx = rv
                            break
                        }
                    }

                    if (base_idx != .) {
                        WATT_plot = WATT_pmat[,2] :- WATT_pmat[base_idx,2]
                        BS_plot   = BS_star :- BS_star[,base_idx]
                    }
                }
                else if ("`rolling'" == "detrend") {
			* --- r = -1 and r = -2 are anchored at 0; they are NOT estimation or inference targets--> excluded from the sup-t critical value
                    for (rv=1; rv<=n_vr; rv++) {
                        if (WATT_pmat[rv,1] == -1 | WATT_pmat[rv,1] == -2) {
                            WATT_plot[rv]  = 0
                            BS_plot[., rv] = J(rows(BS_plot), 1, 0)
                        }
                    }
                }

            * --- Append plotting/reporting estimand as third column
                WATT_pmat = WATT_pmat, WATT_plot
            }

            local alpha  = (100 - `level') / 100

            mata: n_vr_sc = rows(WATT_pmat)
            mata: st_numscalar("n_vr_sc", n_vr_sc)
            local n_vr = n_vr_sc

            * point estimates and pointwise standard errors
            forvalues col = 1/`n_vr' {
                mata: st_numscalar("rv_sc", WATT_pmat[`col', 1])
                mata: st_numscalar("theta_raw_sc",  WATT_pmat[`col', 2])
                mata: st_numscalar("theta_plot_sc", WATT_pmat[`col', 3])
                mata: st_numscalar("se_plot_sc", sqrt(variance(BS_plot[., `col'])))

                local rv_`col'         = rv_sc
                local theta_raw_`col'  = theta_raw_sc
                local theta_plot_`col' = theta_plot_sc
                local se_plot_`col'    = se_plot_sc
            }

            if "`ci_type'" == "simultaneous" {
                mata {
                    se_vec  = sqrt(diagonal(variance(BS_plot)))  // column vector
                    reps_bs = rows(BS_plot)
                    n_cols  = cols(BS_plot)
                    T_star  = J(reps_bs, 1, .)
                    for (b = 1; b <= reps_bs; b++) {
                        z_b = J(1, n_cols, 0)
                        for (j = 1; j <= n_cols; j++) {
                            if (se_vec[j] > 0 & se_vec[j] < .) {
                                z_b[j] = abs(BS_plot[b,j] / se_vec[j])
                            }
                        }
                        T_star[b] = max(z_b)
                    }
                    T_sort = sort(T_star, 1)
                    nT     = rows(T_sort)
                    q_idx  = min((nT, max((1, ceil((1-`alpha') * nT)))))
                    st_numscalar("c_sup_sc", T_sort[q_idx])
                }
                local c_sup = c_sup_sc

                if "`c_sup'" == "" | "`c_sup'" == "." {
                    di as error "Failed to compute simultaneous critical value."
                    exit 198
                }

                forvalues col = 1/`n_vr' {
                    local lo_ci_plot_`col' = `theta_plot_`col'' - `c_sup' * `se_plot_`col''
                    local hi_ci_plot_`col' = `theta_plot_`col'' + `c_sup' * `se_plot_`col''
                }
            }

                   preserve
            qui use "`WATT_point'", clear
            qui gen double se       = .
            qui gen double lower_ci = .
            qui gen double upper_ci = .
            qui gen double watt_plot     = .
            qui gen double se_plot       = .
            qui gen double lower_ci_plot = .
            qui gen double upper_ci_plot = .
            qui gen double base_rminus1  = .
            qui gen double watt_norm     = .
            qui gen double lower_ci_norm = .
            qui gen double upper_ci_norm = .

            quietly forval col = 1/`n_vr' {
				qui replace watt_plot     = `theta_plot_`col''  if ryear == `rv_`col''
				qui replace se_plot       = `se_plot_`col''     if ryear == `rv_`col''
				qui replace lower_ci_plot = `lo_ci_plot_`col''  if ryear == `rv_`col''
				qui replace upper_ci_plot = `hi_ci_plot_`col''  if ryear == `rv_`col''
			}

			if "`rolling'" == "detrend" {
				qui replace watt_plot     = 0 if inlist(ryear,-1,-2)
				qui replace lower_ci_plot = 0 if inlist(ryear,-1,-2)
				qui replace upper_ci_plot = 0 if inlist(ryear,-1,-2)
				qui replace se_plot       = 0 if inlist(ryear,-1,-2)
			}

            qui replace se       = se_plot
            qui replace lower_ci = lower_ci_plot
            qui replace upper_ci = upper_ci_plot

        if "`rolling'" == "demean" {
            qui su watt if ryear == -1, meanonly
            local base = r(mean)
            qui replace base_rminus1  = `base'
            qui replace watt_norm     = watt_plot
            qui replace lower_ci_norm = lower_ci_plot
            qui replace upper_ci_norm = upper_ci_plot

            * force exact zero at r = -1 for display after normalized inference
            qui replace watt_norm     = 0 if ryear == -1
            qui replace lower_ci_norm = 0 if ryear == -1
            qui replace upper_ci_norm = 0 if ryear == -1

            * finalized results for demean: keep normalized series
            qui replace watt     = watt_norm
            qui replace lower_ci = lower_ci_norm
            qui replace upper_ci = upper_ci_norm
        }
        else {
            * finalized results for non-demean cases: keep plotted series
            qui replace watt     = watt_plot
            qui replace lower_ci = lower_ci_plot
            qui replace upper_ci = upper_ci_plot
        }

        * finalized standard errors
        qui replace se = se_plot

        sort ryear

        format watt se lower_ci upper_ci %9.3f

		di as txt "-> WATT(r) with simultaneous `level'% confidence bands  (bootstrap, B=`reps')"
        di as txt "------------------------------------------------------------"

        * temporary variables for cleaner reporting only
        capture drop low_ci up_ci
        qui gen low_ci = lower_ci
        qui gen up_ci  = upper_ci

        format low_ci up_ci %9.3f

        list ryear watt se low_ci up_ci N_cohort N_units , noobs



* --- Graph
        if "`graph'" != "" {
                if "`title'" == "" local title "lwdid: `method' (`rolling')"
                local base_r = 0
                qui su ryear, meanonly
                local xmin = r(min)
                local xmax = r(max)
                local xrange = `xmax' - `xmin'
                local xstep = cond(`xrange'>40, 10, 5)

                * y-axis range based on finalized confidence intervals
                qui su upper_ci if !missing(upper_ci), meanonly
                local yhi = r(max)
                qui su lower_ci if !missing(lower_ci), meanonly
                local ylo = r(min)

                if (`yhi' - `ylo') < 0.2 {
                    local ymid = (`yhi' + `ylo') / 2
                    local yhi  = `ymid' + 0.1
                    local ylo  = `ymid' - 0.1
                }

                * 5% margin instead of fixed 0.1 to avoid over-expanding y axis
                local ymargin = max(0.02, (`yhi' - `ylo') * 0.05)
                local yhi = `yhi' + `ymargin'
                local ylo = `ylo' - `ymargin'
                local raw_step = (`yhi' - `ylo') / 6
                local ystep = 1
                foreach s in 0.01 0.02 0.05 0.1 0.2 0.25 0.5 1 2 5 10 {
                    if `raw_step' <= `s' {
                        local ystep = `s'
                        continue, break
                    }
                }
                if `ystep' < 0.1 local ystep = 0.1
                local ymin = floor(`ylo'/`ystep')*`ystep'
                local ymax = ceil(`yhi'/`ystep')*`ystep'

                * mono scheme detection
                local mono = strpos("`scheme'","mono")
                if `mono' {
                    local col_pre black%50
                    local col_post black%70
                    local mcol_pre black%50
                    local mcol_post black%70
                }
                else {
                    local col_pre navy
                    local col_post cranberry
                    local mcol_pre navy%80
                    local mcol_post cranberry
                }

                local yttl "WATT(r)"
                if "`rolling'" == "demean" {
                    local yttl "WATT(r), normalized at r = -1"
                }
                else if "`rolling'" == "detrend" {
                    local yttl "WATT(r), relative to pre-treatment trend"
                }

                twoway ///
                    (rcap lower_ci upper_ci ryear if ryear < 0, ///
                        lwidth(0.3) lcolor(`col_pre'%50)) ///
                    (rcap lower_ci upper_ci ryear if ryear >= 0, ///
                        lwidth(0.3) lcolor(`col_post'%60)) ///
                    (line watt ryear, ///
                        lcolor(gs8) lwidth(thin)) ///
                    (scatter watt ryear if ryear < 0, ///
                        mcolor(`mcol_pre') msymbol(circle) msize(medlarge)) ///
                    (scatter watt ryear if ryear >= 0, ///
                        mcolor(`mcol_post') msymbol(circle) msize(medlarge)) ///
                    , ///
                    yline(0, lcolor(gs10)) ///
                    xline(0, lcolor(gs10) lpattern(dash)) ///
                    xtitle("Time to Treatment (r)") ///
                    ytitle("`yttl'") ///
                    title(`"`title'"') ///
                    xlabel(`xmin'(`xstep')`xmax' 0, labsize(small)) ///
                    ylabel(`ymin'(`ystep')`ymax', format(%5.2f)) ///
                    legend(off) ///
                    scheme(`scheme') ///
                    `gopts'
        }
		
* --- save
 if "`save'" != "" {
            qui keep ryear ///
                watt se low_ci up_ci ///
                N_cohort N_units

            qui order ryear ///
                watt se low_ci up_ci ///
                N_cohort N_units
            
            * reset display format before saving
            format watt se low_ci up_ci %9.0g

            qui save "`save'", replace
        }		
			restore
						
						
			quietly {
				mata: mata drop IF_mat IF_r IF_col cell_g cell_t cl_vec unit_vec Wmat WATT_pmat BS_star BS_plot WATT_plot
				mata: mata drop n_vr n_vr_sc Nobs_m n_cells_m
				capture mata: mata drop T_star T_sort se_vec
			}
		}
	end


**# RI
		capture mata: mata drop lwdid_ri_inline()

		mata:
		real scalar lwdid_ri_inline(
			real scalar reps,
			real scalar b0,
			string scalar rhs
		)
		{
			real colvector fp, Y, D, Dp, res
			real matrix X, Z, bhat
			real scalar n_fp, r, i, j, tmp
			string rowvector xvars

			fp   = selectindex(st_data(., "firstpost") :== 1)
			n_fp = rows(fp)

			Y = st_data(fp, "ydot_postavg")
			D = st_data(fp, "d_")

			xvars = tokens(rhs)
			if (length(xvars) > 0) {
				X = J(n_fp,1,1), st_data(fp, xvars)
			}
			else {
				X = J(n_fp,1,1)
			}

			res = J(reps,1,.)

			for (r = 1; r <= reps; r++) {
				Dp = D

				for (i = n_fp; i >= 2; i--) {
					j     = ceil(i * runiform(1,1))
					tmp   = Dp[i]
					Dp[i] = Dp[j]
					Dp[j] = tmp
				}

				Z    = X, Dp
				bhat = invsym(quadcross(Z,Z)) * quadcross(Z,Y)
				res[r] = bhat[rows(bhat),1]
			}

			return( (sum(abs(res) :>= abs(b0)) + 1) / (rows(res) + 1) )
		}
		end


