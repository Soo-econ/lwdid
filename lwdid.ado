*
*! version 3.6: 8 November 2026
*! lwdid - Lee & Wooldridge rolling DID estimator (unified: small-N and large-N panels)
*! authors: Soo Jeong Lee, Jeffrey M. Wooldridge
*! contact: soojeong.lee@siu.edu, wooldri1@msu.edu
*! lwdid GitHub repository: https://github.com/Soo-econ/lwdid


**# [1] MAIN PROGRAM: lwdid 
***>> Dispatches to the appropriate subroutine: lwdid_small_single/lwdid_small_staggered/lwdid_large
program define lwdid, eclass sortpreserve
		version 16.0

		syntax varlist(min=1 numeric fv) [if] [in], ///
			IVAR(name) TVAR(name) GVAR(name)     ///
			ROLLING(name)                        ///
			[METHOD(name)                        ///   
			 Small                               ///    
			 REPS(integer 999)                   ///
			 LEVEL(integer 95)                   ///
			 CLUSTER(name)                       ///
			 SEED(integer 12345)                 ///
			 VCE(string)                         ///
			 TABLE(string)                       ///
			 GRAPH                               ///
			 CB                                  ///
			 SCHEME(string)						 ///
			 GOPTS(string asis)                  ///
			 SAVE(string)                        ///
			 GID(string)                         ///
			 RI                                  ///
			 RIREPS(integer 999)                 ///
			 RISEED(string)                      ///
			 ATTGT                               ///
			 PRE(integer -1)                     ///
			 NEVER                               ///
			 YDOT                                ///
			 CELLIF                              ///
			]

		marksample touse, novarlist
		
		* Validate that gvar() is numeric.
		confirm numeric variable `gvar'

		* Negative cohort dates are not supported (0/missing = never treated).
		capture assert `gvar' >= 0 if `touse' & !missing(`gvar')
		if _rc {
			di as err "gvar() cannot be negative: use positive first-treatment periods (0 or missing = never-treated)."
			exit 198
		}

		* Use an internal cohort variable; preserve the user's gvar().
		local __gvar_original `gvar'
		tempvar __gvar_internal
		quietly gen double `__gvar_internal' = `gvar' if `touse'
		quietly replace `__gvar_internal' = 0 if missing(`__gvar_internal') & `touse'
		local gvar `__gvar_internal'
		
		capture assert `gvar' == floor(`gvar') if `touse' & `gvar' > 0
		if _rc {
			di as err "gvar() must be integer-valued."
			exit 198
		}
		
		*-- parse varlist
		local y    : word 1 of `varlist'
		local xlist: list varlist - y       
		* Expand factor-variable controls once, before dispatch. Materialized
		* columns also feed the Mata influence functions (which need real vars).
		* fvrevar respects i., c., # and ## and removes base categories.
		if "`xlist'" != "" {
			quietly fvrevar `xlist' if `touse'
			local __expanded_x `r(varlist)'
			* fvrevar returns physical variables usable by Stata and Mata.
		}

		*-- Covariates must be time-invariant: warn if any varies within units.
		if "`__expanded_x'" != "" {
			foreach __v of local __expanded_x {
				tempvar __sdv
				quietly bysort `ivar': egen double `__sdv' = sd(`__v') if `touse'
				quietly count if `touse' & `__sdv' > 1e-10 & !missing(`__sdv')
				if r(N) > 0 {
					di as txt "Warning: covariate {bf:`__v'} varies within units. " ///
					   "lwdid assumes time-invariant covariates; consider using " ///
					   "time-invariant pre-treatment (baseline) values."
				}
				drop `__sdv'
			}
		}
		local __dispatch_varlist `y' `__expanded_x'
		* Forward validated original options and the expanded numeric controls.
		local __options_pos = strpos(`"`0'"', ",")
		local __original_opts = substr(`"`0'"', `__options_pos', .)
		* Forward the internal cohort variable to each estimation routine.
		local __original_opts : subinstr local __original_opts "gvar(`__gvar_original')" "gvar(`__gvar_internal')", all
		local __original_opts : subinstr local __original_opts "GVAR(`__gvar_original')" "gvar(`__gvar_internal')", all
		local __dispatch_call `"`__dispatch_varlist' `if' `in' `__original_opts'"'

			
		*-- rolling() check
		*   The same six transformations are available in both small-N and large-N modes.
		local rolling = lower("`rolling'")
		local ok = inlist("`rolling'","demean","detrend","demeanq","detrendq","demeanm","detrendm")
		if !`ok' {
			di as err "rolling() must be: demean | detrend | demeanq | detrendq | demeanm | detrendm"
			exit 198
		}

	
		*-- validate pre() option
		if `pre' != -1 {
			if `pre' < 1 {
				di as err "pre() must be a positive integer (number of pre-treatment periods to use)."
				exit 198
			}
			local rolling_lo = lower("`rolling'")
			if "`rolling_lo'" == "detrend" & `pre' < 2 {
				di as err "pre() must be at least 2 when rolling(detrend) is specified (need >= 2 periods to estimate slope)."
				exit 198
			}
			if "`rolling_lo'" == "detrendq" & `pre' < 5 {
				di as err "pre() must be at least 5 when rolling(detrendq) is specified (need within-quarter time variation to identify the trend)."
				exit 198
			}
			if "`rolling_lo'" == "detrendm" & `pre' < 13 {
				di as err "pre() must be at least 13 when rolling(detrendm) is specified (need within-month time variation to identify the trend)."
				exit 198
			}
		}

	*-- DISPATCH (small option)
		if "`small'" != "" {

                * Small-N uses additive OLS; method() is a large-N selector.
                if "`method'" != "" {
                    di as err "Please omit the method() option when using small; small-N uses additive OLS."
                    exit 198
                }

			    * cellif is a large-N inference option
				if "`cellif'" != "" {
					di as err "cellif is available only in large-N mode."
					exit 198
				}
			    * cluster() is not supported in small-N mode
				if "`cluster'" != "" {
					di as err "cluster() is available only in large-N mode."
					exit 198
				}
				
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
					lwdid_small_single `__dispatch_call'
					* Store the control list as typed by the user (e.g. i.region c.x1#c.x2).
					ereturn local controls "`xlist'"
					exit
				}
				*with staggered adoption
				else {
					lwdid_small_staggered `__dispatch_call'
                    * Preserve original, user-specified control notation.
                    ereturn local controls "`xlist'"
					exit
				}
			}
				*with Large-N
		else {
				lwdid_large `__dispatch_call'
				* Store the control list as typed by the user (e.g. i.region c.x1#c.x2).
				ereturn local controls "`xlist'"
			}
			
			
end

**# [2] lwdid_small_single
**>> Subroutine for the small-N common timing case (single treated unit)
		
program define lwdid_small_single, eclass
		version 16.0

		syntax varlist(min=1 numeric fv) [if] [in], ///
			IVAR(name) TVAR(name)  GVAR(name)     ///
			ROLLING(name)                        ///
			[METHOD(name)                        ///   
			 Small                               ///    
			 REPS(integer 999)                   ///
			 LEVEL(integer 95)                   ///
			 CLUSTER(name)                       ///
			 SEED(integer 12345)                 ///
			 VCE(string)                         ///
			 TABLE(string)                       ///
			 GRAPH                               ///
			 POINTWISE                           ///
			 CB                                  ///
			 SCHEME(string)						 ///
			 GOPTS(string asis)                  ///
			 SAVE(string)                        ///
			 TITLE(string)                       ///
			 GID(string)                         ///
			 RI                                  ///
			 RIREPS(integer 999)                 ///
			 RISEED(string)                      ///
			 ATTGT                               ///
			 YDOT                                ///
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
		qui if missing(`riseed') local riseed = 12345

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

		* treatment indicator is created only inside preserve below
		* to avoid modifying or overwriting variables in the user's data.

		* post indicator
		qui gen byte `post_' = (`timevar' >= `gyear') if `touse'
		qui replace `post_' = . if !`touse'
		di as txt "------------------------------------------------------------"
		di as txt "gvar(): single cohort detected -> common timing"
		di as txt "lwdid [small-N mode]  rolling=`rolling'  method=ols"
		di as txt "------------------------------------------------------------"



	* --- Keep analysis sample & create core flags
		tempvar tindex yhat ydot cs ctrlSum ydot_tr
		tempname b V
		
		preserve
		qui keep if `touse'
		qui capture drop d_
		qui gen byte d_ = (`gvar' > 0) if `touse'
		qui replace d_ = 0 if missing(d_)
		qui drop if missing(`y', `id', `post_', d_)

	* ---  (2) Build time index: yearly or year-quarter, plus nice tq label
		quietly su `timevar', meanonly
        local __starttime = r(min)
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


    * ---  For each unit, construct yhat for ALL periods (based on selected rolling)
		qui gen double `yhat' = .

        * rolling(demean): use the mean of the available nonmissing pre-treatment
        * outcomes for each unit directly. This supports unbalanced panels and
        * allows a unit to contribute with only one observed pre-treatment outcome.
        if "`rolling'"=="demean" {
            tempvar __premean
            quietly bysort `id': egen double `__premean' = ///
                mean(cond(`post_'==0, `y', .))
            quietly replace `yhat' = `__premean' if !missing(`__premean')
        }
        else {
            qui levelsof `id', local(IDlist)
            tempvar __fit

            quietly foreach ii of local IDlist {
                capture drop `__fit'
                if "`rolling'"=="detrend" {
                    regress `y' c.`tindex' if `id'==`ii' & `post_'==0
                    predict double `__fit' if `id'==`ii', xb
                }
                else if "`rolling'"=="demeanq" {
                    regress `y' i.`__qvar' if `id'==`ii' & `post_'==0
                    predict double `__fit' if `id'==`ii', xb
                }
                else if "`rolling'"=="detrendq" {
                    regress `y' c.`tindex' i.`__qvar' if `id'==`ii' & `post_'==0
                    predict double `__fit' if `id'==`ii', xb
                }
                else if "`rolling'"=="demeanm" {
                    regress `y' i.`__mvar' if `id'==`ii' & `post_'==0
                    predict double `__fit' if `id'==`ii', xb
                }
                else if "`rolling'"=="detrendm" {
                    regress `y' c.`tindex' i.`__mvar' if `id'==`ii' & `post_'==0
                    predict double `__fit' if `id'==`ii', xb
                }
                qui replace `yhat' = `__fit' if `id'==`ii'
                capture quietly drop `__fit'
            }
        }



	* ---  Residualized outcome for ALL periods: ydot = y - yhat
		* Always also create internal ydot for the subsequent analysis
			quietly gen double `ydot' = `y' - `yhat'
			qui label var `ydot' "y - yhat (residualized via pre-period fit)"


	* ---  Collapse to first post cross-section & Aggregated single treatment effect
		qui capture drop ydot_postavg firstpost
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

		* --- additive covariates only (no D x X interactions; time-invariant X) ---
		local RHS

		local K1 0
		if "`xlist'" != "" {
			local K1 : word count `xlist'
		}

		quietly count if `tindex'==1 & d_==1
		local nt = r(N)
		quietly count if `tindex'==1 & d_==0
		local nc = r(N)

		* HC-type VCE requires both groups to contain more than one unit.
		* In a two-group regression a singleton group has leverage one.
		if "`vce'" != "" & (`nt' < 2 | `nc' < 2) {
			di as err "vce(`vce') requires at least 2 treated and 2 control units (treated N=`nt', control N=`nc')."
			exit 2001
		}

		* Parameters: constant + d_ + K controls = K+2. Need N_tr + N_unt > K+2.
		if "`xlist'" != "" {
			if (`nt' + `nc') > (`K1' + 2) & `nt' >= 1 & `nc' >= 1 {
				local RHS "`xlist'"
			}
			else {
				di as err "Too many controls for small-N estimation: K=`K1', treated N=`nt', controls N=`nc'."
				di as err "Need N_tr + N_unt > K+2. Reduce covariates."
				exit 2001
			}
		}
		else {
			local RHS ""
		}

			
		* --- Main regression & optional randomization inference

				* --- [1] Single ATT regression ---
				di as txt _n "[1] Single ATT"
				if "`vce'" != "" {
					di as txt "    Variance adjustment: `vce'"
					regress ydot_postavg d_ `RHS' if firstpost, vce(`vce')
				}
				else {
					di as txt "    Variance adjustment: none"
					regress ydot_postavg d_ `RHS' if firstpost
				}

                * Warn only if otherwise eligible common-timing units have missing X.
                if "`xlist'" != "" {
                    tempvar __xeligible_ct
                    quietly gen byte `__xeligible_ct' = firstpost & !missing(ydot_postavg)
                    quietly count if `__xeligible_ct'
                    local __prex_ct = r(N)
                    quietly markout `__xeligible_ct' `xlist'
                    quietly count if `__xeligible_ct'
                    local __dropx_ct = `__prex_ct' - r(N)
                    if `__dropx_ct' > 0 {
                        di as txt "Warning: `__dropx_ct' unit(s) excluded from the common-timing" ///
                            " small-N regression due to missing covariates."
                        di as txt "         Estimates may differ because the regression sample changed."
                    }
                }

                * Save the exact OLS estimation sample, N and residual df
                * before period-specific regressions overwrite e().
                tempvar __overall_es __overall_tag
                local __overall_df_r = e(df_r)
                local __overall_N_reg = e(N)
                * Pre-check counts units before missing X are dropped; confirm the
                * actual regression still has residual degrees of freedom.
                if "`xlist'" != "" & (missing(real("`__overall_df_r'")) | real("`__overall_df_r'") < 1) {
                    di as err "No residual degrees of freedom in the small-N regression (N=`__overall_N_reg')."
                    di as err "Covariates with missing values reduced the sample; reduce covariates."
                    exit 2001
                }
                quietly generate byte `__overall_es' = e(sample)
                quietly egen byte `__overall_tag' = tag(`id') if `__overall_es'
                quietly count if `__overall_tag' == 1 & d_ == 1
                local __overall_nt = r(N)
                quietly count if `__overall_tag' == 1 & d_ == 0
                local __overall_nc = r(N)
                local __overall_n = `__overall_nt' + `__overall_nc'
                if "`vce'" != "" & (`__overall_nt' < 2 | `__overall_nc' < 2) {
                    di as err "vce(`vce') needs at least 2 treated and 2 control units in the actual regression sample."
                    exit 2001
                }

				* keep these for e()
				matrix `b' = e(b)
				matrix `V' = e(V)
				scalar __att_overall = _b[d_]
				scalar __se_overall  = _se[d_]

				* optional export
				if `__do_export' & `__has_or2' {
					if "`vce'" != "" {
						outreg2 using "`__xls'", append label nose ///
							ctitle("`rolling' (OLS) with `vce'") addstat(N, e(N))
					}
					else {
						outreg2 using "`__xls'", replace label nose ///
							ctitle("`rolling' (OLS)") addstat(N, e(N))
					}
				}

				* --- Randomization inference ---
				if "`ri'" != "" {
					di as txt _n "    Randomization inference (RI)"
					di as txt "    Replications: `rireps'"

					quietly set rng mt64
					if "`riseed'" != "" {
						quietly set seed `riseed'
						di as txt "    Seed: `riseed'"
					}

					local reps = `rireps'

					scalar __b0 = __att_overall
					* Permute exactly the units used by the original OLS.
					* RHS contains materialized numeric covariates from fvrevar.
					mata: st_numscalar("__p_ri", ///
						lwdid_ri_inline(`reps', st_numscalar("__b0"), "`RHS'", "`__overall_es'"))

					di as txt "    RI p-value: " as res %9.3f __p_ri
				}

		* keep these for e()
		matrix `b' = e(b)
		matrix `V' = e(V)
		scalar __att_overall = _b[d_]
		scalar __se_overall  = _se[d_]


		* ---[2] Period-by-period post effects table

		tempname __pph
		tempfile __ppfile
		qui postfile `__pph' int tindex str20 period double beta se tstat pval N ///
            long N_tr N_unt N_total using `__ppfile', replace

		* first row = overall 'average' from the single-effect above
		local b_avg = _b[d_]
		local s_avg = _se[d_]
		local t_avg = `b_avg'/`s_avg'
		local p_avg = 2*ttail(e(df_r), abs(`t_avg'))
		local N_avg = e(N)
		post `__pph' (.) ("average") (`b_avg') (`s_avg') (`t_avg') (`p_avg') (`N_avg') (.) (.) (.)

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
                local __Nreg = e(N)
                * Count distinct units in the exact regression estimation sample.
                tempvar __sampletag
                quietly egen byte `__sampletag' = tag(`id') if e(sample)
                quietly count if `__sampletag' == 1 & d_ == 1
                local __ntreated = r(N)
                quietly count if `__sampletag' == 1 & d_ == 0
                local __ncontrol = r(N)
                local __ntotal = `__ntreated' + `__ncontrol'
                drop `__sampletag'
                if "`vce'" != "" & (`__ntreated' < 2 | `__ncontrol' < 2) {
                    * Keep the point estimate; its HC-type standard error is not trustworthy.
                    local se1 = .
                    local t = .
                    local p = .
                }
				post `__pph' (`tt') ("`lab'") (`coef') (`se1') (`t') (`p') (`__Nreg') ///
                    (`__ntreated') (`__ncontrol') (`__ntotal')
			}
			else {
				post `__pph' (`tt') ("`lab'") (.) (.) (.) (.) (.) (.) (.) (.)
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

            * Collect period-specific estimates and unique sample counts before
            * converting the numeric time index to a display string.
            tempname __small_attgt
            preserve
                quietly keep if period != "average"
                quietly count
                if r(N) > 0 {
                    quietly gen double cohort = `gyear'
                    quietly gen double time = tindex + `__starttime' - 1
                    quietly gen double ryear = time - cohort
                    quietly rename beta att
                    quietly mkmat cohort time ryear att se N_tr N_unt N_total, matrix(`__small_attgt')
                }
            restore

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
			qui scalar cr = invttail(e(df_r), (100-`level')/200)
			qui gen double ci_lw = beta - cr*se
			qui gen double ci_up = beta + cr*se
			order period tindex beta se ci_lw ci_up tstat pval N
			
            * Default: show up to 10 rows; attgt: show all rows.
            {
                di as txt _n "[2] Period-by-period ATTs"
                * First row is the "average" row; show it plus up to 10 post periods.
                local __maxrows = 11
                quietly count
                local __nrows = r(N)
                if "`attgt'" != "" local __maxrows = `__nrows'

                format beta se ci_lw ci_up %9.0g
                format tstat %7.2f
                format pval %6.3f

                di as txt "{hline 72}"
                di as txt _col(3) "period" ///
                    _col(15) "ATT" ///
                    _col(25) "Std. err." ///
                    _col(36) "t" ///
                    _col(43) "P>|t|" ///
                    _col(53) "[`level'% conf. int.]"
                di as txt "{hline 72}"
                local __shown = min(`__nrows', `__maxrows')
                forvalues __i = 1/`__shown' {
                    di as txt _col(2) %8s period[`__i'] ///
                        as res _col(13) %9.0g beta[`__i'] ///
                        as res _col(23) %9.0g se[`__i'] ///
                        as res _col(33) %7.2f tstat[`__i'] ///
                        as res _col(42) %6.3f pval[`__i'] ///
                        as res _col(51) %9.0g ci_lw[`__i'] ///
                        as res _col(62) %9.0g ci_up[`__i']
                }
                di as txt "{hline 72}"
                di as txt "Note: Confidence intervals, t statistics, and p-values are pointwise t-based."
                if `__nrows' > `__maxrows' {
                    local __nper = `__nrows' - 1
                    di as error "Note: Only the first 10 of `__nper' post-treatment periods are displayed."
                    di as txt "Use attgt to display all rows, or matrix list e(attgt) for the full post-period matrix."
                }
            }

			if "`save'" != "" {

			qui keep period beta se ci_lw ci_up tstat pval N
			qui rename beta att
			order period att se ci_lw ci_up tstat pval N

			* reset display format before saving
			format att se ci_lw ci_up %9.0g
			format tstat %7.2f
			format pval %6.3f
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

            * Period-specific results are already collected into __small_attgt.
            frame change default

			capture frame drop __lwdid_pp
		}

	* --- Control mean trajectory & optional graph
		* Always compute control-group mean over time (d==0)
			qui capture drop y_dot_cont y_dot_tr
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

			local yfmt "%9.2f"
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
					pos(6) ring(1) col(2) ///
					  size(small) ///
				region(fcolor(none) lcolor(none))) ///
				graphregion(color(white)) ///
				plotregion(color(white)) ///
				`gsch' ///
				`gopts' ///
				`aheight'
		}
			
			
   * Store only the original overall treatment coefficient and its variance.
        * The original regression and its inference above remain unchanged.
        tempname __small_b __small_V
        matrix `__small_b' = `b'[1,"d_"]
        matrix `__small_V' = `V'["d_","d_"]
        matrix colnames `__small_b' = Post_avg
        matrix colnames `__small_V' = Post_avg
        matrix rownames `__small_V' = Post_avg
        ereturn post `__small_b' `__small_V', esample(`__overall_es') obs(`__overall_N_reg') dof(`__overall_df_r') depname(`y')
        ereturn scalar N_tr = `__overall_nt'
        ereturn scalar N_unt = `__overall_nc'
        ereturn scalar N_total = `__overall_n'
        ereturn local cmd        "lwdid"
        ereturn local depvar     "`y'"
        ereturn local rolling    "`rolling'"
        ereturn local controls   "`controls'"
        ereturn local gid  "`gid'"

        ereturn scalar K         = `K'
        ereturn scalar tpost1    = `tpost1'
		if "`ri'" != "" ereturn scalar p_ri = __p_ri
		

        * Common-timing small-N: post ATT-by-period and exact unit counts.
        capture confirm matrix `__small_attgt'
        if !_rc {
            matrix colnames `__small_attgt' = cohort time ryear att se N_tr N_unt N_total
            ereturn matrix attgt = `__small_attgt'
        }

	restore

    * Small-N common timing: display the aggregated single-effect summary
    * after any requested period-by-period ATT table and graph output.
    * Display-only; estimation results and e() matrices are unchanged.
    di as res _n "*--- lwdid (rolling: `rolling'): Aggregated single effect = " ///
        %6.3f _b[Post_avg] " (SE = " %5.3f _se[Post_avg] ")"
end


**# [3] lwdid_small_staggered
**>> Subroutine for the small-N staggered adoption (multiple treated units)
program define lwdid_small_staggered, eclass
    version 16.0

    syntax varlist(min=1 numeric fv) [if] [in], ///
        IVAR(name) TVAR(name) GVAR(name)     ///
        ROLLING(name)                        ///
        [METHOD(name)                        ///
         Small                               ///
         REPS(integer 999)                   ///
         LEVEL(integer 95)                   ///
         CLUSTER(name)                       ///
         SEED(integer 12345)                 ///
         VCE(string)                         ///
         TABLE(string)                       ///
         GRAPH                               ///
         POINTWISE                           ///
         CB                                  ///
         SCHEME(string)                      ///
         GOPTS(string asis)                  ///
         SAVE(string)                        ///
         GID(string)                         ///
         RI                                  ///
         RIREPS(integer 999)                 ///
         RISEED(string)                      ///
         ATTGT                               ///
         YDOT                                ///
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

        * Apply the same VCE to all small-N staggered OLS regressions.
        local __small_vceopt ""
        if "`vce'" != "" local __small_vceopt ", vce(`vce')"

		* Built-in RI supports common timing only; see help for staggered ritest.
		if "`ri'" != "" {
			di as err "ri is available only for small-N common timing; for staggered adoption, see the ritest example in help lwdid."
			exit 198
		}

		*-- controls enter additively (time-invariant X; no D x X interactions)


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

		*-- additive controls: need N_treat + N_never > K+2 in the single-effect regression
		if "`xlist'" != "" {
			local K1 : word count `xlist'
			tempvar __xidtag
			quietly egen byte `__xidtag' = tag(`id') if `touse'
			quietly count if `__xidtag' & `gvar' > 0 & !missing(`gvar')
			local __xnt = r(N)
			quietly count if `__xidtag' & `gvar' == 0
			local __xnc = r(N)
			drop `__xidtag'
			if (`__xnt' + `__xnc') <= (`K1' + 2) | `__xnt' < 1 | `__xnc' < 1 {
				di as err "Too many controls for small-N estimation: K=`K1', treated N=`__xnt', never-treated N=`__xnc'."
				di as err "Need N_tr + N_unt > K+2. Reduce covariates."
				exit 2001
			}
		}
        * vce(): a single-unit group has leverage 1, so HC-type SEs are undefined.
        local __vskip 0
        if "`vce'" != "" {
            tempvar __vidtag
            quietly egen byte `__vidtag' = tag(`id') if `touse'
            quietly count if `__vidtag' & `gvar' == 0
            local __vnc = r(N)
            if `__vnc' < 2 {
                di as err "vce(`vce') requires at least 2 never-treated units."
                exit 2001
            }
            quietly levelsof `gvar' if `gvar' > 0 & `touse', local(__vG)
            foreach __vg of local __vG {
                quietly count if `__vidtag' & `gvar' == `__vg'
            if r(N) < 2 {
                    di as err "vce(`vce') requires at least 2 units in each cohort (cohort `__vg' has 1); omit vce()."
                    exit 2001
                }
            }
            drop `__vidtag'
        }

		di as txt "------------------------------------------------------------"
		di as txt "gvar(): `n_cohort' cohorts detected -> Staggered adoption"
		di as txt "lwdid [small-N mode] rolling=`rolling'"
		di as txt "Comparison group: never-treated units only"
		di as txt "------------------------------------------------------------"

		*-- graph support for small-N staggered designs
		*-- The single-effect estimator and regression table below are unchanged.
		tempfile __graph_cells __cell_results __cohort_results
        tempname __gpost __cellpost __cellmatrix __cohortpost __cohortmatrix
        quietly postfile `__cellpost' double cohort time ryear att se N_tr N_unt N_total df_r using "`__cell_results'", replace
        quietly postfile `__cohortpost' double cohort att se N_tr N_unt N_total df_r using "`__cohort_results'", replace

  quietly {
        preserve
            keep if `touse'
            capture drop d_
            gen byte d_ = (`gvar' != 0) if `touse'
            qui replace d_ = 0 if missing(d_)
            label var d_ "ever-treated indicator"
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

            *-- Graph cells are constructed AFTER the overall single-effect OLS,
            *-- using its actual estimation units and its common additive slope.
            *-- The graph remains based on never-treated controls only.
            levelsof `timevar', local(Tlist)
            local __graph_ydots

            *-- cohort weights: computed after the loop from treated UNITS in the
            *-- final (7.18) estimation sample, so w_g = N_g / N_treat holds exactly
            *-- in unbalanced panels (see below).
            local __colist

            *-- build ydot_bar (equation 7.18) without cohort-named variables
            *-- Avoid cohort-specific local names that require nested macro expansion.
            tempvar ydot_bar contsum
            gen double `ydot_bar' = .
            gen double `contsum' = 0

            local i = 0
            foreach g of local Glist {
                local ++i

                tempvar yhat ydotg ybar_tr ybar_co
                gen double `yhat' = .

                * rolling(demean): for each cohort g, use each unit's available
                * nonmissing outcomes from periods t < g. This supports unbalanced
                * panels and permits a unit with only one observed pre-treatment
                * outcome. Units with no observed pre-treatment outcome for cohort g
                * retain missing yhat and are naturally excluded from that comparison.
                if "`rolling'" == "demean" {
                    tempvar __premean_g
                    bysort `id': egen double `__premean_g' = ///
                        mean(cond(`timevar' < `g', `y', .))
                    replace `yhat' = `__premean_g' if !missing(`__premean_g')
                }
                else {
                    levelsof `id', local(IDlist)
                    tempvar __fit

                    foreach idval of local IDlist {
                        capture quietly {
                            if "`rolling'" == "detrend" {
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
                            capture drop `__fit'
                            predict double `__fit' if `id'==`idval', xb
                            replace `yhat' = `__fit' if `id'==`idval'
                            capture drop `__fit'
                        }
                        else {
                            capture drop `__fit'
                        }
                    }
                }

                gen double `ydotg' = `y' - `yhat'
                label var `ydotg' "(rolling=`rolling') Residualized outcome cohort g=`g'"

                *-- Small-N staggered cell-specific regression, never-treated units only.
                *-- Uses the cohort-g pre-treatment transformation for BOTH groups.
                *-- Each cell regression uses the requested VCE (default: classical OLS).
                tempvar __cell_tr
                gen byte `__cell_tr' = (`gvar' == `g')
                foreach tt of local Tlist {
                    local rr = `tt' - `g'
                    * ATT(g,t) is reported only for post-treatment periods (t >= g).
                    if `rr' < 0 continue
                    capture quietly regress `ydotg' `__cell_tr' `xlist' if ///
                        `timevar' == `tt' & inlist(`gvar', 0, `g') `__small_vceopt'
                    if !_rc {
                        local __reg_n = e(N)
                        local __reg_df = e(df_r)
                        local __b = _b[`__cell_tr']
                        local __se = _se[`__cell_tr']
                        quietly count if e(sample) & `__cell_tr' == 1
                        local __ntr = r(N)
                        quietly count if e(sample) & `__cell_tr' == 0
                        local __nco = r(N)
                        if "`vce'" != "" & (`__ntr' < 2 | `__nco' < 2) {
                            local ++__vskip
                        }
                        else if `__ntr' > 0 & `__nco' > 0 & `__reg_df' > 0 & ///
                            !missing(`__b', `__se') {
                            post `__cellpost' (`g') (`tt') (`rr') (`__b') ///
                                (`__se') (`__ntr') (`__nco') (`__reg_n') (`__reg_df')
                        }
                    }
                }
                drop `__cell_tr'

                *-- Retain this cohort's transformed outcome only when a
                *-- post-estimation event-time graph has been requested.
                if "`graph'" != "" {
                    local __graph_ydots `__graph_ydots' `ydotg'
                }

                *-- treated post average for cohort g
                bysort `id': egen double `ybar_tr' = ///
                    mean(cond(`timevar' >= `g' & d_==1, `ydotg', .))
                replace `ydot_bar' = `ybar_tr' if `gvar' == `g'

                *-- control post average for cohort g (unweighted here;
                *-- cohort weights are applied after the loop)
                tempvar __cog
                bysort `id': egen double `__cog' = ///
                    mean(cond(`timevar' >= `g' & d_==0, `ydotg', .))
                local __colist `__colist' `__cog'

                * Cohort-specific post-average ATT: cohort g vs never-treated,
                * using the same g-specific pre-period transformation for both.
                * One observation per unit; additive X coefficients estimated
                * separately for each cohort. This does NOT alter equation (7.18).
                tempvar __group_post __group_tr __group_row
                gen double `__group_post' = `ybar_tr' if `gvar' == `g'
                replace `__group_post' = `__cog' if `gvar' == 0
                gen byte `__group_tr' = (`gvar' == `g')
                bysort `id' (`timevar'): gen byte `__group_row' = (_n == 1)
                capture quietly regress `__group_post' `__group_tr' `xlist' if ///
                    `__group_row' & inlist(`gvar',0,`g') `__small_vceopt'
                if !_rc {
                    local __gdf = e(df_r)
                    local __gN = e(N)
                    capture local __gb = _b[`__group_tr']
                    if !_rc {
                        capture local __gse = _se[`__group_tr']
                        if !_rc & `__gdf' > 0 & !missing(`__gb', `__gse') & `__gse' > 0 {
                            quietly count if e(sample) & `__group_tr' == 1
                            local __gnt = r(N)
                            quietly count if e(sample) & `__group_tr' == 0
                            local __gnc = r(N)
                            local __gmin = cond("`vce'" != "", 2, 1)
                            if `__gnt' >= `__gmin' & `__gnc' >= `__gmin' {
                                post `__cohortpost' (`g') (`__gb') (`__gse') ///
                                    (`__gnt') (`__gnc') (`__gN') (`__gdf')
                            }
                        }
                    }
                }
                drop `__group_post' `__group_tr' `__group_row'

                drop `yhat' `ybar_tr'
                if "`graph'" == "" drop `ydotg'
            }

            postclose `__cellpost'
            postclose `__cohortpost'

            *-- Unit-level estimation sample for (7.18): one row per unit with
            *-- nonmissing ydot_bar (and nonmissing controls, if any).
            tempvar __u1 __esamp
            sort `id' `timevar'
            by `id': gen byte `__u1' = (_n == 1)
            gen byte `__esamp' = `__u1' & d_==1 & !missing(`ydot_bar')
            if "`xlist'" != "" markout `__esamp' `xlist'

            *-- Cohort weights w_g = N_g / N_treat, counting each treated unit
            *-- once and only if it enters the final regression.
            quietly count if `__esamp'
            local __ntreat_w = r(N)
            if `__ntreat_w' == 0 {
                noisily di as err "No treated unit has a nonmissing transformed outcome" ///
                    " (and covariates) for the aggregated single-effect regression."
                exit 2001
            }
            local __ng : word count `Glist'
            matrix w = J(`__ng', 1, 0)
            local i = 0
            foreach g of local Glist {
                local ++i
                quietly count if `__esamp' & `gvar' == `g'
                scalar w_g = r(N) / `__ntreat_w'
                matrix w[`i',1] = scalar(w_g)
                local __cog : word `i' of `__colist'
                if scalar(w_g) > 0 {
                    replace `contsum' = `contsum' + scalar(w_g) * `__cog'
                }
            }

            qui replace `ydot_bar' = `contsum' if d_==0
            drop `contsum' `__colist'

            replace `__esamp' = 1 if `__u1' & d_==0 & !missing(`ydot_bar')
            if "`xlist'" != "" markout `__esamp' `xlist'

            *-- Warn only when nonmissing transformed outcomes were dropped
            *-- from the single-effect regression because of missing X.
            if "`xlist'" != "" {
                tempvar __xeligible
                gen byte `__xeligible' = `__u1' & !missing(`ydot_bar')
                quietly count if `__xeligible'
                local __prex = r(N)
                markout `__xeligible' `xlist'
                quietly count if `__xeligible'
                local __xdropped = `__prex' - r(N)
                if `__xdropped' > 0 {
                    noisily di as txt "Warning: `__xdropped' unit(s) excluded from the aggregated" ///
                        " small-N regression due to missing covariates."
                    noisily di as txt "         Cohort weights use eligible treated units;" ///
                        " estimates can differ due to the changed sample."
                }
            }
    }

		di as txt _n "{bf:[1] Aggregated single treatment effect (Lee & Wooldridge: equation 7.18)}"

		regress `ydot_bar' d_ `xlist' if `__esamp' `__small_vceopt'

		* Pre-check counts units before missing X are dropped; confirm the
		* actual regression still has residual degrees of freedom.
		if "`xlist'" != "" & (missing(e(df_r)) | e(df_r) < 1) {
			di as err "No residual degrees of freedom in the small-N regression (N=`e(N)')."
			di as err "Covariates with missing values reduced the sample; reduce covariates."
			exit 2001
		}

        * Count distinct treated and never-treated units in the actual
        * original overall regression sample (not summed cell counts).
        tempvar __overall_es __overall_tag
        quietly generate byte `__overall_es' = e(sample)
        quietly egen byte `__overall_tag' = tag(`id') if `__overall_es'
        quietly count if `__overall_tag' == 1 & d_ == 1
        local __overall_nt = r(N)
        quietly count if `__overall_tag' == 1 & d_ == 0
        local __overall_nc = r(N)
        local __overall_n = `__overall_nt' + `__overall_nc'
        if "`vce'" != "" & (`__overall_nt' < 2 | `__overall_nc' < 2) {
            di as err "vce(`vce') requires at least 2 treated and 2 never-treated units in the actual regression sample."
            exit 2001
        }

        *-- Event-time paths from the actual single-effect estimation units.
        *-- Adjust the cohort-specific transformed outcomes using the COMMON
        *-- additive X coefficients from the pooled aggregated regression.
        *-- The adjustment uses no treatment-by-covariate interaction.
        if "`graph'" != "" {
            tempvar __graph_unit __graph_adj __graph_xb_unit __graph_xb_all __graph_xpart
            bysort `id': egen byte `__graph_unit' = max(`__overall_es')

            * Factor-variable-safe additive adjustment. Predict the fitted value
            * ONLY on the unit-level estimation row, then carry it to all time
            * observations of the same unit (X is assumed time-invariant).
            * The intercept and D contribution are removed, leaving X'bhat.
            * predict, xb handles i.var, c.var#c.var and other fv notation.
            if "`xlist'" != "" {
                quietly predict double `__graph_xb_unit' if `__overall_es', xb
                bysort `id': egen double `__graph_xb_all' = max(`__graph_xb_unit')
                quietly gen double `__graph_xpart' = `__graph_xb_all' - ///
                    _b[_cons] - _b[d_] * d_ if `__graph_unit'
            }

            quietly postfile `__gpost' double cohort ryear tr_mean co_mean ///
                long Ng N_tr N_co using "`__graph_cells'", replace

            local __gidx = 0
            foreach g of local Glist {
                local ++__gidx
                local __ydotg : word `__gidx' of `__graph_ydots'

                quietly gen double `__graph_adj' = `__ydotg' if `__graph_unit'
                if "`xlist'" != "" {
                    quietly replace `__graph_adj' = `__graph_adj' - ///
                        `__graph_xpart' if `__graph_unit' & !missing(`__graph_adj')
                }

                * Cohort sizes reflect the actual treated units in (7.18).
                quietly count if `__overall_es' & `gvar' == `g'
                local Ng = r(N)
                if `Ng' > 0 {
                    foreach tt of local Tlist {
                        local rr = `tt' - `g'
                        quietly summarize `__graph_adj' if `__graph_unit' & ///
                            `timevar' == `tt' & `gvar' == `g', meanonly
                        local trm = r(mean)
                        local ntr = r(N)
                        quietly summarize `__graph_adj' if `__graph_unit' & ///
                            `timevar' == `tt' & `gvar' == 0, meanonly
                        local com = r(mean)
                        local nco = r(N)
                        if (`ntr' > 0 & `nco' > 0 & ///
                            !missing(`trm', `com')) {
                            post `__gpost' (`g') (`rr') (`trm') (`com') ///
                                (`Ng') (`ntr') (`nco')
                        }
                    }
                }
                drop `__graph_adj'
            }
            postclose `__gpost'
        }

		* Preserve regression results before ereturn post clears e().
		matrix b = e(b)
		matrix V = e(V)
		local N    = e(N)
		local df_r = e(df_r)

		* Repost only the original treatment coefficient and its variance.
		* Preserve the original regression degrees of freedom for t inference.
        tempname __small_b __small_V
        matrix `__small_b' = b[1,"d_"]
        matrix `__small_V' = V["d_","d_"]
        matrix colnames `__small_b' = Post_avg
        matrix colnames `__small_V' = Post_avg
        matrix rownames `__small_V' = Post_avg
        ereturn post `__small_b' `__small_V', obs(`N') dof(`df_r') depname(`y')
        ereturn scalar N_tr = `__overall_nt'
        ereturn scalar N_unt = `__overall_nc'
        ereturn scalar N_total = `__overall_n'

		ereturn local cmd     "lwdid"
		ereturn local depvar  "`y'"
		ereturn local rolling "`rolling'"

        * Always show the original pooled overall ATT before the cohort table.
        di as res _n "*--- lwdid (rolling: `rolling'): Aggregated single effect = " ///
            %6.3f _b[Post_avg] " (SE = " %5.3f _se[Post_avg] ")"



    restore

    * Store cohort-specific post-average ATTs, preserving overall e(b)/e(V).
    preserve
        quietly use "`__cohort_results'", clear
        quietly count
        local __ncohortrows = r(N)
        if `__ncohortrows' > 0 {
            quietly sort cohort
            quietly mkmat cohort att se N_tr N_unt N_total, matrix(`__cohortmatrix')
        }
    restore
    if `__ncohortrows' > 0 {
        matrix colnames `__cohortmatrix' = cohort att se N_tr N_unt N_total
        ereturn matrix cohort_att = `__cohortmatrix'
    }

    di as txt _n "[2] Cohort-specific single effects (post-treatment average)"
    if `__ncohortrows' == 0 {
        di as txt "No estimable cohort-specific single effects."
    }
    else {
        preserve
        quietly use "`__cohort_results'", clear
        quietly sort cohort
        quietly gen double tstat = att/se if se>0 & !missing(se)
        quietly gen double pval = 2*ttail(df_r,abs(tstat)) if df_r>0 & !missing(tstat)
        quietly gen double __crit = invttail(df_r,(100-`level')/200) if df_r>0
        quietly gen double ci_lw = att - __crit*se if !missing(__crit,se)
        quietly gen double ci_up = att + __crit*se if !missing(__crit,se)
        di as txt "{hline 89}"
        di as txt _col(3) "cohort" _col(12) "ATT(avg)" _col(25) "Std. err." ///
            _col(39) "t" _col(48) "P>|t|" _col(62) "[`level'% conf. int.]"
        di as txt "{hline 89}"
        forvalues __j = 1/`__ncohortrows' {
            di as txt _col(2) %8.0g cohort[`__j'] ///
                as res _col(11) %9.0g att[`__j'] ///
                _col(24) %9.0g se[`__j'] ///
                _col(37) %7.2f tstat[`__j'] ///
                _col(47) %6.3f pval[`__j'] ///
                _col(61) %9.0g ci_lw[`__j'] ///
                _col(74) %9.0g ci_up[`__j']
        }
        di as txt "{hline 89}"
        di as txt "Note: Each cohort vs. never-treated (separate additive OLS)."
        restore
    }

    * Store cohort-time estimates in one e(attgt) matrix. Do not alter e(b)/e(V).
    preserve
        quietly use "`__cell_results'", clear
        quietly count
        local __ncells = r(N)
        if `__ncells' > 0 {
            quietly sort cohort time
            quietly mkmat cohort time ryear att se N_tr N_unt N_total, matrix(`__cellmatrix')
        }
    restore
    if `__ncells' > 0 {
        matrix colnames `__cellmatrix' = cohort time ryear att se N_tr N_unt N_total
        ereturn matrix attgt = `__cellmatrix'
    }
    * [3] ATT(g,t) table: displayed only when attgt is specified.
    * e(attgt) always stores all cells; save() works with or without attgt.
    if "`attgt'" != "" {
        di as txt _n "[3] Cohort-by-period ATT(g,t) estimates"
        if `__ncells' == 0 {
            di as txt "No estimable post-treatment cohort-time cells."
        }
    }

    if `__ncells' > 0 & ("`attgt'" != "" | "`save'" != "") {
        * Display and optional export only: do not change the overall e(b)/e(V).
        preserve
        quietly use "`__cell_results'", clear
        quietly sort cohort time
        quietly gen double tstat = att/se if se>0 & !missing(se)
        quietly gen double df = df_r
        quietly gen double pval = 2*ttail(df,abs(tstat)) if df>0 & !missing(tstat)
        quietly gen double __crit = invttail(df,(100-`level')/200) if df>0
        quietly gen double ci_lw = att - __crit*se if !missing(__crit,se)
        quietly gen double ci_up = att + __crit*se if !missing(__crit,se)
        quietly drop __crit

        if "`attgt'" != "" {
            di as txt "{hline 89}"
            di as txt _col(3) "cohort" _col(12) "period" _col(21) "ATT" ///
                _col(34) "Std. err." _col(47) "t" _col(56) "P>|t|" ///
                _col(66) "[`level'% conf. int.]"
            di as txt "{hline 89}"
            forvalues __j = 1/`__ncells' {
                di as txt _col(2) %8.0g cohort[`__j'] ///
                   _col(11) %8.0g time[`__j'] ///
                   as res _col(20) %9.0g att[`__j'] ///
                   _col(33) %9.0g se[`__j'] ///
                   _col(46) %7.2f tstat[`__j'] ///
                   _col(55) %6.3f pval[`__j'] ///
                   _col(65) %9.0g ci_lw[`__j'] ///
                   _col(77) %9.0g ci_up[`__j']
            }
            di as txt "{hline 89}"
            di as txt "Note: Pointwise cell-specific t-based inference; never-treated controls only."
            if `__vskip' > 0 di as txt "Warning: `__vskip' cohort-time cell(s) omitted under vce(`vce') due to singleton groups."
        }

        if "`save'" != "" {
            quietly keep cohort time ryear att se tstat pval ci_lw ci_up N_tr N_unt N_total
            quietly order cohort time ryear att se tstat pval ci_lw ci_up N_tr N_unt N_total
            quietly save "`save'", replace
            di as txt "Cohort-time results saved to " as res "`save'"
        }
        restore
    }
    if `__vskip' > 0 & "`attgt'" == "" {
        di as txt "Warning: `__vskip' ATT(g,t) cell(s) dropped under vce(`vce') (only 1 treated or control unit)."
    }

    * Exposure-time graph: drawn after the [1]-[3] tables so they are not interrupted.
    if "`graph'" != "" {
        preserve
		* Exposure-time graph: cohort-size-weighted residualized paths
		* Runs after the original equation (7.18) regression table.
		* No confidence intervals or pre-trends test are reported.
		if "`graph'" != "" {
			use "`__graph_cells'", clear
			quietly count
			if r(N) == 0 {
				di as err "No valid cohort-by-exposure-time cells are available for the graph."
				restore
				exit 2000
			}

        * In an unbalanced panel, the number of treated observations actually
        * contributing at (g,r) can be below the overall cohort size Ng.
        * Weight each available cohort-period mean by its observed treated
        * count, rather than giving missing treated observations weight.
        * These event-time weights need not equal the fixed cohort weights
        * used in the aggregated single-effect regression (7.18).
        gen double __tr_num = N_tr * tr_mean
        gen double __co_num = N_tr * co_mean
        collapse (sum) __tr_num __co_num weight_sum=N_tr ///
            (count) N_cohorts=cohort, by(ryear)

        gen double y_dot_tr   = __tr_num / weight_sum
        gen double y_dot_cont = __co_num / weight_sum
        gen double gap        = y_dot_tr - y_dot_cont
        sort ryear

        label var ryear      "Exposure time"
        if "`xlist'" == "" {
            label var y_dot_tr   "Observed-treated-count-weighted treated residualized outcome"
            label var y_dot_cont "Observed-treated-count-weighted never-treated residualized outcome"
            label var gap        "Weighted treated-minus-control residualized gap"
        }
        else {
            label var y_dot_tr   "Covariate-adjusted treated residualized outcome"
            label var y_dot_cont "Covariate-adjusted never-treated residualized outcome"
            label var gap        "Covariate-adjusted weighted treated-minus-control gap"
        }
        label var N_cohorts  "Number of contributing treated cohorts"

        quietly summarize ryear, meanonly
        local xmin = floor(r(min))
        local xmax = ceil(r(max))
        local xspan = `xmax' - `xmin'
        if `xspan' <= 12       local xstep = 1
        else if `xspan' <= 30  local xstep = 2
        else if `xspan' <= 60  local xstep = 5
        else                   local xstep = 10

        local mono = 0
        if strpos("`scheme'","mono") local mono = 1
        if `mono' {
            local col_tr black%80
            local col_ct black%60
            local pat_tr solid
            local pat_ct dash
            local gsch "scheme(s1mono)"
        }
        else {
            local col_tr cranberry%90
            local col_ct navy%70
            local pat_tr solid
            local pat_ct dash
            local gsch ""
            if "`scheme'" != "" local gsch "scheme(`scheme')"
        }

        local __graph_ytitle "Residualized outcome"
        if "`xlist'" != "" local __graph_ytitle "Covariate-adjusted residualized outcome"

        twoway ///
            (connected y_dot_tr ryear, msymbol(O) msize(small) ///
                lpattern(`pat_tr') lcolor(`col_tr') mcolor(`col_tr') lwidth(medthick)) ///
            (connected y_dot_cont ryear, msymbol(D) msize(small) ///
                lpattern(`pat_ct') lcolor(`col_ct') mcolor(`col_ct') lwidth(medthick)), ///
            xline(0, lcolor(gs8) lpattern(dash)) ///
            yline(0, lcolor(gs12) lpattern(solid)) ///
            xlabel(`xmin'(`xstep')`xmax', nogrid) ///
            ylabel(, nogrid) ///
            xtitle("Exposure time") ///
            ytitle("`__graph_ytitle'") ///
            legend(order(1 "Eventually treated cohorts (weighted)" 2 "Never-treated units (weighted)") ///
                pos(6) ring(1) col(2) size(small) ///
                region(fcolor(none) lcolor(none))) ///
            graphregion(color(white)) plotregion(color(white)) ///
            `gsch' `gopts'

    }
        restore
    }

end

**# [4] lwdid_large
**>> Subroutine for the LARGE-N (Common&staggered)
program define lwdid_large, eclass
    version 16.0

    syntax varlist(min=1 numeric fv) [if] [in], ///
        IVAR(name) TVAR(name)  GVAR(name)     ///
        ROLLING(name)                        ///
        [METHOD(name)                        ///   
         Small                               ///    
         REPS(integer 999)                   ///
         LEVEL(integer 95)                   ///
         CLUSTER(name)                       ///
         SEED(integer 12345)                 ///
         VCE(string)                         ///
         TABLE(string)                       ///
         GRAPH                               ///
         POINTWISE                           ///
         CB                                  ///
		 SCHEME(string)						 ///
         GOPTS(string asis)                  ///
         SAVE(string)                        ///
         TITLE(string)                       ///
         GID(string)                         ///
         RI                                  ///
         RIREPS(integer 999)                 ///
         RISEED(string)                      ///
         ATTGT                               ///
         PRE(integer -1)                     ///
         NEVER                               ///
         YDOT                                ///
         CELLIF                              ///
        ]

		marksample touse, novarlist
		local y    : word 1 of `varlist'
		local __cellif = ("`cellif'" != "")
		local xlist: list varlist - y
		local rolling = lower("`rolling'")
		local __is_demean  = inlist("`rolling'","demean","demeanq","demeanm")
		local __is_detrend = inlist("`rolling'","detrend","detrendq","detrendm")

		*-- seasonal date validation and internal season index for large-N
		local __tfmt : format `tvar'
		local __gfmt : format `gvar'
		tempvar __qvar __mvar
		if inlist("`rolling'","demeanq","detrendq") {
			if strpos("`__tfmt'","%tq")==0 | strpos("`__gfmt'","%tq")==0 {
				di as err "rolling(`rolling') requires tvar() and gvar() to be Stata quarterly date variables (format %tq)."
				exit 198
			}
			quietly gen byte `__qvar' = quarter(dofq(`tvar')) if `touse'
		}
		if inlist("`rolling'","demeanm","detrendm") {
			if strpos("`__tfmt'","%tm")==0 | strpos("`__gfmt'","%tm")==0 {
				di as err "rolling(`rolling') requires tvar() and gvar() to be Stata monthly date variables (format %tm)."
				exit 198
			}
			quietly gen byte `__mvar' = month(dofm(`tvar')) if `touse'
		}

		local method  = lower("`method'") 
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
			
			*--- Large-N inference
			*    Default display: pointwise confidence intervals.
			*    cb option: simultaneous confidence bands across ALL displayed event times
		local ci_type "simultaneous"
		local interval_display "pointwise"
		if "`cb'" != "" {
			local interval_display "simultaneous"
		}
			
			
			* --- Internal compact numeric id for large-N computations.
			*     Always regroup ivar(), even when it is already numeric, to avoid
			*     storing long/high-dimensional original identifiers in local macros.
			local id_orig `ivar'
			quietly count if `touse' & missing(`ivar')
			if r(N) > 0 {
				di as err "ivar() contains missing values in the estimation sample."
				exit 459
			}
			tempvar __id __xttime
			quietly egen long `__id' = group(`ivar') if `touse'
		* --- Create compact numeric IDs and handle excluded observations for xtset.
			quietly replace `__id' = _N + _n if !`touse'
			quietly gen double `__xttime' = `tvar' if `touse'
			quietly replace `__xttime' = _n if !`touse'
			local id `__id'
				
		* --- Graph-scheme option
			if "`scheme'" == "" {
					local scheme s2color
				}

		* --- LARGE-N PATH  (wild bootstrap, staggered, influence-function)
			di as txt "------------------------------------------------------------"
			di as txt " lwdid [large-N mode]  rolling=`rolling'  method=`method'"
			if `pre' > 0 {
				di as txt " pre(`pre'): using last `pre' pre-treatment period(s) for outcome transformation"
			}
			if "`never'" != "" {
				di as txt " never: comparing each cohort g to never-treated only"
			}
			di as txt "------------------------------------------------------------"


			* --- cluster (default = internal id).
			*     If cluster() is supplied, create a compact internal cluster id as well.
			if "`cluster'" == "" {
				local cluster_var `id'
			}
			else {
				capture confirm variable `cluster'
				if _rc {
					di as err "cluster() variable `cluster' not found."
					exit 111
				}
				quietly count if `touse' & missing(`cluster')
				if r(N) > 0 {
					di as err "cluster() contains missing values in the estimation sample."
					exit 459
				}
				tempvar __cluster
				quietly egen long `__cluster' = group(`cluster') if `touse'
				local cluster_var `__cluster'
			}
			
			* Compact numeric cluster identifier (1..G), built once and used for
			* the cell-level IF cluster sums and the cluster wild bootstrap.
			tempvar __if_cluster_id
			quietly egen long `__if_cluster_id' = group(`cluster_var') if `touse'

			* --- Save user's current xtset settings
			local __lwdid_had_xtset 0
			local __lwdid_old_panel ""
			local __lwdid_old_time  ""
			local __lwdid_old_delta ""

			capture quietly xtset
			if !_rc {
				local __lwdid_old_panel "`r(panelvar)'"
				local __lwdid_old_time  "`r(timevar)'"
				local __lwdid_old_delta "`r(tdelta)'"
				if "`__lwdid_old_panel'" != "" local __lwdid_had_xtset 1
			}

			* --- Temporary internal xtset for lwdid computations
			quietly xtset `id' `__xttime'

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
				tempvar dvar_tmp
				local dvar_`g' `dvar_tmp'
				qui gen byte `dvar_`g'' = (`gvar_clean' == `g') if `touse'
				qui replace `dvar_`g'' = 0 if missing(`dvar_`g'')
				label var `dvar_`g'' "First treated in `g'"
				}

			if "`xlist'" != "" {
				qui {
						foreach v in `xlist' {
							foreach g of local cohorts {
								su `v' if `dvar_`g'' == 1 & `touse', meanonly
								local m_`v'_`g' = r(mean)
								tempvar xvar_tmp
								local xvar_`v'_`g' `xvar_tmp'
								gen double `xvar_`v'_`g'' = `v' - `m_`v'_`g'' if `touse'
							}
						}
					}
				}

		* ---  time dummies
			foreach t of numlist `tmin'/`tmax' {
				local __ti = `t' - `tmin' + 1
				tempvar fvar_tmp
				local fvar_`__ti' `fvar_tmp'
				qui gen byte `fvar_`__ti'' = (`tvar' == `t') if `touse'
				qui replace `fvar_`__ti'' = 0 if missing(`fvar_`__ti'')
			}


		* --- Stage 1: Residualized outcome y{g}d for each cohort
				* >> using the reduced-form representation (computationally efficient)
		
				quietly {
				tempvar yobs
				qui gen byte `yobs' = (`touse' & !missing(`y'))

				foreach g of local cohorts {
					tempvar __y`g'd
					local yvarname_`g' `__y`g'd'
					qui gen double `yvarname_`g'' = . if `touse'
				if `__is_demean' {

					* Demeaning family: demean / demeanq / demeanm
					if "`rolling'" == "demean" {
						tempvar Sy_pre n_pre

						* Unit-specific mean over available nonmissing pre-treatment outcomes.
						if `pre' > 0 {
							local g_pre_lo = `g' - `pre'
							bys `id': egen double `Sy_pre' = total(cond(`tvar' >= `g_pre_lo' & `tvar' < `g' & `yobs', `y', 0)) if `touse'
							bys `id': egen double `n_pre'  = total(cond(`tvar' >= `g_pre_lo' & `tvar' < `g' & `yobs', 1, 0)) if `touse'
						}
						else {
							bys `id': egen double `Sy_pre' = total(cond(`tvar' < `g' & `yobs', `y', 0)) if `touse'
							bys `id': egen double `n_pre'  = total(cond(`tvar' < `g' & `yobs', 1, 0)) if `touse'
						}

						replace `yvarname_`g'' = `y' - (`Sy_pre'/`n_pre') ///
							if `touse' & `yobs' & `n_pre' > 0 ///
							& !missing(`Sy_pre', `n_pre')

						drop `Sy_pre' `n_pre'
					}
					else {
						* Seasonal demeaning.  This is the reduced-form equivalent of
						* regress y i.quarter (or i.month) on the pre-treatment sample:
						* subtract the unit-by-season pre-treatment mean directly.
						local __svar `__qvar'
						if "`rolling'" == "demeanm" local __svar `__mvar'
						tempvar SyS nS

						if `pre' > 0 {
							local g_pre_lo = `g' - `pre'
							bys `id' `__svar': egen double `SyS' = total(cond(`tvar' >= `g_pre_lo' & `tvar' < `g' & `yobs', `y', 0)) if `touse'
							bys `id' `__svar': egen double `nS'  = total(cond(`tvar' >= `g_pre_lo' & `tvar' < `g' & `yobs', 1, 0)) if `touse'
						}
						else {
							bys `id' `__svar': egen double `SyS' = total(cond(`tvar' < `g' & `yobs', `y', 0)) if `touse'
							bys `id' `__svar': egen double `nS'  = total(cond(`tvar' < `g' & `yobs', 1, 0)) if `touse'
						}

						replace `yvarname_`g'' = `y' - (`SyS'/`nS') ///
							if `touse' & `yobs' & `nS' > 0 ///
							& !missing(`SyS', `nS')

						drop `SyS' `nS'
					}
				}
				else if `__is_detrend' {

					* Detrending family: detrend / detrendq / detrendm
					if "`rolling'" == "detrend" {
						tempvar SyP StP SttP StyP nP denomP bP aP fitP

						* Sum pre-treatment observations directly (valid for signed outcomes).
						local __precondition "`tvar' < `g' & `yobs'"
						if `pre' > 0 {
							local g_pre_lo = `g' - `pre'
							local __precondition "`tvar' >= `g_pre_lo' & `tvar' < `g' & `yobs'"
						}
						bys `id': egen double `SyP'  = total(cond(`__precondition', `y', 0)) if `touse'
						bys `id': egen double `StP'  = total(cond(`__precondition', `tvar', 0)) if `touse'
						bys `id': egen double `SttP' = total(cond(`__precondition', `tvar'^2, 0)) if `touse'
						bys `id': egen double `StyP' = total(cond(`__precondition', `tvar'*`y', 0)) if `touse'
						bys `id': egen double `nP'   = total(cond(`__precondition', 1, 0)) if `touse'

						gen double `denomP' = `nP' * `SttP' - (`StP')^2 if `touse'
						gen double `bP' = .
						replace `bP' = (`nP' * `StyP' - `StP' * `SyP') / `denomP' ///
							if `touse' & `nP' > 1 & `denomP' != 0

						gen double `aP' = .
						replace `aP' = (`SyP' - `bP' * `StP') / `nP' ///
							if `touse' & `nP' > 1 & !missing(`bP')

						gen double `fitP' = .
						replace `fitP' = `aP' + `bP' * `tvar' ///
							if `touse' & !missing(`aP', `bP')

						replace `yvarname_`g'' = `y' - `fitP' ///
							if `touse' & `yobs' & !missing(`fitP')

						drop `SyP' `StP' `SttP' `StyP' `nP' `denomP' `bP' `aP' `fitP'
					}
					else {
						* Seasonal detrending.  Use the FWL reduced form for
						* y = season FE + b*t on the pre-treatment sample.
						* This avoids one regression per unit/cohort.
						local __svar `__qvar'
						if "`rolling'" == "detrendm" local __svar `__mvar'

						tempvar preS SyS StS SttS StyS nS ssxyS ssxxS tagS ///
							numB denB bP ybarS tbarS fitP

						gen byte `preS' = 0 if `touse'
						if `pre' > 0 {
							local g_pre_lo = `g' - `pre'
							replace `preS' = (`tvar' >= `g_pre_lo' & `tvar' < `g' & `yobs') if `touse'
						}
						else {
							replace `preS' = (`tvar' < `g' & `yobs') if `touse'
						}

						bys `id' `__svar': egen double `SyS'  = total(cond(`preS', `y',          0)) if `touse'
						bys `id' `__svar': egen double `StS'  = total(cond(`preS', `tvar',       0)) if `touse'
						bys `id' `__svar': egen double `SttS' = total(cond(`preS', `tvar'^2,     0)) if `touse'
						bys `id' `__svar': egen double `StyS' = total(cond(`preS', `tvar'*`y',   0)) if `touse'
						bys `id' `__svar': egen double `nS'   = total(cond(`preS', 1,            0)) if `touse'

						gen double `ssxyS' = .
						replace `ssxyS' = `StyS' - (`StS'*`SyS'/`nS') if `touse' & `nS' > 0
						gen double `ssxxS' = .
						replace `ssxxS' = `SttS' - (`StS'^2/`nS') if `touse' & `nS' > 0

						bys `id' `__svar': gen byte `tagS' = (_n==1) if `touse'
						bys `id': egen double `numB' = total(cond(`tagS'==1 & `nS'>0, `ssxyS', 0)) if `touse'
						bys `id': egen double `denB' = total(cond(`tagS'==1 & `nS'>0, `ssxxS', 0)) if `touse'

						gen double `bP' = .
						replace `bP' = `numB'/`denB' if `touse' & `denB' > 0 & !missing(`denB')

						gen double `ybarS' = .
						replace `ybarS' = `SyS'/`nS' if `touse' & `nS' > 0
						gen double `tbarS' = .
						replace `tbarS' = `StS'/`nS' if `touse' & `nS' > 0

						gen double `fitP' = .
						replace `fitP' = (`ybarS' - `bP'*`tbarS') + `bP'*`tvar' ///
							if `touse' & `nS' > 0 & !missing(`bP', `ybarS', `tbarS')

						replace `yvarname_`g'' = `y' - `fitP' ///
							if `touse' & `yobs' & !missing(`fitP')

						drop `preS' `SyS' `StS' `SttS' `StyS' `nS' `ssxyS' `ssxxS' ///
							`tagS' `numB' `denB' `bP' `ybarS' `tbarS' `fitP'
					}
				}
					label var `yvarname_`g'' "(rolling=`rolling') Residualized outcome cohort g=`g' "

					* --- Anchor period used internally for estimation/output
					* If ydot is requested, saved y`g'd variables are set to missing at the anchor below.
					
					if `__is_demean' {
						* anchor = r = -1 only (the period just before treatment)
						replace `yvarname_`g'' = 0 if `touse' & !missing(`yvarname_`g'') ///
							& (`tvar' - `g' == -1)
					}
					else if `__is_detrend' {
						* anchor = r = -1 and r = -2 (two periods just before treatment)
						replace `yvarname_`g'' = 0 if `touse' & !missing(`yvarname_`g'') ///
							& inlist(`tvar' - `g', -2, -1)
					}
				}
			 } 
			 
	 
	* --- Stage 2: ATT(g,t) point estimates + Influence Functions
			quietly {
					preserve
					keep if `touse' & `gvar' > 0
					keep if `touse' & `gvar' > 0
					capture drop one Ng
					bys `gvar' `id': gen byte one = (_n==1)
					bys `gvar': egen Ng = total(one)
					keep `gvar' Ng
					duplicates drop
					rename `gvar' cohort
					tempfile COHORTSIZE
					save "`COHORTSIZE'", replace
					restore
				}

				qui tempfile ATTfile DROPfile
				qui tempname pf __drop_pf
				qui postfile `pf' int cohort int time int ryear double att double se double tstat double pval double lower_ci double upper_ci long Ngt long N_control using "`ATTfile'", replace
				* reason: 1 insufficient treated/control; 2 propensity failure;
				* 3 outcome regression failure; 4 treatment unidentified;
				* 5 outcome IF rank failure; 6 propensity IF rank failure.
				qui postfile `__drop_pf' int cohort int time int ryear byte reason using "`DROPfile'", replace
				local __n_dropped 0

				quietly tempvar cont
				qui gen byte `cont' = 0

				* Shared work variables for cell-level estimation.
				tempvar __es __ps __w
				qui gen byte   `__es' = 0
				qui gen byte   `__ps' = 0
				qui gen double `__w'  = .

				* Sort by time and ID for efficient cell-specific estimation.
				sort `tvar' `id'
				mata: __lw_tv = st_data(., "`tvar'")
				forval t = `tmin'/`tmax' {
					local __ti = `t' - `tmin' + 1
					mata: __lw_s = selectindex(__lw_tv :== `t')
					mata: st_local("__lw_nrow", strofreal(rows(__lw_s)))
					local __lo_`__ti' ""
					local __hi_`__ti' ""
					if `__lw_nrow' > 0 {
						mata: st_local("__lo_`__ti'", strofreal(__lw_s[1], "%20.0f"))
						mata: st_local("__hi_`__ti'", strofreal(__lw_s[rows(__lw_s)], "%20.0f"))
					}
				}
				capture mata: mata drop __lw_tv
				capture mata: mata drop __lw_s

				
				* Store cell-level IFs as cluster sums for SEs and bootstrap inference.
				* Preallocate matrices for all potential (g,t) cells.

				qui su `__if_cluster_id' if `touse', meanonly
				local __n_cl = r(max)
				local __n_coh : word count `cohorts'
				local __max_cells = `__n_coh' * (`tmax' - `tmin' + 1)
				mata: __lw_C = J(`__n_cl', `__max_cells', 0)
				* __lw_Z[c, k] = number of treated observations of cell k in cluster c
				* (composition term of the estimated weights); cell_b = ATT(g,t).
				mata: __lw_Z = J(`__n_cl', `__max_cells', 0)
				mata: cell_b = J(1, 0, .)
				mata: cell_g = J(1, 0, .)
				mata: cell_t = J(1, 0, .)
				* (event time, compact unit ID, treated flag) for SUCCESSFUL cells only.
				* These records allow unique-unit counts without adding cell counts.
				mata: __lw_watt_members = J(0, 3, .)
				local cell_count = 0

				quietly foreach g of local cohorts {
					local yvar `yvarname_`g''
					local dvar_g `dvar_`g''

					local current_x_list ""
					foreach v in `xlist' {
						local current_x_list `current_x_list' `xvar_`v'_`g''
					}

			forval t = `tmin'/`tmax' {
				local r = `t' - `g'

				* Rows for calendar time t; skip the cell if period t has no rows.
				local __ti = `t' - `tmin' + 1
				if "`__lo_`__ti''" == "" continue
				local __in "in `__lo_`__ti''/`__hi_`__ti''"

				if `t' < `g' {
					* pre-period:
					* never: never-treated + cohort g only
					* default: never-treated + cohort g + cohorts treated after g
					if "`never'" != "" {
						quietly replace `cont' = (`gvar' == 0 | `gvar' == `g') if `touse' `__in'
					}
					else {
						quietly replace `cont' = (`gvar' == 0 | `gvar' == `g' | `gvar' > `g') if `touse' `__in'
					}
				}
				else {
					* post-period:
					* never: never-treated + cohort g only
					* default: never-treated + cohort g + cohorts treated after g (and still untreated at t)
					if "`never'" != "" {
						quietly replace `cont' = (`gvar' == 0 | `gvar' == `g') if `touse' `__in'
					}
					else {
						quietly replace `cont' = (`gvar' == 0 | `gvar' == `g' | (`gvar' > `g' & `gvar' > `t')) if `touse' `__in'
					}
				}

				quietly count if `touse' & `fvar_`__ti'' & `cont' & `dvar_g' == 1 `__in'
				local __nt = r(N)
				quietly count if `touse' & `fvar_`__ti'' & `cont' & `dvar_g' == 0 `__in'
				local __nc = r(N)
				if (`__nt' == 0 & `__nc' == 0) continue
				if (`__nt' == 0 | `__nc' == 0) {
					post `__drop_pf' (`g') (`t') (`r') (1)
					local ++__n_dropped
					continue
				}

				* Variables created inside a cell; both are dropped at the end of
				* the cell (capture drop also handles ones that were not created).
				tempvar p_hat_gt uhat


				* Propensity score estimation and ATT weights
				if inlist("`method'", "ipw", "ipwra") {

					cap qui logit `dvar_g' `xlist' if `touse' & `fvar_`__ti'' & `cont' `__in', nolog
					if _rc != 0 {
						post `__drop_pf' (`g') (`t') (`r') (2)
						local ++__n_dropped
						continue
					}

					* Record exactly which logit coefficients Stata omitted.
					* The IF score/Hessian must use that same effective design.
					tempname __ps_omit
					_ms_omit_info e(b)
					matrix `__ps_omit' = r(omit)
					if e(converged) != 1 {
						post `__drop_pf' (`g') (`t') (`r') (2)
						local ++__n_dropped
						continue
					}
					qui replace `__ps' = e(sample) `__in'

					qui predict double `p_hat_gt' if `__ps' `__in', pr
					* ATT weights: 1 for treated, p/(1-p) for controls; missing outside psamp.
					qui replace `__w' = cond(`__ps', cond(`dvar_g' == 1, 1, `p_hat_gt' / (1 - `p_hat_gt')), .) `__in'
				}


				* VCE for cohort-time ATT(g,t) (RA only; see note below)
				* Default: heteroskedasticity-robust
				* cluster(): cluster-robust at user-specified level
				local __vce_opt ", vce(robust)"
				if "`cluster'" != "" {
					local __vce_opt ", vce(cluster `cluster_var')"
				}
				* IPW/IPWRA: the cell SE is the IF-based cluster sandwich SE computed
				* below, so no regression VCE is needed (point estimate and e(sample)
				* do not depend on vce()).
				if inlist("`method'", "ipw", "ipwra") local __vce_opt ""


				* Point estimation: ATT(g,t)
                local __reg_rc 0
				if "`method'" == "ra" {
					if "`xlist'" == "" {
						capture quietly regress `yvar' `dvar_g' if `touse' & `fvar_`__ti'' & `cont' `__in' `__vce_opt'
                        local __reg_rc = _rc
					}
					else {
						local __lwdid_int_list ""
						foreach __cx of local current_x_list {
							local __lwdid_int_list `__lwdid_int_list' c.`dvar_g'#c.`__cx'
						}
						capture quietly reg `yvar' `dvar_g' `current_x_list' `__lwdid_int_list' ///
							if `touse' & `fvar_`__ti'' & `cont' `__in' `__vce_opt'
                        local __reg_rc = _rc
					}
				}
				else if "`method'" == "ipw" {
					capture quietly reg `yvar' `dvar_g' [aw=`__w'] if `__ps' `__in' `__vce_opt'
                        local __reg_rc = _rc
				}
				else if "`method'" == "ipwra" {
					if "`xlist'" == "" {
						capture quietly reg `yvar' `dvar_g' [aw=`__w'] if `__ps' `__in' `__vce_opt'
                        local __reg_rc = _rc
					}
					else {
						local __lwdid_int_list ""
						foreach __cx of local current_x_list {
							local __lwdid_int_list `__lwdid_int_list' c.`dvar_g'#c.`__cx'
						}
						capture quietly reg `yvar' `dvar_g' `current_x_list' `__lwdid_int_list' ///
							[aw=`__w'] if `__ps' `__in' `__vce_opt'
                        local __reg_rc = _rc
					}
				}
				else if "`method'" == "psmatch" {
					capture quietly teffects psmatch (`yvar') (`dvar_g' `xlist') ///
						if `touse' & `fvar_`__ti'' & `cont' `__in', atet
                    local __reg_rc = _rc
				}

				if `__reg_rc' != 0 {
					post `__drop_pf' (`g') (`t') (`r') (3)
					local ++__n_dropped
					capture drop `p_hat_gt'
					continue
				}

				* Align outcome IF with Stata's actual estimated design, including
				* automatically omitted collinear covariates/interactions.
				tempname __or_omit
				if inlist("`method'", "ra", "ipwra") {
					_ms_omit_info e(b)
					matrix `__or_omit' = r(omit)
					if `__or_omit'[1,1] != 0 {
						post `__drop_pf' (`g') (`t') (`r') (4)
						local ++__n_dropped
						capture drop `p_hat_gt'
						continue
					}
				}
				qui replace `__es' = e(sample) `__in'

				* Store point estimate and cell-specific treated size
				qui tempname b_att se_att tstat_att pval_att crit_att lo_att hi_att
				if inlist("`method'", "ra", "ipw", "ipwra") {
					scalar `b_att'  = _b[`dvar_g']
					scalar `se_att' = _se[`dvar_g']
				}
				else {
					scalar `b_att'  = _b[ATET:r1vs0.`dvar_g']
					scalar `se_att' = _se[ATET:r1vs0.`dvar_g']
				}

				* CI calculations are deferred until the cell IF has been computed.

				qui count if `__es' & `dvar_g' == 1 `__in'
				local Ngt = r(N)
				qui count if `__es' & `dvar_g' == 0 `__in'
				local N_control = r(N)
				* Post only after the IF for this cell has passed all diagnostics.

				* ------------------------------------------------------------
				* Influence function, computed in Mata on the estimation-sample
				* rows of this cell (idx = row numbers within the time-t block).
				* ------------------------------------------------------------
				if inlist("`method'", "ra", "ipw", "ipwra") {

					mata: __lw_rng = (`__lo_`__ti''::`__hi_`__ti'')
					mata: idx = __lw_rng[selectindex(st_data(__lw_rng, "`__es'") :== 1)]

					if "`method'" == "ra" {

					qui predict double `uhat' if `__es' `__in', resid

					* Z = (1, D, Xc, D*Xc), built in Mata
					mata: D = st_data(idx, "`dvar_g'")
					if "`xlist'" != "" {
						mata: Xc = st_data(idx, tokens("`current_x_list'"))
						mata: Z  = J(rows(D),1,1), D, Xc, D :* Xc
					}
					else {
						mata: Z  = J(rows(D),1,1), D
					}
					mata: u   = st_data(idx, "`uhat'")
					* e(b) order: treatment, X, D#X, constant; Z order: constant, treatment, X, D#X
					mata: om = st_matrix("`__or_omit'")
					mata: keep = selectindex((om[1,cols(om)], om[|1,1 \ 1,cols(om)-1|]) :== 0)
					* Treatment is original Z column 2; omission may change its position.
					mata: att_pos = selectindex(keep :== 2)
					mata: Z = Z[,keep]
					mata: st_local("__lwdid_if_rank_ok", strofreal(length(att_pos)==1 & rank(Z)==cols(Z)))
					if !`__lwdid_if_rank_ok' {
						post `__drop_pf' (`g') (`t') (`r') (5)
						local ++__n_dropped
						capture drop `p_hat_gt'
						capture drop `uhat'
						continue
					}
					mata: ZZinv = invsym(quadcross(Z,Z))
					* treatment column of (Z*ZZinv):*u
					mata: __lw_ifs = (Z * ZZinv[., att_pos]) :* u

					* Add sampling variation of the treated covariate mean (as in
					* teffects atet): D_i*(Xc_i - xbar1)'*delta/n1, where delta =
					* coefficients on D#Xc, i.e., elements K+2 to 2K+1 of e(b).
					if "`xlist'" != "" {
						mata: __lw_b  = st_matrix("e(b)")
						mata: __lw_K  = cols(Xc)
						mata: __lw_dl = __lw_b[|1, __lw_K+2 \ 1, 2*__lw_K+1|]
						mata: __lw_dl = __lw_dl'
						mata: __lw_n1 = quadsum(D)
						mata: __lw_xm = quadcross(D, Xc) / __lw_n1
						mata: __lw_ifs = __lw_ifs :+ D :* ((Xc :- __lw_xm) * __lw_dl) / __lw_n1
					}
				}
				else if "`method'" == "ipw" {

					* IF of the normalized IPW ATT, sum(D*Y)/sum(D) - sum(q*Y)/sum(q),
					* q = (1-D)*p/(1-p), including the logit-score correction.
					* Contributions are already scaled: do not divide by N.
					mata: y = st_data(idx, "`yvar'")
					mata: d = st_data(idx, "`dvar_g'")
					mata: p = st_data(idx, "`p_hat_gt'")
					if "`xlist'" != "" {
						mata: X = st_data(idx, tokens("`xlist'"))
					}
					else {
						mata: X = J(rows(y),0,.)
					}
					mata: X = J(rows(y),1,1), X
					mata: omps = st_matrix("`__ps_omit'")
					if "`xlist'" == "" {
						mata: keep_ps = 1
					}
					else {
						mata: keep_ps = selectindex((omps[1,cols(omps)], omps[|1,1 \ 1,cols(omps)-1|]) :== 0)
					}
					mata: X = X[,keep_ps]
					mata: st_local("__lwdid_ps_rank_ok", strofreal(rank(X)==cols(X)))
					if !`__lwdid_ps_rank_ok' {
						post `__drop_pf' (`g') (`t') (`r') (6)
						local ++__n_dropped
						capture drop `p_hat_gt'
						continue
					}
					mata: q = (1 :- d) :* p :/ (1 :- p)
					mata: nd = sum(d)
					mata: nq = sum(q)
					mata: st_local("__lwdid_ipw_denom_ok", strofreal(nd>0 & nq>0 & nd<. & nq<.))
					if !`__lwdid_ipw_denom_ok' {
						post `__drop_pf' (`g') (`t') (`r') (2)
						local ++__n_dropped
						capture drop `p_hat_gt'
						continue
					}
					mata: mu1 = sum(d :* y)/nd
					mata: mu0 = sum(q :* y)/nq
					* Derivative d(mu0)/d(gamma) = sum(q*(Y-mu0)*X)/sum(q).
					mata: G = quadcross(X, q :* (y :- mu0))' / nq
					* Logit MLE IF contribution: score_i * [sum p(1-p)XX']^-1.
					mata: A = quadcross(X, X :* (p :* (1 :- p)))
					mata: Ainv = invsym(A)
					mata: IFg = (X :* (d :- p)) * Ainv
					mata: plug = d :* (y :- mu1) / nd :- q :* (y :- mu0) / nq
					mata: __lw_ifs = plug :- IFg * G'
				}
				else if "`method'" == "ipwra" {

					qui predict double `uhat' if `__es' `__in', resid

					if "`xlist'" == "" {

						mata: Z    = st_data(idx, "`dvar_g'")
						mata: Z    = J(rows(Z),1,1), Z
						* Align with outcome e(b), including an omitted intercept if any.
						mata: om0 = st_matrix("`__or_omit'")
						mata: keep0 = selectindex((om0[1,cols(om0)], om0[|1,1 \ 1,cols(om0)-1|]) :== 0)
						mata: att_pos = selectindex(keep0 :== 2)
						mata: Z = Z[,keep0]
						mata: u    = st_data(idx, "`uhat'")
						mata: w    = st_data(idx, "`__w'")
						mata: Qw   = quadcross(Z, Z :* w)
						mata: st_local("__lwdid_or_rank_ok", strofreal(length(att_pos)==1 & rank(Z)==cols(Z)))
						if !`__lwdid_or_rank_ok' {
							post `__drop_pf' (`g') (`t') (`r') (5)
							local ++__n_dropped
							capture drop `p_hat_gt'
							capture drop `uhat'
							continue
						}
						mata: Qwinv = invsym(Qw)
						mata: M    = Z :* (w :* u)
						* treatment column of M*Qwinv
						mata: __lw_ifs = M * Qwinv[., att_pos]
					}
					else {

						* Z = (1, D, X, D*X) built in Mata (no per-cell interaction variables)
						mata: D     = st_data(idx, "`dvar_g'")
						mata: Xc    = st_data(idx, tokens("`current_x_list'"))
						mata: Z     = J(rows(D),1,1), D, Xc, D :* Xc
						mata: X     = st_data(idx, tokens("`xlist'"))
						mata: X     = J(rows(X),1,1), X
						mata: omor = st_matrix("`__or_omit'")
						mata: keep_or = selectindex((omor[1,cols(omor)], omor[|1,1 \ 1,cols(omor)-1|]) :== 0)
						* Track treatment location after omitting collinear outcome columns.
						mata: att_pos = selectindex(keep_or :== 2)
						mata: Z = Z[,keep_or]
						mata: omps = st_matrix("`__ps_omit'")
						mata: keep_ps = selectindex((omps[1,cols(omps)], omps[|1,1 \ 1,cols(omps)-1|]) :== 0)
						mata: X = X[,keep_ps]
						mata: u     = st_data(idx, "`uhat'")
						mata: w     = st_data(idx, "`__w'")
						mata: p     = st_data(idx, "`p_hat_gt'")
						mata: d     = D
						mata: Qw    = quadcross(Z, Z :* w)
						mata: st_local("__lwdid_or_rank_ok", strofreal(length(att_pos)==1 & rank(Z)==cols(Z)))
						if !`__lwdid_or_rank_ok' {
							post `__drop_pf' (`g') (`t') (`r') (5)
							local ++__n_dropped
							capture drop `p_hat_gt'
							capture drop `uhat'
							continue
						}
						mata: st_local("__lwdid_ps_rank_ok", strofreal(rank(X)==cols(X)))
						if !`__lwdid_ps_rank_ok' {
							post `__drop_pf' (`g') (`t') (`r') (6)
							local ++__n_dropped
							capture drop `p_hat_gt'
							capture drop `uhat'
							continue
						}
						mata: Qwinv = invsym(Qw)
						mata: M     = Z :* (w :* u)
						mata: A     = quadcross(X, X :* (p :* (1 :- p)))
						mata: Ainv  = invsym(A)
						mata: S     = X :* (d :- p)
						mata: IFg   = S * Ainv
						mata: H     = quadcross(Z :* ((1 :- d) :* w :* u), X)
						* treatment column of (M + IFg*H')*Qwinv
						mata: __lw_ifs = (M + IFg * H') * Qwinv[., att_pos]

						* Add sampling variation of the treated covariate mean (as in RA).
						* e(b) order: D, Xc (K), D#Xc (K), _cons.
						mata: __lw_b  = st_matrix("e(b)")
						mata: __lw_K  = cols(Xc)
						mata: __lw_dl = __lw_b[|1, __lw_K+2 \ 1, 2*__lw_K+1|]
						mata: __lw_dl = __lw_dl'
						mata: __lw_n1 = quadsum(D)
						mata: __lw_xm = quadcross(D, Xc) / __lw_n1
						mata: __lw_ifs = __lw_ifs :+ D :* ((Xc :- __lw_xm) * __lw_dl) / __lw_n1
					}
				}

				* Aggregate cell IFs by cluster for SEs and WATT bootstrap.
				mata: __lw_ifs  = editmissing(__lw_ifs, 0)
				mata: __lw_cl   = st_data(idx, "`__if_cluster_id'")
				mata: __lw_o    = order(__lw_cl, 1)
				mata: __lw_cls  = __lw_cl[__lw_o]
				mata: __lw_info = panelsetup(__lw_cls, 1)
				mata: __lw_S    = panelsum(__lw_ifs[__lw_o], __lw_info)
				mata: __lw_G    = rows(__lw_S)

				* IPW/IPWRA: report cell-specific IF-based sandwich SE, using
				* the same estimated-propensity-score-adjusted cell IF as WATT.
				* Cluster sums of finite-sample IF contributions give the
				* IF-based cluster sandwich variance, without further division by N.
				if inlist("`method'", "ipw", "ipwra") {
					mata: st_numscalar("`se_att'", (__lw_G>1 ? sqrt(quadcross(__lw_S,__lw_S)) : .))
					mata: st_local("__df_att", strofreal(__lw_G - 1))
				}
				else {
					local __df_att = e(df_r)
				}

				scalar `tstat_att' = cond(`se_att' > 0 & `se_att' < ., `b_att' / `se_att', .)
				if (`tstat_att' < .) {
					if (`__df_att' < . & `__df_att' > 0) scalar `pval_att' = 2*ttail(`__df_att', abs(`tstat_att'))
					else scalar `pval_att' = 2*normal(-abs(`tstat_att'))
				}
				else scalar `pval_att' = .

				local __alpha_att = (100 - `level') / 100
				if (`__df_att' < . & `__df_att' > 0) scalar `crit_att' = invttail(`__df_att', `__alpha_att'/2)
				else scalar `crit_att' = invnormal(1 - `__alpha_att'/2)
				scalar `lo_att' = cond(`se_att' < ., `b_att' - `crit_att' * `se_att', .)
				scalar `hi_att' = cond(`se_att' < ., `b_att' + `crit_att' * `se_att', .)

				* Store this cell's cluster sums
				local cell_count = `cell_count' + 1
				mata: __lw_C[__lw_cls[__lw_info[., 1]], `cell_count'] = __lw_S
				mata: __lw_Z[__lw_cls[__lw_info[., 1]], `cell_count'] = panelsum(st_data(idx, "`dvar_g'")[__lw_o], __lw_info)
				mata: cell_b = cell_b, st_numscalar("`b_att'")
				mata: cell_g = cell_g, `g'
				mata: cell_t = cell_t, `t'

						capture drop `p_hat_gt'
						capture drop `uhat'
					}
					else {
						capture drop `p_hat_gt'
					}

					* Record original unit IDs actually in THIS cell regression sample.
					* Only reached for cells that passed all estimation/IF checks.
					mata: __lw_mm = st_data((`__lo_`__ti''::`__hi_`__ti''), ("`__es'", "`id'", "`dvar_g'"))
					mata: __lw_mm = select(__lw_mm, __lw_mm[,1] :== 1)
					mata: __lw_watt_members = __lw_watt_members \ (J(rows(__lw_mm),1,`r'), __lw_mm[,2], __lw_mm[,3])
					post `pf' (`g') (`t') (`r') (`b_att') (`se_att') (`tstat_att') (`pval_att') (`lo_att') (`hi_att') (`Ngt') (`N_control')
				}
						}
							postclose `pf'
							postclose `__drop_pf'

						* --- Tempvars are cleaned automatically; avoid wildcard drops that may touch user variables.


				* --- attgt: optional cell-specific ATT(g,t) table, before WATT aggregation
				if "`attgt'" != "" {
					preserve
					use "`ATTfile'", clear
					qui drop if missing(att)
					qui keep if ryear >= 0
					qui sort cohort time

					* Printed large-N results use seven digits after the decimal point.
					format att se lower_ci upper_ci tstat pval %12.7f

					di as txt _n " {bf:Group-time ATT(g,t) estimates (post-period only)}"

					quietly count
					local __N_attgt = r(N)
					quietly levelsof cohort, local(__attgt_cohorts)

					foreach __g of local __attgt_cohorts {
						di as txt "{hline 102}"
						di as txt "g=`__g'"
						di as txt _col(6)  "year" ///
							_col(16) "ATT" ///
							_col(31) "Std. err." ///
							_col(46) "t" ///
							_col(60) "P>|t|" ///
							_col(74) "[`level'% conf. interval]"
						di as txt "{hline 102}"

						forvalues __i = 1/`__N_attgt' {
							if cohort[`__i'] == `__g' {
								di as txt _col(5)  %6.0f time[`__i'] ///
									as res _col(14) %9.0g att[`__i'] ///
									as res _col(29) %9.0g se[`__i'] ///
									as res _col(44) %8.2f tstat[`__i'] ///
									as res _col(58) %7.3f pval[`__i'] ///
									as res _col(73) %9.0g lower_ci[`__i'] ///
									as res _col(87) %9.0g upper_ci[`__i']
							}
						}
					}
					di as txt "{hline 102}"
					di as txt ""
					restore
				}

			
            if `cell_count' == 0 & inlist("`method'", "ra", "ipw", "ipwra") {
                di as err "lwdid: no estimable ATT(g,t) cells remain."
                exit 459
            }

	* --- Stage 2b: WATT(r) aggregation
			quietly {
			tempfile WATT_point WATT_weights AVG_point AVG_weights
			preserve
			use "`ATTfile'", clear
			if `__is_demean' {
				drop if missing(att)
			}
	else if `__is_detrend' {
				replace att = 0 if inlist(ryear,-1,-2)
				drop if missing(att) & !inlist(ryear,-1,-2)
			}
			tempvar tag
			egen byte `tag' = tag(ryear cohort)
			bys ryear: egen N_cohort = total(`tag')
			bys ryear: egen N_units  = total(Ngt)
			bys ryear: egen double N_unt = total(N_control)
			gen double N_tr = N_units
			gen double N_total = N_tr + N_unt
			gen double weight = Ngt / N_units
			duplicates drop ryear cohort, force
			keep ryear cohort time Ngt N_control N_tr N_unt N_total N_units weight
			save "`WATT_weights'", replace
			restore

			preserve
			use "`ATTfile'", clear
			if `__is_demean' {
				drop if missing(att)
			}
	else if `__is_detrend' {
				replace att = 0 if inlist(ryear,-1,-2)
				drop if missing(att) & !inlist(ryear,-1,-2)
			}
			tempvar tag
			qui egen byte `tag' = tag(ryear cohort)
			bys ryear: egen N_cohort = total(`tag')
			bys ryear: egen N_units  = total(Ngt)
			bys ryear: egen double N_unt = total(N_control)
			gen double N_tr = N_units
			gen double N_total = N_tr + N_unt
			gen double weight = Ngt / N_units
			gen double w_att  = weight * att
			duplicates drop ryear cohort, force
			collapse (sum) watt = w_att (sum) weight (first) N_cohort N_units N_tr N_unt N_total, by(ryear)
			sort ryear
			qui save "`WATT_point'", replace
			restore


			* --- Single pre- and post-treatment averages.
			*     Each ATT(g,t) receives weight Ngt divided by the total Ngt
			*     across all group-time cells included in that average.
			preserve
			use "`ATTfile'", clear
			drop if missing(att)
			gen byte avg_type = .
			if `__is_demean' {
				replace avg_type = 1 if ryear < -1
			}
			else if `__is_detrend' {
				replace avg_type = 1 if ryear < -2
			}
			replace avg_type = 2 if ryear >= 0
			keep if inlist(avg_type,1,2)
			bys avg_type: egen double N_units = total(Ngt)
			bys avg_type: egen double N_unt = total(N_control)
			gen double N_tr = N_units
			gen double N_total = N_tr + N_unt
			bys avg_type: gen long N_cells = _N
			gen double weight = Ngt / N_units
			gen double w_att  = weight * att
				* Save cell weights before collapsing. Do not use nested preserve here.
				qui save "`AVG_weights'", replace
			collapse (sum) watt = w_att (sum) weight (first) N_cells N_units N_tr N_unt N_total, by(avg_type)
			gen str8 effect = cond(avg_type==1, "Pre_avg", "Post_avg")
			sort avg_type
			qui save "`AVG_point'", replace
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

                
                preserve
                use "`AVG_point'", clear
                mkmat avg_type watt, matrix(AVG_pmat)
                restore
                mata: AVG_pmat = st_matrix("AVG_pmat")

                preserve
                use "`AVG_weights'", clear
                mkmat avg_type cohort time weight, matrix(Wavgmat)
                restore
                mata: Wavgmat = st_matrix("Wavgmat")

                * Number of estimation-sample observations per cluster
                * (used to center the IF at the observation level).
                mata: __lw_clv = sort(st_data(., "`__if_cluster_id'", "`touse'"), 1)
                mata: __lw_nc  = panelsetup(__lw_clv, 1)
                mata: __lw_nc  = __lw_nc[., 2] :- __lw_nc[., 1] :+ 1
                local n_clusters = `__n_cl'
            }

            mata {
                n_vr      = rows(WATT_pmat)
                Nobs_m    = `Nobs'
                n_cells_m = cols(cell_g)

                // cluster-level cell IFs (n_clusters x n_cells)
                IFc    = __lw_C[|1,1 \ rows(__lw_C), n_cells_m|]
                __lw_C = J(0, 0, .)
                Zc     = __lw_Z[|1,1 \ rows(__lw_Z), n_cells_m|]
                __lw_Z = J(0, 0, .)
                cellif_m = `__cellif'

                // --- cell -> event-time weights; aggregated IF = IFc * Wr
                Wr = J(n_cells_m, n_vr, 0)
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

                        Wr[k, rv] = w_k
                    }
                }
                IF_r = IFc * Wr

                * --- Composition adjustment for estimated weights (omitted with cellif).
                if (!cellif_m) {
                    for (rv=1; rv<=n_vr; rv++) {
                        for (k=1; k<=n_cells_m; k++) {
                            if (Wr[k, rv] == 0) continue
                            nk = colsum(Zc[., k])
                            if (nk <= 0 | cell_b[1, k] >= .) continue
                            IF_r[., rv] = IF_r[., rv] + Zc[., k] * (Wr[k, rv] * (cell_b[1, k] - WATT_pmat[rv, 2]) / nk)
                        }
                    }
                }

                * --- Pre_avg and Post_avg weights
                n_avg = rows(AVG_pmat)
                Wa = J(n_cells_m, n_avg, 0)

                for (av=1; av<=n_avg; av++) {
                    avg_val = AVG_pmat[av, 1]

                    for (k=1; k<=n_cells_m; k++) {
                        g_k = cell_g[1, k]
                        t_k = cell_t[1, k]

                        w_k = 0
                        for (m=1; m<=rows(Wavgmat); m++) {
                            if (Wavgmat[m,1] == avg_val & Wavgmat[m,2] == g_k & Wavgmat[m,3] == t_k) {
                                w_k = Wavgmat[m,4]
                                break
                            }
                        }

                        Wa[k, av] = w_k
                    }
                }
                IF_avg = IFc * Wa
                if (!cellif_m) {
                    for (av=1; av<=n_avg; av++) {
                        for (k=1; k<=n_cells_m; k++) {
                            if (Wa[k, av] == 0) continue
                            nk = colsum(Zc[., k])
                            if (nk <= 0 | cell_b[1, k] >= .) continue
                            IF_avg[., av] = IF_avg[., av] + Zc[., k] * (Wa[k, av] * (cell_b[1, k] - AVG_pmat[av, 2]) / nk)
                        }
                    }
                }

                * --- Center cluster IF sums, equivalent to observation-level mean centering.
                IF_r = IF_r :- __lw_nc * (colsum(IF_r) / Nobs_m)
                if (n_avg > 0) {
                    IF_avg = IF_avg :- __lw_nc * (colsum(IF_avg) / Nobs_m)
                }

                reps_m = `reps'
                n_cl   = `n_clusters'

                * --- Generate bootstrap fluctuations using cluster-level IF sums
                * and one Rademacher multiplier per cluster.
                BS_star = J(reps_m, n_vr, .)
                BS_avg  = J(reps_m, n_avg, .)

                for (rep=1; rep<=reps_m; rep++) {
                    xi_cl = ((runiform(n_cl,1):>0.5) :- 0.5) * 2

                    BS_star[rep,.] = xi_cl' * IF_r
                    if (n_avg > 0) {
                        BS_avg[rep,.] = xi_cl' * IF_avg
                    }
                }

            * --- What gets plotted / reported depends on rolling()
                WATT_plot = WATT_pmat[,2]
                BS_plot   = BS_star

                if (`__is_demean') {
				* --- r = -1 is anchored at 0
                    for (rv=1; rv<=n_vr; rv++) {
                        if (WATT_pmat[rv,1] == -1) {
                            WATT_plot[rv]  = 0
                            BS_plot[., rv] = J(rows(BS_plot), 1, 0)
                        }
                    }
                }
                else if (`__is_detrend') {
			* --- r = -1 and r = -2 are anchored at 0
                    for (rv=1; rv<=n_vr; rv++) {
                        if (WATT_pmat[rv,1] == -1 | WATT_pmat[rv,1] == -2) {
                            WATT_plot[rv]  = 0
                            BS_plot[., rv] = J(rows(BS_plot), 1, 0)
                        }
                    }
                }

            * --- Append plotting/reporting estimand as third column
                WATT_pmat = WATT_pmat, WATT_plot

			* --- Bootstrap matrix used for simultaneous event-study bands.
                BS_all = BS_plot
            }

            local alpha  = (100 - `level') / 100
            local z_point = invnormal(1 - `alpha'/2)

            mata: n_vr_sc = rows(WATT_pmat)
            mata: st_numscalar("n_vr_sc", n_vr_sc)
            local n_vr = n_vr_sc

			* point estimates, standard errors, and pointwise confidence intervals
            forvalues col = 1/`n_vr' {
                mata: st_numscalar("rv_sc", WATT_pmat[`col', 1])
                mata: st_numscalar("theta_raw_sc",  WATT_pmat[`col', 2])
                mata: st_numscalar("theta_plot_sc", WATT_pmat[`col', 3])
                mata: st_numscalar("se_plot_sc", sqrt(variance(BS_plot[., `col'])))

                local rv_`col'         = rv_sc
                local theta_raw_`col'  = theta_raw_sc
                local theta_plot_`col' = theta_plot_sc
                local se_plot_`col'    = se_plot_sc
				local lo_ci_point_`col' = `theta_plot_`col'' - `z_point' * `se_plot_`col''
				local hi_ci_point_`col' = `theta_plot_`col'' + `z_point' * `se_plot_`col''
            }

            mata: n_avg_sc = rows(AVG_pmat)
            mata: st_numscalar("n_avg_sc", n_avg_sc)
            local n_avg = n_avg_sc
            forvalues col = 1/`n_avg' {
                mata: st_numscalar("avg_type_sc", AVG_pmat[`col', 1])
                mata: st_numscalar("avg_theta_sc", AVG_pmat[`col', 2])
                mata: st_numscalar("avg_se_sc", sqrt(variance(BS_avg[., `col'])))

                local avg_type_`col'  = avg_type_sc
                local avg_theta_`col' = avg_theta_sc
                local avg_se_`col'    = avg_se_sc
                local avg_label_`col' = cond(`avg_type_`col''==1, "Pre_avg", "Post_avg")
            }

            if "`ci_type'" == "simultaneous" {
                * Simultaneous critical value over the full displayed event-study path
                * (all pre- and post-treatment event times; anchor columns are zero).
                mata: st_numscalar("c_sup_sc", lwdid_csup_all(BS_all, `alpha'))
                local c_sup = c_sup_sc

                if "`c_sup'" == "" | "`c_sup'" == "." {
                    di as error "Failed to compute simultaneous critical value."
                    exit 198
                }

				forvalues col = 1/`n_vr' {
					local lo_band_`col' = `theta_plot_`col'' - `c_sup' * `se_plot_`col''
					local hi_band_`col' = `theta_plot_`col'' + `c_sup' * `se_plot_`col''
				}
                forvalues col = 1/`n_avg' {
                    local avg_lo_`col' = `avg_theta_`col'' - `z_point' * `avg_se_`col''
                    local avg_hi_`col' = `avg_theta_`col'' + `z_point' * `avg_se_`col''
                }
            }

            * Cache first five excluded cells before the outer preserve.
            if `__n_dropped' > 0 {
                preserve
                quietly use "`DROPfile'", clear
                quietly sort cohort time
                local __drop_show = min(_N, 5)
                forvalues __k = 1/`__drop_show' {
                    local __drop_cohort_`__k' = cohort[`__k']
                    local __drop_time_`__k'   = time[`__k']
                    local __drop_ryear_`__k'  = ryear[`__k']
                    local __drop_reason_`__k' = reason[`__k']
                }
                restore
            }

                   preserve
            qui use "`WATT_point'", clear
			qui gen double se         = .
			qui gen double lower_ci   = .
			qui gen double upper_ci   = .
			qui gen double lower_band = .
			qui gen double upper_band = .
			qui gen double watt_plot       = .
			qui gen double se_plot         = .
			qui gen double lower_ci_plot   = .
			qui gen double upper_ci_plot   = .
			qui gen double lower_band_plot = .
			qui gen double upper_band_plot = .
			qui gen double base_rminus1    = .
			qui gen double watt_norm       = .
			qui gen double lower_ci_norm   = .
			qui gen double upper_ci_norm   = .
			qui gen double lower_band_norm = .
			qui gen double upper_band_norm = .

            quietly forval col = 1/`n_vr' {
					qui replace watt_plot       = `theta_plot_`col''   if ryear == `rv_`col''
					qui replace se_plot         = `se_plot_`col''      if ryear == `rv_`col''
					qui replace lower_ci_plot   = `lo_ci_point_`col''  if ryear == `rv_`col''
					qui replace upper_ci_plot   = `hi_ci_point_`col''  if ryear == `rv_`col''
					qui replace lower_band_plot = `lo_band_`col''      if ryear == `rv_`col''
					qui replace upper_band_plot = `hi_band_`col''      if ryear == `rv_`col''
			}

			if `__is_detrend' {
					qui replace watt_plot       = 0 if inlist(ryear,-1,-2)
					qui replace lower_ci_plot   = 0 if inlist(ryear,-1,-2)
					qui replace upper_ci_plot   = 0 if inlist(ryear,-1,-2)
					qui replace lower_band_plot = 0 if inlist(ryear,-1,-2)
					qui replace upper_band_plot = 0 if inlist(ryear,-1,-2)
					qui replace se_plot         = 0 if inlist(ryear,-1,-2)
			}

			qui replace se         = se_plot
			qui replace lower_ci   = lower_ci_plot
			qui replace upper_ci   = upper_ci_plot
			qui replace lower_band = lower_band_plot
			qui replace upper_band = upper_band_plot

        if `__is_demean' {
            qui su watt if ryear == -1, meanonly
            local base = r(mean)
			qui replace base_rminus1    = `base'
			qui replace watt_norm       = watt_plot
			qui replace lower_ci_norm   = lower_ci_plot
			qui replace upper_ci_norm   = upper_ci_plot
			qui replace lower_band_norm = lower_band_plot
			qui replace upper_band_norm = upper_band_plot

            * r = -1 is an anchor only: report it as exactly zero with no band
			qui replace watt_norm       = 0 if ryear == -1
			qui replace lower_ci_norm   = 0 if ryear == -1
			qui replace upper_ci_norm   = 0 if ryear == -1
			qui replace lower_band_norm = 0 if ryear == -1
			qui replace upper_band_norm = 0 if ryear == -1

            * finalized results for demean: keep anchored series without subtracting WATT(-1)
			qui replace watt       = watt_norm
			qui replace lower_ci   = lower_ci_norm
			qui replace upper_ci   = upper_ci_norm
			qui replace lower_band = lower_band_norm
			qui replace upper_band = upper_band_norm
        }
        else {
            * finalized results for non-demean cases: keep plotted series
			qui replace watt       = watt_plot
			qui replace lower_ci   = lower_ci_plot
			qui replace upper_ci   = upper_ci_plot
			qui replace lower_band = lower_band_plot
			qui replace upper_band = upper_band_plot
        }

        * finalized standard errors
        qui replace se = se_plot

        * Add Pre_avg and Post_avg rows 
        qui gen long N_cells = N_cohort
        qui gen str12 effect = string(ryear, "%9.0g")
        qui gen byte is_avg = 0
        qui gen double sort_order = 1000 + ryear

        tempfile __WATT_event
        qui save "`__WATT_event'", replace

        qui clear
        qui set obs `n_avg'
        qui gen str12 effect = ""
        qui gen double ryear = .
        qui gen double watt = .
        qui gen double weight = .
        qui gen long N_cells = .
        qui gen double N_units = .
		qui gen double N_tr = .
		qui gen double N_unt = .
		qui gen double N_total = .
		qui gen double N_cohort = .
        qui gen double se = .
        qui gen double lower_ci = .
        qui gen double upper_ci = .
	        qui gen double lower_band = .
	        qui gen double upper_band = .
        qui gen byte is_avg = 1
        qui gen double sort_order = .

        forvalues col = 1/`n_avg' {
            qui replace effect     = "`avg_label_`col''" in `col'
            qui replace watt       = `avg_theta_`col'' in `col'
            qui replace se         = `avg_se_`col'' in `col'
            qui replace lower_ci   = `avg_lo_`col'' in `col'
            qui replace upper_ci   = `avg_hi_`col'' in `col'
            qui replace sort_order = cond("`avg_label_`col''"=="Pre_avg", -2, -1) in `col'
        }
	        * Bring in the N of contributing cells/units for Pre_avg and Post_avg.
	        qui gen byte avg_type = cond(effect=="Pre_avg",1,cond(effect=="Post_avg",2,.))
	        qui merge m:1 avg_type using "`AVG_point'", nogen update replace keepusing(N_cells N_units N_tr N_unt N_total)
	        qui drop avg_type

        qui append using "`__WATT_event'"
        qui replace N_cells = N_cohort if missing(N_cells) & !missing(N_cohort)

        * t statistics and pointwise large-N normal p-values
        capture drop t_stat p_value
        qui gen double t_stat  = watt / se if se > 0 & !missing(watt, se)
        qui gen double p_value = 2*normal(-abs(t_stat)) if !missing(t_stat)

        sort sort_order

	        * Printed large-N results use seven digits after the decimal point.
	        format watt se lower_ci upper_ci lower_band upper_band t_stat p_value %12.7f

		* --- Display exposure-time estimates first.
		di as txt _n "{bf: Aggregated WATT(r) estimates}"
	        * Save pointwise confidence intervals and simultaneous confidence bands once each.
	        capture drop low_ci up_ci low_band up_band
	        qui gen double low_ci   = lower_ci
	        qui gen double up_ci    = upper_ci
	        qui gen double low_band = lower_band
	        qui gen double up_band  = upper_band

	        format low_ci up_ci low_band up_band %12.7f
	        capture drop __Ncells_disp
	        qui gen str12 __Ncells_disp = string(N_cells, "%9.0g") if is_avg==0
	        qui replace __Ncells_disp = "" if is_avg==0 & missing(N_cells)
			
			* Display -1 & -2 (with detrending) periods as missing in the table only.
			tempvar anchor watt_disp se_disp t_disp p_disp lo_disp hi_disp

			qui gen byte `anchor' = is_avg == 0 & ( ///
				(`__is_demean'  & ryear == -1) | ///
				(`__is_detrend' & inlist(ryear, -2, -1)) ///
			)

			qui gen double `watt_disp' = watt
			qui gen double `se_disp'   = se
			qui gen double `t_disp'    = t_stat
			qui gen double `p_disp'    = p_value
			if "`interval_display'" == "pointwise" {
				qui gen double `lo_disp' = low_ci
				qui gen double `hi_disp' = up_ci
				local interval_header "[`level'% conf. interval]"
			}
			else {
				qui gen double `lo_disp' = low_band
				qui gen double `hi_disp' = up_band
				local interval_header "[`level'% conf. band]"
			}

			qui replace `watt_disp' = . if `anchor'
			qui replace `se_disp'   = . if `anchor'
			qui replace `t_disp'    = . if `anchor'
			qui replace `p_disp'    = . if `anchor'
			qui replace `lo_disp'   = . if `anchor'
			qui replace `hi_disp'   = . if `anchor'

		di as txt "{hline 93}"
		di as txt _col(2)  "r" ///
			_col(13) "WATT" ///
			_col(25) "Std. err." ///
			_col(38) "t" ///
			_col(50) "P>|t|" ///
			_col(60) "`interval_header'" ///
			_col(85) "N cohorts"
		di as txt "{hline 93}"

		quietly count
		local __N_display = r(N)

		forvalues __i = 1/`__N_display' {
			if is_avg[`__i'] == 0 {
		di as txt _col(2)  %7s effect[`__i'] ///
			as res _col(11) %9.0g `watt_disp'[`__i'] ///
			as res _col(23) %9.0g `se_disp'[`__i'] ///
			as res _col(35) %8.2f `t_disp'[`__i'] ///
			as res _col(47) %7.3f  `p_disp'[`__i'] ///
			as res _col(59) %9.0g `lo_disp'[`__i'] ///
			as res _col(71) %9.0g `hi_disp'[`__i'] ///
			as txt _col(85) %7s __Ncells_disp[`__i']
			}
		}

		di as txt "{hline 93}"
		* --- Display the weighted pre- and post-treatment averages separately, below WATT(r).
		di as txt _n "{bf: Cohort-size-weighted averages of ATT(g,t)}"

		di as txt "{hline 82}"
		di as txt _col(2)  "effect" ///
			_col(13) "ATT" ///
			_col(25) "Std. err." ///
			_col(38) "t" ///
			_col(50) "P>|t|" ///
			_col(60) "[`level'% conf. interval]"
		di as txt "{hline 82}"

		forvalues __i = 1/`__N_display' {
			if is_avg[`__i'] == 1 {
				di as txt _col(2)  %8s effect[`__i'] ///
					as res _col(11) %9.0g watt[`__i'] ///
					as res _col(23) %9.0g se[`__i'] ///
					as res _col(35) %8.2f t_stat[`__i'] ///
					as res _col(47) %7.3f  p_value[`__i'] ///
					as res _col(59) %9.0g low_ci[`__i'] ///
					as res _col(71) %9.0g up_ci[`__i']
			}
		}

		di as txt "{hline 82}"

		if "`interval_display'" == "simultaneous" {
			di as txt "Current run: {bf:simultaneous `level'% confidence bands} over all estimated event times (pre + post)."
			di as txt "Simultaneous critical value: " as res %7.4f `c_sup'
		}
		else {
			di as txt "Current run: {bf:pointwise `level'% confidence intervals}."
		}
		if `__cellif' di as txt "Inference: cohort-size weights held fixed in the IF (cellif; as in Lee and Wooldridge 2026a)."
		else          di as txt "Inference: IF includes the composition term of the estimated cohort-size weights (use {bf:cellif} for the fixed-weight IF)." 
		* --- Concise exclusion warning AFTER the WATT and average ATT tables.
		if `__n_dropped' > 0 {
			di as err _n "Warning: Some ATT(g,t) estimates were excluded."
			di as txt _n "  Cohort     Time     Relative time    Reason code"
			di as txt "  {hline 74}"
            forvalues __k = 1/`__drop_show' {
                di as txt "  " %6.0f `__drop_cohort_`__k'' _col(14) %6.0f `__drop_time_`__k'' ///
                    _col(29) %6.0f `__drop_ryear_`__k'' _col(42) %2.0f `__drop_reason_`__k''
            }
			di as txt "  {hline 74}"
			di as err "  Total excluded ATT(g,t): `__n_dropped'"
			di as txt "  See e(ATTgt_dropped) for the complete list."
			di as txt _n "  Reason codes in e(ATTgt_dropped):"
			di as txt "    1 = Insufficient treated/control observations"
			di as txt "    2 = Propensity model failed/nonconverged"
			di as txt "    3 = Outcome regression failed"
			di as txt "    4 = Treatment not identified"
			di as txt "    5 = Outcome IF singular"
			di as txt "    6 = Propensity IF singular"
			di as txt _n "  Note: WATT(r) and average ATT estimates are reweighted"
			di as txt "        over successfully estimated ATT(g,t) cells."
		}
		* --- Graph
        if "`graph'" != "" {
                if "`title'" == "" local title "lwdid: `method' (`rolling')"
                local base_r = 0
                qui su ryear, meanonly
                local xmin = r(min)
                local xmax = r(max)
                local xrange = `xmax' - `xmin'
                local xstep = cond(`xrange'>40, 10, 5)

	                * y-axis range based on the interval type displayed in the graph
	                local graph_lo "low_band"
	                local graph_hi "up_band"
	                if "`interval_display'" == "pointwise" {
	                    local graph_lo "low_ci"
	                    local graph_hi "up_ci"
	                }
	                qui su `graph_hi' if !missing(`graph_hi') & !missing(ryear), meanonly
                local yhi = r(max)
	                qui su `graph_lo' if !missing(`graph_lo') & !missing(ryear), meanonly
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
                * robust y-axis ticks for both small and large effect scales
                local y_max_ticks = 6
                local y_min_step  = 0.1
                local yrange = `yhi' - `ylo'

                if (`yrange' <= 0) {
                    local ystep = `y_min_step'
                    local ymin = `ylo' - 5*`ystep'
                    local ymax = `yhi' + 5*`ystep'
                }
                else {
                    local raw  = max(`yrange'/`y_max_ticks', `y_min_step')
                    local k    = floor(log10(`raw'))
                    local base = 10^`k'
                    local frac = `raw'/`base'
                    local nice = cond(`frac'<=1, 1, cond(`frac'<=2, 2, cond(`frac'<=5, 5, 10)))
                    local ystep = max(`nice'*`base', `y_min_step')
                    local ymin = `ystep'*floor(`ylo'/`ystep')
                    local ymax = `ystep'*ceil(`yhi'/`ystep')

                    local nticks = floor((`ymax' - `ymin')/`ystep') + 1
                    local guard = 0
                    while (`nticks' > `y_max_ticks' & `guard' < 10) {
                        local ystep = `ystep'*2
                        local ymin = `ystep'*floor(`ylo'/`ystep')
                        local ymax = `ystep'*ceil(`yhi'/`ystep')
                        local nticks = floor((`ymax' - `ymin')/`ystep') + 1
                        local guard = `guard' + 1
                    }
                }

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
                if `__is_demean' {
                    local yttl "WATT(r)"
                }
                else if `__is_detrend' {
                    local yttl "WATT(r)"
                }

                twoway ///
	                    (rcap `graph_lo' `graph_hi' ryear if !missing(ryear) & ryear < 0, ///
                        lwidth(0.3) lcolor(`col_pre'%50)) ///
	                    (rcap `graph_lo' `graph_hi' ryear if !missing(ryear) & ryear >= 0, ///
                        lwidth(0.3) lcolor(`col_post'%60)) ///
                    (line watt ryear if !missing(ryear), ///
                        lcolor(gs8) lwidth(thin)) ///
                    (scatter watt ryear if !missing(ryear) & ryear < 0, ///
                        mcolor(`mcol_pre') msymbol(circle) msize(medlarge)) ///
                    (scatter watt ryear if !missing(ryear) & ryear >= 0, ///
                        mcolor(`mcol_post') msymbol(circle) msize(medlarge)) ///
                    , ///
                    yline(0, lcolor(gs10)) ///
                    xline(0, lcolor(gs10) lpattern(dash)) ///
                    xtitle("Time to Treatment (r)") ///
                    ytitle("`yttl'") ///
                    title(`"`title'"') ///
                    xlabel(`xmin'(`xstep')`xmax' 0, labsize(small) nogrid) ///
                    ylabel(`ymin'(`ystep')`ymax', format(%9.2f) nogrid) ///
                    legend(off) ///
                    graphregion(color(white)) plotregion(color(white)) ///
                    scheme(`scheme') ///
                    `gopts'
        }
		
* --- save
		if "`save'" != "" {

			if "`cb'" != "" {
				qui keep effect ryear ///
					watt se t_stat p_value ///
					low_ci up_ci low_band up_band ///
					N_cells N_units N_tr N_unt N_total

				qui order effect ryear ///
					watt se t_stat p_value ///
					low_ci up_ci low_band up_band ///
					N_cells N_units N_tr N_unt N_total

				format watt se t_stat low_ci up_ci low_band up_band %10.0g
				format p_value %5.3f

				label var low_ci   "Lower bound: pointwise `level'% confidence interval"
				label var up_ci    "Upper bound: pointwise `level'% confidence interval"
				label var low_band "Lower bound: simultaneous `level'% confidence band"
				label var up_band  "Upper bound: simultaneous `level'% confidence band"
			}
			else {
				qui keep effect ryear ///
					watt se t_stat p_value ///
					low_ci up_ci ///
					N_cells N_units N_tr N_unt N_total

				qui order effect ryear ///
					watt se t_stat p_value ///
					low_ci up_ci ///
					N_cells N_units N_tr N_unt N_total

				format watt se t_stat low_ci up_ci %10.0g
				format p_value %5.3f

				label var low_ci "Lower bound: pointwise `level'% confidence interval"
				label var up_ci  "Upper bound: pointwise `level'% confidence interval"
			}

			label var N_cells "Number of group-time cells contributing to estimate"
			label var N_units "Treated unit contributions (legacy name; sum over ATT cells)"
			label var N_tr "Treated unit contributions across ATT(g,t) cells"
			label var N_unt "Untreated (never-/not-yet-treated) unit contributions across ATT(g,t) cells"
			label var N_total "Total treated + untreated unit contributions across ATT cells"

			qui save "`save'", replace

			* Save companion ATT(g,t) dataset ONLY when attgt is requested.
			* Counts use the ACTUAL cell estimation e(sample), not potential-eligibility counts.
			if "`attgt'" != "" {
			qui use "`ATTfile'", clear
			qui keep if !missing(att)
			qui gen long N_tr = Ngt
			qui rename N_control N_unt
			qui gen long N_total = N_tr + N_unt
			qui keep cohort time ryear att se tstat pval lower_ci upper_ci N_tr N_unt N_total
			qui order cohort time ryear att se tstat pval lower_ci upper_ci N_tr N_unt N_total
			label var N_tr "Treated units actually used in cell ATT(g,t)"
			label var N_unt "Untreated comparison units (never-/not-yet-treated) used in cell ATT(g,t)"
			label var N_total "Total units actually used in cell ATT(g,t)"
			local __gt_save = "`save'"
			if regexm(lower("`__gt_save'"), "[.]dta$") {
				local __gt_save = substr("`__gt_save'", 1, length("`__gt_save'")-4)
			}
			qui save "`__gt_save'_attgt.dta", replace
			}
		}	
			restore

			* --- Optionally save residualized outcomes y{g}d (only if ydot requested)
			*     After all ATT/WATT calculations, set normalization periods to missing.
			if "`ydot'" != "" {
				quietly {
					foreach g of local cohorts {
						capture confirm new variable y`g'd
						if _rc {
							noisily di as err "ydot requested, but variable y`g'd already exists. Drop or rename it before rerunning with ydot."
							exit 110
						}
						gen double y`g'd = `yvarname_`g''
						label var y`g'd "(rolling=`rolling') Residualized outcome cohort g=`g' "
						if `__is_demean' {
							replace y`g'd = . if `touse' & !missing(y`g'd) ///
								& (`tvar' - `g' == -1)
						}
						else if `__is_detrend' {
							replace y`g'd = . if `touse' & !missing(y`g'd) ///
								& inlist(`tvar' - `g', -2, -1)
						}
					}
				}
			}

			* --- e(attgt): one row per estimated ATT(g,t) cell (incl. t < g);
			*     columns: cohort time ryear att se N_tr N_unt N_total.
			*     Excluded cells are in e(attgt_dropped); e(b) holds only the aggregates.
			tempname __lwdid_gt __lwdid_watt __lwdid_wattvar __lwdid_wattmeta
			quietly {
				preserve
				use "`ATTfile'", clear
				keep if !missing(att)
				sort cohort time
				count
				local __lwdid_ngt = r(N)
				if `__lwdid_ngt' > 0 {
					gen double N_total = Ngt + N_control
					mkmat cohort time ryear att se Ngt N_control N_total, matrix(`__lwdid_gt')
					matrix colnames `__lwdid_gt' = cohort time ryear att se N_tr N_unt N_total
				}
				restore
			}
            tempname __lwdid_dropped
            if `__n_dropped' > 0 {
                quietly {
                    preserve
                    use "`DROPfile'", clear
                    sort cohort time
                    mkmat cohort time ryear reason, matrix(`__lwdid_dropped')
                    matrix colnames `__lwdid_dropped' = cohort time ryear reason
                    restore
                }
            }
            * WATT: distinct contributing treated / comparison unit IDs per r.
            * A unit is counted once per r and role even if used in several cells.
            quietly {
                preserve
                use "`WATT_point'", clear
                sort ryear
                mkmat ryear N_cohort, matrix(`__lwdid_wattmeta')
                restore
            }
            mata {
                __lw_wm = st_matrix("`__lwdid_wattmeta'")
                __lw_wu = J(rows(__lw_wm), 2, 0)
                if (rows(__lw_watt_members)) {
                    __lw_watt_members = uniqrows(sort(__lw_watt_members, (1,2,3)))
                    for (__lw_k=1; __lw_k<=rows(__lw_wm); __lw_k++) {
                        __lw_sub = select(__lw_watt_members, __lw_watt_members[,1] :== __lw_wm[__lw_k,1])
                        if (rows(__lw_sub)) {
                            __lw_wu[__lw_k,1] = sum(__lw_sub[,3] :== 1)
                            __lw_wu[__lw_k,2] = sum(__lw_sub[,3] :== 0)
                        }
                    }
                }
                st_matrix("`__lwdid_wattmeta'", (__lw_wm, __lw_wu))
            }
            mata: st_matrix("`__lwdid_wattvar'", sqrt(diagonal(variance(BS_plot))))
            mata: st_matrix("`__lwdid_watt'", (WATT_pmat[,1], WATT_pmat[,3], st_matrix("`__lwdid_wattvar'"), st_matrix("`__lwdid_wattmeta'")[,2..4]))
            matrix colnames `__lwdid_watt' = ryear watt se N_cohorts N_unique_tr N_unique_unt

			* --- Standard e-class coefficients: average ATT and event-study WATT.
			* Wild-bootstrap covariance includes cross-covariances, not only SEs.
			* Zero-normalized event times have zero variance (as in plotted output).
			tempname __lwdid_b __lwdid_V
			mata: st_matrix("`__lwdid_b'", (AVG_pmat[,2] \ WATT_pmat[,3])')
			mata: st_matrix("`__lwdid_V'", variance((BS_avg, BS_plot)))
			local __lwdid_names
			forvalues __j = 1/`n_avg' {
				local __lwdid_names `__lwdid_names' `avg_label_`__j''
			}
			forvalues __j = 1/`n_vr' {
				local __r = `rv_`__j''
				if `__r' < 0 local __nm Rm`=abs(`__r')'
				else local __nm Rp`__r'
				local __lwdid_names `__lwdid_names' `__nm'
			}
			matrix colnames `__lwdid_b' = `__lwdid_names'
			matrix colnames `__lwdid_V' = `__lwdid_names'
			matrix rownames `__lwdid_V' = `__lwdid_names'
			ereturn post `__lwdid_b' `__lwdid_V', depname(`y')
			if `__lwdid_ngt' > 0 ereturn matrix attgt = `__lwdid_gt'

            if `__n_dropped' > 0 ereturn matrix attgt_dropped = `__lwdid_dropped'
            ereturn scalar N_dropped = `__n_dropped'
			ereturn matrix watt = `__lwdid_watt'
			ereturn scalar N_attgt = `__lwdid_ngt'
			ereturn local cmd "lwdid"
			ereturn local depvar "`y'"
			ereturn local rolling "`rolling'"
			ereturn local method "`method'"
			ereturn local ifagg = cond(`__cellif', "cell", "composition")
			ereturn local controls "`xlist'"
			ereturn scalar reps = `reps'
			forvalues __j = 1/`n_avg' {
				if "`avg_label_`__j''" == "Post_avg" {
					ereturn scalar att = `avg_theta_`__j''
					ereturn scalar se_att = `avg_se_`__j''
				}
			}

			* --- Store inference metadata
			ereturn scalar cb_crit = `c_sup'
			ereturn local ci_default "pointwise"
			ereturn local ci_display "`interval_display'"

			* --- Restore user's original xtset settings
			if `__lwdid_had_xtset' {
				if "`__lwdid_old_time'" != "" {
					if "`__lwdid_old_delta'" != "" {
						capture quietly xtset `__lwdid_old_panel' `__lwdid_old_time', delta(`__lwdid_old_delta')
						if _rc {
							quietly xtset `__lwdid_old_panel' `__lwdid_old_time'
						}
					}
					else {
						quietly xtset `__lwdid_old_panel' `__lwdid_old_time'
					}
				}
				else {
					quietly xtset `__lwdid_old_panel'
				}
			}
			else {
				capture quietly xtset, clear
			}
						
			quietly {
				* Free Mata objects created during estimation
				foreach __m in IF_r IF_avg IFc Wr Wa __lw_C __lw_Z Zc cell_b cellif_m nk __lw_clv __lw_nc ///
					cell_g cell_t Wmat Wavgmat WATT_pmat AVG_pmat ///
					BS_star BS_avg BS_plot BS_all WATT_plot ///
					n_vr n_vr_sc n_avg n_avg_sc Nobs_m n_cells_m ///
					T_star T_sort se_vec ///
					__lw_rng __lw_ifs __lw_cl __lw_o __lw_cls __lw_info __lw_S __lw_G ///
					__lw_b __lw_K __lw_dl __lw_n1 __lw_xm ///
                    __lw_watt_members __lw_mm __lw_wm __lw_wu __lw_sub __lw_k {
					capture mata: mata drop `__m'
				}
			}
		}
	end


* --- Mata helper: simultaneous critical value over all event-time columns
		capture mata: mata drop lwdid_csup_all()
		mata:
		real scalar lwdid_csup_all(real matrix BS, real scalar alpha)
		{
			real colvector se, T, Ts
			real rowvector z
			real scalar b, j, B, K, q

			B  = rows(BS)
			K  = cols(BS)
			se = sqrt(diagonal(variance(BS)))
			T  = J(B, 1, .)

			for (b=1; b<=B; b++) {
				z = J(1, K, 0)
				for (j=1; j<=K; j++) {
					if (se[j] > 0 & se[j] < .) z[j] = abs(BS[b,j] / se[j])
				}
				T[b] = max(z)
			}

			Ts = sort(T, 1)
			q  = min((B, max((1, ceil((1-alpha)*B)))))
			return(Ts[q])
		}
		end


**# RI
		capture mata: mata drop lwdid_ri_inline()

		mata:
		real scalar lwdid_ri_inline(
			real scalar reps,
			real scalar b0,
			string scalar rhs,
			string scalar sample_var
		)
		{
			real colvector fp, Y, D, Dp, res
			real matrix X, Z, bhat
			real scalar n_fp, r, i, j, tmp
			string rowvector xvars

			// The permutation set must match the main OLS e(sample),
			// including exclusions from missing additive controls.
			fp   = selectindex(st_data(., sample_var) :== 1)
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


