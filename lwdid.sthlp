{smcl}
{* *! version 2.0 04MAR2026 Soo Jeong Lee and Jeffrey M. Wooldridge *}
{viewerjumpto "Syntax" "lwdid##syntax"}{...}
{viewerjumpto "Description" "lwdid##description"}{...}
{viewerjumpto "Options" "lwdid##options"}{...}
{viewerjumpto "Examples" "lwdid##examples"}{...}
{viewerjumpto "Saved results" "lwdid##results"}{...}
{viewerjumpto "Citation" "lwdid##citation"}{...}
{viewerjumpto "Authors" "lwdid##author"}{...}

{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{bf:lwdid} {hline 2}}Transformation-based rolling DID estimator (Lee & Wooldridge, 2023, 2025){p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:lwdid} {it:varlist},
{cmd:ivar(}{it:id}{cmd:)} 
{cmd:tvar(}{it:year [quarter]}{cmd:)} 
{cmd:gvar(}{it:first_treated_year}{cmd:)}
{cmd:rolling(}{it:demean | detrend | demeanq | detrendq}{cmd:)}
{break}
[{cmd:method(}{it:ra | ipw | ipwra}{cmd:)} {cmd:small} {cmd:vce(}{it:vartype}{cmd:)}
{cmd:table(}{it:"filename"}{cmd:)} {cmd:gid(}{it:string}{cmd:)} 
{cmd:graph} {cmd:gopts(}{it:string}{cmd:)}{break}
{cmd:ri} {cmd:rireps(}{it:#}{cmd:)} {cmd:riseed(}{it:#}{cmd:)}]

	  {synoptset 30 tabbed}{...}
{marker options_table}{...}
{synopthdr:Arguments and Options}
{synoptline}

{syntab:Main variables}
{synopt:{it:varlist} Outcome variable followed by optional covariates
(Covariates are included only when both treated and control groups satisfy
N1 > K+1 and N0 > K+1.){p_end}

{syntab:Required options}
{synopt:{opt ivar(varname)}}Panel identifier (numeric or string){p_end}
{synopt:{opt tvar(varlist)}}Time variable(s): 
specify {it:year} for yearly data or {it:year quarter} for quarterly data 
(quarterly specification is supported when the {cmd:small} option is used){p_end}
{synopt:{opt gvar(varname)}}Treatment cohort variable indicating the first
period each unit becomes treated. Units that are never treated should be coded as 0{p_end}

{synopt:{opt rolling(type)}}Type of transformation applied to residualize {it:yvar}:{break}
  {bf:demean} (remove pre-period mean), {bf:detrend} (remove linear pre-trend), {break}  {bf:demeanq} (demeaning+deseaonalizing), or {bf:detrendq} (detrending+deseasonalizing){p_end}

{syntab:Optional options}
{synopt:{opt saving("filename")}} Save the estimation results to disk (Stata .dta). {p_end}
{synopt:{opt graph}}Plot treated vs. control mean of residualized outcomes over time{p_end}
{synopt:{opt gopts(string)}}Additional {cmd:twoway} graph options (effective only with {cmd:graph}){p_end}

{syntab:Large-N options}
{synopt:{opt method(ra | ipw | ipwra)}}Estimation method for the large-N
implementation: regression adjustment (RA), inverse probability weighting (IPW),
or the doubly robust IPWRA estimator. This option is required unless the
{cmd:small} option is specified{p_end}
{synopt:{opt reps(#)}}Number of bootstrap repetitions used for inference
in the large-N implementation (default = 999){p_end}
{synopt:{opt vce(vartype)}}Variance estimator for regression (e.g., {bf:robust}, {bf:cluster id}, or {bf:hc3}){p_end}

{syntab:Small-N options (with {cmd:small})}
{synopt:{opt small}}Use the small-sample inference procedure described in
Lee and Wooldridge (2025). When this option is specified, the command
switches to the small-N implementation{p_end}
{synopt:{opt gid(id)}}Choose which treated unit's residualized series to plot (if omitted, plots treated-group average){p_end}
{synopt:{opt ri}}Perform fully reproducible randomization inference (RI) {p_end}
{synopt:{opt rireps(#)}}Number of randomization repetitions (default = 1000){p_end}
{synopt:{opt riseed(#)}}Set seed for RI reproducibility; if omitted, a random seed is drawn automatically{p_end}


{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:lwdid} implements the transformation-based rolling Difference-in-Differences
estimators developed in Lee and Wooldridge (2023, 2025).
The command provides a unified implementation that accommodates both
{it:large-N} and {it:small-N} panel settings. For panels with staggered
treatment adoption and heterogeneous treatment effects, the transformation
approach of Lee and Wooldridge (2023) delivers consistent estimates of
average treatment effects. When the cross-sectional dimension is small
({it:small-N}), conventional large-N asymptotic approximations may be unreliable.
In such settings, specifying the {cmd:small} option applies the exact
small-sample inference procedures developed in Lee and Wooldridge (2025).

The treatment cohort variable specified in {cmd:gvar()} should contain the first
period in which each unit becomes treated. Units that are never treated should
have value of 0. Based on this variable, {cmd:lwdid} automatically detects whether
the design involves a single treatment cohort (common timing) or
multiple cohorts (staggered adoption), and applies the appropriate
estimation procedure.

{pstd}
The procedure transforms outcomes within each unit to remove pre-treatment
means, trends, or seasonal components, producing residualized outcomes.
These transformed outcomes can then be used to estimate both a single
treatment effect and period-specific treatment effects through simple
cross-sectional regressions in each post-treatment period.

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{cmd:ivar(}{it:id}{cmd:)} — Panel identifier variable (numeric or string).  
String IDs are internally grouped automatically.

{phang}
{cmd:tvar(}{it:year [quarter]}{cmd:)} — Time variable(s).
Specify {it:year} for yearly data, or {it:year quarter} for quarterly data.
The quarterly specification is currently supported when the {cmd:small}
option is used.

{phang}
{cmd:gvar(}{it:first_treated_year}{cmd:)} — Treatment cohort variable indicating
the first period in which each unit becomes treated. Units that are never
treated should be coded as 0.

{phang}
{cmd:rolling(}{it:type}{cmd:)} — Specifies the transformation applied to residualize the outcome variable prior to estimation:  
  - {bf: demean}   removes unit-specific pre-treatment mean
  - {bf: detrend}  removes unit-specific pre-treatment linear trend 
  - {bf: demeanq}  removes unit-specific pre-treatment mean and quarter effects (requires {cmd:tvar(year quarter)})  
  - {bf: detrendq} removes unit-specific pre-treatment trend and quarter effects (requires {cmd:tvar(year quarter)})

{dlgtab:Optional}

{phang}
{cmd:vce(}{it:vartype}{cmd:)} — Variance estimator for regressions (e.g., {bf:hc3}, {bf:robust}, or {bf:cluster id}).

{phang}
{cmd:saving(}{it:"filename"}{cmd:)} — Save the estimation results to disk as a
Stata dataset (.dta).

{phang}
{cmd:graph} — Plot the treated and control mean of the transformed outcome over time.

{phang}
{cmd:gopts(}{it:string}{cmd:)} — Additional {cmd:twoway} graph options passed directly to the plotting command. 
This option is effective only when the {cmd:graph} option is specified.  
For example:  
{cmd:gopts(ytitle("Residualized outcome") xtitle("Year") legend(pos(1) ring(0)) title("Control vs Treated"))}

{dlgtab:Large-N estimation}

{phang}
{cmd:method(}{it:ra | ipw | ipwra}{cmd:)} — Estimation method for the large-N
implementation: regression adjustment (RA), inverse probability weighting
(IPW), or the doubly robust IPWRA estimator. This option is required unless
the {cmd:small} option is specified.

{phang}
{cmd:reps(}{it:#}{cmd:)} — Number of bootstrap repetitions used for inference
in the large-N implementation (default = 999).

{phang}
{cmd:vce(}{it:vartype}{cmd:)} — Variance estimator for regressions
(e.g., {bf:robust}, {bf:cluster id}, or {bf:hc3}).

dlgtab:Small-N inference (with {cmd:small})}

{phang}
{cmd:small} — Use the small-sample inference procedure described in
Lee and Wooldridge (2025). When specified, the command switches to the
small-N implementation.

{phang}
{cmd:gid(}{it:id}{cmd:)} — Choose which treated unit's residualized series
to display in the graph. If omitted, the treated-group average is plotted.

{phang}
{cmd:ri} — Performs randomization inference (RI) for treatment effects.

{phang}
{cmd:rireps(}{it:#}{cmd:)} — Number of randomization repetitions (default = 1000).

{phang}
{cmd:riseed(}{it:#}{cmd:)} — Random seed used for RI reproducibility.
If omitted, a seed is drawn automatically.

{marker results}{...}

{pstd}
{cmd:lwdid} with {cmd:small} option saves the following results in {cmd:e()}:
{synoptset 25 tabbed}{...}
{synopt:{cmd:e(b)}}Vector of estimated ATT_t values{p_end}
{synopt:{cmd:e(V)}}Variance-covariance matrix of ATT_t{p_end}
{synopt:{cmd:e(sample)}}Estimation sample indicator{p_end}
{synopt:{cmd:e(cmd)}}"lwdid"{p_end}
{synopt:{cmd:e(options)}}Command options used{p_end}
{synopt:{cmd:e(version)}}Command version{p_end}


{marker examples}{...}
{title:Examples}

{dlgtab:Large-N implementation}

{pstd}
Example 1:Large-N estimation (RA, demean transformation)

{phang2}{cmd:. use example_panel, clear}{p_end}
{phang2}{cmd:. lwdid y, ivar(id) tvar(year) gvar(first_treat) rolling(demean) method(ra) graph}{p_end}
{pstd}

Example 2: Large-N estimation (IPWRA, detrend transformation)
{pstd}
This example estimates treatment effects using the IPWRA estimator with
the detrend transformation. The {cmd:graph} option plots the treated and
control means of the transformed outcome over time. The option
{cmd:saving(mydata)} saves the period-by-period WATT estimates to
{cmd:mydata.dta}. The {cmd:gopts()} option customizes the graph by setting
axis labels and adding a title.

{phang2}{cmd:. lwdid y x1 x2 x3, ivar(id) tvar(year) gvar(first_treat) rolling(detrend) method(ipwra) graph saving(mydata) gopts(ytitle("Residualized average outcome") xtitle("Year") ///
                      title("The Effects of Walmart Opening")) }{p_end}

{dlgtab:Small-N implementation}

{pstd}
Example 3: Small-N estimation (Quarterly data with detrending)
{pstd}
This example illustrates the small-N implementation using quarterly data
with the  {cmd:detrendq} transformation. When {cmd:detrendq} is used,
the time variable must be specified as {cmd:tvar(year quarter)} to
account for quarter-specific effects.

The treatment cohort variable specified
in {cmd:gvar(first_treat)} indicates the first period each unit becomes
treated. Based on this variable, {cmd:lwdid} automatically detects whether
the design involves a single treatment cohort (common treatment timing)
or multiple cohorts (staggered adoption), and applies the corresponding
estimation procedure described in Lee and Wooldridge (2025).

{phang2}{cmd:. lwdid y, ivar(id) tvar(year quarter) gvar(first_treat) ///
rolling(detrendq) graph}{p_end}

{pstd}
Example 4: You can modify graph appearance directly through the {cmd:gopts()} option.

{phang2}{cmd:. lwdid y d, ivar(id) tvar(year) post(post) rolling(detrend) graph gopts(ytitle("Residualized outcome") xtitle("Year") legend(pos(1) ring(0)))}{p_end}

{marker citation}{...}
{title:Citation}

{pstd}
Lee, Soo Jeong, and Jeffrey M. Wooldridge (2025), 
"Simple Approaches to Inference with Difference-in-Differences Estimators with Small Cross-Sectional Sample Sizes,"  
Working Paper, Available at {browse "https://dx.doi.org/10.2139/ssrn.5325686":SSRN 5325686}

{pstd}
Lee, Soo Jeong, and Jeffrey M. Wooldridge (2023), 
"A Simple Transformation Approach to Difference-in-Differences Estimation for Panel Data," 
Working Paper, Available at {browse "https://dx.doi.org/10.2139/ssrn.4516518":SSSRN 4516518}.

{title:Author}

{pstd}
{bf:Soo Jeong Lee}, Southern Illinois University Carbondale,   
{browse "mailto:soojeong.lee@siu.edu":soojeong.lee@siu.edu}

{pstd}
{bf:Jeffrey M. Wooldridge}, Michigan State University,   
{browse "mailto:wooldri1@msu.edu":wooldri1@msu.edu}



