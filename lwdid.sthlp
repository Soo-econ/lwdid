{smcl}
{* *! version 1.0 15oct2025 Soo Jeong Lee and Jeffrey M. Wooldridge *}
{viewerjumpto "Syntax" "lwdid##syntax"}{...}
{viewerjumpto "Description" "lwdid##description"}{...}
{viewerjumpto "Options" "lwdid##options"}{...}
{viewerjumpto "Examples" "lwdid##examples"}{...}
{viewerjumpto "Saved results" "lwdid##results"}{...}
{viewerjumpto "Citation" "lwdid##citation"}{...}
{viewerjumpto "Authors" "lwdid##author"}{...}

{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{bf:lwdid} {hline 2}}Rolling Approach (Lee & Wooldridge, 2025){p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:lwdid} {it:yvar dvar}, {cmd:ivar(}{it:id}{cmd:)} {cmd:tvar(}{it:year [quarter]}{cmd:)} {cmd:post(}{it:postvar}{cmd:)}
     {cmd:rolling(}{it: demean | detrend | demeanq | detrendq }{cmd:)}
{break}     [{cmd:vce(}{it:vartype}{cmd:)} {cmd:table(}{it:"filename"}{cmd:)} {cmd:gid(}{it:string}{cmd:)} {cmd:graph}
      {cmd:gopts(}{it:string}{cmd:)} {cmd:controls(}{it:varlist}{cmd:)}
      {break}{cmd:ri} {cmd:rireps(}{it:#}{cmd:)} {cmd:riseed(}{it:#}{cmd:)}]


	  {synoptset 30 tabbed}{...}
{marker options_table}{...}
{synopthdr:Arguments and Options}
{synoptline}

{syntab:Main variables}
{synopt:{it:yvar}}Outcome variable{p_end}
{synopt:{it:dvar}}Binary treatment indicator{p_end}

{syntab:Required options}
{synopt:{opt ivar(varname)}}Panel identifier (numeric or string){p_end}
{synopt:{opt tvar(varlist)}}Time variable(s): year  or year-quarter{p_end}
{synopt:{opt post(varname)}}Post-period indicator (1 = post, 0 = pre){p_end}
{synopt:{opt rolling(type)}}Type of transformation applied to residualize {it:yvar}:{break}
  {bf:demean} (remove pre-period mean), {bf:detrend} (remove linear pre-trend), {break}  {bf:demeanq} (demeaning+deseaonalizing), or {bf:detrendq} (detrending+deseasonalizing){p_end}

{syntab:Optional options}
{synopt:{opt vce(vartype)}}Variance estimator (e.g., {bf:robust}, {bf:cluster id}, {bf:hc3}){p_end}
{synopt:{opt controls(varlist)}}Including covariates only if both treated and control groups satisfy N1 > K+1 and N0 > K+1{p_end}
{synopt:{opt table("filename")}}Save results to Excel (.xls) and CSV files ({it:name}_byperiod.csv){p_end}
{synopt:{opt graph}}Plot treated vs. control mean of residualized outcomes over time{p_end}
{synopt:{opt gopts(string)}}Additional {cmd:twoway} graph options (effective only with {cmd:graph}){p_end}
{synopt:{opt gid(id)}}Choose which treated unit's residualized series to plot (if omitted, plots treated-group average){p_end}
{synopt:{opt ri}}Perform fully reproducible randomization inference (RI) {p_end}
{synopt:{opt rireps(#)}}Number of randomization repetitions (default = 1000){p_end}
{synopt:{opt riseed(#)}}Set seed for RI reproducibility; if omitted, a random seed is drawn automatically{p_end}

{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:lwdid} implements the rolling approach developed in Lee and Wooldridge (2025). 
The command estimates period-specific treatment effects (ATT_t) after transforming outcomes to remove 
pre-treatment means, linear trends, or seasonal components. 
It is designed for {it:small-N} panel settings with common treatment timing. 
A forthcoming version will extend the method to accommodate staggered treatment adoption and {it:large-N} panels (see, e.g., Lee and Wooldridge (2023))

{pstd}
The procedure transforms outcomes within each unit to eliminate pre-treatment averages or trends,
creating "residualized" outcomes suitable for both a single treatment effect estimation and cross-sectional ATTs estimation in each post-treatment period.
The command can also implement Randomization Inference (RI) for small-sample inference robustness.

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{cmd:ivar(}{it:id}{cmd:)} — Panel identifier variable (numeric or string).  
String IDs are internally grouped automatically.

{phang}
{cmd:tvar(}{it:year [quarter]}{cmd:)} — Time variable(s).  
For yearly data, specify only {it:year}.  
For quarterly data, specify both {it:year quarter}, where {it:quarter} must be 1–4 or labeled.

{phang}
{cmd:post(}{it:postvar}{cmd:)} — Post-period indicator (1 = post, 0 = pre).

{phang}
{cmd:rolling(}{it:type}{cmd:)} — Specifies the type of transformation applied to residualize the outcome variable:  
  - {bf: demean}   : Removes unit-specific mean over pre-treatment periods  
  - {bf: detrend}  : Removes unit-specific linear time trend  
  - {bf: demeanq}  : As demean, plus quarter fixed effects (requires {cmd:tvar(year quarter)})  
  - {bf: detrendq} : As detrend, plus quarter fixed effects (requires {cmd:tvar(year quarter)})

{dlgtab:Optional}

{phang}
{cmd:vce(}{it:vartype}{cmd:)} — Variance estimator for regressions (e.g., {bf:hc3}, {bf:robust}, or {bf:cluster id}).

{phang}
{cmd:table(}{it:"filename"}{cmd:)} — Save results to files in the current working directory:
  - Single-effect (overall post average) → Excel file via {cmd:outreg2}, saved as `"name.xls"`.
  - Period-by-period ATT estimates → CSV file saved as `"name_byperiod.csv"`.

{phang}
{cmd:graph} — Plot the residualized outcome: control-group mean vs. treated units over time.

{phang}
{cmd:gopts(}{it:string}{cmd:)} — Additional {cmd:twoway} graph options passed directly to the plot command  
(effective only when the {cmd:graph} option is specified).  
For example:  
{cmd:gopts(ytitle("Residualized outcome") xtitle("Year") legend(pos(1) ring(0)) title("Control vs Treated"))}

{phang}
{cmd:controls(}{it:varlist}{cmd:)} — List of additional covariates to include in the regressions.  
These variables are used as adjustment covariates when estimating ATT_t.  
As discussed in the paper, controls are included in the cross-sectional regressions only when  
{it:N₁ > K+1} and {it:N₀ > K+1}, where {it:K} is the number of control variables, and {it:N₁ (N₀)}  
is the number of treated (control) units.

{phang}
{cmd:gid(}{it:id}{cmd:)} — Choose what to plot on the treated side in the graph.  
  - If omitted: plots the treated-group average of the residualized outcome over time.  
  - If provided: plots the series for that specific treated unit (numeric or string ID),  
    e.g. {cmd:gid(101)} or {cmd:gid("CA")}.  
  - The specified value must correspond to a treated unit ({cmd:d==1}); otherwise, an error is thrown.

{dlgtab:Randomization Inference}

{phang}
{cmd:ri} — Performs a manual randomization inference (RI) test (instead of using {cmd:ritest})  
for the estimated ATT_t and reports the simulated p-value.  
This version ensures fully reproducible RI results when a random seed is specified.

{phang}
{cmd:rireps(}{it:#}{cmd:)} — Number of randomization repetitions (default = 1000).

{phang}
{cmd:riseed(}{it:#}{cmd:)} — Specifies the random seed used for reproducibility in the RI procedure.  
Although technically parsed as a string, users should enter a numeric value (e.g., {cmd:riseed(12345)}).  
When omitted, a random seed is drawn automatically during execution.
{marker results}{...}
{title:Saved results}

{pstd}
{cmd:lwdid} saves the following results in {cmd:e()}:

{synoptset 25 tabbed}{...}
{synopt:{cmd:e(b)}}Vector of estimated ATT_t values{p_end}
{synopt:{cmd:e(V)}}Variance-covariance matrix of ATT_t{p_end}
{synopt:{cmd:e(sample)}}Estimation sample indicator{p_end}
{synopt:{cmd:e(cmd)}}"lwdid"{p_end}
{synopt:{cmd:e(options)}}Command options used{p_end}
{synopt:{cmd:e(version)}}Command version{p_end}

{marker examples}{...}
{title:Examples}

{pstd}
Example 1: Demeaned transformation and ATT_t plot

{phang2}{cmd:. use example_panel, clear}{p_end}
{phang2}{cmd:. lwdid y d, ivar(id) tvar(year) post(post) rolling(demean) graph}{p_end}

{pstd}
Example 2: Detrended transformation and RI

{phang2}{cmd:. lwdid y d, ivar(id) tvar(year) post(post) rolling(detrend) ri }{p_end}

{pstd}
Example 3: Quarterly data with save and table options

{phang2}{cmd:. lwdid y d, ivar(id) tvar(year quarter) post(post) rolling(detrendq) table(results)}{p_end}


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
Lee, Soo Jeong, and Jeffrey M. Wooldridge (2023), "A Simple Transformation Approach to Difference-in-Differences Estimation for Panel Data," Working Paper, Available at {browse "https://dx.doi.org/10.2139/ssrn.4516518":SSSRN 4516518}.

{title:Author}

{pstd}
{bf:Soo Jeong Lee}, Southern Illinois University Carbondale,   
{browse "mailto:soojeong.lee@siu.edu":soojeong.lee@siu.edu}

{pstd}
{bf:Jeffrey M. Wooldridge}, Michigan State University,   
{browse "mailto:wooldri1@msu.edu":wooldri1@msu.edu}


