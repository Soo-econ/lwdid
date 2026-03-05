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
{cmd:lwdid} {it:varlist} [{it:if}] [{it:in}],
{cmd:ivar(}{it:varname}{cmd:)}
{cmd:tvar(}{it:varlist}{cmd:)}
{cmd:gvar(}{it:varname}{cmd:)}
{cmd:rolling(}{it:type}{cmd:)}
[{it:options}]

{marker Arguments and Options}{...}
{title:Arguments and Options}

{synopthdr:Options}
{synoptline}

{syntab:Main variables}

{synopt:{it:varlist}} Outcome variable followed by optional covariates(x-variables). Covariates are included only when both treated and control groups satisfy
N1 > K+1 and N0 > K+1.{p_end}

{syntab:Required options}
{synopt:{opt ivar(varname)}}Panel identifier (numeric or string){p_end}

{synopt:{opt tvar(varlist)}}Time variable(s): {cmd:twoway({it:year}} for yearly data or  {cmd:twoway({it:year quarter}} for quarterly data-only with {cmd:small} option{p_end}

{synopt:{opt gvar(varname)}}Treatment cohort variable indicating the first treated period for each unit. Units never treated should be coded as 0{p_end}

{synopt:{opt rolling(type)}}Outcome transformation applied to {it:yvar}:{break}
{bf:demean}    remove pre-treatment mean{break}
{bf:detrend}   remove pre-treatment linear trend{break}
{bf:demeanq}   demean + quarter effects (requires {cmd:tvar(year quarter)}){break}
{bf:detrendq}  detrend + quarter effects (requires {cmd:tvar(year quarter)}){p_end}

{syntab:Large-N Required option}
{synopt:{opt method(ra|ipw|ipwra)}}Estimation method for the large-N
implementation: regression adjustment (RA), inverse probability weighting (IPW),
or the doubly robust IPWRA estimator. This option is required unless the
{cmd:small} option is specified{p_end}

{syntab:Small-N Required option}
{synopt:{opt small}}Use the small-sample inference procedure described in
Lee and Wooldridge (2025). When this option is specified, the command
switches to the small-N implementation{p_end}

{synoptline}

{syntab:Optional options}
{synopt:{opt saving(filename)}} Save the estimation results to disk (filename.dta).{p_end}

{synopt:{opt graph}}Display graphical results. For the large-N implementation,
plots weighted ATT estimates by relative time. For the small-N implementation,
plots treated vs. control means of residualized outcomes over time{p_end}

{synopt:{opt gopts(string)}}Additional {cmd:twoway} graph options (effective only with {cmd:graph}){p_end}

{syntab:Large-N options}
{synopt:{opt reps(#)}}Number of bootstrap repetitions used for inference in the large-N implementation (default = 999){p_end}

{syntab:Small-N options}
{synopt:{opt gid(id)}}Select certain treated unit to plot (default: treated-group average){p_end}

{synopt:{opt ri}}Perform randomization inference (RI) {p_end}

{synopt:{opt rireps(#)}}Number of randomization repetitions (default = 1000){p_end}

{synopt:{opt riseed(#)}}Set seed for RI reproducibility; if omitted, a random seed is drawn automatically{p_end}

{synopt:{opt vce(vartype)}}Variance estimator for regression (e.g., {bf:robust}, {bf:cluster id}, or {bf:hc3}){p_end}

{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:lwdid} implements the transformation-based rolling Difference-in-Differences estimators 
developed in Lee and Wooldridge (2023, 2025). 
The command provides a unified implementation for panel data settings with either 
a {it:large-N} or a {it:small-N} cross-sectional dimension, allowing treatment effects
to be estimated under both the staggered treatment adoption and the common timing case. 

{pstd}
By default, {cmd:lwdid} implements the transformation-based Difference-in-Differences
approach proposed in Lee and Wooldridge (2023). This approach is designed for
panels with a large cross-sectional dimension and allows for heterogeneous treatment effects 
and unit-specific heterogeneous linear trends.

{pstd}
When the cross-sectional dimension is small ({it:small-N}), conventional
large-N asymptotic approximations may be unreliable. In such cases, specifying
the {cmd:small} option applies the exact small-sample inference procedures
developed in Lee and Wooldridge (2025).

{pstd}
A convenient feature of {cmd:lwdid} is that, based on the treatment cohort
variable specified in {cmd:gvar()}, the command automatically detects
whether the design involves a single treatment cohort (common timing) or
multiple cohorts (staggered adoption) and applies the appropriate
estimation procedure.

{pstd}
The treatment cohort variable specified in {cmd:gvar()} should contain the first
period in which each unit becomes treated. Units that are never treated should
have value of 0. Based on this variable, {cmd:lwdid} automatically detects whether
the design involves a single treatment cohort (common timing) or
multiple cohorts (staggered adoption), and applies the appropriate
estimation procedure.

{pstd}
The central idea is to transform outcomes within each unit to remove
pre-treatment means, trends, or seasonal components, yielding residualized
outcomes. These transformed outcomes allow treatment effects to be estimated
using simple cross-sectional regressions in each post-treatment period,
facilitating both overall and period-specific ATT estimation.


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

{marker author}{...}
{title:Author}

{pstd}
{bf:Soo Jeong Lee}, Southern Illinois University Carbondale,   
{browse "mailto:soojeong.lee@siu.edu":soojeong.lee@siu.edu}

{pstd}
{bf:Jeffrey M. Wooldridge}, Michigan State University,   
{browse "mailto:wooldri1@msu.edu":wooldri1@msu.edu}



