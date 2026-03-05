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

{marker argsopts}{...}
{title:Arguments and Options}

{synoptset 28 tabbed}{...}
{synopthdr:Options}
{synoptline}

{syntab:Main variables}
{synopt:{it:varlist}}Outcome variable followed by optional covariates (x-variables). Covariates are included only if both treated and control groups satisfy N1 > K+1 and N0 > K+1.{p_end}

{syntab:Required options}
{synopt:{opt ivar(varname)}}Panel identifier (numeric or string).{p_end}

{synopt:{opt tvar(varlist)}}Time variable(s): {it:year} or {it:year quarter}. Quarterly specification is supported when {cmd:small} is used.{p_end}

{synopt:{opt gvar(varname)}}Treatment cohort variable (first treated period). Never-treated units should be coded as 0.{p_end}

{synopt:{opt rolling(type)}}Outcome transformation applied to {it:yvar}:{break}
{space 2}{bf:demean}   removes pre-treatment mean{break}
{space 2}{bf:detrend}  removes pre-treatment linear trend{break}
{space 2}{bf:demeanq}  removes pre-treatment mean and quarter effects (requires {cmd:tvar(year quarter)}){break}
{space 2}{bf:detrendq} remove pre-treatment trend and quarter effects (requires {cmd:tvar(year quarter)}){p_end}

{syntab:Required (depends on implementation)}
{synopt:{opt method(ra|ipw|ipwra)}}{it:Large-N only.} Estimation method for the large-N implementation. Required unless {cmd:small} is specified.{p_end}

{synopt:{opt small}}{it:Small-N only.} Switches to the small-sample inference procedure (Lee & Wooldridge, 2025).{p_end}

{syntab:Optional options}
{synopt:{opt saving(filename)}}Save estimation results to disk (filename.dta).{p_end}

{synopt:{opt graph}}Display graphical results.{break}
Large-N: weighted ATT estimates by relative time.{break}
Small-N: treated vs. control means of residualized outcomes over time.{p_end}

{synopt:{opt gopts(string)}}Additional {cmd:twoway} graph options (only with {cmd:graph}).{p_end}

{synopt:{opt vce(vartype)}}Variance estimator for regression (e.g., {bf:robust}, {bf:cluster id}, {bf:hc3}).{p_end}

{synopt:{opt reps(#)}}{it:Large-N only.} Bootstrap repetitions for large-N inference (default = 999).{p_end}

{synopt:{opt gid(id)}}Select treated unit to plot (default: treated-group average).{p_end}

{synopt:{opt ri}}{it:Small-N only.}Perform randomization inference (RI).{p_end}

{synopt:{opt rireps(#)}}{it:Small-N only.}Number of RI repetitions (default = 1000).{p_end}

{synopt:{opt riseed(#)}}{it:Small-N only.}Seed for RI reproducibility (default: randomly drawn).{p_end}

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
The central idea is to transform outcomes within each unit to remove
pre-treatment means, trends, or seasonal components, yielding residualized
outcomes. These transformed outcomes allow treatment effects to be estimated
using simple cross-sectional regressions in each post-treatment period,
facilitating both overall and period-specific ATT estimation.

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

{marker examples}{...}
{title:Examples}

{dlgtab:{bf:Example 1:} Large-N estimation (RA, demean transformation)}

{phang}
{cmd:. lwdid y, ivar(id) tvar(year) gvar(first_treat) rolling(demean) method(ra) graph}

{bf:Example 2:} Large-N estimation (IPWRA, detrend transformation)

{cmd:. lwdid y x1 x2 x3, ivar(id) tvar(year) gvar(first_treat) rolling(detrend) method(ipwra) graph saving(mydata) gopts(ytitle("Residualized average outcome") xtitle("Year") title("The Effects of Walmart Opening"))}

{pstd}
This example estimates treatment effects using the IPWRA estimator with the
{cmd:detrend} transformation. The option {cmd:saving(mydata)} saves the
estimates to {cmd:mydata.dta}. The {cmd:gopts()} option customizes the graph.

{bf:Example 3:} Small-N estimation (Quarterly data with detrending)

{cmd:. lwdid y, small ivar(id) tvar(year quarter) gvar(first_treat) rolling(detrendq) graph}

{pstd}
With the {cmd:small} option, {cmd:lwdid} implements the small-N inference
procedure. Here the example uses quarterly data with the {cmd:detrendq}
transformation. When {cmd:detrendq} is used, the time variable must be
specified as {cmd:tvar(year quarter)}.

{pstd}
Based on the treatment cohort variable {cmd:gvar(first_treat)}, {cmd:lwdid}
automatically detects whether the design involves a single treatment cohort
(common timing) or multiple cohorts (staggered adoption), and applies the
corresponding estimation procedure described in Lee and Wooldridge (2025).

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
{title:Authors}

    Soo Jeong Lee
    Southern Illinois University Carbondale
    {browse "mailto:soojeong.lee@siu.edu":soojeong.lee@siu.edu}

    Jeffrey M. Wooldridge
    Michigan State University
    {browse "mailto:wooldri1@msu.edu":wooldri1@msu.edu}






