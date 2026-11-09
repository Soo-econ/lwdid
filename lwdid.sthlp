{smcl}
{* *! version 3.6, 8 NOVEMBER 2026 Soo Jeong Lee and Jeffrey M. Wooldridge *}
{viewerjumpto "Syntax" "lwdid##syntax"}{...}
{viewerjumpto "Options" "lwdid##options"}{...}
{viewerjumpto "Saved results" "lwdid##results"}{...}
{viewerjumpto "Description" "lwdid##description"}{...}
{viewerjumpto "Examples" "lwdid##examples"}{...}
{viewerjumpto "Citation" "lwdid##citation"}{...}
{viewerjumpto "Authors" "lwdid##author"}{...}

{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{bf:lwdid} {hline 2}}Transformation-based rolling DID estimator (Lee & Wooldridge, 2026a, 2026b){p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:lwdid} {it:varlist} [{it:if}] [{it:in}],
{cmd:ivar(}{it:varname}{cmd:)}
{cmd:tvar(}{it:varname}{cmd:)}
{cmd:gvar(}{it:varname}{cmd:)}
{cmd:rolling(}{it:type}{cmd:)}
[{it:options}]

{marker argsopts}{...}
{title:Arguments and Options}

{synoptset 28 tabbed}{...}
{synopthdr:Options}
{synoptline}

{syntab:Main variables}
{synopt:{it:varlist}}Outcome variable followed by optional covariates (x-variables).{p_end}

{syntab:Required options}
{synopt:{opt ivar(varname)}}Panel identifier (numeric or string).{p_end}

{synopt:{opt tvar(varname)}}Time variable, supplied as a single numeric time index. {p_end}

{synopt:{opt gvar(varname)}}Treatment cohort variable indicating the first-treated period for each unit. {cmd:gvar()} should be constant across all observations within each treated unit. Never-treated units should be coded as 0 or missing. {cmd:gvar()} must be measured on the same scale as {cmd:tvar()}.{p_end}

{synopt:{opt rolling(type)}}Unit-specific outcome transformation for {it:yvar}:{break}
{space 2}{bf:demean}   removes the pre-treatment mean{break} 
{space 2}{bf:detrend}  removes the pre-treatment linear trend{break} 
{space 2}{bf:demeanq}  removes the pre-treatment mean and quarterly seasonal effects{break} 
{space 2}{bf:detrendq} removes the pre-treatment linear trend and quarterly seasonal effects{break} 
{space 2}{bf:demeanm}  removes the pre-treatment mean and monthly seasonal effects{break} 
{space 2}{bf:detrendm} removes the pre-treatment linear trend and monthly seasonal effects{p_end}

{syntab:Required (depends on implementation)}
{synopt:{opt method(ra|ipw|ipwra)}}Required for the large-N implementation:{break}
{cmd:ra} for regression adjustment{break}
{cmd:ipw} for inverse probability weighting{break}
{cmd:ipwra} for inverse probability weighted regression adjustment (doubly robust){break}
{cmd:ipw} and {cmd:ipwra} require covariates.{p_end}

{synopt:{opt small}}Required for the small-N implementation (few treated or control units). Additive covariates are supported for both common timing and staggered adoption. In staggered designs, the pooled ATT uses common covariate slopes, while cohort and cohort-time regressions use separate slopes. The overall regression requires N_tr + N_unt > K+2 (K covariates).{p_end}

{syntab:Optional options}

{synopt:{opt pre(#)}}{it:Large-N only.} Specifies the number of pre-treatment periods used to construct the outcome transformation.{break}
By default, all available pre-treatment periods are used.{break}
For example, {cmd:pre(1)} uses only the last pre-treatment period, while
{cmd:pre(3)} uses the last three pre-treatment periods.{p_end}

{synopt:{opt never}}{it:Large-N only.} Uses only never-treated units as controls; by default, {cmd:lwdid} uses both never-treated and not-yet-treated units.{p_end}

{synopt:{opt attgt}}Displays cell-specific estimates: ATT(g,t) in large-N mode, cohort-time (staggered) or period-specific (common timing) estimates in small-N mode. In small-N staggered mode, cohort-specific averages appear by default; {cmd:attgt} displays all cohort-time estimates. Full results remain in {cmd:e(attgt)}.{p_end}

{synopt:{opt cluster(varname)}}{it:Large-N only.} Specifies the clustering level for ATT(g,t) standard errors and multiplier-bootstrap inference. By default, clustering is at the unit level.{p_end}

{synopt:{opt ydot}}{it:Large-N only.}  Saves the cohort-specific residualized outcome variables used to construct the ATT(g,t) estimates. By default, these variables are not saved to avoid conflicts with existing user variables.{p_end}

{synopt:{opt save(filename)}}Saves estimation results as a .dta file. In large-N mode, specifying {cmd:attgt} additionally saves a companion {cmd:_attgt.dta} file containing cell-level results. In small-N staggered mode, {cmd:save()} saves all estimable cohort-time results even when {cmd:attgt} is omitted; {cmd:e(attgt)} is always stored.{p_end}

{synopt:{opt graph}}Displays graphical results.{break}
Large-N: plots weighted ATT estimates by relative time with pointwise confidence intervals by default. {break}
Small-N: plots average residualized outcomes of treated and control units over time.{p_end}

{synopt:{opt cb}}{it:Large-N only.} Uses simultaneous confidence bands over all estimated relative times r (pre- and post-treatment) in both the event-study table and graph.{p_end}

{synopt:{opt gopts(string)}}Additional {cmd:twoway} graph options (only with {cmd:graph}).{p_end}

{synopt:{opt gid(id)}}Select treated unit to plot (default: treated-group average).{p_end}

{synopt:{opt vce(vartype)}} Small-N only. Specifies the variance estimator (e.g., robust or hc3). HC3 requires sufficient observations and may be unavailable with only one treated unit.{p_end}

{synopt:{opt cellif}}{it:Large-N only.} Holds the cohort-size weights fixed in the influence functions, reproducing the inference reported in Lee and Wooldridge (2026a). Point estimates are unchanged.{p_end}

{synopt:{opt reps(#)}}{it:Large-N only.} Number of bootstrap replications for large-N inference (default = 999).{p_end}

{synopt:{opt seed(#)}}{it:Large-N only.} Random-number seed for the bootstrap, for reproducible standard errors and confidence bands (default =12345).{p_end}

{synopt:{opt level(#)}}Confidence level, in percent (default = 95).{p_end}

{synopt:{opt ri}}{it:Small-N common timing only.} Built-in randomization inference (RI). For staggered adoption, see the external {cmd:ritest} example below.{p_end}

{synopt:{opt rireps(#)}}{it:Small-N common timing only.} Number of RI repetitions (default = 999).{p_end}

{synopt:{opt riseed(#)}}{it:Small-N common timing only.} Seed for RI reproducibility (default =12345).{p_end}

{synoptline}

{marker results}{...}
{title:Saved results}

{pstd}
The following saved results apply to large-N estimation; type {cmd:ereturn list} to display them. Small-N results are described below:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(att)}}post-treatment average ATT ({cmd:Post_avg}){p_end}
{synopt:{cmd:e(se_att)}}standard error of {cmd:e(att)}{p_end}
{synopt:{cmd:e(N_attgt)}}number of estimated ATT(g,t) cells{p_end}
{synopt:{cmd:e(N_dropped)}}number of excluded ATT(g,t) cells{p_end}
{synopt:{cmd:e(reps)}}number of bootstrap replications{p_end}
{synopt:{cmd:e(cb_crit)}}critical value for simultaneous confidence bands{p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:lwdid}{p_end}
{synopt:{cmd:e(depvar)}}outcome variable{p_end}
{synopt:{cmd:e(controls)}}covariates, as specified{p_end}
{synopt:{cmd:e(method)}}{cmd:ra}, {cmd:ipw}, or {cmd:ipwra}{p_end}
{synopt:{cmd:e(rolling)}}transformation type{p_end}
{synopt:{cmd:e(ci_default)}}default interval type ({cmd:pointwise}){p_end}
{synopt:{cmd:e(ci_display)}}{cmd:pointwise} or {cmd:simultaneous}{p_end}
{synopt:{cmd:e(ifagg)}}{cmd:composition} (default) or {cmd:cell} (with {cmd:cellif}){p_end}

{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}{cmd:Pre_avg}, {cmd:Post_avg}, and event-study WATT(r), named {cmd:Rm}{it:#} (r < 0) and {cmd:Rp}{it:#} (r >= 0){p_end}
{synopt:{cmd:e(V)}}bootstrap variance-covariance matrix of {cmd:e(b)}{p_end}
{synopt:{cmd:e(watt)}}event-study results: {cmd:ryear watt se N_cohorts N_unique_tr N_unique_unt}{p_end}
{synopt:{cmd:e(attgt)}}cell-level results (small-N staggered: post-treatment only): {cmd:cohort time ryear att se N_tr N_unt N_total}{p_end}
{synopt:{cmd:e(attgt_dropped)}}excluded cells, if any: {cmd:cohort time ryear reason}{p_end}

{pstd}
{cmd:e(attgt)} reports ATT(g,t), its standard error, and the units used in each
cell: {cmd:N_tr} treated units of cohort g, {cmd:N_unt} untreated comparison
units (never- and not-yet-treated; never-treated only with {cmd:never}), and
{cmd:N_total} = {cmd:N_tr} + {cmd:N_unt}.

{pstd}
{cmd:e(watt)} reports WATT(r), its standard error, the number of contributing
cohorts, and the numbers of distinct treated and comparison units for each r
({cmd:N_unique_tr}, {cmd:N_unique_unt}). Because a unit can appear in several
cells, and as both treated and comparison unit, these counts should not be
summed across r. Use {cmd:matrix list} to display either matrix.

{pstd}
Reason codes in {cmd:e(attgt_dropped)}:
1 = insufficient treated or control observations;
2 = propensity model failed or did not converge;
3 = outcome regression failed;
4 = treatment not identified;
5 = outcome influence function singular;
6 = propensity influence function singular.


{pstd}
In large-N mode, {cmd:e(watt)} stores the event-time WATT estimates and
sample counts, while {cmd:e(attgt)} stores cell-level ATT(g,t) estimates.
The {cmd:e(b)} matrix contains the reported aggregate and event-study
coefficients, accessible through {cmd:_b[]} and {cmd:_se[]}. Coefficient
names can be listed with, {cmd:. ereturn display, coeflegend}, and any coefficient can be saved for later use, for example,

{space 4}{cmd:. scalar post_att = _b[Post_avg]}

{p2col 5 24 28 2: Small-N saved results}{p_end}
{synopt:{cmd:e(b)}, {cmd:e(V)}}overall single-effect ATT ({cmd:Post_avg}) and its variance{p_end}
{synopt:{cmd:e(N_tr)}}distinct treated units in the overall regression{p_end}
{synopt:{cmd:e(N_unt)}}distinct untreated comparison units in the overall regression{p_end}
{synopt:{cmd:e(N_total)}}{cmd:N_tr} + {cmd:N_unt}{p_end}
{synopt:{cmd:e(cohort_att)}}cohort-specific average ATT results in staggered small-N designs{p_end}
{synopt:{cmd:e(attgt)}}period-specific (common timing) or cohort-time (staggered) estimates: {cmd:cohort time ryear att se N_tr N_unt N_total}{p_end}
{synopt:{cmd:e(p_ri)}}randomization-inference p-value, if {cmd:ri} is specified{p_end}

{pstd}
In small-N mode, {cmd:e(b)} and {cmd:e(V)} always refer to the overall
single-effect regression. The estimates in {cmd:e(attgt)} are supplementary:
with common timing, they are period-specific regressions; with staggered
adoption, they are separate regressions for each cohort and post-treatment
period ({cmd:ryear >= 0}), using never-treated units as controls. All use
t-based inference with the specified {cmd:vce()} (classical OLS by default).

{pstd}
For small-N staggered designs, {cmd:attgt} displays all estimable cohort-time cells.
Without {cmd:attgt}, they remain in {cmd:e(attgt)} and can be saved with {cmd:save()}.

{marker description}{...}
{title:Description}

{pstd}
{cmd:lwdid} implements the transformation-based rolling Difference-in-Differences estimators 
developed in Lee and Wooldridge (2026a, 2026b). 
The command provides a unified implementation for panel data settings with either 
a {it:large-N} or a {it:small-N} cross-sectional dimension, allowing treatment effects
to be estimated under both the staggered treatment adoption and the common timing case. 

{pstd}
The central idea is to transform outcome variables within each unit to remove
pre-treatment means, trends, or seasonal components. This within-unit transformation
converts panel data into a series of cross-sectional datasets where the residualized outcomes serve as the dependent variable.
Treatment effects are estimated from cross-sectional regressions of transformed
outcomes. The aggregate small-N single-effect regression is retained separately
from any optional cohort-time regressions.

{pstd}
By default, {cmd:lwdid} uses the large-N procedure of Lee and Wooldridge (2026a),
which is designed for panels with a large cross-sectional dimension and allows for
heterogeneous treatment effects and unit-specific heterogeneous linear trends.

{pstd}
When the cross-sectional dimension is small ({it:small-N}), conventional
large-N inference may be unreliable. In such settings, specifying
the {cmd:small} option invokes the small-sample inference procedures
developed in Lee and Wooldridge (2026b). Cell-level small-N staggered
regressions reported through {cmd:attgt} use conventional t-based inference.

{pstd}
In the large-N procedure, ATT(g,t) is estimated for each cohort and period and
aggregated, with treated cohort-size weights, into WATT(r) and pre- and
post-treatment averages. Cells that cannot be estimated are excluded; up to five
are listed in a warning, and all are stored in {cmd:e(attgt_dropped)}.

{pstd}
Inference uses a multiplier bootstrap of the influence functions, with one
multiplier per cluster (unit by default; see {cmd:cluster()}). By default, the
influence functions account for the estimated weights (Callaway and Sant'Anna
2021, Corollary 2), i.e., treated contributions are centered at the reported
aggregate rather than at each ATT(g,t); {cmd:cellif} holds the weights fixed.

{pstd}
Based on the treatment cohort variable specified in {cmd:gvar()}, {cmd:lwdid}
automatically detects whether the design involves a single treatment cohort
(common timing) or multiple cohorts (staggered adoption) and applies the
appropriate estimation procedure. 

{pstd}
Details of the implementation are given in Hur, Lee and Wooldridge (2026); 
see also the {browse "https://github.com/Soo-econ/lwdid":lwdid GitHub repository}
for additional examples and documentation.


{marker examples}{...}
{title:Examples}

{bf:Example 1:} Large-N estimation (RA, demean transformation)

{space 4}{cmd:. lwdid y, ivar(id) tvar(year) gvar(first_treat) rolling(demean) method(ra) graph}


{bf:Example 2:} Large-N estimation (IPWRA, detrend transformation)

{space 4}{cmd:. lwdid y x1 x2, ivar(id) tvar(year) gvar(first_treat) ///}
{space 8}{cmd:rolling(detrend) method(ipwra) save(myresult) graph ///}
{space 8}{cmd:gopts(ytitle("Residualized average outcome") xtitle("Year") ///}
{space 8}{cmd:title("The Effects of Walmart Opening"))}

{pstd}
This example estimates treatment effects using the IPWRA estimator with the
detrend transformation. The option save(myresult) saves the estimates to
myresult.dta. The gopts() option customizes the graph.


{bf:Example 3:} Large-N staggered adoption: cohort-time and event-time results

{space 4}{cmd:. lwdid y, ivar(id) tvar(year) gvar(first_treat) rolling(detrend) method(ra)}
{space 4}{cmd:. matrix list e(attgt)}
{space 4}{cmd:. matrix list e(watt)}

{pstd}
{cmd:e(attgt)} contains the cohort-time estimates ATT(g,t), and {cmd:e(watt)}
contains the event-time estimates WATT(r), with their standard errors and unit counts.

{bf:Example 4:} Small-N estimation (Quarterly data with detrending and seasonal adjustment)

{space 4}{cmd:. gen qdate = yq(year, quarter)}
{space 4}{cmd:. format qdate %tq}
{space 4}{cmd:. gen gq = yq(first_year, first_quarter)}
{space 4}{cmd:. format gq %tq}
{space 4}{cmd:. }
{space 4}{cmd:. lwdid y, small ivar(id) tvar(qdate) gvar(gq) rolling(detrendq) attgt graph}

{pstd}
With the {cmd:small} option, {cmd:lwdid} implements the small-N inference
procedure. Here the example uses quarterly data with the {cmd:detrendq}
transformation. For quarterly seasonal adjustment, {cmd:tvar()} should be
a single Stata quarterly date variable created by {cmd:yq()}, and {cmd:gvar()}
should be on the same scale.

{bf:Example 5:} Small-N estimation (Monthly data with detrending and seasonal adjustment)

{space 4}{cmd:. gen mdate = ym(year, month)}
{space 4}{cmd:. format mdate %tm}
{space 4}{cmd:. gen gm = ym(first_year, first_month)}
{space 4}{cmd:. format gm %tm}
{space 4}{cmd:. lwdid y, small ivar(id) tvar(mdate) gvar(gm) rolling(detrendm) graph}

{pstd}
This example uses monthly data with the {cmd:detrendm} transformation. For
monthly seasonal adjustment, {cmd:tvar()} should be a single Stata monthly
date variable created by {cmd:ym()}, and {cmd:gvar()} should be on the same
scale.


{bf:Example 6:} Small-N staggered adoption: permutation test with {cmd:ritest}

{space 4}{cmd:. ssc install ritest}
{space 4}{cmd:. ritest first_treat _b[Post_avg], cluster(id) reps(999) seed(12345) nodots: ///}
{space 8}{cmd:lwdid y x, small ivar(id) tvar(year) gvar(first_treat) rolling(demean)}

{pstd}
The built-in {cmd:ri} option is available only with common timing. With
staggered adoption, {cmd:ritest} (Heß 2017) can permute cohort assignments
across units ({cmd:cluster(id)}) and re-estimate {cmd:_b[Post_avg]}. This is
valid only if such permutations are justified by the assignment mechanism;
parallel trends alone is not sufficient. Do not specify {cmd:ri}, {cmd:graph},
{cmd:attgt}, or {cmd:save()} inside {cmd:ritest}.

{marker citation}{...}
{title:Citation}

{pstd}
Lee, S. J., and Wooldridge, J. M. (2026a), 
"A Simple Transformation Approach to Difference-in-Differences Estimation for Panel Data," forthcoming at
{browse "https://doi.org/10.1080/07350015.2026.2683047":{it:Journal of Business & Economic Statistics}}.

{pstd}
Lee, S. J., and Wooldridge, J. M. (2026b), 
"Simple Approaches to Inference with Difference-in-Differences Estimators with Small Cross-Sectional Sample Sizes,"  
Working Paper, Available at {browse "https://dx.doi.org/10.2139/ssrn.5325686":SSRN 5325686}

{pstd}
Hur, E. K., Lee, S. J., and Wooldridge, J. M. (2026), 
"Rolling Difference-in-Differences Estimation for Small and Large Panels,"  
Working Paper, Available at {browse "https://dx.doi.org/10.2139/ssrn.6502558":SSRN 6502558}

{marker references}{...}
{title:References}

{pstd}
Callaway, B., and Sant'Anna, P. H. C. (2021),
"Difference-in-Differences with Multiple Time Periods,"
{it:Journal of Econometrics} 225(2), 200-230.

{pstd}
Heß, S. (2017), "Randomization Inference with Stata: A Guide and Software,"
{it:Stata Journal} 17(3), 630-651.


{marker author}{...}
{title:Authors}

    Soo Jeong Lee
    Southern Illinois University
    {browse "mailto:soojeong.lee@siu.edu":soojeong.lee@siu.edu}

    Jeffrey M. Wooldridge
    Michigan State University
    {browse "mailto:wooldri1@msu.edu":wooldri1@msu.edu}











