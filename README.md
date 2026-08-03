# lwdid

`lwdid` is a Stata package that implements the **rolling difference-in-differences methods** developed in **Lee and Wooldridge ([2026a](https://doi.org/10.1080/07350015.2026.2683047), [2026b](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5325686))**. These methods transform the original panel-data problem into a sequence of
cross-sectional treatment-effect problems by constructing transformed outcomes
that remove unit-specific pre-treatment means or trends.

A companion manuscript describing the method and software is available at SSRN: [Hur, Lee, and Wooldridge (2026)](https://dx.doi.org/10.2139/ssrn.6502558). The manuscript provides details on the implementation of `lwdid` and illustrates its use through several empirical examples. Users interested in applying the command are encouraged to consult the manuscript for a more complete discussion.

This page provides a concise overview of the command syntax, along with selected examples from the companion manuscript that users can adapt directly for their own applications.

# Contents

* [Important Update Note](#important-update-note)
* [Citations](#citations)
* [Installation](#installation)
* [Syntax](#syntax)
* [Examples: Small-N](#small-n-case)
* [Examples: Large-N](#large-n-case)
* [Contact and Updates](#contact-and-updates)

<br>

# Citations

`lwdid` is a user-written Stata command freely available for academic and research use.

If you use `lwdid`, or develop software based on the original `lwdid` Stata module, please acknowledge the original module and cite the relevant paper(s) below.

**[1] Large-N procedure**
Lee, S. J. and Wooldridge, J. M. (2026a). 
A Simple Transformation Approach to Difference-in-Differences Estimation for Panel Data, 
*Journal of Business and Economic Statistics*,  1–27. [https://doi.org/10.1080/07350015.2026.2683047](https://doi.org/10.1080/07350015.2026.2683047).
> Working Paper, Available at [SSRN 4516518](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4516518).

```text
Lee, S. J. and Wooldridge, J. M. (2026a). A Simple Transformation Approach to Difference-in-Differences Estimation for Panel Data. Journal of Business & Economic Statistics, 1–27. https://doi.org/10.1080/07350015.2026.2683047
```

**[2] Small-N procedure** <br>

Lee, S. J. and Wooldridge, J. M. (2026b). 
Simple Approaches to Inference with Difference-in-Differences Estimators with Small Cross-Sectional Sample Sizes.
Working Paper, Available at [SSRN 5325686](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5325686).

```text
Lee, S. J. and Wooldridge, J. M. (2026b). Simple Approaches to Inference with Difference-in-Differences Estimators with Small Cross-Sectional Sample Sizes. Working Paper. Available at SSRN 5325686
```

**[3] Original Stata module**<br>

Lee, S. J. and Wooldridge, J. M. (2026c). `lwdid`: Stata module to implement rolling difference-in-differences estimator for small-(N) and large-(N) panel data, *Statistical Software Components* S459672, Boston College Department of Economics, revised 1 June 2026.

```text
Lee, S. J. and Wooldridge, J. M. (2026c). lwdid: Stata module to implement rolling difference-in-differences estimator for small-N and large-N panel data. Statistical Software Components S459672, Boston College Department of Economics, revised 1 June 2026.
```

**[4] lwdid: command guidance** <br>

Hur, E. K., Lee, S. J. and Wooldridge, J. M. (2026). Rolling Difference-in-Differences Estimation for Small and Large Panels.
Working Paper, Available at [SSRN 6502558](https://dx.doi.org/10.2139/ssrn.6502558).

```text
Hur, E. K., Lee, S. J. and Wooldridge, J. M. (2026). Rolling Difference-in-Differences Estimation for Small and Large Panels. Working Paper. Available at SSRN 6502558
```

<a id="important-update-note"></a>
## 🚨 Important Update Note

> **Version 2.4.2 (June 15, 2026) of `lwdid.ado` is now available on both GitHub and SSC.**
> 
> To install or update the package, type: 
> ```stata
> ssc install lwdid, replace
> ```
> 
> **This version 2.4.2 include:**
> 
> - Fixed bugs; a potential naming conflict by no longer saving internally generated residualized outcome variables.
> - For the large-N case, users can use the `ydot` option to save the residualized/transformed outcomes for additional analyses.
> - The small-N `graph` option now displays y-axis labels with at least two decimal places.
>   
> **Previous SSC update — June 1, 2026**
> - new Large-\(N\) estimation options;
> - `pre(#)`: allows users to choose how many pre-treatment periods are used in the outcome transformation. For example, `pre(1)` uses only the period immediately prior to the
intervention, while `pre(3)` uses the three periods immediately prior to
the intervention.
> - `attgt`: reports group-time specific \(ATT(g,t)\) estimates. By default, `lwdid` reports event-time aggregated \(WATT(r)\) estimates only;
> - `never`: uses only never-treated units as the comparison group. By default, `lwdid` uses both never-treated and not-yet-treated units.
>   
> **Previous SSC update — April 27, 2026**
> - improved handling of large datasets and large unit identifiers, including cases with values such as `1,000,000`;
> - revised result-table presentation for cleaner reporting;
> - updated event-study plotting and result-table reporting for the Large-\(N\) `rolling(demean)` procedure;
> - improved storage of additional results through the `save()` option.
> - 
<br>



# Installation

You can install the package directly from SSC by running:
```
ssc install lwdid, replace
``` 
> To ensure that you are using the most recent version, run `which lwdid` to verify the installed version,
>
>or reinstall it using `ssc install lwdid, replace`.

 To view the SSC package description, version, and metadata, run:
```
ssc describe lwdid
```

To download the `.dta` files used to replicate the results in the paper and the [examples](#small-n-case) below, run:
```
net get lwdid
```

To view the help file and code description, please run:
```
help lwdid
```

> 🚨 To install the latest version directly from GitHub, run the following command in Stata:
```
 net install lwdid, from("https://raw.githubusercontent.com/Soo-econ/lwdid/main/") replace
 ```
 

<br>

# Syntax
```
lwdid yvar [covariates], ivar(varname) tvar(varname) gvar(varname) rolling(type)
          [method(estimator) small graph gopts(string) scheme(string)
           gid(string) table("filename") save(name)
           ri rireps(#) riseed(#)]
```

To view the full syntax and option descriptions directly in Stata:
```
help lwdid
```

## Option Details 

### **Main Variables**

| Argument | Description                |
|--------------|----------------------------|
| `yvar`       | Outcome variable           |
| `covariates` | Optional control variables  |

---

### **Required Options**

| Option | Description |
|--------|--------------|
| `ivar(varname)` | Panel identifier (numeric or string) |
| `tvar(varname)` | Time variable. The format of the variable determines the data frequency:<br><br>• Annual: numeric year variable (e.g., `year`)<br>• Quarterly: Stata quarterly date (e.g., `%tq`)<br>• Monthly: Stata monthly date (e.g., `%tm`) |
| `gvar(varname)` | Treatment cohort variable (first treated period). Never‐treated units should be coded as 0.|
| `rolling(type)` | Unit-specific Transformation type used to residualize the outcome variable. <br><br> **Available types:** <br> • `demean` – removes the unit-specific pre-treatment mean <br> • `detrend` – removes the unit-specific linear pre-treatment trend <br> • `demeanq` – demeaning with quarter-of-year effects removed <br> • `detrendq` – detrending with quarter-of-year effects removed <br> • `demeanm` – demeaning with month-of-year effects removed <br> • `detrendm` – detrending with month-of-year effects removed <br><br> *Note:* `demeanq`, `detrendq`, `demeanm`, and `detrendm` are currently available only in the small-`N` implementation. These options require `tvar()` and `gvar()` to be defined on the same Stata date scale. <br> See [Example 2: Using Quarterly (or Monthly) Data](#example-2-using-quarterly-data) for details.|


### **Required Estimation Options (depends on implementation)**
| Option | Description |
|--------|--------------|
| `small` | **Small-N Only** <br> Requests the small-sample inference procedure designed for settings with few treated (or control) units (default: large-N inference).|
|`method(ra\|ipw\|ipwra)`| <br>**Large-N Only: required unless `small` is specified**<br>•`method(ra)` –Regression Adjustment estimator <br>•`method(ipw)` – Inverse Probability Weighting estimator. <br>•`method(ipwra)` – Doubly‑robust IPW‑RA estimator.<br><br>|

---

## **Optional Options**

### Estimation Options for Large-\(N\) Panels

| Option | Description |
|-------|-------------|
| `pre(#)` | **Large-\(N\) only.** Specifies the number of pre-treatment periods used for averaging or detrending. The default uses all available pre-treatment periods; `pre(1)` uses only the last pre-treatment period for demeaning, while detrending requires at least `pre(2)`. |
| `never` | **Large-\(N\) only.** Uses only never-treated units as controls; the default also includes not-yet-treated units. |
| `attgt` | **Large-\(N\) only.** Reports group-time specific \(ATT(g,t)\) estimates before aggregation to \(W\!ATT(r)\), with robust standard errors and confidence intervals. |
| `ydot` | **Large-(N) only.** Saves the internally generated residualized/transformed outcome variables for additional analyses. By default, these variables are not saved. |


### Graph Options

| Option | Description |
|-------|-------------|
| `graph` | Produces a graph of estimated treatment effects or residualized outcome paths |
| `gopts(string)` | Additional graph options passed directly to the `twoway` command |
| `scheme(string)` | Graph scheme (e.g., `scheme(s1mono)` for black-and-white figures suitable for publication)) |
| `gid(#\|string)` |  **Available option for small-N:** <br> Displays the treated path for a specific unit instead of the average treated trajectory <br> For example, if Michigan and Illinois are treated states and the identifier for Illinois is 9, specifying `gid(9)` plots the treated path for Illinois only.|

### Other Options

| Option | Description |
|------|-------------|
| `save(name)` | Saves estimation results under the specified name |
| `reps(#)` |  **Large-N Only**. Number of wild bootstrap repetitions for large-N inference (default = 999) |
| `vce(vartype)` | **Small-N Only:** Variance estimator for the regression (e.g., `var(robust)`, `var(cluster id)`, `var(hc3)`) |
| `ri` |  **Small-N Only**. Requests randomization inference (RI) <br>•`rireps(#)` –Number of RI repetitions (default = 999).<br>•`riseed(#)` – Seed for reproducible RI p-value (default: randomly drawn; results may differ across runs if not specified)|

---

<br>

# Examples

This section illustrates basic usage of `lwdid` and shows how to replicate the main results reported in the accompanying papers.

# Small-N case

First, load the dataset:
```stata
use lw_smoking, clear
```

> ⚠️ **Important:** The `small` option must be specified for small-N settings.

<br>

## \[Example 1\] Detrended transformation and graph
```
lwdid lcigsale, small ivar(state) tvar(year) gvar(first_year) rolling(detrend) graph
```

By default, `lwdid` reports both:

* the **overall (single) treatment effect**, and
* **period-by-period ATT estimates**.

<div align="center">
<img src="subfolder/ex1.png">
</div>
<br>

> To avoid lengthy output for monthly or quarterly data, period-by-period ATTs are not displayed when there are more than 15 rows. To view the full results, specify the `save()` option; the results will be saved in a `.dta` file.

<br>

Here, `gvar(first_year)` specifies the first treatment year for each unit (never-treated units should be coded as `0`).  

By specifying `rolling(detrend)`, the command performs a __unit-specific detrending tranformation__, allowing each unit to follow its own heterogeneous linear trend and thereby relaxing the parallel trends (PT) assumption.

Alternatively, you can use `rolling(demean)` to remove **unit-specific means** (when the PT assumption hold), instead of **unit-specific trends**.


__The default `graph` option produces:__
<div align="center">
<img src="subfolder/ex2.png">
</div>

<br>

## [Example 2] Customizing graphs using `gopts()` and `scheme()`

You can __easily customize the graph__ through the `gopts()` option, which allows full customization of titles, axes, legends, and other graphical elements. For example:
```
lwdid lcigsale, small ivar(state) tvar(year) gvar(first_year) rolling(detrend) ///
      graph gopts(ytitle("Residualized average outcome") xtitle("Year") ///
      legend(pos(1) ring(0)) ///
      title("The Effects of California’s Tobacco Control Program"))
```
__This produces:__

<div align="center">
<img src="subfolder/ex3.png">
</div>

To generate a **black-and-white version** suitable for journal submissions, you can use the `scheme(s1mono)` option:

```
lwdid lcigsale, small ivar(state) tvar(year) gvar(first_year) rolling(detrend) scheme(s1mono) graph gopts(ytitle("Residualized average outcome") xtitle("Year") ///
      legend(pos(1) ring(0)) ///
      title("The Effects of California's Tobacco Control Program"))
```
<div align="center">
<img src="subfolder/ex4.png">
</div>

## \[Example 3\] Using Quarterly (or Monthly) data 

For quarterly data, `rolling(demeanq)` or `rolling(detrendq)` can be used to account for seasonality.  
For monthly data, `rolling(demeanm)` or `rolling(detrendm)` can be used to remove month-of-year effects.

To apply these transformations, `tvar()` must be specified as a **single Stata date variable**:

- quarterly data: a Stata quarterly date variable generated by `yq()` and formatted as `%tq`
- monthly data: a Stata monthly date variable generated by `ym()` and formatted as `%tm`

> ⚠️ **Important:** In all cases, `tvar()` and `gvar()` must use the **same Stata date scale**.  

For example, with quarterly data:

```stata
* Generating `tvar()`
gen tq = yq(year, q)
format tq %tq

* Generating `gvar()`
gen first_treat_q = yq(first_year, first_quarter)
replace first_treat_q = 0 if missing(first_treat_q)
format first_treat_q %tq
```

Then estimation can be carried out as follows:

```stata
lwdid y, small ivar(state) tvar(tq) gvar(first_treat_q) rolling(detrendq) graph
```

In practice, never-treated control units may be coded as `0` or missing. Under Stata’s internal date representation, such values are mapped to the base period (for example, `1960q1` for quarterly dates or `1960m1` for monthly dates), and this convention works properly in the data. For example, the data may look as follows:

| id         | tq     | first_treat_q |
|------------|--------|---------------|
| California | 2008q1 | 2008q3        |
| California | 2008q2 | 2008q3        |
| California | 2008q3 | 2008q3        |
| Illinois   | 2008q1 | 1960q1        |
| Illinois   | 2008q2 | 1960q1        |
| Illinois   | 2008q3 | 1960q1        |

Here, California is first treated in `2008q3`, while Illinois is a never-treated control unit. Since `first_treat_q` is coded as `0` before formatting, Stata displays it as the base quarterly date, `1960q1`.

For monthly data, similarly:

```stata
gen tm = ym(year, month)
format tm %tm

gen first_treat_m = ym(first_year, first_month)
replace first_treat_m = 0 if missing(first_treat_m)
format first_treat_m %tm
```

Then estimation can be carried out as follows:

```stata
lwdid y, small ivar(state) tvar(tm) gvar(first_treat_m) rolling(detrendm) graph
```

<br>

### [Example 4] Randomization Inference (RI) p-values

When using the `ri` option, the command performs a manual randomization inference procedure.
```
lwdid lcigsale, small ivar(state) tvar(year) gvar(first_year) rolling(detrend) ri 
```
By default, it runs `rireps(999)` replications and does __not__ set a random seed. Therefore, the reported RI p-value will vary slightly across runs.

The RI results are displayed immediately after the single-ATT regression output (see below):

<div align="center">
<img src="subfolder/ex5.png">
</div>

<br>

To reproduce the same RI p-value, specify a seed using `riseed()`. Simillarly you can specify the number of replications using `rireps()`.

```
lwdid lcigsale, small ivar(state) tvar(year) gvar(first_year) rolling(detrend) ri riseed(349139) rireps(999)
```
This specification produces __identical RI p-values__ each time the command is executed.

<br>

# Large-N case

First, load the dataset:
```
use lw_walmart, clear
```
> ⚠️ **Important:** For the large-N implementation, the `method()` option must be specified.  
>> Available options: <br>
>> `ra` (Regression Adjustment), `ipw` (Inverse Probability Weighting), `ipwra` (Doubly-robust IPWRA).

<br>

### [Example 1] Detrending transformation with the IPWRA estimator

Then run:
```
lwdid log_wholesale_emp x1 x2 x3, ivar(cid) tvar(year) gvar(first_year) ///
rolling(detrend) method(ipwra) graph 
```

Here, `ivar(cid)` identifies states. the dataset includes __multiple treatment groups__, corresponding to a **staggered adoption setting**. If there is only a single treated group, `lwdid` automatically detects this and applies the common timing procedure.

In this example, we use `rolling(detrend)` because the __parallel trends (PT)__ assumption appears to be violated, and apply the the doubly-robust estimator via `method(ipwra)`.

This command produces the default graph plotting **weighted ATTs by relative time `r` (time to treatment)**:

* `r = 0`: the immediate treatment effect
* `r ≥ 1`: post-treatment (dynamic) effects
* `r < 0`: pre-treatment periods

<div align="center">
<img src="subfolder/ex6.png">
</div>

<br>

But, similarly, you can customize your graphs with `gopts` and `scheme` options: more over you can save the results, and apply your own graph stayle to plot the event-study graphs.

<br>

### [Example 2] Customizing graphs and saving results

```stata
lwdid log_wholesale_emp x1 x2 x3, ivar(cid) tvar(year) gvar(first_year) ///
      rolling(detrend) method(ipwra) save(myresult) graph scheme(s2color) ///
      gopts( ytitle("WATT") xtitle("Time to Treatment(r)") title("The Effects of Walmart opening on ln(Wholesale Emp)") ///
				ylabel(-.4(.1).25, format(%3.1f) angle(horizontal)) ///
             graphregion(color(white)) plotregion(color(white)) )		
```
In this example:

* `gopts()` adds custom titles and axis labels.
* `scheme()` changes the graph color scheme.
> For a black-and-white version suitable for journal submissions, you can use scheme(s2mono).
> `graphregion(color(white)) plotregion(color(white))` replace the default gray background in `scheme(s2color)` with a clean white background, which is often preferred for presentations and publication-quality figures.

<div align="center">
<img src="subfolder/ex7.png">
</div>

By default, the graph reports **simultaneous confidence bands** for the event-study path. 

These bands are designed to provide uniform coverage across the reported event times $r$, rather than pointwise coverage for each $r$ separately.

The option `save(myresult)` also creates a new dataset, `myresult.dta`, in the working directory. This file contains:

* `ryear`: event time, $r=t-g$;
* `watt`: the weighted average treatment effect on the treated, $W\!ATT(r)$;
* `se`: the standard error for $WATT(r)$;
* `t_stat`: the pointwise test statistic, computed as `watt/se`;
* `p_value`: the pointwise normal-approximation p-value;
* `low_ci` and `up_ci`: the lower and upper bounds of the simultaneous confidence band;
* `N_cohort`: the number of treated cohorts contributing to each event-time estimate;
* `N_units`: the number of units included in the $W\!ATT(r)$ estimation sample.

Note that `t_stat` and `p_value` are **pointwise** quantities, while `low_ci` and `up_ci` correspond to the **simultaneous confidence band**.

<br>
<div align="center">
<img src="subfolder/ex8.png">
</div>



<br>

## Contact and Updates

Minor updates related to graph plotting (particularly for **Small-N: Staggered adoption** cases) are planned for upcoming releases — please stay tuned.


__For questions or suggestions, feel free to reach out to the authors:__

**[Soo Jeong Lee](https://sites.google.com/view/sjlee-econ/home)** ([soojeong.lee@siu.edu](mailto:soojeong.lee@siu.edu)) , Southern Illinois University Carbondale

**Jeffrey M. Wooldridge** ([wooldri1@msu.edu](mailto:wooldri1@msu.edu)), Michigan State University

<br><br><br><br> <br>


