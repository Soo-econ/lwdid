
# lwdid

A Stata package that implements the **Rolling Difference-in-Differences Estimator** proposed by Lee and Wooldridge ([2025a](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4516518), [2025b](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5325686)). It provides fast and flexible estimation for staggered and common timing treatment settings in panel data, and offers a unified implementation that covers both standard large-N asymptotic settings and cases with small cross-sectional sample sizes, where conventional large-N inference may not be reliable.

`lwdid` is a user-written Stata command freely available for academic and research use. 
> A manuscript describing the `lwdid` method and software in detail is currently in preparation: <br> Soo Jeong Lee and Jeffrey M. Wooldridge (2025c), “*lwdid: Rolling Difference-in-Differences Estimation for Panel data,*” working paper.

<br>

# Citation

If you use `lwdid` in your research, please cite the following papers:

**Large-N Procedure**<br>

Soo Jeong Lee, and Jeffrey M. Wooldridge (2025a),
"_A Simple Transformation Approach to Difference-in-Differences Estimation for Panel Data,"
Working Paper_, Available at [SSRN 4516518](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4516518).
> Revise and resubmit at *Journal of Business & Economic Statistics*.

**Small-N Procedure** <br>
Soo Jeong Lee, and Jeffrey M. Wooldridge (2025b),
"_Simple Approaches to Inference with Difference-in-Differences Estimators with Small Cross-Sectional Sample Sizes,"
Working Paper_, Available at [SSRN 5325686](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5325686).
> Under review

<br>

# Contents
- [Citation](#citation)
- [Installation](#installation)
- [Syntax](#syntax)
- [Examples: Small-N](#small-n-case)
- [Examples: Large-N](#large-n-case)
- [Contact and Updates](#contact-and-updates)

<br>

# Installation

To install the latest version from GitHub, type the following command in Stata:

```stata
net describe lwdid, from(https://raw.githubusercontent.com/Soo-econ/lwdid/main/)
```

 You should see a window like this:
 
![The Screenshot of results window](subfolder/ex1.png)

Then, install the package:
```
net install lwdid, replace
```
The `replace` option overwrites any previously installed version of **lwdid**.

To downlowd the accompanying example datasets for the manuscript replication or examples in the  [Examples](#examples) section, run:
```
net get lwdid, replace
```

Alternatively, you can install and download everything directly from GitHub in one step: 
```
net install lwdid, from(https://raw.githubusercontent.com/Soo-econ/lwdid/main/)
net get lwdid, from(https://raw.githubusercontent.com/Soo-econ/lwdid/main/)
```
If installation from within Stata fails, the files can be downloaded manually from [my GitHub page](https://github.com/Soo-econ/lwdid.git)


<br>

# Syntax
```
lwdid yvar [covariates], ivar(varname) tvar(varlist) gvar(varname) rolling(type)
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
| `tvar(varlist)` | Time variable(s). Either `tvar(year)` or `tvar(year quarter)` for quarterly data |
| `gvar(varname)` | Treatment cohort variable (first treated period). Never‐treated units should be coded as 0.|
| `rolling(type)` | Transformation type applied to residualize `yvar`. <br><br>**Available types:**<br>• `demean` – remove unit-specific pre-period mean <br>• `detrend` – remove unit-specific linear pre-trend <br>• `demeanq` – demeaning + deseasonalizing *(requires quarterly data)* <br>• `detrendq` – detrending + deseasonalizing *(requires quarterly data)* <br><br>**Note:** `demeanq` and `detrendq` require `tvar(year quarter)` format. |

### **Required Estimation Options (depends on implementation)**
| Option | Description |
|--------|--------------|
| `small` | **Small-N Only** <br> Requests the small-sample inference procedure designed for settings with few treated (or control) units (default: large-N inference).|
|`method(ra\|ipw\|ipwra)`| <br>**Large-N Only: required unless `small` is specified**<br>•`method(ra)` –Regression Adjustment estimator <br>•`method(ipw)` – Inverse Probability Weighting estimator. <br>•`method(ipwra)` – Doubly‑robust IPW‑RA estimator.<br><br>|

---

## **Optional Options**

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
  
Here, `gvar(first_year)` specifies the first treatment year for each unit (never-treated units should be coded as `0`).  
In this example, only California is treated in 1989, so this is a **single-treatment case**. If multiple units are treated in different years, `lwdid` automatically detects this and applies the **small-N staggered adoption** procedure.

By specifying `rolling(detrend)`, the command performs a __unit-specific detrending tranformation__, allowing each unit to follow its own heterogeneous linear trend and thereby relaxing the parallel trends (PT) assumption.

Alternatively, you can use

```stata
rolling(demean)
```

to remove **unit-specific means** (when the PT assumption hold), instead of **unit-specific trends**.

## \[Example 2\] Using Quarterly data 

For year–quarter data, the options `rolling(demeanq)` or `rolling(detrendq)` can be used to account for __seaonality__. 

To apply these transformations, the time variable must include __two identifiers__ in the form `tvar(year quarter)`.

For example:
```
lwdid y, ivar(state) tvar(year q) gvar(first_year_q) rolling(detrendq) graph
```
This plots the residualized treated and control series over year-quarter time.

>  ⚠️ **Important:** When using quarterly time variables (`tvar(year q)`),
the treatment timing variable in `gvar()` must be defined as __a Stata quarterly date__ (format `%tq`).

For example, if the treatment occurs in a specific quarter:
```stata
gen first_treat_q = yq(first_year, first_quarter)
format first_treat_q %tq
```
When the detrendq option is used, both time variables must be specified in `tvar()`, namely `year` and `q`. 
Internally, the command uses `tvar(year q)` together with `gvar(first_year)` to identify post-treatment periods:
```stata
post=(tvar(year q) >= gvar(first_year))
```
For example:
|id| year  |q| first_year   | post |
|----| ------|--- | ------ | ---- |
|California| 2008 | 2 | 2008q3 | 0    |
|California| 2008 | 3 | 2008q3 | 1    |
|California| 2008 | 4 | 2008q3 | 1    |


<br>

## [Example 3] Customizing graphs using `gopts()` and `scheme()`

The default graph (from Example 1 using the `graph` option) produces a plot comparing **residualized ourcomes of treated and control units**:

![](subfolder/ex2.png)

You can __easily customize the graph__ through the `gopts()` option, which allows full customization of titles, axes, legends, and other graphical elements. For example:
```
lwdid lcigsale, small ivar(state) tvar(year) gvar(first_year) rolling(detrend) ///
      graph gopts(ytitle("Residualized average outcome") xtitle("Year") ///
      legend(pos(1) ring(0)) ///
      title("The Effects of California’s Tobacco Control Program"))
```
__This produces:__

![](subfolder/ex3.png)

To generate a **black-and-white version** suitable for journal submissions, you can use the `scheme(s1mono)` option:

![](subfolder/ex4.png)

<br>

### [Example 4] Randomization Inference (RI) p-values

When using the `ri` option, the command performs a manual randomization inference procedure.
```
lwdid lcigsale, ivar(state) tvar(year) gvar(first_year) rolling(detrend) ri 
```
By default, it runs `rireps(1000)` replications and does __not__ set a random seed -- meaning that the computed RI p-value will vary slightly each time you run the command.


The calculated RI p-value will be reported at __the end of the results__ (see below):

<div align="center">
<img src="subfolder/ex5.png">
</div>

<br>

If you want fully reproducible results, you can specify the seed using `riseed()`. Simillarly you can specify the number of replications using `rireps()`.

```
lwdid lcigsale, small ivar(state) tvar(year) gvar(first_year) rolling(detrend) ri riseed(123) rireps(2000)
```
This specification produces __identical RI p-values__ each time the command is executed, with 2,000 randomization replications and a fixed seed of 123.

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
lwdid log_wholesale_emp x1 x2 x3, ivar(fips) tvar(year) gvar(first_year) ///
rolling(detrend) method(ipwra) graph 
```

Here, `ivar(fips)` identifies states. the dataset includes __multiple treatment groups__, corresponding to a **staggered adoption setting**. If there is only a single treated group, `lwdid` automatically detects this and applies the common timing procedure.

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
lwdid log_wholesale_emp x1 x2 x3, ivar(fips) tvar(year) gvar(first_year)  ///
      rolling(detrend) method(ipwra) save(myresult) graph scheme(s1color) ///
      gopts(ytitle("WATT") xtitle("Time to Treatment(r)") title("The Effects of Walmart opening on ln(Wholesale Emp)"))  
```
In this example:

* `gopts()` adds custom titles and axis labels.
* `scheme(s1color)` changes the graph color scheme.
> For a black-and-white version suitable for journal submissions, you can use `scheme(s1mono)`.

<div align="center">
<img src="subfolder/ex7.png">
</div>

The option `save(myresult)` also creates a new dataset `myresult.dta` in your working directory.
This file contains:

* weighted ATT estimates across relative time,
* corresponding **standard errors** and **confidence intervals** (computed via wild bootstrap),
* `N_cohort`: the number of treated cohorts used compute the estimates.
* `N_units`: the total number of units included in the WATT estimation sample.

<br>
<div align="center">
<img src="subfolder/ex8.png">
</div>


<br>

## Contact and Updates

The `lwdid` package has been updated to **version 2.0**, implementing the procedures described in the accompanying papers and replicating their main results.

Minor updates related to graph plotting (particularly for **Small-N: Staggered adoption** cases) are planned for upcoming releases — please stay tuned.


__For questions or suggestions, feel free to reach out to the authors:__

**[Soo Jeong Lee](https://sites.google.com/view/sjlee-econ/home)** ([soojeong.lee@siu.edu](mailto:soojeong.lee@siu.edu)) , Southern Illinois University Carbondale

**Jeffrey M. Wooldridge** ([wooldri1@msu.edu](mailto:wooldri1@msu.edu)), Michigan State University

<br><br><br><br> <br>


