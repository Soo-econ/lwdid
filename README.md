
# lwdid

A Stata package that implements the **Rolling Difference-in-Differences Estimator** proposed by Lee and Wooldridge ([2025a](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4516518), [2025b](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5325686)). It provides fast and flexible estimation for staggered and common timing treatment settings in panel data, and offers a unified implementation that covers both standard large-N asymptotic settings and cases with small cross-sectional sample sizes, where conventional large-N inference may not be reliable.

`lwdid` is a user-written Stata command freely available for academic and research use. 
> A manuscript describing the `lwdid` method and software in detail is currently in preparation: <br> Soo Jeong Lee and Jeffrey M. Wooldridge (2025c), “*lwdid: Rolling Difference-in-Differences Estimation for Small and Large Panels,*” working paper.

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


# Contents
- [Citation](#citation)
- [Installation](#installation)
- [Syntax](#syntax)
- [Examples: Large-N](#large-n-case)
- [Examples: Small-N](#small-n-case)
- [Contact and Updates](#contact-and-updates)

<br>

# Installation

To install the latest version from GitHub, type the following command in Stata:

```stata
net describe lwdid, from(https://raw.githubusercontent.com/Soo-econ/lwdid/main/)
```

 You should see a window like this:
 
![The Screenshot of results window](https://raw.githubusercontent.com/Soo-econ/lwdid/main/subfolder/ex1.png)

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
| `ri` |  **Small-N Only**. Requests randomization inference (RI) <br>•`rireps(#)` –Number of RI repetitions (default = 999).<br>•`riseed(#)` – Seed for reproducible RI (default: randomly drawn; results may differ across runs if not specified)|

---

<br>

# Examples
# Large-N case

# Small-N case
First, load the dataset:
```
use smoking, clear
```

## \[Example 1\] Detrended transformation and graph
```
lwdid lcigsale, ivar(state) tvar(year) gvar(first_year) post(post) rolling(detrend) graph
```
By specifying `rolling(detrend)`, the command performs a __unit-specific detrending tranformation__, allowing each unit to follow its own heterogeneous linear trend and thereby relaxing the parallel trends (PT) assumption.

Alternatively, you can use

```stata
rolling(demean)
```

to remove **unit-specific means** (when PT assumption hold) instead of **unit-specific trends**.

For year–quarter data, the options `rolling(demeanq)` or `rolling(detrendq)` can be used to account for __seaonality__. 

To apply these transformations, the time variable must include __two identifiers__ in the form `tvar(year quarter)`.

For example:
```
lwdid y, ivar(state) tvar(year q) gvar(first_year) rolling(detrendq) graph
```
This plots the residualized treated and control series over year-quarter time.


By default, `lwdid` reports both:

* the **overall (single) treatment effect**, and
* **period-by-period ATT estimates**.

<br>


## [Example 2] Customizing graphs using `gopts()` 

The default graph (from Example 1 using the `graph` option) produces a plot comparing **residualized ourcomes of treated and control units:

![](https://raw.githubusercontent.com/Soo-econ/lwdid/main/subfolder/ex2.png)

You can __easily customize the graph__ through the `gopts()` option, which allows full customization of titles, axes, legends, and other graphical elements. For example:
```
lwdid lcigsale, ivar(state) tvar(year) gvar(first_year) rolling(detrend) ///
      graph gopts(ytitle("Residualized average outcome") xtitle("Year") ///
                  legend(pos(1) ring(0)) ///
                  title("The Effects of California’s Tobacco Control Program"))
```
__This produces:__

![](https://raw.githubusercontent.com/Soo-econ/lwdid/main/subfolder/ex3.png)

To generate a **black-and-white version** suitable for journal submissions, you can use the `scheme(s1mono)` option:
![](https://raw.githubusercontent.com/Soo-econ/lwdid/main/subfolder/ex2.png)


<br>

### [Example 3] Randomization Inference (RI) p-values

When using the `ri` option, the command performs a manual randomization inference procedure.
```
lwdid lcigsale, ivar(state) tvar(year) gvar(first_year) rolling(detrend) ri 
```
By default, it runs `rireps(1000)` replications and does __not__ set a random seed -- meaning that the computed RI p-value will vary slightly each time you run the command.


The calculated RI p-value will be reported at __the end of the results__ (see below).

![](https://raw.githubusercontent.com/Soo-econ/lwdid/main/subfolder/ex4.png)


<br>

If you want fully reproducible results, you can specify the seed using `riseed()`. Simillarly you can specify the number of replications using `rireps()`.

```
lwdid lcigsale d,ivar(state) tvar(year) post(post) rolling(detrend) ri riseed(123) rireps(2000)
```
This setup will generate __identical RI p-values__ each time the command is executed, with 2,000 randomization replications and a fixed seed of 123.

<br>


## Contact and Updates

The `lwdid` package is now __largely updated__. The current version can replicate the main results reported in the accompanying paper.

Minor updates related to graph plotting (particularly for **Small-N: Staggered adoption** cases) are planned for upcoming releases — please stay tuned.

For questions or suggestions, feel free to reach out to the authors:

**Soo Jeong Lee** ([soojeong.lee@siu.edu](mailto:soojeong.lee@siu.edu)) , Southern Illinois University Carbondale

**Jeffrey M. Wooldridge** ([wooldri1@msu.edu](mailto:wooldri1@msu.edu)), Michigan State University

<br><br><br><br>


