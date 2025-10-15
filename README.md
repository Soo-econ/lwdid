
# lwdid
Stata package implementing the Rolling Estimation Method proposed in
[Lee and Wooldridge (2025)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5325686).

The current version (v1.0) of `lwdid` is designed for panel data with small-N Cross-sectional sample sizes and common treatment timing. A forthcoming update will extend it to staggered intervention designs, and  large-N panels ([Lee and Wooldridge (2023)](https://dx.doi.org/10.2139/ssrn.4516518)

The command estimates a single post-treatment effect and period-by-period ATTs by transforming unit-level outcomes to remove pre-treatment means or trends (seaonality), following the rolling estimation framework.
<br>

## Contents
- [Installation](#installation)
- [Citation](#citation)
- [Syntax](#syntax)
- [Examples](#examples)
- [Contact and Updates](#Contact and Updates)

<br>

## Installation

To install the latest version from GitHub, type the following command in Stata:

```stata
net describe lwdid, from(https://raw.githubusercontent.com/Soo-econ/lwdid/main/)
```

 You should see a window like this:
 

![The Screenshot of results window](https://raw.githubusercontent.com/Soo-econ/lwdid/main/subfolder/ex1.png)

Then, install the package:
```
net install lwdid
```

To downlowd the accompanying example dataset (SMOKING.DTA) for the manuscript replication or examples in the [Example](##Example) section, run:
```
net get lwdid
```

Alternatively, you can install and download everything directly from GitHub in one step: 
```
net install lwdid, from(https://raw.githubusercontent.com/Soo-econ/lwdid/main/)
net get lwdid, from(https://raw.githubusercontent.com/Soo-econ/lwdid/main/)
```
If installation from within Stata fails, the files can be downloaded manually from 
- https://raw.githubusercontent.com/Soo-econ/lwdid/main

<br>

## Citation

lwdid is a user-written Stata command freely available for academic and research use.
If you use this command in your work, please cite:

Lee, Soo Jeong, and Jeffrey M. Wooldridge (2025),
"_Simple Approaches to Inference with Difference-in-Differences Estimators with Small Cross-Sectional Sample Sizes,"
Working Paper_, Available at [SSRN 5325686](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5325686).

<br>

## Syntax
```
lwdid yvar dvar, ivar(id) tvar(year [quarter]) post(post) rolling(type)
       [vce(vartype) controls(varlist) table("filename") graph gopts(string)
        save(string) gid(id) ri rireps(#) riseed(#)]
```

To view the full syntax and options, type:
```
help lwdid
```

### Option Details 

#### **Main Variables**

| Argument | Description                |
|-----------|----------------------------|
| `yvar`    | Outcome variable           |
| `dvar`    | Binary treatment indicator |

---

#### **Required Options**

| Option | Description |
|--------|--------------|
| `ivar(varname)` | Panel identifier (numeric or string) |
| `tvar(varlist)` | Time variable(s): `year` or `year quarter` |
| `post(varname)` | Post-period indicator (1 = post, 0 = pre) |
| `rolling(type)` | Transformation type applied to residualize `yvar`. <br><br>**Available types:**<br>• `demean` – remove unit-specific pre-period mean <br>• `detrend` – remove unit-specific linear pre-trend <br>• `demeanq` – demeaning + deseasonalizing *(requires quarterly data)* <br>• `detrendq` – detrending + deseasonalizing *(requires quarterly data)* <br><br>**Note:** `demeanq` and `detrendq` require `tvar(year quarter)` format. |

---

#### **Optional Options**

| Option | Description |
|--------|--------------|
| `controls(varlist)` | Additional covariates included in the cross-sectional regressions |
| `vce(vartype)` | Variance estimator (`robust`, `cluster id`, etc.) |
| `graph` | Displays estimated treatment effects over post-treatment periods |
| `gopts(string)` | Graph options passed directly to the `twoway` command |
| `table("filename")` | Saves period-by-period ATT estimates to a table (e.g., `"att_results.txt"`) |
| `save(string)` | Saves estimation results under a specified name |
| `gid(id)` | Group identifier variable (used for labeling or subgroup estimation) |
| `ri` | Requests Randomization Inference (RI) p-values |
| `rireps(#)` | Number of RI repetitions (default = 1000) |
| `riseed(#)` | Seed for reproducible RI results |

---

<br>

## Examples

First, load the dataset:
```
use smoking, clear
```

<br>

### [Example 1] Detrended transformation and graph
```
lwdid lcigsale d,ivar(state) tvar(year) post(post) rolling(detrend) graph
```
By selecting `rolling(detrend)`, the command performs a __unit-specific detrending tranformation__, allowing each unit to follow its own rend.

You can easily change this to `demean`, or with year-quarter data, choose `demeanq` or `detrendq` to account for __seaonality__.

By default, `lwdid` reports both the __overall (single)__ treatment effect and the __period-by-period__ ATT estimates.


<br>


### [Example 2] Customizing graphs using `gopts()` 

The default graph (from Example 1 using the `graph` option) produces a plot comparing residualized treated and control units outcomes:

![](https://raw.githubusercontent.com/Soo-econ/lwdid/main/subfolder/ex2.png)

You can __easily customize the graph__ through the `gopts()` option, allowing full customization of titles, axes, legends, and colors. For example:
```
lwdid lcigsale d, ivar(state) tvar(year) post(post) rolling(detrend) ///
      graph gopts(ytitle("Residualized average outcome") xtitle("Year") ///
                  legend(pos(1) ring(0)) ///
                  title("The Effects of California’s Tobacco Control Program"))

```
__This produces:__

![](https://raw.githubusercontent.com/Soo-econ/lwdid/main/subfolder/ex3.png)


<br>

### [Example 3] Year-quarter data

As mentioned earlier, using `rolling(demeanq)` or `rolling(detrendq)` for quarterly data automatically remove __seaonal effects__, and seperately estimates treatment effects across year-quarter periods in the result table.

To apply these transformations, the time variable must include __two identifiers__ in the form `tvar(year quarter)`.

For example:

```
lwdid y d, ivar(id) tvar(year qt) post(post) rolling(detrendq) graph
```
This specification performs __a quarter-specific detrending transformation__ and plots the residualized treated and control series over year-quarter time.



### [Example 4] Randomization Inference (RI) p-values

When using the `ri` option, the command performs a manual randomization inference procedure.
```
lwdid lcigsale d,ivar(state) tvar(year) post(post) rolling(detrend) ri 
```
By default, it runs `rireps(1000)` replications and does __not__ set a random seed -- meaning that the computed RI p-value will vary slightly each time you run the command.


The calculated RI p-value will be reported at __the end of the results__ (see below).

![](https://raw.githubusercontent.com/Soo-econ/lwdid/main/subfolder/ex4.png)


If you want fully reproducible results, you can specify the seed using `riseed()`. Simillarly you can specify the number of replications using `rireps()`.

```
lwdid lcigsale d,ivar(state) tvar(year) post(post) rolling(detrend) ri riseed(123) rireps(2000)
```
This setup will generate __identical RI p-values__ each time the command is executed, with 2,000 randomization replications and a fixed seed of 123.

<br>

## Contact and Updates

This package is __actively being updated__ with new features and extensions
(e.g., staggered treatment timing and large-N implementation).

Stay tuned for upcoming releases and documentation updates.

For questions or suggestions, feel free to reach out to the authors:

**Soo Jeong Lee** ([soojeong.lee@siu.edu](mailto:soojeong.lee@siu.edu))  
**Jeffrey M. Wooldridge** ([wooldri1@msu.edu](mailto:wooldri1@msu.edu))

<br><br>


