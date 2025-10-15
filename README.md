
# lwdid
Stata package implementing the Rolling Estimation Method proposed in
[Lee and Wooldridge (2025)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5325686).

The current version (v1.0) of `lwdid` is designed for panel data with small-N Cross-sectional sample sizes and common treatment timing. A forthcoming update will extend it to staggered intervention designs, and  large-N panels ([Lee and Wooldridge (2023)](https://dx.doi.org/10.2139/ssrn.4516518)

The command estimates a single post-treatment effect and period-by-period ATTs by transforming unit-level outcomes to remove pre-treatment means or trends (seaonality), following the rolling estimation framework.

## Contents
- [Installation](#installation)
- [Citation](#citation)
- [Syntax](#syntax)
- [Example](#example)


## Installation

To install the latest version from GitHub, type the following command in Stata:

```stata
net describe lwdid, from(https://raw.githubusercontent.com/Soo-econ/lwdid/main/)
```

 You should see a window like this:
 

![The Screenshot of results window](https://raw.githubusercontent.com/Soo-econ/lwdid/main/subfolder/fig1.png)

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

## Citation

lwdid is a user-written Stata command freely available for academic and research use.
If you use this command in your work, please cite:

Lee, Soo Jeong, and Jeffrey M. Wooldridge (2025),
"Simple Approaches to Inference with Difference-in-Differences Estimators with Small Cross-Sectional Sample Sizes,"
Working Paper, Available at [SSRN 5325686](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5325686).

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


## Example















