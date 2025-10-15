
# lwdid
Stata package to implement he Rolling Estimation Method proposed in Lee and Wooldridge (2025).
[``Simple Approaches to Inference with Difference-in-Differences Estimators with Small Cross-Sectional Sample Sizes'']([https://pages.github.com/](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5325686)).


The current version 1.0 of lwdid is for small-N (Cross-sectional sample size) panel data with common treatment timing. 
A forthcoming version will extend the method to large-N panels and staggered intervention designs.
The command estimates period-by-period ATTs by transforming unit-level outcomes to remove pre-treatment means or trends.

## Contents
- [Install](##install)
- [Citation](##citation)
- [Syntax](##syntax)
- [Example](##Example)
--

## Install

To obtain the latest version through github, from the main command window in Stata, please run:
```
net describe lwdid, from(https://raw.githubusercontent.com/Soo-econ/lwdid/main/)
```

It will then show you
![The Screenshot of results window](https://raw.githubusercontent.com/Soo-econ/lwdid/main/subfolder/fig1.png)

Then, install, by run:
```
net install lwdid
```
To downlowd SMOKING,DTA dataset and replicate our empirical work in manuscript (or running example codes in [Example](##Example) section, run:
```
net get lwdid
```

To install directly from GitHub in Stata:
```
net install lwdid, from(https://raw.githubusercontent.com/Soo-econ/lwdid/main/)
net get lwdid, from(https://raw.githubusercontent.com/Soo-econ/lwdid/main/)
```
If the download from within Stata fails, you can manually download the files directly from https://raw.githubusercontent.com/Soo-econ/lwdid/main

## Citation

lwdid is a user-written Stata command freely available for academic and research use.
If you use this command in your work, please cite:

Lee, Soo Jeong, and Jeffrey M. Wooldridge (2025),
"Simple Approaches to Inference with Difference-in-Differences Estimators with Small Cross-Sectional Sample Sizes,"
Working Paper, Available at SSRN 5325686.

## Syntax
If you want to learn more about lwdid Syntax:
```
help lwdid
```
```
lwdid yvar dvar, ivar(id) tvar(year [quarter]) post(post) rolling(type)
       [vce(vartype) controls(varlist) table("filename") graph gopts(string)
        save(string) gid(id) ri rireps(#) riseed(#)]
```

### Options
__Main Variable__ 

| Argument | Description                |
| -------- | -------------------------- |
| `yvar`   | Outcome variable           |
| `dvar`   | Binary treatment indicator |



__Required Options__

| Option          | Description                                                         |
| --------------- | --------------------------------------------------------------------|
| `ivar(varname)` | Panel identifier (numeric or string)                                |
| `tvar(varlist)` | Time variable(s): `year` or `year quarter`                          |
| `post(varname)` | Post-period indicator (1 = post, 0 = pre)                           |
| `rolling(type)` | Transformation applied to residualize `yvar`:  4 possble types      |
                     [1] _demean_   : remove unit-specific pre-period average)          |
                     [2] _detreand_ : remove unit-specific linear pre-trend)            |
                     [3] _demeanq_ : unit-specific demeaning + deseasonalizing          |  
                     [4] _detreandq_ : unit-specific detrending + deseasonalizing       |
                     Note:[3]&[4] requires year-quarter data, with_tvar(year quarter)_  |
__Optional__

| Option                 | Description                                                 |
| ---------------------- | ----------------------------------------------------------- |
| `ivar(id)`             | Panel identifier                                            |
| `tvar(year [quarter])` | Time variable (supports quarterly data)                     |
| `post(post)`           | Post-treatment indicator variable                           |
| `rolling()`            | Specify transformation type: `demean`, `detrend`, etc.      |
| `controls(varlist)`    | Optional covariates included in cross-sectional regressions |
| `graph`                | Displays estimated treatment effects over time              |
| `vce(vartype)`         | Variance estimator (`robust`, `cluster id`, etc.)           |



## Example















