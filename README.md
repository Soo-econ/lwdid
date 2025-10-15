
# lwdid
Stata command implementing the Rolling Estimation Method proposed in Lee and Wooldridge (2025)
[``Simple Approaches to Inference with Difference-in-Differences Estimators with Small Cross-Sectional Sample Sizes'']([https://pages.github.com/](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5325686)).

The current version 1.0 of lwdid is for small-N panel data with common treatment timing. A forthcoming version will extend the method to large-N panels and staggered intervention designs.
The command estimates period-by-period ATTs by transforming unit-level outcomes to remove pre-treatment means or trends.

- [Install](##install)
- [Citation](##citation)

---

## Install

```
net describe lwdid, from(https://raw.githubusercontent.com/Soo-econ/lwdid/main/)
```
This will show you
![The Screenshot of results window](https://raw.githubusercontent.com/Soo-econ/lwdid/main/subfolder/fig1.png)
Then, install:
```
net install lwdid
net get lwdid
```

or directly
```
net install lwdid, from(https://raw.githubusercontent.com/Soo-econ/lwdid/main/)
net get lwdid, from(https://raw.githubusercontent.com/Soo-econ/lwdid/main/)
```
If you want to learn more about lwdid Syntax:
```
help lwdid
```

If the download from within Stata fails, you can manually download the files directly from https://raw.githubusercontent.com/Soo-econ/lwdid/main
  

## Citation

lwdid is a user-written Stata command freely available for academic and research use.
If you use this command in your work, please cite:

Lee, Soo Jeong, and Jeffrey M. Wooldridge (2025),
"Simple Approaches to Inference with Difference-in-Differences Estimators with Small Cross-Sectional Sample Sizes,"
Working Paper, Available at SSRN 5325686.























