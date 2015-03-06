# cv.lm
Cross Validation in Linear Model, and the subset models with various criteria along with ridge regression, PLS and PCR

The function takes the parameter `dataSet`, `x.var`, `y.var`, `step=FALSE`, `criteria=NULL`, `split=12`, where `x.var` is a character vector of x variables, `y.var` is a chacter vector with one element representing predictor variable.

Set `step` to true if you want to have stepwise subset model where you have to specify the `criteria`. The `criteria` can take one of the values among "AIC", "BIC", "Cp", "R2adj", "forward" and  "backward" as a character. `split` is the number of split you want in your dataset during cross-validation. Currently, the split will split the dataset into consecutive segments. I will later add some more alternatives.

The function is based on another function called `makeForumla` to create the linear model formula. Before running the function `leaps` and `mixlm` packages should be installed.

Added ability to compute `RMSEP` and `R2pred` for `linearRidge` from `ridge` package with automatic selection of ridge parameter, `plsr` and `pcr` from `pls` package.

<span style="color:red">The results from this function for PLS and PCR are not same as `plsr` function with cross-validation enabeled. I am still configuring the problem.</span>