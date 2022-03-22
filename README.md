## Introduction
A tree-based algorithm to implement heterogeneous causal inference with observational data. 

## Environment setup
**Required packages** will be installed and loaded by `./Functions/library.R`.

Load **basic functions** that achieve construction of tree:
```
setwd("..\\Functions")
functions <- list.files(getwd())
functions <- paste(paste(getwd(), "\\", sep = ""), functions, sep = "")
for (f in functions) {
  source(i)
}

## Usage
```
Use `Treebuild()` to construct regular version of GCT, and use `honest.Treebuild()` to construct honest version of GCT. 
Make inference of subgroup-level treatment effects by using function `Treebuild.GCTwithVar()` or `honest.Treebuild()`, which provides confidence interval by regular and honest version of GCT respectively.

Here is an example showing how to use `Treebuild()` and `honest.Treebuild()`:
```R
gct <- Treebuild(data = data, form.outcome = "Y ~ .", form.pscore = "A ~ .", 
                 type.outcome = "continuous", train.frac = 0.8, 
                 modelFitBefore = TRUE, outcome.mthd = "GBM", pscore.mthd = "GBM", 
                 cp = qchisq(0.95, 1), minsplit = 60, minsize = 30)

opt_gct <- gct$optree
ste_nonhonest <- opt_gct$frame$yval[which(opt_gct$frame[,1] == "leaf")]
fancyRpartPlot(opt_gct, type = 4, palettes = "Greens",
               sub = "optimal tree")

honest_gct <- honest.Treebuild(data = data, form.outcome = "Y ~ .", form.pscore = "A ~ .", 
                               type.outcome = "continuous", train.frac = 0.8,
                               modelFitBefore = TRUE, est.modelFitBefore = TRUE,
                               outcome.mthd = "GBM", pscore.mthd = "GBM", cp = qchisq(0.95, 1), 
                               minsplit = 60, minsize = 30)

opt_honest_gct <- honest_gct$optree
ste_honest <- opt_honest_gct$frame$yval[which(opt__honest_gct$frame[,1] == "leaf")]
CI <- cbind(opt_honest_gct$frame$ci.lower[which(opt__honest_gct$frame[,1] == "leaf")],
            opt_honest_gct$frame$ci.upper[which(opt__honest_gct$frame[,1] == "leaf")]) 
fancyRpartPlot(opt_honest_gct, type = 4, palettes = "Greens", 
               sub = "honest optimal tree")
```
## Data Analysis
A data analysis example is offered in `./Data analysis/`. We build a generic causal tree to find the heterogeneity of racial disparities between American Whites and American Africans. The data we use is from Treatment Episode Data Set: Admissions 2015 (TEDS-A-2015), a program that collects the information of admissions to substance abuse treatment occurring in 2015. TEDS-A-2015 is available on the website of Substance Abuse & Mental Health Data Archive (SAMHDA).
