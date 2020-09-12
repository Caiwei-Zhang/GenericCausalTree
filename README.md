# GenericCausalTree
Generic Causal Tree: Inference on Subgroup-Level Treatment Effects in Observational Studies

## Brief Intro
An implement for tree-based algorithm to do heterogeneous treatment effects estimation and inference with observational data. 

## Required Packages
All the required packages will be installed and loaded by the file `./functions/library.R`.

## Usage
Build dependent environment for constructing Generic Casual Tree:
```
setwd("../functions")
fileName <- list.files(getwd())
fileList <- paste(getwd(), fileName, sep = “/”)
for (file in fileList) {
  source(file)
}
```

Use `Treebuild()` to construct the regular version of GCT, and use `honest.Treebuild()` to construct the honest version Tree. 

If you want to make inference on the treatment effects, please use function `honest.Treebuild()`, which provides confidence interval in the *frame* of the *optree*.
An example of how to use function `Treebuild()` and `honest.Treebuild()` is given below:

```
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
A data analysis example is offered in `./Data analysis/`. We build a generic causal tree to find the heterogeneity of racial disparities between American Whites and American Africans. The data we use is from Treatment Episode Data Set: Admissions 2015 (TEDS-A-2015), a program that collects the information of admissions to substance abuse treatment occurring in 2015. TEDS-A-2015 is available on the website of [Substance Abuse & Mental Health Data Archive (SAMHDA)](https://www.datafiles.samhsa.gov/).
