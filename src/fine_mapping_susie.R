# remotes::install_github("stephenslab/susieR")
library(dplyr)
library(susieR)


data(N3finemapping)
attach(N3finemapping)

b <- true_coef[,1]
plot(b, pch=16, ylab='effect size')

sumstats <- univariate_regression(X, Y[,1])
z_scores <- sumstats$betahat / sumstats$sebetahat
susie_plot(z_scores, y = "z", b=b)
