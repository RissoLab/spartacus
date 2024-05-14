# **SPA**tially **R**esolved **T**r**A**nscriptomics **C**lustering **U**sing **S**urroundings (SpaRTaCUS) <img src="logo2.jpeg" align="right" width="115" />

</br>



## Installation instructions

Install the development version from
[GitHub](https://https://github.com/RissoLab/spartacus) with:

``` r
remotes::install_github("RissoLab/spartacus")
```

## Run the model

Let `x` be the spatial experiment matrix containing the expression of `nrow(x)` genes measured over `ncol(x)` spots. The spatial coordinates of the spots are stored in the matrix `coordinates`. By default, the sparsity of the estimation is set to 'n.neighbors' = 20. You can run SpaRTaCUS and search for K gene clusters and R spot clusters, execute the provided code:

``` r
library(spartacus)
spartacus(x = x, coordinates = coordinates, K = K, R = R) 
```

## Starting points

The estimation algorithm is automatically run using random starting points. However, it is also possible to start from a previous output of the `spartacus` function:

``` r 
output1 <- spartacus(x = x, coordinates = coordinates, K = K, R = R)
output2 <- spartacus(x = x, coordinates = coordinates, input.values = output1)
```

## Convergence

By default, the estimation procedure is run for at most `max.iter` iterations, but it is previously stopped if a certain convergence criterion is reached (`conv.criterion`). If `conv.criterion = NULL`, there is no ending condition and the  procedure is run for `max.iter` iterations, otherwise it is stopped when the increment of the classification log-likelihood is smaller than a certain threshold `conv.criterion$epsilon` for `conv.criterion$iterations` times in a row. 
