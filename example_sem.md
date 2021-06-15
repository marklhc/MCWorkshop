Simulation Example on Structural Equation Modeling
================
Mark Lai
April 29, 2019, last updated by Winnie Tse on June 14, 2021

  - [Simulation Example on Structural Equation Modeling (SEM) Using the
    `SimDesign`
    Package](#simulation-example-on-structural-equation-modeling-sem-using-the-simdesign-package)
      - [Simulate Multivariate Data](#simulate-multivariate-data)
      - [Workflow for Using `SimDesign`](#workflow-for-using-simdesign)
          - [Step 1: Create data frame of design
            factors](#step-1-create-data-frame-of-design-factors)
          - [Step 2: Define a function for generating
            data](#step-2-define-a-function-for-generating-data)
          - [Step 3: Define a function for analyzing the simulated
            data](#step-3-define-a-function-for-analyzing-the-simulated-data)
          - [Step 4: Define a function for evaluating the sample
            results](#step-4-define-a-function-for-evaluating-the-sample-results)
          - [Step 5: Test Run the
            simulation](#step-5-test-run-the-simulation)
          - [Step 6: Full simulation](#step-6-full-simulation)
      - [Summarize Simulation Results](#summarize-simulation-results)
          - [ANOVA](#anova)
      - [Exercise](#exercise)

# Simulation Example on Structural Equation Modeling (SEM) Using the `SimDesign` Package

Recently, the `SimDesign` package was developed so that designing and
running simulation studies can be more structured and organized. The
package also provides some great features such as parallel computing,
fail-safe stopping, gathering of error or warning messages, among
others. Check out the paper by Sigal & Chalmers (2016) as well as the
package vignettes
(<https://cran.r-project.org/web/packages/SimDesign/index.html>) for
more information. Below is a quick hands-on example for using the
package.

``` r
# Load required packages
library(tidyverse)
theme_set(theme_classic() +
            theme(panel.grid.major.y = element_line(color = "grey92")))
library(SimDesign)
library(mnormt)
library(lavaan)
```

## Simulate Multivariate Data

In SEM, when multivariate normality is assumed, one can either generate
data directly using matrix algebra, or generate the latent variables
first before generating the observed variables. The first method is
faster, but the second method is more general and can be applied to
situations like categorical data or multilevel data. Therefore, in this
note we’ll use the second method.

Let’s do a latent growth model (LGM) similar to the one in the note
“Simulating Multilevel Data.” Here is the model in a latent growth
model representation:   
![\\begin{bmatrix}&#10; y\_{0i} \\\\&#10; y\_{1i} \\\\&#10; y\_{2i}
\\\\&#10; y\_{3i}&#10; \\end{bmatrix} = &#10;
\\boldsymbol{\\mathbf{\\Lambda }}&#10; \\begin{bmatrix}&#10; \\eta\_{1i}
\\\\&#10; \\eta\_{2i}&#10; \\end{bmatrix} + &#10; \\begin{bmatrix}&#10;
e\_{0i} \\\\&#10; e\_{1i} \\\\&#10; e\_{2i} \\\\&#10; e\_{3i}&#10;
\\end{bmatrix}](https://latex.codecogs.com/png.latex?%5Cbegin%7Bbmatrix%7D%0A%20%20%20%20y_%7B0i%7D%20%5C%5C%0A%20%20%20%20y_%7B1i%7D%20%5C%5C%0A%20%20%20%20y_%7B2i%7D%20%5C%5C%0A%20%20%20%20y_%7B3i%7D%0A%20%20%5Cend%7Bbmatrix%7D%20%3D%20%0A%20%20%5Cboldsymbol%7B%5Cmathbf%7B%5CLambda%20%7D%7D%0A%20%20%5Cbegin%7Bbmatrix%7D%0A%20%20%20%20%5Ceta_%7B1i%7D%20%5C%5C%0A%20%20%20%20%5Ceta_%7B2i%7D%0A%20%20%5Cend%7Bbmatrix%7D%20%2B%20%0A%20%20%5Cbegin%7Bbmatrix%7D%0A%20%20%20%20e_%7B0i%7D%20%5C%5C%0A%20%20%20%20e_%7B1i%7D%20%5C%5C%0A%20%20%20%20e_%7B2i%7D%20%5C%5C%0A%20%20%20%20e_%7B3i%7D%0A%20%20%5Cend%7Bbmatrix%7D
"\\begin{bmatrix}
    y_{0i} \\\\
    y_{1i} \\\\
    y_{2i} \\\\
    y_{3i}
  \\end{bmatrix} = 
  \\boldsymbol{\\mathbf{\\Lambda }}
  \\begin{bmatrix}
    \\eta_{1i} \\\\
    \\eta_{2i}
  \\end{bmatrix} + 
  \\begin{bmatrix}
    e_{0i} \\\\
    e_{1i} \\\\
    e_{2i} \\\\
    e_{3i}
  \\end{bmatrix}")  
where ![y\_{0i}, \\ldots,
y\_{3i}](https://latex.codecogs.com/png.latex?y_%7B0i%7D%2C%20%5Cldots%2C%20y_%7B3i%7D
"y_{0i}, \\ldots, y_{3i}") are the outcome values for person
![i](https://latex.codecogs.com/png.latex?i "i") from time 0 to time 3,
![\\eta\_{1i}](https://latex.codecogs.com/png.latex?%5Ceta_%7B1i%7D
"\\eta_{1i}") is the specific intercept for person
![i](https://latex.codecogs.com/png.latex?i "i"),
![\\eta\_{2i}](https://latex.codecogs.com/png.latex?%5Ceta_%7B2i%7D
"\\eta_{2i}") is the specific slope for person
![i](https://latex.codecogs.com/png.latex?i "i"), and ![e\_{0i},
\\ldots,
e\_{3i}](https://latex.codecogs.com/png.latex?e_%7B0i%7D%2C%20%5Cldots%2C%20e_%7B3i%7D
"e_{0i}, \\ldots, e_{3i}") are the within-person level error term. The
distributional assumptions are

  
![&#10; \\begin{aligned}&#10; \\begin{bmatrix}&#10; \\eta\_{1i}
\\\\&#10; \\eta\_{2i}&#10; \\end{bmatrix} & \\sim &#10;
\\mathcal{N}\\left(\\begin{bmatrix}&#10; \\alpha\_1 \\\\&#10;
\\alpha\_2&#10; \\end{bmatrix}, &#10; \\begin{bmatrix}&#10; \\phi\_{11}
& \\phi\_{21} \\\\&#10; \\phi\_{21} & \\phi\_{22}&#10;
\\end{bmatrix}\\right) \\\\&#10; \\begin{bmatrix}&#10; e\_{0i} \\\\&#10;
e\_{1i} \\\\&#10; e\_{2i} \\\\&#10; e\_{3i}&#10; \\end{bmatrix} &
\\sim&#10; \\mathcal{N}\\left(\\begin{bmatrix}&#10; 0 \\\\&#10; 0
\\\\&#10; 0 \\\\&#10; 0&#10; \\end{bmatrix}, &#10; \\begin{bmatrix}&#10;
\\theta\_{11} & 0 & 0 & 0 \\\\&#10; 0 & \\theta\_{22} & 0 & 0
\\\\&#10; 0 & 0 & \\theta\_{33} & 0 \\\\&#10; 0 & 0 & 0 & \\theta\_{44}
\\\\&#10; \\end{bmatrix}\\right)&#10;
\\end{aligned}&#10;](https://latex.codecogs.com/png.latex?%0A%20%20%5Cbegin%7Baligned%7D%0A%20%20%5Cbegin%7Bbmatrix%7D%0A%20%20%20%20%5Ceta_%7B1i%7D%20%5C%5C%0A%20%20%20%20%5Ceta_%7B2i%7D%0A%20%20%5Cend%7Bbmatrix%7D%20%26%20%5Csim%20%0A%20%20%5Cmathcal%7BN%7D%5Cleft%28%5Cbegin%7Bbmatrix%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Calpha_1%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Calpha_2%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Bbmatrix%7D%2C%20%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cbegin%7Bbmatrix%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cphi_%7B11%7D%20%26%20%5Cphi_%7B21%7D%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cphi_%7B21%7D%20%26%20%5Cphi_%7B22%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Bbmatrix%7D%5Cright%29%20%5C%5C%0A%20%20%5Cbegin%7Bbmatrix%7D%0A%20%20%20%20e_%7B0i%7D%20%5C%5C%0A%20%20%20%20e_%7B1i%7D%20%5C%5C%0A%20%20%20%20e_%7B2i%7D%20%5C%5C%0A%20%20%20%20e_%7B3i%7D%0A%20%20%5Cend%7Bbmatrix%7D%20%26%20%5Csim%0A%20%20%5Cmathcal%7BN%7D%5Cleft%28%5Cbegin%7Bbmatrix%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%200%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%200%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%200%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%200%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Bbmatrix%7D%2C%20%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cbegin%7Bbmatrix%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Ctheta_%7B11%7D%20%26%200%20%26%200%20%26%200%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%200%20%26%20%5Ctheta_%7B22%7D%20%26%200%20%26%200%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%200%20%26%200%20%26%20%5Ctheta_%7B33%7D%20%26%200%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%200%20%26%200%20%26%200%20%26%20%5Ctheta_%7B44%7D%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Bbmatrix%7D%5Cright%29%0A%20%20%5Cend%7Baligned%7D%0A
"
  \\begin{aligned}
  \\begin{bmatrix}
    \\eta_{1i} \\\\
    \\eta_{2i}
  \\end{bmatrix} & \\sim 
  \\mathcal{N}\\left(\\begin{bmatrix}
                \\alpha_1 \\\\
                \\alpha_2
              \\end{bmatrix}, 
              \\begin{bmatrix}
                \\phi_{11} & \\phi_{21} \\\\
                \\phi_{21} & \\phi_{22}
              \\end{bmatrix}\\right) \\\\
  \\begin{bmatrix}
    e_{0i} \\\\
    e_{1i} \\\\
    e_{2i} \\\\
    e_{3i}
  \\end{bmatrix} & \\sim
  \\mathcal{N}\\left(\\begin{bmatrix}
                0 \\\\
                0 \\\\
                0 \\\\
                0
              \\end{bmatrix}, 
              \\begin{bmatrix}
                \\theta_{11} & 0 & 0 & 0 \\\\
                0 & \\theta_{22} & 0 & 0 \\\\
                0 & 0 & \\theta_{33} & 0 \\\\
                0 & 0 & 0 & \\theta_{44} \\\\
              \\end{bmatrix}\\right)
  \\end{aligned}
")  

A path diagram is shown below:

``` r
growth_model <- "i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4
                 s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4
                 i ~~ 1 * i
                 s ~~ 0.2 * s + 0.1 * i
                 y1 ~~ 0.5 * y1
                 y2 ~~ 0.5 * y2
                 y3 ~~ 0.5 * y3
                 y4 ~~ 0.5 * y4
                 i ~ 1 * 1
                 s ~ 0.5 * 1"
library(semPlot)
semPaths(semPlotModel_lavaanModel(growth_model))
```

![](example_sem_files/figure-gfm/semPaths-growth-1.png)<!-- -->

## Workflow for Using `SimDesign`

### Step 1: Create data frame of design factors

In a methodological experiment with Monte Carlo simulation, one usually
generates millions of data sets across tens or hundreds of carefully
chosen conditions. As an example, here is a small scale simulation study
on LGM. The two goals are: (a) to understand the bias on the mean of
slopes and its standard error estimates, and (b) to illustrate the
difference between estimated and empirical standard error.

For simplicity, I’ll only choose three **designed factors** (i.e.,
manipulated independent variables), namely sample size, variance of the
random slope in the data generating model, and the mean of the slopes.
The design factors are summarized here:

  - Sample size (*N*): 50, 100, 200
  - Variance of slopes
    (![\\phi\_{22}](https://latex.codecogs.com/png.latex?%5Cphi_%7B22%7D
    "\\phi_{22}")): 0.1, 0.5 (i.e., 1/10 and 1/2 of the intercept
    variance)
  - Mean of slopes
    (![\\alpha\_2](https://latex.codecogs.com/png.latex?%5Calpha_2
    "\\alpha_2")): 1, 0.5

Therefore, it’s a 3
![\\times](https://latex.codecogs.com/png.latex?%5Ctimes "\\times") 2
![\\times](https://latex.codecogs.com/png.latex?%5Ctimes "\\times") 2
factorial design.

``` r
# Design factors:
designfactor <- createDesign(
  N = c(50, 100, 200), 
  phi22 = c(0.1, 0.5), 
  alpha2 = c(0, 0.5)
)
designfactor
```

    ># # A tibble: 12 x 3
    >#        N phi22 alpha2
    >#    <dbl> <dbl>  <dbl>
    >#  1    50   0.1    0  
    >#  2   100   0.1    0  
    >#  3   200   0.1    0  
    >#  4    50   0.5    0  
    >#  5   100   0.5    0  
    >#  6   200   0.5    0  
    >#  7    50   0.1    0.5
    >#  8   100   0.1    0.5
    >#  9   200   0.1    0.5
    ># 10    50   0.5    0.5
    ># 11   100   0.5    0.5
    ># 12   200   0.5    0.5

#### Fixed Values for the Study

In a simulation study, we cannot manipulate every possible variables. So
while we have designed on three design factors before, each with
multiple levels, there are values that we want to set to some constant
values, such as the mean and the variance of intercepts, etc. In
`SimDesign` these will be passed to different functions using a named
list called `fixed_objects`. Let’s create such a list:

``` r
fixed_obj <- list(phi11 = 1, 
                  Lambda = cbind(1, seq_len(4)), 
                  Theta = diag(0.5, nrow = 4))
```

### Step 2: Define a function for generating data

Create function for generating the data:

``` r
gen_lgm_data <- function(condition, fixed_objects = NULL) {
  N <- condition$N
  phi22 <- condition$phi22
  alpha2 <- condition$alpha2
  alpha <- c(1, alpha2)
  Phi <- matrix(c(fixed_objects$phi11, phi22 / 2,
                  phi22 / 2, phi22), nrow = 2)
  Lambda <- fixed_objects$Lambda
  Theta <- fixed_objects$Theta
  # Generate latent factor scores
  eta <- rmnorm(N, mean = alpha, varcov = Phi)
  # Generate residuals:
  e <- rmnorm(N, varcov = Theta)
  # Compute outcome scores
  y <- tcrossprod(eta, Lambda) + e
  colnames(y) <- paste0("y", 1:4)
  # Make it a data frame
  as.data.frame(y)
}
```

We can test that the function works for, say, condition 1 (i.e., row 1
of `designfactor`):

``` r
(test_data <- gen_lgm_data(designfactor[1, ], fixed_objects = fixed_obj))
```

    >#            y1          y2         y3          y4
    ># 1   0.8347251  0.38217213 -1.4155241 -0.41825775
    ># 2  -0.1424525 -2.38517074 -3.7443364 -3.40481824
    ># 3   0.3797961  1.20308378 -1.4583647 -0.97223312
    ># 4   0.5916468  1.29397939 -0.5166293  2.56479120
    ># 5   2.0020448  1.55798267  1.9176765  2.23577907
    ># 6   0.5332657  1.60585857  1.8757878  1.75464607
    ># 7   0.7454691  1.38998092  0.7432844  1.32394529
    ># 8   0.8061445  1.06064797  2.2033998  0.73335182
    ># 9  -0.4369852 -2.19125016 -0.8819522 -0.75603159
    ># 10  1.4510327  1.06650862  0.9763176 -0.16133342
    ># 11  1.2962571  1.50806920  0.3051535  0.03997832
    ># 12  3.2295427  3.33753630  2.4134247  3.38483689
    ># 13  2.1213320  1.80285328  3.4001704  3.35713400
    ># 14  2.0206412  1.07579671  1.6443448  1.59151284
    ># 15  2.0971434  1.77102291  2.1268183  1.15008252
    ># 16  1.2634529  0.82169633  1.6064522  1.42548834
    ># 17  0.9595218  2.27576483  2.8485264  1.98048267
    ># 18  2.3072508  2.96029323  1.3214246  1.40699840
    ># 19  0.9051683  0.84363518  0.3818444 -0.18729063
    ># 20  3.6980739  2.96498891  2.5519515  5.68663702
    ># 21  0.8015359  2.27374453  3.2517967  2.58267941
    ># 22  0.9734395  0.01437083  0.2802065 -1.11441149
    ># 23  0.6424167  1.86418297  0.1401472 -0.08316836
    ># 24  2.2127877  1.23128527  2.8367679  1.76033822
    ># 25 -1.4450568 -0.36299300 -0.3771272 -0.29563176
    ># 26  1.2651490  1.54138803  2.3497551  1.91271419
    ># 27  2.3547232  1.88436299  2.6966139  3.18166597
    ># 28  0.6100655  1.53344261  1.6709941  3.05642719
    ># 29  1.5458041 -1.20180279 -0.6420674 -1.00322995
    ># 30  2.8078939  1.65552217  3.2377611  4.10374016
    ># 31 -1.7920763  0.03233447  0.7198937  0.11673607
    ># 32  3.1931196  2.31931558  3.2943432  3.45514100
    ># 33  0.7427201  2.90433796  1.3814859  1.24387292
    ># 34  0.2463625  1.79780005  4.1162760  3.04780754
    ># 35 -0.4028881 -2.93738693 -1.3082866 -2.03025938
    ># 36  2.3116239  1.37546912  4.9092816  4.60624001
    ># 37  0.1093947 -1.65197220 -1.8409591 -1.27570642
    ># 38  0.4622808 -0.33853207 -0.7094562 -1.74262134
    ># 39  1.5742422  2.89724522  3.4178988  2.88181432
    ># 40  0.7211275  2.16263475  1.7268297  3.24536779
    ># 41 -1.1731955 -1.29762878 -0.1765694 -1.07734763
    ># 42  0.4202583  2.38155831  3.0002366  2.51968933
    ># 43  1.7608830  3.54144445  2.9000664  3.23403951
    ># 44  0.5401468  1.52338036  0.5101402  0.95034944
    ># 45 -0.9035784 -0.09535578 -0.1897216 -0.03142138
    ># 46  2.8860946  3.21722105  2.1093320  4.58258239
    ># 47  1.2463991 -0.31083396  1.2084012  0.50516282
    ># 48  2.2858351  2.95323388  1.9999386  2.49125960
    ># 49  1.2110286  0.85491597  3.3022137  1.29476807
    ># 50 -0.9917841 -2.13292363 -2.3112845 -1.32547799

If you want to check whether the simulated data is correct, generate
with a large sample size, and check the means and covariances:

``` r
large_test_df <- gen_lgm_data(tibble(N = 1e5, phi22 = 0.1, alpha2 = 0.5), 
                              fixed_objects = fixed_obj)
colMeans(large_test_df)
```

    >#       y1       y2       y3       y4 
    ># 1.497476 1.994736 2.492501 2.990415

``` r
cov(large_test_df)
```

    >#          y1       y2       y3       y4
    ># y1 1.707542 1.352779 1.502971 1.650655
    ># y2 1.352779 2.103302 1.846582 2.096394
    ># y3 1.502971 1.846582 2.694122 2.543183
    ># y4 1.650655 2.096394 2.543183 3.490395

#### Note: Other methods for generating SEM data

Many SEM software or packages have capability in generating data with
input of an SEM model. For example, in R, you can call Mplus using the
`MplusAutomation` package and use their `MONTECARLO` routine. In R, you
can generate SEM data using the `lavaan` package with the
`simulateData()` function, like the following example:

``` r
# Using a previously defined SEM model:
lavaan::simulateData(growth_model) %>% 
  head() # shows only the first six cases
```

    >#            y1       y2        y3        y4
    ># 1 -0.70726387 1.092568 1.6066896 1.6273346
    ># 2  2.79605629 3.805439 3.8631354 4.6831369
    ># 3  0.02203814 2.441585 3.9447603 4.7847719
    ># 4  0.78487534 2.895133 2.5946433 3.9875209
    ># 5  0.81832362 1.287995 0.9569316 0.5092388
    ># 6  3.77554598 3.462055 2.8186255 5.0123033

Personally, however, I prefer directly simulating data in R because

  - it forces you to specify everything in the model in the way you
    want. Especially in Mplus there are a lot of hidden default settings
    that may mess up with your simulation;
  - it makes the process of generating data more transparent;
  - it helps you learn the math behind the model;
  - it is more flexible as you can specify any distributional
    assumptions or models not supported by the SEM packages.

### Step 3: Define a function for analyzing the simulated data

In R, for running SEM models, the most common options are `lavaan`,
`OpenMx`, and Mplus (via `MplusAutomation`). When possible, I’ll stick
to `lavaan` to avoid jumping between programs, so let’s analyze the
simulated data twice, first with the true model and second with a
misspecified model where the random slope term is omitted (i.e., the
variance of `s` is constrained to zero).

``` r
analyze_lgm <- function(condition, dat, fixed_objects = NULL) {
  m1 <- 'i =~ 1 * y1 + 1 * y2 + 1 * y3 + 1 * y4
         s =~ 0 * y1 + 1 * y2 + 2 * y3 + 3 * y4
         i ~~ s'
  m2 <- 'i =~ 1 * y1 + 1 * y2 + 1 * y3 + 1 * y4
         s =~ 0 * y1 + 1 * y2 + 2 * y3 + 3 * y4
         s ~~ 0 * i + 0 * s'
  # Run model 1
  m1_fit <- growth(m1, data = dat)
  # Run model 2
  m2_fit <- growth(m2, data = dat)
  # Extract parameter estimates and standard errors
  ret <- c(coef(m1_fit)["s~1"], 
           sqrt(vcov(m1_fit)["s~1", "s~1"]),
           coef(m2_fit)["s~1"], 
           sqrt(vcov(m2_fit)["s~1", "s~1"]))
  names(ret) <- c("m1_est", "m1_se", "m2_est", "m2_se")
  ret
}
```

Test the analysis function:

``` r
analyze_lgm(designfactor[1, ], dat = test_data)
```

    >#     m1_est      m1_se     m2_est      m2_se 
    ># 0.07560676 0.06462994 0.08168169 0.05821217

Behind the scene, `SimDesign` saves the analysis results `ret` to a
matrix row by row per simulated data set. As an illustration, below
shows how the result matrix looks like.

``` r
set.seed(123)
rep <- 3
test_ret <- NULL
for (i in 1:rep) {
  test_data <- gen_lgm_data(designfactor[1, ], fixed_objects = fixed_obj)
  test_ret <- rbind(test_ret, 
                    analyze_lgm(designfactor[1, ], dat = test_data))
}
test_ret
```

    >#            m1_est      m1_se       m2_est      m2_se
    ># [1,] -0.004638175 0.05925965 -0.002177909 0.05417852
    ># [2,] -0.022546013 0.06004321 -0.010528237 0.05512071
    ># [3,] -0.039045161 0.06807877 -0.026733546 0.05703484

### Step 4: Define a function for evaluating the sample results

``` r
# Helper function for computing relative SE bias
rse_bias <- function(est_se, est) {
  est_se <- as.matrix(est_se)
  est <- as.matrix(est)
  est_se <- colMeans(est_se)
  emp_sd <- apply(est, 2L, sd)
  est_se / emp_sd - 1
}
evaluate_lgm <- function(condition, results, fixed_objects = NULL) {
  alpha2 <- condition$alpha2
  c(bias = bias(results[ , c("m1_est", "m2_est")], parameter = alpha2), 
    std_bias = bias(results[ , c("m1_est", "m2_est")], parameter = alpha2, 
                    type = "standardized"), 
    rmse = RMSE(results[ , c("m1_est", "m2_est")], parameter = alpha2), 
    rse_bias = rse_bias(results[ , c("m1_se", "m2_se")], 
                        results[ , c("m1_est", "m2_est")])
  )
}
```

Test the evaluation function:

``` r
evaluate_lgm(designfactor[1, ], test_ret)
```

    >#     bias.m1_est     bias.m2_est std_bias.m1_est std_bias.m2_est     rmse.m1_est 
    >#     -0.02207645     -0.01314656     -1.28289559     -1.05295041      0.02616843 
    >#     rmse.m2_est  rse_bias.m1_se  rse_bias.m2_se 
    >#      0.01663600      2.62967580      3.44074265

### Step 5: Test Run the simulation

Trial run with 2 replications

``` r
sim_trial <- runSimulation(designfactor, 
                           replications = 2, 
                           generate = gen_lgm_data, 
                           analyse = analyze_lgm, 
                           summarise = evaluate_lgm, 
                           fixed_objects = fixed_obj)
```

    ># 
    ># Design row: 1/12;   Started: Mon Jun 14 13:58:38 2021;   Total elapsed time: 0.00s 
    ># 
    ># Design row: 2/12;   Started: Mon Jun 14 13:58:38 2021;   Total elapsed time: 0.16s 
    ># 
    ># Design row: 3/12;   Started: Mon Jun 14 13:58:38 2021;   Total elapsed time: 0.28s 
    ># 
    ># Design row: 4/12;   Started: Mon Jun 14 13:58:39 2021;   Total elapsed time: 0.42s 
    ># 
    ># Design row: 5/12;   Started: Mon Jun 14 13:58:39 2021;   Total elapsed time: 0.54s 
    ># 
    ># Design row: 6/12;   Started: Mon Jun 14 13:58:39 2021;   Total elapsed time: 0.67s 
    ># 
    ># Design row: 7/12;   Started: Mon Jun 14 13:58:39 2021;   Total elapsed time: 0.80s 
    ># 
    ># Design row: 8/12;   Started: Mon Jun 14 13:58:39 2021;   Total elapsed time: 0.92s 
    ># 
    ># Design row: 9/12;   Started: Mon Jun 14 13:58:39 2021;   Total elapsed time: 1.05s 
    ># 
    ># Design row: 10/12;   Started: Mon Jun 14 13:58:39 2021;   Total elapsed time: 1.17s 
    ># 
    ># Design row: 11/12;   Started: Mon Jun 14 13:58:39 2021;   Total elapsed time: 1.30s 
    ># 
    ># Design row: 12/12;   Started: Mon Jun 14 13:58:40 2021;   Total elapsed time: 1.43s

    ># 
    ># Simulation complete. Total execution time: 1.57s

Hopefully it runs fine. If there’s any error, we have to go back and
check each component. The errors from `SimDesign` provides some hints.

### Step 6: Full simulation

Now we’re ready to run 500 replications. The `runSimulation()` function
has a number of handy arguments. Here I will use `parallel = TRUE` and
set `ncores` so that I use two cores.

An important thing to note is that when `parallel = TRUE`, one needs to
also export the packages by specifying all the packages needed for the
simulation via the `packages` argument. Here we will set `packages =
c("mnormt", "lavaan")`.

``` r
sim_result <- runSimulation(designfactor, 
                            replications = 500,
                            generate = gen_lgm_data,
                            analyse = analyze_lgm,
                            summarise = evaluate_lgm,
                            fixed_objects = fixed_obj,
                            parallel = TRUE,
                            ncores = min(parallel::detectCores() - 1, 2),
                            packages = c("mnormt", "lavaan"), 
                            save_results = TRUE)
```

## Summarize Simulation Results

A handy feature of `SimDesign` is that it saves the package and other
information used to simulation the data, as shown below

``` r
summary(sim_result)
```

    ># $sessionInfo
    ># R version 4.0.5 (2021-03-31)
    ># Platform: x86_64-pc-linux-gnu (64-bit)
    ># Running under: Ubuntu 20.04.2 LTS
    ># 
    ># Matrix products: default
    ># BLAS/LAPACK: /opt/OpenBLAS/lib/libopenblas-r0.3.13.so
    ># 
    ># locale:
    >#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    >#  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    >#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    >#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    >#  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ># [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ># 
    ># attached base packages:
    ># [1] stats     graphics  grDevices utils     datasets  methods   base     
    ># 
    ># other attached packages:
    >#  [1] semPlot_1.1.2   lavaan_0.6-8    mnormt_2.0.2    SimDesign_2.3  
    >#  [5] forcats_0.5.1   stringr_1.4.0   dplyr_1.0.5     purrr_0.3.4    
    >#  [9] readr_1.4.0     tidyr_1.1.3     tibble_3.1.0    tidyverse_1.3.1
    ># [13] ggplot2_3.3.3  
    ># 
    ># loaded via a namespace (and not attached):
    >#   [1] readxl_1.3.1        backports_1.2.1     Hmisc_4.5-0        
    >#   [4] systemfonts_1.0.1   plyr_1.8.6          igraph_1.2.6       
    >#   [7] splines_4.0.5       digest_0.6.27       foreach_1.5.1      
    >#  [10] htmltools_0.5.1.1   matrixcalc_1.0-3    rsconnect_0.8.17   
    >#  [13] fansi_0.4.2         magrittr_2.0.1      Rsolnp_1.16        
    >#  [16] checkmate_2.0.0     lisrelToR_0.1.4     cluster_2.1.1      
    >#  [19] openxlsx_4.2.3      modelr_0.1.8        svglite_2.0.0      
    >#  [22] jpeg_0.1-8.1        sem_3.1-11          colorspace_2.0-0   
    >#  [25] rvest_1.0.0         haven_2.4.0         xfun_0.22          
    >#  [28] crayon_1.4.1        jsonlite_1.7.2      lme4_1.1-26        
    >#  [31] regsem_1.6.2        survival_3.2-10     iterators_1.0.13   
    >#  [34] glue_1.4.2          kableExtra_1.3.4    gtable_0.3.0       
    >#  [37] webshot_0.5.2       mi_1.0              abind_1.4-5        
    >#  [40] scales_1.1.1        DBI_1.1.1           Rcpp_1.0.6         
    >#  [43] viridisLite_0.4.0   xtable_1.8-4        htmlTable_2.1.0    
    >#  [46] tmvnsim_1.0-2       foreign_0.8-81      Formula_1.2-4      
    >#  [49] stats4_4.0.5        truncnorm_1.0-8     htmlwidgets_1.5.3  
    >#  [52] httr_1.4.2          RColorBrewer_1.1-2  ellipsis_0.3.1     
    >#  [55] pkgconfig_2.0.3     XML_3.99-0.6        farver_2.1.0       
    >#  [58] nnet_7.3-15         sass_0.3.1          kutils_1.70        
    >#  [61] dbplyr_2.1.1        utf8_1.2.1          tidyselect_1.1.0   
    >#  [64] labeling_0.4.2      rlang_0.4.10        reshape2_1.4.4     
    >#  [67] munsell_0.5.0       cellranger_1.1.0    tools_4.0.5        
    >#  [70] cli_2.4.0           generics_0.1.0      broom_0.7.6        
    >#  [73] fdrtool_1.2.16      evaluate_0.14       arm_1.11-2         
    >#  [76] yaml_2.2.1          tables_0.9.6        knitr_1.32         
    >#  [79] fs_1.5.0            zip_2.1.1           glasso_1.11        
    >#  [82] pbapply_1.4-3       nlme_3.1-152        xml2_1.3.2         
    >#  [85] compiler_4.0.5      rstudioapi_0.13     png_0.1-7          
    >#  [88] reprex_2.0.0        statmod_1.4.35      modelsummary_0.6.6 
    >#  [91] bslib_0.2.4         pbivnorm_0.6.0      stringi_1.5.3      
    >#  [94] highr_0.8           qgraph_1.6.9        rockchalk_1.8.144  
    >#  [97] lattice_0.20-41     Matrix_1.3-2        psych_2.1.3        
    ># [100] nloptr_1.2.2.2      vctrs_0.3.7         pillar_1.6.0       
    ># [103] lifecycle_1.0.0     jquerylib_0.1.3     OpenMx_2.19.1      
    ># [106] data.table_1.14.0   corpcor_1.6.9       R6_2.5.0           
    ># [109] latticeExtra_0.6-29 bookdown_0.21       gridExtra_2.3      
    ># [112] codetools_0.2-18    gtools_3.8.2        boot_1.3-27        
    ># [115] MASS_7.3-53.1       assertthat_0.2.1    withr_2.4.1        
    ># [118] parallel_4.0.5      hms_1.0.0           grid_4.0.5         
    ># [121] rpart_4.1-15        coda_0.19-4         minqa_1.2.4        
    ># [124] rmarkdown_2.7       carData_3.0-4       lubridate_1.7.10   
    ># [127] base64enc_0.1-3    
    ># 
    ># $packages
    >#   packages versions
    ># 1   mnormt    2.0.2
    ># 2   lavaan    0.6.8
    ># 
    ># $save_info
    ># save_results_dirname 
    >#    "sim-lgm-results" 
    ># 
    ># $ncores
    ># [1] 2
    ># 
    ># $number_of_conditions
    ># [1] 12
    ># 
    ># $date_completed
    ># [1] "Sat Apr 17 21:00:24 2021"
    ># 
    ># $total_elapsed_time
    ># [1] "09m 0.26s"

### ANOVA

First, because the model (m1 vs m2) is a within-condition factor, let’s
restructure the data

``` r
sim_result_long <- sim_result %>%
  # Add condition ID
  rownames_to_column("con_id") %>%
  pivot_longer(
    bias.m1_est:rse_bias.m2_se,
    names_to = c(".value", "model"),
    names_pattern =
      "(bias|std_bias|rmse|rse_bias)\\.(m1|m2)_.*",
    names_ptypes =
      list(model = factor(levels = c("m1", "m2")))
  )
sim_result_long
```

    ># # A tibble: 24 x 14
    >#    con_id     N phi22 alpha2 REPLICATIONS SIM_TIME COMPLETED       SEED WARNINGS
    >#    <chr>  <dbl> <dbl>  <dbl>        <int>    <dbl> <chr>          <int>    <int>
    >#  1 1         50   0.1      0          500     42.4 Sat Apr 17 2… 1.10e9       86
    >#  2 1         50   0.1      0          500     42.4 Sat Apr 17 2… 1.10e9       86
    >#  3 2        100   0.1      0          500     46.9 Sat Apr 17 2… 1.57e9       19
    >#  4 2        100   0.1      0          500     46.9 Sat Apr 17 2… 1.57e9       19
    >#  5 3        200   0.1      0          500     45.1 Sat Apr 17 2… 1.70e9        2
    >#  6 3        200   0.1      0          500     45.1 Sat Apr 17 2… 1.70e9        2
    >#  7 4         50   0.5      0          500     48.3 Sat Apr 17 2… 1.25e9      106
    >#  8 4         50   0.5      0          500     48.3 Sat Apr 17 2… 1.25e9      106
    >#  9 5        100   0.5      0          500     41.9 Sat Apr 17 2… 5.05e8       29
    ># 10 5        100   0.5      0          500     41.9 Sat Apr 17 2… 5.05e8       29
    ># # … with 14 more rows, and 5 more variables: model <fct>, bias <dbl>,
    ># #   std_bias <dbl>, rmse <dbl>, rse_bias <dbl>

``` r
# Function for computing eta-squared
aov_etasq <- function(id, dv, data, between, within) {
  trans_data <- data
  trans_data[c(between, within)] <- 
    lapply(trans_data[c(between, within)], as.factor)
  form <- paste0(dv, " ~ ", paste(c(between, within), collapse = " * "), 
                 " + Error(", id, "/", paste(within, collapse = " * "), ")")
  form <- as.formula(form)
  aov_mixed <- aov(form, data = trans_data)
  sum_aov <- summary(aov_mixed)
  tab <- do.call(rbind, unname(unlist(sum_aov, recursive = FALSE)))
  tab$`Eta Sq` <- tab$`Sum Sq` / sum(tab$`Sum Sq`)
  tab
}
aov_etasq("con_id", dv = "bias", data = sim_result_long, 
          between = c("N", "phi22", "alpha2"), within = "model") %>%
  knitr::kable(digits = 5L)
```

|                      | Df | Sum Sq | Mean Sq |  Eta Sq |
| :------------------- | -: | -----: | ------: | ------: |
| N                    |  2 |  1e-05 |   1e-05 | 0.07284 |
| phi22                |  1 |  2e-05 |   2e-05 | 0.12221 |
| alpha2               |  1 |  0e+00 |   0e+00 | 0.00141 |
| N:phi22              |  2 |  2e-05 |   1e-05 | 0.13182 |
| N:alpha2             |  2 |  7e-05 |   3e-05 | 0.36616 |
| phi22:alpha2         |  1 |  0e+00 |   0e+00 | 0.01255 |
| N:phi22:alpha2       |  2 |  5e-05 |   3e-05 | 0.27960 |
| model                |  1 |  0e+00 |   0e+00 | 0.00050 |
| N:model              |  2 |  0e+00 |   0e+00 | 0.00255 |
| phi22:model          |  1 |  0e+00 |   0e+00 | 0.00009 |
| alpha2:model         |  1 |  0e+00 |   0e+00 | 0.00014 |
| N:phi22:model        |  2 |  0e+00 |   0e+00 | 0.00068 |
| N:alpha2:model       |  2 |  0e+00 |   0e+00 | 0.00442 |
| phi22:alpha2:model   |  1 |  0e+00 |   0e+00 | 0.00068 |
| N:phi22:alpha2:model |  2 |  0e+00 |   0e+00 | 0.00434 |

With relatively small number of conditions, one can present the results
in a table (and it’s handy in R):

``` r
sim_result_long %>% 
  arrange(alpha2, phi22, N, phi22) %>%
  select(`$\\alpha_{2}$` = alpha2,
         `$\\phi_{22}$` = phi22,
         N,
         model,
         `Standardized Bias` = std_bias, 
         `Relative SE Error` = rse_bias) %>% 
  knitr::kable(digits = 3L)
```

| ![\\alpha\_{2}](https://latex.codecogs.com/png.latex?%5Calpha_%7B2%7D "\\alpha_{2}") | ![\\phi\_{22}](https://latex.codecogs.com/png.latex?%5Cphi_%7B22%7D "\\phi_{22}") |   N | model | Standardized Bias | Relative SE Error |
| -----------------------------------------------------------------------------------: | --------------------------------------------------------------------------------: | --: | :---- | ----------------: | ----------------: |
|                                                                                  0.0 |                                                                               0.1 |  50 | m1    |             0.019 |           \-0.033 |
|                                                                                  0.0 |                                                                               0.1 |  50 | m2    |             0.023 |           \-0.154 |
|                                                                                  0.0 |                                                                               0.1 | 100 | m1    |             0.002 |           \-0.011 |
|                                                                                  0.0 |                                                                               0.1 | 100 | m2    |             0.006 |           \-0.128 |
|                                                                                  0.0 |                                                                               0.1 | 200 | m1    |             0.041 |             0.002 |
|                                                                                  0.0 |                                                                               0.1 | 200 | m2    |             0.035 |           \-0.125 |
|                                                                                  0.0 |                                                                               0.5 |  50 | m1    |             0.005 |             0.028 |
|                                                                                  0.0 |                                                                               0.5 |  50 | m2    |             0.023 |           \-0.215 |
|                                                                                  0.0 |                                                                               0.5 | 100 | m1    |           \-0.021 |           \-0.052 |
|                                                                                  0.0 |                                                                               0.5 | 100 | m2    |           \-0.027 |           \-0.289 |
|                                                                                  0.0 |                                                                               0.5 | 200 | m1    |           \-0.010 |             0.000 |
|                                                                                  0.0 |                                                                               0.5 | 200 | m2    |           \-0.016 |           \-0.258 |
|                                                                                  0.5 |                                                                               0.1 |  50 | m1    |             0.022 |           \-0.039 |
|                                                                                  0.5 |                                                                               0.1 |  50 | m2    |             0.028 |           \-0.151 |
|                                                                                  0.5 |                                                                               0.1 | 100 | m1    |             0.083 |             0.019 |
|                                                                                  0.5 |                                                                               0.1 | 100 | m2    |             0.088 |           \-0.110 |
|                                                                                  0.5 |                                                                               0.1 | 200 | m1    |           \-0.044 |           \-0.037 |
|                                                                                  0.5 |                                                                               0.1 | 200 | m2    |           \-0.039 |           \-0.166 |
|                                                                                  0.5 |                                                                               0.5 |  50 | m1    |           \-0.062 |           \-0.065 |
|                                                                                  0.5 |                                                                               0.5 |  50 | m2    |           \-0.066 |           \-0.307 |
|                                                                                  0.5 |                                                                               0.5 | 100 | m1    |             0.026 |           \-0.038 |
|                                                                                  0.5 |                                                                               0.5 | 100 | m2    |             0.033 |           \-0.270 |
|                                                                                  0.5 |                                                                               0.5 | 200 | m1    |             0.028 |           \-0.070 |
|                                                                                  0.5 |                                                                               0.5 | 200 | m2    |             0.021 |           \-0.297 |

It is, however, recommended you try to plot the results, both for
exploratory purpose and for better presentation of the results.

## Exercise

1.  From the simulation results, evaluate the relative efficiency
    (`rse_bias`) of the estimated average slope (i.e.,
    ![\\alpha\_2](https://latex.codecogs.com/png.latex?%5Calpha_2
    "\\alpha_2")) under model 2 relative to that under model 1.
