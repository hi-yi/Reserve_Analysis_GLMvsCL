GLM\_Reserve
================

``` r
library(ChainLadder)
```

    ## Warning: package 'ChainLadder' was built under R version 3.6.3

    ## 
    ## Welcome to ChainLadder version 0.2.11
    ## 
    ## Type vignette('ChainLadder', package='ChainLadder') to access
    ## the overall package documentation.
    ## 
    ## See demo(package='ChainLadder') for a list of demos.
    ## 
    ## More information is available on the ChainLadder project web-site:
    ## https://github.com/mages/ChainLadder
    ## 
    ## To suppress this message use:
    ## suppressPackageStartupMessages(library(ChainLadder))

``` r
# Make triangle with incremental payments
incremental_triangle = as.matrix(RAA)
for (i in 2:10){
  incremental_triangle[,i] = RAA[,i] - RAA[,i-1]
}  
# Prepare the data as dataframe (X:AT, Dev, Y:Loss)
df_orig = as.data.frame(incremental_triangle)
colnames(df_orig) = c("AY","Dev","Loss")
df_orig = na.omit(df_orig)
df_orig$origin = factor(df_orig$AY)
df_orig$Dev = factor(df_orig$Dev)
# impute negative value as 1 (avoid error when log)
df_orig[which(df_orig$Loss<0),3] = 1
print(incremental_triangle)
```

    ##       dev
    ## origin    1    2    3    4    5    6    7   8   9  10
    ##   1981 5012 3257 2638  898 1734 2642 1828 599  54 172
    ##   1982  106 4179 1111 5270 3116 1817 -103 673 535  NA
    ##   1983 3410 5582 4881 2268 2594 3479  649 603  NA  NA
    ##   1984 5655 5900 4211 5500 2159 2658  984  NA  NA  NA
    ##   1985 1092 8473 6271 6333 3786  225   NA  NA  NA  NA
    ##   1986 1513 4932 5257 1233 2917   NA   NA  NA  NA  NA
    ##   1987  557 3463 6926 1368   NA   NA   NA  NA  NA  NA
    ##   1988 1351 5596 6165   NA   NA   NA   NA  NA  NA  NA
    ##   1989 3133 2262   NA   NA   NA   NA   NA  NA  NA  NA
    ##   1990 2063   NA   NA   NA   NA   NA   NA  NA  NA  NA

### Prepared the triangle data into tabular format

#### X variables : origin(AY or UY), dev(development year)

#### Y variable : incremental loss

``` r
# Run GLM model with poisson distribution assumption & log link
model = glm(Loss ~ AY+Dev,data=df_orig, family=poisson(link='log'))
summary(model)
```

    ## 
    ## Call:
    ## glm(formula = Loss ~ AY + Dev, family = poisson(link = "log"), 
    ##     data = df_orig)
    ## 
    ## Deviance Residuals: 
    ##    Min      1Q  Median      3Q     Max  
    ## -61.40  -22.60    0.00   14.47   53.58  
    ## 
    ## Coefficients:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  7.653960   0.009860 776.296  < 2e-16 ***
    ## AY1982      -0.104634   0.010634  -9.840  < 2e-16 ***
    ## AY1983       0.245808   0.009832  25.000  < 2e-16 ***
    ## AY1984       0.421234   0.009569  44.018  < 2e-16 ***
    ## AY1985       0.430239   0.009661  44.531  < 2e-16 ***
    ## AY1986       0.035944   0.010922   3.291 0.000999 ***
    ## AY1987      -0.058181   0.011788  -4.935 8.00e-07 ***
    ## AY1988       0.244326   0.011685  20.909  < 2e-16 ***
    ## AY1989      -0.159131   0.015875 -10.024  < 2e-16 ***
    ## AY1990      -0.022044   0.024123  -0.914 0.360832    
    ## Dev2         0.692826   0.008290  83.574  < 2e-16 ***
    ## Dev3         0.626028   0.008595  72.835  < 2e-16 ***
    ## Dev4         0.276947   0.009618  28.796  < 2e-16 ***
    ## Dev5         0.060559   0.010550   5.740 9.47e-09 ***
    ## Dev6        -0.195820   0.011995 -16.325  < 2e-16 ***
    ## Dev7        -1.052591   0.018531 -56.801  < 2e-16 ***
    ## Dev8        -1.274260   0.024364 -52.300  < 2e-16 ***
    ## Dev9        -1.917732   0.042067 -45.588  < 2e-16 ***
    ## Dev10       -2.506466   0.076884 -32.601  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 92002  on 54  degrees of freedom
    ## Residual deviance: 37568  on 36  degrees of freedom
    ## AIC: 38117
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
coef = model$coefficients
```

### Most of the variables have been proven to be significant.

### Using these coefficients, we can fill the lower-down triangle to predict future loss.

``` r
fitted = matrix(coef[1],10,10)
for (i in 2:10){
  fitted[i,] = fitted[i,] + coef[i]
}
for (j in 2:10){
  fitted[,j] = fitted[,j] + coef[9+j]
}
fitted = exp(fitted)

# Incurred Loss
inc = c()
for (i in 1:10){
  inc = rbind(inc, RAA[i,11-i])
}

# Ultimate Loss
ult_glm = rowSums(fitted)

# INBR
ibnr_glm = ult_glm-inc

glm_result=data.frame("AY"=c(1981:1990),"Incurred Loss"=inc,"Ultimate, Loss"=round(ult_glm,0),"IBNR"=round(ibnr_glm,0))
glm_result
```

    ##      AY Incurred.Loss Ultimate..Loss  IBNR
    ## 1  1981         18834          18834     0
    ## 2  1982         16704          16963   259
    ## 3  1983         23466          24082   616
    ## 4  1984         27067          28700  1633
    ## 5  1985         26180          28960  2780
    ## 6  1986         15852          19523  3671
    ## 7  1987         12314          17769  5455
    ## 8  1988         13112          24047 10935
    ## 9  1989          5395          16063 10668
    ## 10 1990          2063          18423 16360

### We can calculate IBNR using Fitted Triangle from GLM Model

# IBNR Calculation with Incurred Loss Development Method (weighted average LDF assumption)

``` r
ldf = c()
for (i in 1:9){
  ldf = cbind(ldf, sum(RAA[1:(10-i),i+1])/sum(RAA[1:(10-i),i]))
}
ult_cl = c(inc[1])
for (i in 2:10){
  ult_cl = rbind(ult_cl, inc[i]*prod(ldf[(11-i):9]))
}
ibnr_cl = ult_cl - inc
CL=data.frame("AY"=c(1981:1990),"Incurred Loss"=inc,"Ultimate Loss"=round(ult_cl,0),"IBNR"=round(ibnr_cl,0))
rownames(CL)=NULL
print(CL)
```

    ##      AY Incurred.Loss Ultimate.Loss  IBNR
    ## 1  1981         18834         18834     0
    ## 2  1982         16704         16858   154
    ## 3  1983         23466         24083   617
    ## 4  1984         27067         28703  1636
    ## 5  1985         26180         28927  2747
    ## 6  1986         15852         19501  3649
    ## 7  1987         12314         17749  5435
    ## 8  1988         13112         24019 10907
    ## 9  1989          5395         16045 10650
    ## 10 1990          2063         18402 16339

``` r
print(paste("IBNR from GLM:",round(sum(ibnr_glm),0)))
```

    ## [1] "IBNR from GLM: 52378"

``` r
print(paste("IBNR from ILDM:",round(sum(ibnr_cl),0)))
```

    ## [1] "IBNR from ILDM: 52135"

``` r
print(paste("IBNR difference:",round(sum(ibnr_glm-ibnr_cl),0)))
```

    ## [1] "IBNR difference: 242"
