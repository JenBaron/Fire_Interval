---
title: "Fire Intervals"
author: "Jen Baron, j.baron@alumni.ubc.ca, UBC Tree Ring Lab"
date: 'May 25, 2021'
output:
  html_document:
    keep_md: yes
    number_sections: yes
    toc: yes
    toc_float: yes
---


# Load Packages


```r
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(ggplot2)
library(strucchange)
```

```
## Loading required package: zoo
```

```
## 
## Attaching package: 'zoo'
```

```
## The following objects are masked from 'package:base':
## 
##     as.Date, as.Date.numeric
```

```
## Loading required package: sandwich
```

```r
source("opt_breakpoints.R")
```


# Load Data


```r
fire <- read.csv("data/fire_curve.csv")
fire_years <- read.csv("data/fire_years.csv")
str(fire_years)
```

```
## 'data.frame':	25 obs. of  2 variables:
##  $ fire_year: int  1600 1619 1630 1662 1677 1694 1704 1717 1730 1747 ...
##  $ fire     : int  0 1 2 3 4 5 6 7 8 9 ...
```

```r
str(fire)
```

```
## 'data.frame':	416 obs. of  2 variables:
##  $ year: int  1600 1601 1602 1603 1604 1605 1606 1607 1608 1609 ...
##  $ fire: int  0 0 0 0 0 0 0 0 0 0 ...
```

```r
plot(fire)
```

![](Fire_Interval_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

Need to convert to time series before analysis


```r
fire_ts <- ts(fire$fire, start=c(1600), end=c(2015), frequency=1)

plot(fire_ts)
```

![](Fire_Interval_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
hist(fire$fire, breaks=40)
```

![](Fire_Interval_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```r
ggplot(fire, aes(fire)) +
  geom_density(fill = "lightgray") +
  theme_classic()
```

![](Fire_Interval_files/figure-html/unnamed-chunk-4-2.png)<!-- -->

```r
 ggplot(fire, aes(x=year, y=fire)) +
  geom_point(size=0.7, alpha=0.3) +
   geom_smooth(method="loess", col="black") +
  ylab("Number of Fires") +
  theme_classic() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))
```

```
## `geom_smooth()` using formula 'y ~ x'
```

![](Fire_Interval_files/figure-html/unnamed-chunk-4-3.png)<!-- -->

# Changepoint Detection


```r
set.seed(1210)
```



**Generalized fluctuation test & F test**


```r
ocus_f <- efp(fire_ts ~ 1, type = "OLS-CUSUM")
plot(ocus_f)
```

![](Fire_Interval_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
sctest(ocus_f)
```

```
## 
## 	OLS-based CUSUM test
## 
## data:  ocus_f
## S0 = 9.1219, p-value < 2.2e-16
```

```r
f.p <- sctest(ocus_f)$p.value

fs_f <- Fstats(fire_ts ~ 1)
plot(fs_f)
```

![](Fire_Interval_files/figure-html/unnamed-chunk-6-2.png)<!-- -->

**Regression breakpoints**


```r
bp_f <- breakpoints(fire ~ year, data = fire)
plot(bp_f)
```

![](Fire_Interval_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

```r
bpts_sum <- summary(bp_f)
opt_brks <- opt_breakpoints(bpts_sum$RSS["BIC",])
opt_brks <- 2
```


```r
bp.1 <- breakpoints(bp_f, breaks = opt_brks)
breaks.1 <- bp.1$breakpoints

## confidence intervals
ci.1 <- confint(bp_f, breaks= opt_brks)
```

```
## Warning in confint.breakpointsfull(bp_f, breaks = opt_brks): Confidence interval
## 2 cannot be computed: P(argmax V <= 0) = 1271310319616
```

```r
ci.1$confint
```

```
##   2.5 % breakpoints 97.5 %
## 1   157         158    161
## 2    NA         343     NA
```

```r
fire[bp.1$breakpoints,] #get the Year from row number
```

```
##     year fire
## 158 1757    8
## 343 1942   21
```


```r
plot(fire ~ year, data = fire, type = "l")
for (i in 1: opt_brks) {
  abline(v = fire$year[ci.1$confint[i,2]], col = "red")
  abline(v = fire$year[ci.1$confint[i,1]], col = "black", lty = 3)
  abline(v = fire$year[ci.1$confint[i,3]], col = "black", lty = 3)
}
```

![](Fire_Interval_files/figure-html/unnamed-chunk-10-1.png)<!-- -->



```r
est.1 <- fire[ci.1$confint[,2],] #get the Year from row number
low.1 <- fire[ci.1$confint[,1],] #lower CI
up.1 <- fire[ci.1$confint[,3],] #upper CI

est.1 <- est.1 %>% 
  select(year) %>%
  rename(year, "estimate" = "year")
low.1 <- low.1 %>% 
  select(year) %>%
  rename(year, "lower" = "year")
up.1 <- up.1 %>%
  select(year) %>%
  rename(year, "upper" = "year")
est.f <- cbind(est.1, low.1, up.1)
```


```r
bf.1 <- breakfactor(bp.1)
#bf.1<- relevel(bf.1, ref=4)

fm0.f <- lm(data=fire, fire ~ year) #null hypothesis
fm1.f <- lm(data=fire, fire~year*bf.1) #alternative hypothesis
coef(fm1.f)
```

```
##       (Intercept)              year      bf.1segment2      bf.1segment3 
##      -84.02550357        0.05223934      -23.67174445      106.02550357 
## year:bf.1segment2 year:bf.1segment3 
##        0.01441659       -0.05223934
```

```r
summary(fm1.f)
```

```
## 
## Call:
## lm(formula = fire ~ year * bf.1, data = fire)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -1.0171 -0.2162  0.0000  0.2488  1.0503 
## 
## Coefficients:
##                     Estimate Std. Error t value Pr(>|t|)    
## (Intercept)       -8.403e+01  1.084e+00  -77.54   <2e-16 ***
## year               5.224e-02  6.454e-04   80.95   <2e-16 ***
## bf.1segment2      -2.367e+01  1.436e+00  -16.48   <2e-16 ***
## bf.1segment3       1.060e+02  4.209e+00   25.19   <2e-16 ***
## year:bf.1segment2  1.442e-02  8.221e-04   17.54   <2e-16 ***
## year:bf.1segment3 -5.224e-02  2.154e-03  -24.25   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.37 on 410 degrees of freedom
## Multiple R-squared:  0.9977,	Adjusted R-squared:  0.9976 
## F-statistic: 3.498e+04 on 5 and 410 DF,  p-value: < 2.2e-16
```


```r
fm1.f$coefficients
```

```
##       (Intercept)              year      bf.1segment2      bf.1segment3 
##      -84.02550357        0.05223934      -23.67174445      106.02550357 
## year:bf.1segment2 year:bf.1segment3 
##        0.01441659       -0.05223934
```

```r
#intercepts
b1 <- fm1.f[["coefficients"]][["(Intercept)"]]
b2 <- fm1.f[["coefficients"]][["(Intercept)"]]+fm1.f[["coefficients"]][["bf.1segment2"]]
b3 <- fm1.f[["coefficients"]][["(Intercept)"]]+fm1.f[["coefficients"]][["bf.1segment3"]]
# b4 <- fm1.f[["coefficients"]][["(Intercept)"]]+fm1.f[["coefficients"]][["bf.1segment4"]]
# b5 <- fm1.f[["coefficients"]][["(Intercept)"]]+fm1.f[["coefficients"]][["bf.1segment5"]]
#slopes
m1 <- fm1.f[["coefficients"]][["year"]]
m2 <- fm1.f[["coefficients"]][["year"]]+fm1.f[["coefficients"]][["year:bf.1segment2"]]
m3 <- fm1.f[["coefficients"]][["year"]]+fm1.f[["coefficients"]][["year:bf.1segment3"]]
# m4 <- fm1.f[["coefficients"]][["year"]]+fm1.f[["coefficients"]][["year:bf.1segment4"]]
# m5 <- fm1.f[["coefficients"]][["year"]]+fm1.f[["coefficients"]][["year:bf.1segment5"]]

f.coef <- data.frame(metric = "nf", segment = c(1:3), slope = c(m1,m2,m3), intercept = c(b1,b2,b3)) 
f.coef
```

```
##   metric segment        slope intercept
## 1     nf       1 5.223934e-02  -84.0255
## 2     nf       2 6.665593e-02 -107.6972
## 3     nf       3 2.359224e-16   22.0000
```


```r
fit.a <- fitted(fm1.f)
mod.a <- data.frame(year = fire$year, fire = fit.a)

# plot the fitted model
fig.1 <- ggplot() +
  #geom_line(data=fire, aes(x = year, y = fire), colour = "darkgrey", size=0.75) +
  geom_point(data=fire, aes(x=year, y=fire), size=0.6, alpha=0.3) +
  geom_vline(xintercept = est.f$estimate, size=0.6, col = "black", linetype="dashed") +
  #geom_vline(data=fire_years, aes(xintercept = fire_year)) +
  #geom_vline(xintercept = (est.f$lower), linetype = "dashed", size=0.3, alpha=0.7) +
  #geom_vline(xintercept = (est.f$upper), linetype = "dashed", size=0.3, alpha=0.7) +
  geom_line(data=mod.a, aes(x = year, y = fire), alpha=1, col="black", size=0.6) + 
  #geom_abline(slope=m1, intercept=b1, col = "blue") +
  #geom_abline(slope=m2, intercept=b2, col = "blue") +
  #geom_abline(slope=m3, intercept=b3, col = "blue") +
  scale_x_continuous(limits=c(1600,2015), name = "") +
  ylab("Cumulative Number of Fires") +  
  theme_classic() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.title.x = element_blank())

fig.1
```

![](Fire_Interval_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

```r
ggsave("figures/fig_1.png", fig.1, dpi=300, width=6, height=4)
```

```r
est.f
```

```
##     estimate lower upper
## 158     1757  1756  1760
## 343     1942    NA    NA
```

# Assumptions


```r
plot(fm1.f, which=2) #normality
```

![](Fire_Interval_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

```r
plot(fm1.f, which=3) #equal variance
```

![](Fire_Interval_files/figure-html/unnamed-chunk-16-2.png)<!-- -->

```r
plot(fm1.f, which=1) #linearity
```

![](Fire_Interval_files/figure-html/unnamed-chunk-16-3.png)<!-- -->

```r
plot(fm1.f, which=4) #Cook's distance
```

![](Fire_Interval_files/figure-html/unnamed-chunk-16-4.png)<!-- -->
These look good


```r
plot(residuals(fm1.f), type="b")
abline(h=0,lty=3)
```

![](Fire_Interval_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

```r
# ACF function on residuals
acf(resid(fm1.f), main="acf(resid(fm1.f))", plot=T)
```

![](Fire_Interval_files/figure-html/unnamed-chunk-17-2.png)<!-- -->

```r
acf(resid(fm1.f), main="acf(resid(fm1.f))", plot=F)
```

```
## 
## Autocorrelations of series 'resid(fm1.f)', by lag
## 
##      0      1      2      3      4      5      6      7      8      9     10 
##  1.000  0.821  0.663  0.526  0.411  0.317  0.244  0.192  0.144  0.118  0.107 
##     11     12     13     14     15     16     17     18     19     20     21 
##  0.100  0.109  0.139  0.106  0.057  0.012 -0.019 -0.047 -0.054 -0.067 -0.077 
##     22     23     24     25     26 
## -0.068 -0.038 -0.042 -0.026 -0.007
```

```r
acf((fire_ts), main="fire_ts", plot=T)
```

![](Fire_Interval_files/figure-html/unnamed-chunk-17-3.png)<!-- -->

```r
acf((fire_ts), main="fire_ts", plot=F)
```

```
## 
## Autocorrelations of series '(fire_ts)', by lag
## 
##     0     1     2     3     4     5     6     7     8     9    10    11    12 
## 1.000 0.994 0.989 0.983 0.978 0.972 0.967 0.961 0.955 0.950 0.944 0.938 0.933 
##    13    14    15    16    17    18    19    20    21    22    23    24    25 
## 0.927 0.921 0.915 0.909 0.902 0.896 0.890 0.884 0.878 0.872 0.866 0.861 0.855 
##    26 
## 0.849
```

Temporal autocorrelation is high in the residuals - increases mean absolute error and probability of a type I error (false positive).

# Save Outputs


```r
write.csv(mod.a, "outputs/model_fitted.csv")
write.csv(est.f, "outputs/estimates.csv")
write.csv(f.coef, "outputs/coef.csv")
```

# Reproducibility


```r
Sys.time()
```

```
## [1] "2021-05-25 14:10:38 PDT"
```

```r
git2r::repository()
```

```
## Local:    main C:/Users/jenbaron/Documents/UBC/Research/Tree Ring Lab/Manuscripts/Fire Interval/Fire_Interval
## Remote:   main @ origin (https://github.com/JenBaron/Fire_Interval.git)
## Head:     [47ffe07] 2021-05-25: Initial commit
```

```r
sessionInfo()
```

```
## R version 4.0.5 (2021-03-31)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 19042)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_United States.1252 
## [2] LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] strucchange_1.5-2 sandwich_3.0-0    zoo_1.8-9         ggplot2_3.3.3    
## [5] dplyr_1.0.6      
## 
## loaded via a namespace (and not attached):
##  [1] git2r_0.28.0      highr_0.9         pillar_1.6.0      bslib_0.2.4      
##  [5] compiler_4.0.5    jquerylib_0.1.4   tools_4.0.5       digest_0.6.27    
##  [9] nlme_3.1-152      lattice_0.20-41   jsonlite_1.7.2    evaluate_0.14    
## [13] lifecycle_1.0.0   tibble_3.1.1      gtable_0.3.0      mgcv_1.8-34      
## [17] pkgconfig_2.0.3   rlang_0.4.11      Matrix_1.3-3      DBI_1.1.1        
## [21] yaml_2.2.1        xfun_0.22         withr_2.4.2       stringr_1.4.0    
## [25] knitr_1.33        generics_0.1.0    vctrs_0.3.8       sass_0.3.1       
## [29] grid_4.0.5        tidyselect_1.1.1  glue_1.4.2        R6_2.5.0         
## [33] fansi_0.4.2       rmarkdown_2.7     farver_2.1.0      purrr_0.3.4      
## [37] magrittr_2.0.1    splines_4.0.5     scales_1.1.1      ellipsis_0.3.2   
## [41] htmltools_0.5.1.1 assertthat_0.2.1  colorspace_2.0-1  labeling_0.4.2   
## [45] utf8_1.2.1        stringi_1.5.3     munsell_0.5.0     crayon_1.4.1
```
