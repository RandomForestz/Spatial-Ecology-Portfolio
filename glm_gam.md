GLMs or GAMs: Does Allowing Non-linear Relationships Improve Model
Predictive Power?
================

- [The Purpose](#the-purpose)
- [The Subject](#the-subject)
- [The Prep](#the-prep)
  - [Libraries](#libraries)
  - [Occurrence](#occurrence)
  - [Predictors](#predictors)
  - [Pseudo-absences](#pseudo-absences)
  - [Extracting predictors](#extracting-predictors)
- [The GLM](#the-glm)
  - [GLM Comparison](#glm-comparison)
- [The GAM](#the-gam)
  - [Interpreting the GAM](#interpreting-the-gam)
- [GLM vs. GAM](#glm-vs-gam)
- [The Conclusion](#the-conclusion)

## The Purpose

Generalized linear models (GLMs) are a common framework for evaluating
relationships between a response variable and a set of predictors. GLMs
extend linear regression by allowing the response variable to follow a
non-normal distribution (e.g., binomial or Poisson). However, they still
assume a linear relationship between predictors and the response on the
link scale, which may not capture more complex ecological patterns.

Generalized additive models (GAMs) offer greater flexibility by allowing
non-linear relationships between predictors and the response through the
use of smooth functions. This makes GAMs especially useful when
ecological responses are expected to vary along environmental gradients
in non-linear ways.

Here, I evaluate whether the linearity assumption improves model
predictive performance. Specifically, I compare the fit and predictive
accuracy of a GLM and a GAM applied to species occurrence data, using
soil characteristics as environmental predictors. I also see how models
perform individually and how predictors drive the response along the
way.

------------------------------------------------------------------------

## The Subject

Common juniper is a small tree/shrub of the cypress family. It is has
the largest geographical range of any woody plant, practically covering
the cool temperate ecosystems of the Northern Hemisphere.

With such a range, I’m excited to see which soil characteristics are
influencing common juniper distribution. I’m also just a big fan of this
plant.

------------------------------------------------------------------------

## The Prep

### Libraries

I need some spatial and statistical libraries to handle the data prep
and analysis. Here are my preferred:

``` r
library(sf) # vector data
library(terra) # raster data
library(tidyverse) # data manipulation and plotting
library(tmap) # I like tmap maps
library(mgcv) # for GAMs
library(ggeffects) # Response curves
library(gratia) # GAM response curves
library(car) # VIF
```

### Occurrence

I first want to map the occurrence data, but my data is in .csv format
from an iNaturalist download. I’ll take care of that and map!

``` r
wrld <- st_read("D:/Projects/Spatial-Ecology-Portfolio/data/glm/World_Continents.shp")
```

    ## Reading layer `World_Continents' from data source 
    ##   `D:\Projects\Spatial-Ecology-Portfolio\data\glm\World_Continents.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 8 features and 6 fields
    ## Geometry type: MULTIPOLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: -20037510 ymin: -30240970 xmax: 20037510 ymax: 18418390
    ## Projected CRS: WGS 84 / Pseudo-Mercator

``` r
wrld <- wrld %>% 
  filter(CONTINENT != "Antarctica")

occ <- read.csv("D:/Projects/Spatial-Ecology-Portfolio/data/glm/presence.csv") %>% 
  na.omit() # get rid of NA's - sf cant handle missing values for vector creation

occ_sf <- st_as_sf(occ, coords = c("longitude", "latitude"), crs = 4326) # to match world
```

A quick but not so pretty map..

``` r
tm_shape(wrld) +
  tm_polygons(col = "CONTINENT") +
  tm_shape(occ_sf) +
    tm_dots() +
  tm_layout(frame.lwd = 3, legend.position = c("left", "bottom"))
```

![](glm_gam_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

a much prettier map and a bit more insightful, showing observations per
large grid cell.

``` r
grid <- st_make_grid(wrld,
                     cellsize = c(1000000,1000000),
                     square = T) %>% 
  st_intersection(., wrld) %>% 
  st_transform(., crs = st_crs(occ_sf)) %>% 
  st_cast(., "MULTIPOLYGON") %>% 
  st_sf() %>% 
  st_make_valid()

sf_use_s2(FALSE) # so the grid doesn't give me any issues!
```

    ## Spherical geometry (s2) switched off

``` r
grid_occ <- grid %>% 
  mutate(counts = lengths(st_intersects(., occ_sf)))
```

    ## although coordinates are longitude/latitude, st_intersects assumes that they
    ## are planar

``` r
tm_shape(grid_occ) +
  tm_polygons(col = "counts", palette = "viridis", style = "fixed",
              breaks = c(0, 1, 10, 50, 100, 200, 500, 1000, 2000, 3000, Inf)) +
  tm_layout(legend.position = c("left", "bottom"))
```

![](glm_gam_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### Predictors

``` r
rast.list <- list.files(path = "D:/Projects/Spatial-Ecology-Portfolio/data/glm",
                        pattern = "\\.tif$",
                        full.names = T)
preds <- lapply(rast.list, rast)
soils <- rast(preds)
```

### Pseudo-absences

The response variable in this study is binary, indicating species
presence (1) or absence (0). While presence data are commonly available
through field surveys or databases, true absence data are often lacking
or unreliable due to detectability issues and incomplete sampling. To
address this, pseudo-absences were generated — randomly selected
background points assumed to represent absences.

Pseudo-absences help overcome the imbalance in presence-only datasets
and reduce the risk of Type II errors (false negatives). In accordance
with standard practice, the number of pseudo-absences was set at
approximately 2–4 times the number of presences to balance statistical
power and model sensitivity.

``` r
set.seed(111)
abs <- st_sample(st_geometry(wrld), type = "random", size = (nrow(occ_sf)*2)) %>% 
  st_sf() %>% 
  mutate(type = 0) %>% 
  st_transform(., st_crs(occ_sf))

one <- occ_sf %>% 
  select(geometry) %>% 
  mutate(type = 1)
  
data <- rbind(one, abs)
```

### Extracting predictors

Now we need to extract predictor variable values to each presence or
absence.

``` r
data <- st_transform(data, crs = crs(soils)) # change crs for data to match soils
train <- terra::extract(soils, data) # extract soil information to points
train$type <- data$type # add in 0,1 info
train <- train[c(8, 2:7)] # reorder for preferece


train <- st_drop_geometry(train)
train <- na.omit(train) %>% 
  as.data.frame()

pander::pander(head(train)) # view
```

| type | clay | nitrogen | ph  | sand | silt | soc  |
|:----:|:----:|:--------:|:---:|:----:|:----:|:----:|
|  1   | 310  |   635    | 52  | 311  | 379  | 1072 |
|  1   | 120  |   395    | 62  | 517  | 363  | 791  |
|  1   | 258  |   908    | 54  | 295  | 447  | 1653 |
|  1   | 224  |   507    | 66  | 353  | 423  | 896  |
|  1   | 233  |   560    | 58  | 370  | 397  | 1080 |
|  1   | 236  |   795    | 56  | 328  | 436  | 920  |

------------------------------------------------------------------------

## The GLM

``` r
mod_glm <- glm(type ~ clay + nitrogen + ph + sand + silt + soc,
               data = train, 
               family = binomial)

summary(mod_glm)
```

    ## 
    ## Call:
    ## glm(formula = type ~ clay + nitrogen + ph + sand + silt + soc, 
    ##     family = binomial, data = train)
    ## 
    ## Coefficients:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -8.152e+00  1.268e+01  -0.643    0.520    
    ## clay         2.154e-03  1.268e-02   0.170    0.865    
    ## nitrogen     2.340e-03  3.172e-05  73.762   <2e-16 ***
    ## ph          -6.048e-02  8.786e-04 -68.837   <2e-16 ***
    ## sand         1.051e-02  1.268e-02   0.829    0.407    
    ## silt         1.725e-02  1.268e-02   1.361    0.174    
    ## soc         -1.692e-03  1.749e-05 -96.727   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 180236  on 142198  degrees of freedom
    ## Residual deviance: 149102  on 142192  degrees of freedom
    ## AIC: 149116
    ## 
    ## Number of Fisher Scoring iterations: 4

The model summary provides estimates for the intercept and predictor
coefficients, along with their standard errors, z-values, and p-values.
Among the soil variables, nitrogen (positive effect), pH (negative
effect), and SOC (negative effect) are all statistically significant
predictors of presence. In contrast, clay, sand, and silt are
non-significant, and may be contributing unnecessary noise to the model
— which we can check for potential collinearity.

It is important to emphasize that the coefficient estimates reflect the
change in log-odds of presence per unit increase in each predictor, not
a direct change in probability. This is a fundamental feature of
logistic regression models with binary response variables.

The model shows good explanatory power: the null deviance of 180,236 is
reduced to a residual deviance of 149,102, indicating that the
predictors substantially improve model fit. Additionally, the Akaike
Information Criterion (AIC) value of 149,116 provides a useful metric
for model comparison. A lower AIC indicates a better-fitting model, and
this value will be used to compare against more flexible alternatives,
such as a generalized additive model (GAM).

------------------------------------------------------------------------

#### Log-odds to odds ratios

Log-odds can be a difficult thing to interpret. I’m applying the code
below to examine it differently, making our interpretation a bit more
intuitive.

``` r
pander::pander(exp(coef(mod_glm)))
```

| (Intercept) | clay  | nitrogen |   ph   | sand  | silt  |  soc   |
|:-----------:|:-----:|:--------:|:------:|:-----:|:-----:|:------:|
|  0.0002881  | 1.002 |  1.002   | 0.9413 | 1.011 | 1.017 | 0.9983 |

| Predictor   | Odds Ratio | Interpretation                                                        |
|-------------|------------|-----------------------------------------------------------------------|
| (Intercept) | 0.00029    | Baseline odds of presence when all predictors = 0 (not interpretable) |
| Clay        | 1.0022     | Each 1-unit increase in clay increases odds of presence by 0.2%       |
| Nitrogen    | 1.0023     | Each 1-unit increase in nitrogen increases odds by 0.23%              |
| pH          | 0.9413     | Each 1-unit increase in pH decreases odds by about 5.9%               |
| Sand        | 1.0106     | Each 1-unit increase in sand increases odds by 1.1%                   |
| Silt        | 1.0174     | Each 1-unit increase in silt increases odds by 1.7%                   |
| SOC         | 0.9983     | Each 1-unit increase in SOC decreases odds by about 0.17%             |

------------------------------------------------------------------------

#### Checking Multicollinearity

Checking collinearity is extremely important in modeling. Using the
vif() function, we get value to describe multicollinearity in our model.
For reference, 5-10 presents problems, 10 and up is pretty bad.

``` r
pander::pander(vif(mod_glm))
```

| clay  | nitrogen |  ph   | sand  | silt  | soc  |
|:-----:|:--------:|:-----:|:-----:|:-----:|:----:|
| 17639 |  3.005   | 1.436 | 55880 | 30333 | 3.35 |

We see that our non-significant variables from the glm are extremely
high. This might be due to each of these variables being reported in a
0-100% content. I’ll probably drop 2 of them in a model 2.0.

#### Reponse Curves

This is a fun part of modeling. How does each variable influence the
predicted probability of common juniper presence?

``` r
par(mfrow = c(1,3))
plot(ggpredict(mod_glm, terms = "ph"))
```

    ## Data were 'prettified'. Consider using `terms="ph [all]"` to get smooth
    ##   plots.

![](glm_gam_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
plot(ggpredict(mod_glm, terms = "nitrogen"))
```

    ## Data were 'prettified'. Consider using `terms="nitrogen [all]"` to get
    ##   smooth plots.

![](glm_gam_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
plot(ggpredict(mod_glm, terms = "soc"))
```

    ## Data were 'prettified'. Consider using `terms="soc [all]"` to get smooth
    ##   plots.

![](glm_gam_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->

``` r
par(mfrow=c(1,1))
plot(ggpredict(mod_glm, terms = c("ph", "nitrogen", "soc")))
```

    ## Data were 'prettified'. Consider using `terms="ph [all]"` to get smooth
    ##   plots.

![](glm_gam_files/figure-gfm/unnamed-chunk-11-4.png)<!-- -->

#### Model 2.0

``` r
mod_glm2 <- glm(type ~ ph + nitrogen + soc + clay, data = train, family = binomial)

summary(mod_glm2)
```

    ## 
    ## Call:
    ## glm(formula = type ~ ph + nitrogen + soc + clay, family = binomial, 
    ##     data = train)
    ## 
    ## Coefficients:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  3.912e+00  5.978e-02   65.44   <2e-16 ***
    ## ph          -5.557e-02  8.189e-04  -67.86   <2e-16 ***
    ## nitrogen     2.516e-03  3.084e-05   81.58   <2e-16 ***
    ## soc         -1.586e-03  1.730e-05  -91.64   <2e-16 ***
    ## clay        -6.733e-03  9.718e-05  -69.29   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 180236  on 142198  degrees of freedom
    ## Residual deviance: 158224  on 142194  degrees of freedom
    ## AIC: 158234
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
pander::pander(exp(coef(mod_glm2)))
```

| (Intercept) |   ph   | nitrogen |  soc   |  clay  |
|:-----------:|:------:|:--------:|:------:|:------:|
|    49.98    | 0.9459 |  1.003   | 0.9984 | 0.9933 |

| Predictor   | Odds Ratio | Interpretation                                                        |
|-------------|------------|-----------------------------------------------------------------------|
| (Intercept) | 50.089     | Baseline odds of presence when all predictors = 0                     |
| pH          | 0.9459     | Each 1-unit increase in pH decreases odds of presence by ~5.4%        |
| Nitrogen    | 1.0025     | Each 1-unit increase in nitrogen increases odds of presence by ~0.25% |
| SOC         | 0.9984     | Each 1-unit increase in SOC decreases odds of presence by ~0.16%      |
| Clay        | 0.9933     | Each 1-unit increase in clay decreases odds of presence by ~0.67%     |

``` r
pander::pander(vif(mod_glm2))
```

|  ph   | nitrogen |  soc  | clay  |
|:-----:|:--------:|:-----:|:-----:|
| 1.405 |   3.29   | 3.619 | 1.143 |

``` r
par(mfrow = c(1,3))
plot(ggpredict(mod_glm2, terms = "ph"))
```

    ## Data were 'prettified'. Consider using `terms="ph [all]"` to get smooth
    ##   plots.

![](glm_gam_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
plot(ggpredict(mod_glm2, terms = "nitrogen"))
```

    ## Data were 'prettified'. Consider using `terms="nitrogen [all]"` to get
    ##   smooth plots.

![](glm_gam_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

``` r
plot(ggpredict(mod_glm2, terms = "soc"))
```

    ## Data were 'prettified'. Consider using `terms="soc [all]"` to get smooth
    ##   plots.

![](glm_gam_files/figure-gfm/unnamed-chunk-15-3.png)<!-- -->

``` r
plot(ggpredict(mod_glm2, terms = "clay"))
```

    ## Data were 'prettified'. Consider using `terms="clay [all]"` to get
    ##   smooth plots.

![](glm_gam_files/figure-gfm/unnamed-chunk-15-4.png)<!-- -->

``` r
par(mfrow=c(1,1))
plot(ggpredict(mod_glm2, terms = c("ph", "nitrogen", "soc", "clay")))
```

    ## Data were 'prettified'. Consider using `terms="ph [all]"` to get smooth
    ##   plots.

    ## Package {see} needed to plot multiple panels in one integrated figure.
    ##   Please install it by typing `install.packages("see", dependencies =
    ##   TRUE)` into the console.

    ## [[1]]

![](glm_gam_files/figure-gfm/unnamed-chunk-15-5.png)<!-- -->

    ## 
    ## [[2]]

![](glm_gam_files/figure-gfm/unnamed-chunk-15-6.png)<!-- -->

    ## 
    ## [[3]]

![](glm_gam_files/figure-gfm/unnamed-chunk-15-7.png)<!-- -->

### GLM Comparison

``` r
AIC(mod_glm, mod_glm2)
```

    ##          df    AIC
    ## mod_glm   7 149116
    ## mod_glm2  5 158234

``` r
library(pROC)
```

    ## Warning: package 'pROC' was built under R version 4.4.1

    ## Type 'citation("pROC")' for a citation.

    ## 
    ## Attaching package: 'pROC'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     cov, smooth, var

``` r
prob_glm <- predict(mod_glm, type = "response")
prob_glm2 <- predict(mod_glm2, type = "response")

roc_glm  <- roc(train$type, prob_glm)
```

    ## Setting levels: control = 0, case = 1

    ## Setting direction: controls < cases

``` r
roc_glm2 <- roc(train$type, prob_glm2)
```

    ## Setting levels: control = 0, case = 1
    ## Setting direction: controls < cases

``` r
plot(roc_glm, col = "blue", lwd = 2, main = "ROC Curves: GLM vs GLM2")
plot(roc_glm2, col = "red", lwd = 2, add = TRUE)
legend("bottomright", legend = c("mod_glm", "mod_glm2"),
       col = c("blue", "red"), lwd = 2)
```

![](glm_gam_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
roc.test(roc_glm, roc_glm2)
```

    ## 
    ##  DeLong's test for two correlated ROC curves
    ## 
    ## data:  roc_glm and roc_glm2
    ## Z = 52.547, p-value < 2.2e-16
    ## alternative hypothesis: true difference in AUC is not equal to 0
    ## 95 percent confidence interval:
    ##  0.04034700 0.04347342
    ## sample estimates:
    ## AUC of roc1 AUC of roc2 
    ##   0.7809060   0.7389958

Modeling is most exciting when it requires weighing statistical rigor
against predictive performance. In this case, two models offer competing
advantages:

1.  mod_glm provides higher predictive performance, as evidenced by a
    lower AIC, higher AUC, and improved classification accuracy. If the
    goal were solely predictive mapping, this model would clearly
    dominate.

2.  mod_glm2, however, addresses important structural issues. It removes
    highly collinear predictors, includes only statistically significant
    variables, and results in a more interpretable and stable model.
    While it shows slightly reduced predictive power, it offers greater
    confidence in coefficient estimates and the ecological
    interpretation of soil variables.

Though mod_glm performs better numerically, we favor mod_glm2 for its
stronger inferential validity. A responsible modeling approach should
prioritize the integrity of the model structure and the interpretability
of its predictors, especially when informing ecological understanding or
conservation decisions. Predictive accuracy is valuable, but not at the
expense of robustness and transparency.

Let’s see if a GAMs allowance of non-linear relationships improves our
model.

------------------------------------------------------------------------

## The GAM

Now I’ll fit a GAM using the same predictors, examine response smooths,
compare AIC and AUC, and decide whether flexibility improves
performance.

``` r
mod_gam <- gam(type ~ s(ph) + s(nitrogen) + s(soc) + s(clay),
               data = train,
               family = binomial)

summary(mod_gam)
```

    ## 
    ## Family: binomial 
    ## Link function: logit 
    ## 
    ## Formula:
    ## type ~ s(ph) + s(nitrogen) + s(soc) + s(clay)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -1.91343    0.05484  -34.89   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##               edf Ref.df Chi.sq p-value    
    ## s(ph)       8.121  8.267   1707  <2e-16 ***
    ## s(nitrogen) 8.914  8.996   2434  <2e-16 ***
    ## s(soc)      8.932  8.998  11097  <2e-16 ***
    ## s(clay)     7.809  8.210   5818  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.288   Deviance explained = 26.6%
    ## UBRE = -0.069221  Scale est. = 1         n = 142199

The GAM output differs slightly from the GLM summary in structure and
interpretation. While the intercept is not of primary interest, it is
significant and serves as the baseline from which smoothed relationships
are built.

Each smooth term is reported with its estimated degrees of freedom (edf)
and an associated p-value. The edf reflects the flexibility of the
spline: values near 1 suggest a linear relationship, while values
significantly above 1 indicate non-linear patterns. In this case, all
predictors exhibit strongly non-linear responses:

- ph (edf = 8.12),

- nitrogen (edf = 8.91),

- soc (edf = 8.93), and

- clay (edf = 7.81)

These high edf values reveal that the relationship between each soil
variable and the probability of presence is complex and curved, not
simply linear. Moreover, all smooth terms are highly statistically
significant (p \< 2e-16), reinforcing the value of allowing non-linear
structures in this model.

Additionally, the R-squared value of .288 explains ~29% in the binary
response, which isn’t outstanding by any means. The deviance explained
is 26.6%, which seems reasonable for this models response of presence
and absence. Lastly, the UBRE value is -0.069221, which is negative and
small, suggesting a good fit (this is the gam verison of glm AIC).

### Interpreting the GAM

One of the strengths of using GAMs is the ability to visualize partial
effect response curves. These plots display the estimated effect of each
predictor on the response variable (in this case, presence probability)
while holding all other variables constant.

- The x-axis represents the range of observed values for a given
  predictor.

- The y-axis represents the partial effect on the log-odds of presence —
  effectively how that variable alone shifts the predicted probability.

- The black line is the estimated smooth function — a flexible,
  data-driven fit showing how the response changes.

- The grey shaded area represents the 95% confidence interval, providing
  a sense of uncertainty around the estimate. Narrow bands indicate high
  certainty; wider bands suggest more variability or sparse data.

These plots are particularly useful for identifying non-linear
thresholds, plateaus, or unimodal relationships that GLMs would miss due
to their assumption of linearity.

``` r
draw(mod_gam)
```

![](glm_gam_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

#### pH

The relationship between pH and presence probability is relatively flat
across most of the gradient, with a noticeable dip beyond a value of
~80. This suggests that soil pH has little influence on common juniper
presence up to a certain point, after which probability declines. This
pattern may reflect physiological intolerance to higher pH levels,
consistent with known ecological thresholds in many plant species.

#### Nitrogen

There is a strong positive effect of nitrogen up to approximately 800
units, after which the curve levels off, followed by minor dips and
rises as uncertainty increases. This indicates that common juniper
presence increases sharply with nitrogen availability to a threshold,
beyond which additional nitrogen no longer enhances probability of
occurrence — a classic saturation effect.

#### Soil Organic Carbon (SOC)

SOC displays a highly non-linear and fluctuating relationship with
presence probability. An initial positive peak around 500, followed by a
dip and subtle oscillations, suggests that juniper favors low to
moderate organic soils. High SOC values are associated with declines in
presence, indicating avoidance of highly organic, possibly poorly
drained soils.

#### Clay Content

The relationship shows a steady decline in presence probability as clay
content increases. This aligns well with field observations: common
juniper prefers well-drained, sandy or rocky soils. High clay content
likely reflects heavier, wetter soils that may be unsuitable for
establishment or survival.

#### Response Curves using ggeffect

I like ggeeffects and while the partial plots above are informative,
these stand-along plots work great for visualizing each variable. I can
also visualize them together as was done with the glm model above.

``` r
library(ggeffects)
par(mfrow = c(2, 2))
plot(ggpredict(mod_gam, terms = "ph"))
```

![](glm_gam_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
plot(ggpredict(mod_gam, terms = "nitrogen"))
```

![](glm_gam_files/figure-gfm/unnamed-chunk-20-2.png)<!-- -->

``` r
plot(ggpredict(mod_gam, terms = "soc"))
```

![](glm_gam_files/figure-gfm/unnamed-chunk-20-3.png)<!-- -->

``` r
plot(ggpredict(mod_gam, terms = "clay"))
```

![](glm_gam_files/figure-gfm/unnamed-chunk-20-4.png)<!-- -->

``` r
par(mfrow = c(1, 1))
plot(ggpredict(mod_gam, terms = c("ph", "nitrogen", "soc", "clay")))
```

    ## Package {see} needed to plot multiple panels in one integrated figure.
    ##   Please install it by typing `install.packages("see", dependencies =
    ##   TRUE)` into the console.

    ## [[1]]

![](glm_gam_files/figure-gfm/unnamed-chunk-20-5.png)<!-- -->

    ## 
    ## [[2]]

![](glm_gam_files/figure-gfm/unnamed-chunk-20-6.png)<!-- -->

    ## 
    ## [[3]]

![](glm_gam_files/figure-gfm/unnamed-chunk-20-7.png)<!-- -->

------------------------------------------------------------------------

## GLM vs. GAM

``` r
prob_gam <- predict(mod_gam, type = "response")

roc_gam <- roc(train$type, prob_gam)
```

    ## Setting levels: control = 0, case = 1

    ## Setting direction: controls < cases

``` r
plot(roc_glm2, col = "blue", lwd = 2, main = "ROC Curves: GLM vs GAM")
plot(roc_gam, col = "darkgreen", lwd = 2, add = TRUE)
legend("bottomright", legend = c("mod_glm", "mod_gam"),
       col = c("blue", "darkgreen"), lwd = 2)
```

![](glm_gam_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
roc.test(roc_glm2, roc_gam)
```

    ## 
    ##  DeLong's test for two correlated ROC curves
    ## 
    ## data:  roc_glm2 and roc_gam
    ## Z = -80.315, p-value < 2.2e-16
    ## alternative hypothesis: true difference in AUC is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.09133699 -0.08698529
    ## sample estimates:
    ## AUC of roc1 AUC of roc2 
    ##   0.7389958   0.8281569

Just as was done with comparing GLMs, DeLong’s test for correlated ROC
curves shows that that the test was significant that there is a
difference between the model performance, but that GAM was much higher
(73.9% vs. 82.8%). The improvement is real.

------------------------------------------------------------------------

## The Conclusion

Modeling is rarely as straightforward as comparing numbers. While some
models may perform statistically better (as with the GLM vs. GLM2),
selecting the “best” model depends just as much on how well it aligns
with the ecological questions, underlying assumptions, and the way we
intend to interpret and communicate the results.

In this case, the GAM outperformed the GLMs in predictive accuracy. But
more importantly, it better reflects how we understand ecological
systems. The relationships between ecological responses such as species
presence, abundance, or canopy cover and environmental covariates like
climate, soil, or topography are rarely linear. They bend, peak, and
dip. The GAM allows us to embrace that complexity, not ignore it.

This isn’t to say GLMs aren’t useful — they are. They are interpretable,
robust, and efficient. But when ecological realism matters, and when the
response surface may be non-linear, GAMs offer a way to model the world
more like it really is: dynamic, flexible, and full of gradients.
