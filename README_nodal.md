# Risk Factor Modeling of Tumor Nodal Involvement

## Project Overview

This project uses **binomial logistic regression** to investigate factors associated with lymph-node involvement among **53 patients with prostate cancer**.

The objective is to evaluate the association of age, clinical stage, biopsy grade, X-ray findings, and serum acid phosphatase with nodal involvement, and to identify a parsimonious model using **backward selection and likelihood-ratio tests**.

## Dataset

The analysis uses the `nodal` dataset from the R package `boot`; no external data file is required. The outcome (`r`) and predictors are recorded as follows:

| Variable | Description |
|---|---|
| `r` | Nodal involvement: 1 = present, 0 = absent |
| `aged` | Age category: 1 = at least 60 years, 0 = younger than 60 |
| `stage` | More serious tumor size/position assessment (1 versus 0) |
| `grade` | More serious biopsy-based tumor grade (1 versus 0) |
| `xray` | More serious X-ray assessment (1 versus 0) |
| `acid` | Serum acid phosphatase indicator (binary in the dataset) |

The dataset also contains `m`, a constant column of ones, which is not used as a predictor.

Dataset documentation: [R documentation for `boot::nodal`](https://stat.ethz.ch/R-manual/R-devel/library/boot/help/nodal.html).

## Statistical Methods

1. **Exploratory data analysis:** Frequency distributions, missing-value checks, and contingency tables.
2. **Bivariate associations:** Fisher's exact tests for the outcome and candidate predictors, and exploratory associations among predictors.
3. **Logistic regression:** A full model including all five candidate predictors.
4. **Backward model selection:** Sequential elimination of predictors using likelihood-ratio tests to identify a reduced model.
5. **Model interpretation:** Adjusted odds ratios, 95% confidence intervals, and a forest plot.
6. **Apparent discrimination:** Confusion matrix at a 0.5 cutoff, sensitivity, specificity, ROC curve, and AUC, evaluated on the model development data.
7. **Model checks:** Likelihood-ratio comparison with the intercept-only model and an exploratory Hosmer–Lemeshow test.

## Main Results

Of the 53 patients, **20 (37.7%)** had observed nodal involvement. Backward selection retained **`stage`, `xray`, and `acid`** in the final logistic regression model.

The estimated model is:

```text
logit[P(r = 1)] = -3.0518 + 1.6453*stage + 1.9116*xray + 1.6378*acid
```

The adjusted odds ratios for the retained predictors were:

| Predictor | Adjusted OR | 95% CI |
|---|---:|---:|
| `stage` | 5.18 | 1.24–21.66 |
| `xray` | 6.76 | 1.47–31.02 |
| `acid` | 5.14 | 1.17–22.54 |

The selected model showed a statistically significant improvement over the intercept-only model (likelihood-ratio test, **p = 0.00010**).

### In-sample classification performance

| Measure | Result |
|---|---:|
| Correct classifications | 43 / 53 |
| Accuracy | 81.1% |
| Sensitivity | 80.0% |
| Specificity | 81.8% |
| ROC AUC | 0.834 |

These measures were calculated using the **same 53 patients used for model development and variable selection**; they may therefore overestimate performance on independent data.

The fitted model estimated a probability of nodal involvement of approximately **4.5%** when all retained indicators were 0, compared with **89.5%** when all three were 1. These are model-based probabilities for the specified predictor combinations, not causal effects.

## Interpretation and Limitations

The retained predictors were associated with the **odds** of nodal involvement in this sample. Odds ratios are not risk ratios, and the analysis does not establish causal relationships.

With only 53 patients, coefficient estimates and variable selection may be unstable; the confidence intervals are wide. In-sample accuracy and AUC may be optimistic, particularly after backward selection. The Hosmer–Lemeshow test could not form all ten requested groups because of tied predicted probabilities; its non-significant result should not be treated as proof of adequate calibration. Independent validation and a more comprehensive calibration assessment would be needed before considering clinical predictive use.

## Requirements

The script uses R and these packages:

```r
install.packages(c("boot", "ggplot2", "pROC", "ResourceSelection"))
```

## How to Run

1. Download or clone the repository.
2. Open `01_nodal_involvement.R` in RStudio.
3. Install the packages listed above if necessary.
4. Run the script from top to bottom.

The script loads the dataset directly from `boot`, prints the statistical analyses, and displays plots in RStudio. **No files are exported automatically.**

## Tools and Technologies

R | Logistic Regression | Likelihood-Ratio Testing | Backward Selection | Odds Ratios | ROC Analysis
