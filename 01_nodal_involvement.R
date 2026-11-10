# RISK FACTOR MODELING OF TUMOR NODAL INVOLVEMENT
# Logistic regression in 53 prostate cancer patients

library(boot)
library(ggplot2)
library(pROC)

# DATASET
# r: nodal involvement (1 = present, 0 = absent)
# aged: age >= 60 years (1), < 60 years (0)
# stage: tumour size/position assessment (1 = more serious)
# grade: biopsy-based tumour grade (1 = more serious)
# xray: X-ray severity indicator (1 = more serious)
# acid: serum acid phosphatase indicator (binary in this dataset)
# m: constant column of ones, excluded from the analysis
# Objective: identify factors associated with nodal involvement.

nodal = boot::nodal

head(nodal)
dim(nodal)
str(nodal)
summary(nodal)
colSums(is.na(nodal))

# GRAPHIC SETTINGS
colors = c("Absent" = "#2864A0", "Present" = "#E18B45")
theme_set(theme_classic(base_size = 12))

# UNIVARIATE EXPLORATORY ANALYSIS
# Outcome and binary clinical predictors.
table(nodal$r)
prop.table(table(nodal$r))

predictors = c("aged", "stage", "grade", "xray", "acid")
summary(nodal[, predictors])

for (variable in predictors) {print(table(nodal[[variable]]))
  print(prop.table(table(nodal[[variable]])))}

# Bar plot of the proportion with each positive risk indicator.
proportions = colMeans(nodal[, predictors] == 1)
plot_prev = ggplot(data.frame(predictor = names(proportions),
                              proportion = as.numeric(proportions)),
                   aes(x = predictor, y = proportion)) +
  geom_col(fill = "#2864A0", width = 0.65) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(title = "Prevalence of clinical risk indicators",
       x = NULL, y = "Proportion with indicator = 1") +
  theme_classic(base_size = 12)
print(plot_prev)

# BIVARIATE ANALYSIS
for (variable in predictors) {cat("\nRisk factor:", variable, "\n")
  tab = table(nodal$r, nodal[[variable]])
  print(tab)
  print(fisher.test(tab))}

# Associations between the candidate predictors.
round(cor(nodal[, predictors]), 3)
for (i in 1:(length(predictors) - 1)) 
  {for (j in (i + 1):length(predictors)) 
    {cat("\n", predictors[i], "vs", predictors[j], "\n")
    print(fisher.test(table(nodal[[predictors[i]]],
                            nodal[[predictors[j]]])))}}

# LOGISTIC REGRESSION
# Binary logistic model: logit[P(r = 1)] = beta0 + beta1*x1 + ...
# Start with all five predictors.
fitT = glm(r ~ aged + stage + grade + xray + acid,
           family = binomial, data = nodal)
summary(fitT)

# Backward variable selection using likelihood-ratio tests.
# At each step, remove the predictor with the weakest evidence
# of association, based on the likelihood-ratio test.
drop1(fitT, test = "Chisq")
fit1 = update(fitT, . ~ . - aged)
summary(fit1)
drop1(fit1, test = "Chisq")
fit2 = update(fit1, . ~ . - grade)
summary(fit2)
drop1(fit2, test = "Chisq")

# Check successive nested-model comparisons.
anova(fit2, fit1, fitT, test = "Chisq")

# Compare the selected model with the intercept-only model.
fit0 = glm(r ~ 1, family = binomial, data = nodal)
anova(fit0, fit2, test = "Chisq")
AIC(fit0, fitT, fit1, fit2)

# FINAL MODEL
# Logistic regression with the selected predictors.
# Adjusted odds ratios and 95% confidence intervals.
summary(fit2)
coef(fit2)
exp(coef(fit2))             # Adjusted odds ratios
exp(confint.default(fit2))  # 95% Wald confidence intervals for ORs

# Forest plot of adjusted odds ratios (excluding the intercept).
coef_table = summary(fit2)$coefficients
or_data = data.frame(predictor = rownames(coef_table)[-1],
  OR = exp(coef_table[-1, 1]),
  lower = exp(coef_table[-1, 1] - 1.96 * coef_table[-1, 2]),
  upper = exp(coef_table[-1, 1] + 1.96 * coef_table[-1, 2]))

plot_or = ggplot(or_data, aes(x = predictor, y = OR)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey55") +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  color = "#2864A0", linewidth = 0.8) +
  coord_flip() +
  scale_y_log10() +
  labs(title = "Adjusted odds ratios for nodal involvement",
       x = NULL, y = "Odds ratio (log scale), 95% CI") +
  theme_classic(base_size = 12)
print(plot_or)

# IN-SAMPLE CLASSIFICATION
# The same 53 patients are used to fit and evaluate the model.
# These are apparent (training-set) results, not external validation.
predicted_prob = fitted(fit2)
predicted_class = ifelse(predicted_prob >= 0.5, 1, 0)

confusion = table(Observed = factor(nodal$r, levels = c(0, 1)),
                  Predicted = factor(predicted_class, levels = c(0, 1)))
confusion

TP = sum(nodal$r == 1 & predicted_class == 1)
TN = sum(nodal$r == 0 & predicted_class == 0)
FP = sum(nodal$r == 0 & predicted_class == 1)
FN = sum(nodal$r == 1 & predicted_class == 0)

accuracy = (TP + TN) / nrow(nodal)
sensitivity = TP / (TP + FN)
specificity = TN / (TN + FP)
c(accuracy = accuracy, sensitivity = sensitivity, specificity = specificity)

# ROC curve and area under the curve (evaluated in sample).
roc_fit = roc(nodal$r, predicted_prob,
              levels = c(0, 1), direction = "<", quiet = TRUE)
auc(roc_fit)

roc_data = data.frame(FPR = 1 - roc_fit$specificities,
                      TPR = roc_fit$sensitivities)
plot_roc = ggplot(roc_data, aes(x = FPR, y = TPR)) +
  geom_line(color = "#2864A0", linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey55") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = paste0("ROC curve (in sample): AUC = ",
                      round(as.numeric(auc(roc_fit)), 3)),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_classic(base_size = 12)
print(plot_roc)

# MODEL GOODNESS OF FIT
# Hosmer-Lemeshow test.
# The number of groups may be reduced because of tied predictions.
# Results should be interpreted cautiously in this small dataset.
library(ResourceSelection)
hoslem.test(nodal$r, predicted_prob, g = 10)

# Predicted probability for patients with none / all of the retained factors.
new_patients = data.frame(stage = c(0, 1), xray = c(0, 1),
                          acid = c(0, 1))
predict(fit2, newdata = new_patients, type = "response")


# CONCLUSIONS
# The analysis included 53 patients, of whom 20 (37.7%) presented
# nodal involvement.
# The initial logistic regression model included five predictors:
# aged, stage, grade, xray and acid.
# Backward selection using likelihood-ratio tests identified
# stage, xray and acid as the predictors retained in the final model.
# All three predictors showed statistically significant
# associations with nodal involvement in the selected model.

# Adjusted odds ratios:
# stage = 5.18 (95% CI: 1.24-21.66)
# xray  = 6.76 (95% CI: 1.47-31.02)
# acid  = 5.14 (95% CI: 1.17-22.54)

# The selected model showed a statistically significant improvement
# over the intercept-only model (likelihood-ratio test, p < 0.001).

# At a classification threshold of 0.5, the model achieved:
# Accuracy = 81.1%
# Sensitivity = 80.0%
# Specificity = 81.8%
# AUC = 0.834

# These performance measures were calculated on the training data
# and may overestimate performance on new patients.
# The analysis is limited by the small sample size, variable
# selection uncertainty and the absence of external validation.
# The findings indicate associations between the selected clinical
# factors and nodal involvement, but do not establish causality.