# ---------------------------------------------- #
# Effect Size Estimation in Mixed Linear Models
# ---------------------------------------------- #
library(effectsize)
library(lme4)
# ---------------------------------------------- #
# Section 3: Example 
# ---------------------------------------------- #
# ---------------------------------------------- #
# Generating Artificial Data 
# ---------------------------------------------- #
set.seed(123)
n <- 1000
# -
library(mvtnorm)
Sig = matrix(c(1,0.7,-0.5,
               0.7,1,0.1,
               -0.5,0.1,1), nrow = 3, byrow = TRUE)
data = rmvnorm(n, mean = rep(0,nrow(Sig)), sigma=Sig)
X1_raw <- as.numeric(data[,1] > 0.5)
X2_raw <- data[,2]
Z_raw <- as.numeric(cut(data[,3], breaks = c(-Inf,qnorm(seq(1/15,14/15,by =1/15)),Inf))) 
# dependent variable:
Y <- 450 + 20 * X1_raw + X2_raw + 3* Z_raw + rnorm(n, sd = 20)
# independent variables: 
# Binary variable
X1 <- as.factor(X1_raw)
# Continuous variable
X2 <- X2_raw
# Categorical variable with 15 levels
Z <- as.factor(Z_raw)
# ------------------------- #
# Figure 1
# ------------------------- #
# - Histogram of X2 with boxplots of X in the 2 groups of X1
oldpar <- par(mar = c(5,5,2,4), cex.lab = 1.4)
hist(Y, ylim= c(0,300), yaxt="n",main="")
boxplot(Y~X1, horizontal =TRUE, 
        add=TRUE, at= c(220,280),boxwex=20,
        axes=FALSE)
axis(side=2,at = c(0, 50,100,150),labels=c(0, 50,100,150))
axis(side=4,at = c(0, 220,280),labels=c("","0","1"))
par(oldpar)
# ---------------------------------------------- #
# Basic function to compute f2 (to avoid code repetition)
# ---------------------------------------------- #
f2LMM <- function(fit, lmTrue, R) {
  n <- length(fitted(fit))
  bhat <- if (lmTrue) coef(fit) else fixef(fit)
  p <- length(bhat) 
  # Estimated covariance matrix of beta-hat, see Formula (12)
  Vbhat <-vcov(fit)
  # 
  Rb <- R %*% bhat
  # - Cohen's f2 as in Formula (13) with nu=n-p
  f2 <- as.vector((t(Rb) %*% solve(R %*% Vbhat %*% t(R)) %*%  Rb)/(n-p))
  return(f2)
}
# ---------------------------------------------- #
# Section 3.1: Variable of Interest 
# ---------------------------------------------- #
# ------------------------- #
# t-Test for group effect of X1
# ------------------------- #
t.test(Y ~ X1)
# ------------------------- #
# Cohen's d
# ------------------------- #
# -------------#
# from package effectsize
# -------------#
cohens_d(Y ~ X1)$Cohens_d
# --> the unconditional effect size of X1 is medium
# -------------#
# directly from model fit (see Groß & Möller, 2023)
# -------------#
fit1 <-  lm(Y ~ 1 + X1)
sig1 <- summary(fit1)$sigma
coef(fit1)[2]/sig1
# --> the same as before apart from sign
# -------------#
# from Cohen's f2
# -------------#
# --- #
# with package effectsize
# --- #
f2_1 <- cohens_f_squared(model = fit1, partial=FALSE)$Cohens_f2
f2_1
# --- #
# with formula (13) from our paper
# -- #
f2LMM(fit1, lmTrue = TRUE, R = matrix(c(0,1), nrow = 1))
# -- > the same as before
# -------------#
# conversion to Cohen's d
# -------------#
gam1 <- vcov(fit1)[2,2] /sig1^2
w1 <- 0
sqrt((f2_1*(n-2-w1)*gam1))
# --> the same d as before
# ---------------------------------------------- #
# Section 3.2: Additional Fixed Effects 
# ---------------------------------------------- #
# ------------------------- #
# Figure 2
# ------------------------- #
# - Scatter plot of response Y and continuous X2, with points indicating groups of X1
oldpar <- par(mar = c(5,5,2,4), cex.lab = 1.4)
plot(X2, Y, pch = as.character(X1), cex = 0.7)
lines(lowess(X2,Y))
par(oldpar)
# ------------------------- #
# Cohen's f2
# ------------------------- #
fit2 <- lm(Y ~ 1 + X1 + X2)
fit21 <- lm(Y ~ 1 + X2)
# ------------- #
# from package effectsize
# ------------- #
f2_2 <- cohens_f_squared(model = fit2, model2 = fit21, partial=FALSE)$Cohens_f2
f2_2
# ------------- #
# with formula (13)
# ------------- #
f2LMM(fit2, lmTrue = TRUE, R = matrix(c(0,1,0), nrow = 1))
# -- > the same as before
# ------------- #
# generalized Cohen's d (formula (14))
# ------------- #
# - Cohen's d*: effect size of binary grouping variable X1 conditional on X2 
# - Corresponds to regression model with X1 and X2 as independent variables
# where Cohen's d* measures the group effect of X1 given X2 is fixed in the model
gam2 <- vcov(fit2)[2,2] /summary(fit2)$sigma^2
gam2 # see paper
w2 <- 1
sqrt((f2_2*(n-2-w2)*gam2)) # see paper
# --> the effect size now conditioned on X2 is only small, thus the effect of X1 on Y
#     is to some extent explained by X2.
# ---------------------------------------------- #
# Section 3.3: Additional Fixed Effects and Random Effects 
# ---------------------------------------------- #
# ------------------------- #
# Figure 3
# ------------------------- #
# - Mean of Y in the groups of Z where binary variable X1=1
group_means <- by(data.frame(Y,X1), Z, function(x){mean(x[,1])})
group_names <- names(group_means)
out1 = round(as.vector(group_means), 1)
# -
out2 = as.vector(by(data.frame(Y, X1), Z, function(x){sum(x[,2] == 1)}))
# -
oldpar <- par(mar = c(5,5,2,4), cex.lab = 1.4)
plot(out2,out1, xlab = expression('# {X1 = 1}'), 
     ylab='Y mean in Z group', cex=1.2,
     ylim=c(465,505))
points(out2,out1, cex=1.2)
text(out2, out1, labels=group_names, pos=3)
par(oldpar)
# ------------------------- # 
# Computation of f2 in Liner Mixed Model (LMM) according to Section 2.2
# ------------------------- #
# LMM with fixed effects X1, X2, and random effect Z on intercept level (random intercept)
fit3 <- lmer(Y ~ 1 + X1 + X2 + (1|Z))
# - Information from model fit:
f2_3 <- f2LMM(fit3, lmTrue = FALSE, R = matrix(c(0,1,0), nrow = 1))
f2_3 # see paper
# --> Small effect size of X1 conditional on X2 and Z in LMM
# ------------------------- #
# - estimated k factor from mixed model fit
# ------------------------- #
vce <- as.data.frame(VarCorr(fit3))$vcov
k_est <- vce[1]/vce[2] 
k_est # see paper
# ---------------------------------------------- #
# Section 3.4: Coefficient of Determination
# ---------------------------------------------- #
# ------------------------- #
# Coefficient of determination (full model, X1, X2, Z)
# ------------------------- #
# Matrix of Restrictions for F Statistic from Section 3.4
R0 <- matrix(c(0,1,0, 0,0,1), nrow = 2, byrow = TRUE)

f2_0 <- f2LMM(fit3, lmTrue = FALSE, R = R0)
# R_A,B:
R_AB <- (f2_0/(1 + f2_0))
R_AB
# ------------------------- #
# Coefficient of determination (partial model, without X1)
# ------------------------- #
# Fit reduced model without X1
fit31 <- lmer(Y ~ 1 + X2 + (1|Z))
# - Information from fit:
f2_31 <- f2LMM(fit31, lmTrue = FALSE, R = matrix(c(0,1), nrow = 1, byrow = TRUE))
# R_A:
R_A <- (f2_31/(1 + f2_31))
R_A
# ---------------------------- #
# Alternative computation of f2:
# ---------------------------- #
(R_AB - R_A)/(1 - R_AB)
# see paper
# ---------------------------------------- #
# EOF
# ---------------------------------------- #
