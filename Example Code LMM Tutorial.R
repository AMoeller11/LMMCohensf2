# ---------------------------------------------- #
# Effect Size Estimation in Mixed Linear Models
# ---------------------------------------------- #
# ---------------------------------------------- #
# Examples with data sets and code from Brysbaert & Debeer (2025) 
# ---------------------------------------------- #
library(readxl)
library(effectsize)
library(lme4)
library(lmerTest)
library(reshape2)
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

# ------------------------- #
# Data in Table 12
# ------------------------- #

path <- "C:/Example Code/Table 12.xlsx"
# -
Table_12 <- read_excel(path2)
# Convert to long format, disentangle the variables Participant, Stimulus, and measurement
Table_12long <- melt(Table_12, id.vars=c("Participant","Age"), value.name="Rating")
colnames(Table_12long) <- c("Participant", "Age", "Stimulus", "Rating")
# ------------------------- #
# ---------------------------------------------------------------------- #
# LMM with only random effect for Participant on intercept level
# ---------------------------------------------------------------------- #
fit1 <- lmerTest::lmer(Rating ~ Age + (1 | Participant), data=Table_12long,
                       contrasts = list(Age = "contr.sum"))
summary(fit1)

f2_1 <- f2LMM(fit1,lmTrue = FALSE, R = matrix(c(0,1), nrow = 1))
f2_1
# Calculate eta2 with effectsize package
eta2_1 <- eta_squared(fit1, partial=TRUE)$Eta2
eta2_1
# ------------------------------------------------------------------------- #
# LMM with only random effect for Stimulus, both on intercept and slope level
# ------------------------------------------------------------------------- #
fit2 <- lmerTest::lmer(Rating ~ Age + (Age | Stimulus), data=Table_12long,
                             contrasts = list(Age = "contr.sum"))
summary(fit2)
# --
f2_2 <- f2LMM(fit2,lmTrue = FALSE, R = matrix(c(0,1), nrow = 1))
f2_2
# Calculate eta2 with effectsize package
eta2_2 <- eta_squared(fit2, partial=TRUE)$Eta2
eta2_2
# 
# ------------------------------------------------------------------------- #
# LMM with 2 random effects, one for Participant (intercept level only) 
# and one for Stimulus (intercept and slope level)
# ------------------------------------------------------------------------- #
fit3 <- lmerTest::lmer(Rating ~ Age + (1 | Participant) + (Age | Stimulus), data=Table_12long,
                             contrasts = list(Age = "contr.sum"))
summary(fit3)
#
f2_3 <- f2LMM(fit3,lmTrue = FALSE, R = matrix(c(0,1), nrow = 1))
f2_3
# smaller than f2_2
# Calculate eta2 with effectsize package
eta2_3 <- eta_squared(fit3, partial=TRUE)$Eta2
eta2_3
# smaller than eta2_2
# ---------------------------------------- #
# EOF
# ---------------------------------------- #