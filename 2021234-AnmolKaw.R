rm(list = ls())
set.seed(12345)

cat("Question 1\n")

# Load the data
data <- read.csv("/Users/anmolkaw/Desktop/6th_Semester/SI/Assignment/data_csv.csv")
x <- data$x

# Estimate parameters for Normal distribution
norm_params <- c(mean = mean(x), sd = sd(x))

# Estimate parameters for Exponential distribution
exp_rate <- 1 / mean(x)

# Estimate parameters for Log-normal distribution
lognorm_meanlog <- mean(log(x))
lognorm_sdlog <- sd(log(x))

# Estimate parameters for Gamma distribution using method of moments
gamma_shape <- (mean(x)^2) / (var(x))
gamma_rate <- mean(x) / var(x)

# Set up the sequence of values for x (for plotting PDFs)
x_values <- seq(min(x), max(x), length.out = 300)

# Plot histogram of the data
hist(x, probability = TRUE, main = "Histogram with Fitted Distributions", xlab = "Data Values", breaks = 30, col = "gray")

# Add lines for the PDFs of each fitted distribution
lines(x_values, dnorm(x_values, mean = norm_params['mean'], sd = norm_params['sd']), col = "red", lwd = 2)
lines(x_values, dexp(x_values, rate = exp_rate), col = "blue", lwd = 2)
lines(x_values, dlnorm(x_values, meanlog = lognorm_meanlog, sdlog = lognorm_sdlog), col = "orange", lwd = 2)
lines(x_values, dgamma(x_values, shape = gamma_shape, rate = gamma_rate), col = "purple", lwd = 2)

# Add a legend
legend("topright", legend=c("Normal", "Exponential", "Log-normal", "Gamma"), 
       col=c("red", "blue", "orange", "purple"), lty=1, lwd=2)

# Perform Kolmogorov-Smirnov tests
ks_norm <- ks.test(x, "pnorm", mean = norm_params['mean'], sd = norm_params['sd'])
ks_exp <- ks.test(x, "pexp", rate = exp_rate)
ks_lognorm <- ks.test(x, "plnorm", meanlog = lognorm_meanlog, sdlog = lognorm_sdlog)
ks_gamma <- ks.test(x, "pgamma", shape = gamma_shape, rate = gamma_rate)

# Print KS test results
cat("Kolmogorov-Smirnov Test Results:\n")
cat(sprintf("Normal Distribution: D = %f, p-value = %f\n", ks_norm$statistic, ks_norm$p.value))
cat(sprintf("Exponential Distribution: D = %f, p-value = %f\n", ks_exp$statistic, ks_exp$p.value))
cat(sprintf("Log-normal Distribution: D = %f, p-value = %f\n", ks_lognorm$statistic, ks_lognorm$p.value))
cat(sprintf("Gamma Distribution: D = %f, p-value = %f\n", ks_gamma$statistic, ks_gamma$p.value))

# Determine the best fitting distribution
p_values <- c(ks_norm$p.value, ks_exp$p.value, ks_lognorm$p.value, ks_gamma$p.value)
best_fit <- which.max(p_values)
distribution_names <- c("Normal", "Exponential", "Log-normal", "Gamma")
cat(sprintf("\nThe best fitting distribution is: %s\n", distribution_names[best_fit]))

cat("Question-2\n")

# Estimate parameters for Gamma distribution using method of moments
gamma_shape <- (mean(x)^2) / (var(x))
gamma_rate <- mean(x) / var(x)

cat("MoM Estimates:\n")
cat(sprintf("Shape = %f, Rate = %f\n", gamma_shape, gamma_rate))

# Function to calculate log-likelihood for Gamma distribution
# Parameters are transformed using exp to ensure they are positive
gamma_logLik <- function(log_params) {
  params <- exp(log_params)  # Transform back to original scale
  shape <- params[1]
  rate <- params[2]
  likelihoods <- dgamma(x, shape = shape, rate = rate, log = TRUE)
  -sum(likelihoods)  # Return negative log-likelihood
}

# Initial guesses for parameters (log-transformed)
start_params <- log(c(shape = 1, rate = 1))  # Start at reasonable positive values

# Maximum Likelihood Estimation using optim
gamma_mle <- optim(start_params, gamma_logLik, method = "BFGS")

# Transform parameters back to original scale for output
mle_params <- exp(gamma_mle$par)

cat("Maximum Likelihood Estimates:\n")
cat(sprintf("Shape = %f, Rate = %f\n", mle_params[1], mle_params[2]))

cat("Question-3\n")
cat("Method of Moments Estimates (Potentially Consistent, Not Necessarily Unbiased or Efficient):\n")
# Estimate parameters for Gamma distribution using method of moments
gamma_shape <- (mean(x)^2) / (var(x))
gamma_rate <- mean(x) / var(x)
cat(sprintf("Shape = %f, Rate = %f\n", gamma_shape, gamma_rate))

# Function to calculate log-likelihood for Gamma distribution
# Parameters are transformed using exp to ensure they are positive
gamma_logLik <- function(log_params) {
  params <- exp(log_params)  # Transform back to original scale
  shape <- params[1]
  rate <- params[2]
  likelihoods <- dgamma(x, shape = shape, rate = rate, log = TRUE)
  -sum(likelihoods)  # Return negative log-likelihood
}

# Initial guesses for parameters (log-transformed)
start_params <- log(c(shape = 1, rate = 1))  # Start at reasonable positive values

# Maximum Likelihood Estimation using optim
gamma_mle <- optim(start_params, gamma_logLik, method = "BFGS")

# Transform parameters back to original scale for output
mle_params <- exp(gamma_mle$par)

cat("Maximum Likelihood Estimates (Asymptotically Unbiased, Consistent, and Efficient):\n")
cat(sprintf("Shape = %f, Rate = %f\n", mle_params[1], mle_params[2]))

cat("Question-4\n")

# Calculate the sample mean and variance
sample_mean <- mean(x)
sample_variance <- var(x)

# Method of Moments Estimation for Gamma parameters
gamma_shape <- (sample_mean^2) / sample_variance
gamma_rate <- sample_mean / sample_variance

cat("Method of Moments Estimates:\n")
cat(sprintf("Shape (α) = %f\n", gamma_shape))
cat(sprintf("Rate (β) = %f\n", gamma_rate))


cat("Question 5\n")
# Assume alpha is known or estimate alpha using MoM or MLE
estimated_alpha <- (mean(x)^2) / var(x)  # Method of Moments estimation for alpha

# Calculate the sample mean
sample_mean <- mean(x)

# UMVUE for beta
umvue_beta <- estimated_alpha / sample_mean

cat("UMVUE for Rate (β) assuming alpha is known or estimated:\n")
cat(sprintf("UMVUE Beta = %f\n", umvue_beta))



cat("Question 6\n")
# Parameters
alpha_estimated <- (mean(x)^2) / var(x)  # Estimate of alpha using Method of Moments
n <- length(x)
sample_mean <- mean(x)
beta_estimated <- alpha_estimated / sample_mean  # UMVUE for beta

# Standard error for 1/beta
se_one_over_beta <- sqrt(1 / (n * alpha_estimated))

# Confidence levels and corresponding Z-values
conf_levels <- c(0.01, 0.05, 0.1)
z_values <- qnorm(1 - conf_levels / 2)

# Compute confidence intervals
ci <- sapply(z_values, function(z) {
  error_margin <- z * se_one_over_beta
  lower_bound <- 1 / (1/beta_estimated + error_margin)
  upper_bound <- 1 / (1/beta_estimated - error_margin)
  return(c(lower_bound, upper_bound))
})

# Print results
cat("Confidence Intervals for Beta at Different Confidence Levels:\n")
sapply(1:length(conf_levels), function(i) {
  cat(sprintf("Alpha = %.2f: (%.4f, %.4f)\n", conf_levels[i], ci[1, i], ci[2, i]))
})

cat("Question 7\n")
# Sample statistics
sample_mean <- mean(x)
sample_sd <- sd(x)
n <- length(x)

# Assume mu0 (null hypothesis mean value) is the nearest integer to sample mean or a specific reasonable guess
mu0 <- 0.52

# Z-test for large samples
z_stat <- (sample_mean - mu0) / (sample_sd / sqrt(n))
p_value <- 2 * (1 - pnorm(abs(z_stat)))  # Two-tailed test

# Output results
cat(sprintf("Sample Mean: %.4f\n", sample_mean))
cat(sprintf("Hypothesized Mean (mu0): %g\n", mu0))
cat(sprintf("Z-Statistic: %.4f\n", z_stat))
cat(sprintf("P-Value: %.4f\n", p_value))

# Conclusion
if (p_value < 0.05) {
  cat("Reject the null hypothesis: There is sufficient evidence that the mean is not mu0.\n")
} else {
  cat("Fail to reject the null hypothesis: There is not sufficient evidence that the mean is not mu0.\n")
}


cat("Question 8\n")

# Sample statistics
sample_variance <- var(x)
n <- length(x)

# Assume sigma0 squared (null hypothesis variance value) is the nearest integer or specific reasonable estimate
sigma0_squared <- 0.25

# Chi-squared test for variance
chi_squared_stat <- (n - 1) * sample_variance / sigma0_squared
p_value <- 2 * min(pchisq(chi_squared_stat, df = n-1), 1 - pchisq(chi_squared_stat, df = n-1))  # Two-tailed test

# Output results
cat(sprintf("Sample Variance: %.4f\n", sample_variance))
cat(sprintf("Hypothesized Variance (sigma0 squared): %g\n", sigma0_squared))
cat(sprintf("Chi-Squared Statistic: %.4f\n", chi_squared_stat))
cat(sprintf("P-Value: %.4f\n", p_value))

# Conclusion
if (p_value < 0.05) {
  cat("Reject the null hypothesis: There is sufficient evidence that the variance is not sigma0 squared.\n")
} else {
  cat("Fail to reject the null hypothesis: There is not sufficient evidence that the variance is not sigma0 squared.\n")
}


cat("Question 9\n")

# Parameters estimated using Method of Moments
alpha_est <- (mean(x)^2) / var(x)  # Shape parameter
beta_est <- mean(x) / var(x)       # Rate parameter

# Define the number of bins
k <- 10  # Typically 10 bins are used, but this might need adjustment

# Create histogram bins for observed data
observed <- hist(x, breaks=k, plot=FALSE)$counts

# Define breaks for expected probability calculations
breaks <- hist(x, breaks=k, plot=FALSE)$breaks

# Calculate expected counts using the fitted Gamma parameters
expected <- diff(pgamma(breaks, shape=alpha_est, rate=beta_est)) * length(x)

# Chi-squared Goodness-of-Fit Test
chi_squared_stat <- sum((observed - expected)^2 / expected)
p_value <- 1 - pchisq(chi_squared_stat, df=k-1-2)  # df = number of bins - 1 - number of estimated parameters

# Output the results
cat(sprintf("Chi-Squared Statistic: %.4f\n", chi_squared_stat))
cat(sprintf("P-Value: %.4f\n", p_value))

# Conclusion based on p-value
if (p_value < 0.05) {
  cat("Reject the null hypothesis: The data does not follow the specified Gamma distribution.\n")
} else {
  cat("Fail to reject the null hypothesis: The data follows the specified Gamma distribution.\n")
}






