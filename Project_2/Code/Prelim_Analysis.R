library(pwr)

prelim <- read.csv("../DataRaw/PrelimData.csv")

# Correlations from preliminary data (n=30)
round(cor(prelim, use = "complete.obs"), 4)

cor.test(prelim$IL_6,  prelim$CVLT_CNG3)
cor.test(prelim$IL_6,  prelim$CORT_CNG3)
cor.test(prelim$MCP_1, prelim$CVLT_CNG3)
cor.test(prelim$MCP_1, prelim$CORT_CNG3)
cor.test(prelim$IL_6,  prelim$MCP_1)

sapply(prelim, function(x) shapiro.test(x)$p.value)

# Parameters
N <- 175
alpha1 <- 0.05
alpha2 <- 0.0125

# Aim 1: pwr.r.test
obs_r <- c(IL6_Mem = 0.259, IL6_Cort = 0.599, MCP1_Mem = 0.318, MCP1_Cort = 0.685)

sapply(obs_r, function(r) pwr.r.test(n = N, r = r, sig.level = alpha1)$power)
sapply(obs_r, function(r) pwr.r.test(n = N, r = r, sig.level = alpha2)$power)

pwr.r.test(n = N, power = 0.80, sig.level = alpha1)$r
pwr.r.test(n = N, power = 0.80, sig.level = alpha2)$r

r_range <- seq(0.10, 0.50, by = 0.05)
aim1_table <- data.frame(
  r = r_range,
  power_05   = sapply(r_range, function(r) pwr.r.test(n=N, r=r, sig.level=alpha1)$power),
  power_0125 = sapply(r_range, function(r) pwr.r.test(n=N, r=r, sig.level=alpha2)$power)
)
aim1_table
write.csv(aim1_table, "../Results/Aim1_PowerTable.csv", row.names = FALSE)

# Aim 2: pwr.f2.test (interaction term)
# Model: outcome ~ cytokine + SUVR + cytokine:SUVR + 7 covariates (p=10)
u <- 1
v <- N - 10 - 1

pwr.f2.test(u = u, v = v, power = 0.80, sig.level = alpha1)$f2
pwr.f2.test(u = u, v = v, power = 0.80, sig.level = alpha2)$f2

f2_range <- c(0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.15)
aim2_table <- data.frame(
  f2 = f2_range,
  delta_R2   = round(f2_range * 0.70, 4),
  power_05   = sapply(f2_range, function(f2) pwr.f2.test(u=u, v=v, f2=f2, sig.level=alpha1)$power),
  power_0125 = sapply(f2_range, function(f2) pwr.f2.test(u=u, v=v, f2=f2, sig.level=alpha2)$power)
)
aim2_table
write.csv(aim2_table, "../Results/Aim2_PowerTable.csv", row.names = FALSE)

# Power curve figure
r_seq  <- seq(0.05, 0.55, by = 0.005)
f2_seq <- seq(0.005, 0.20, by = 0.002)

pdf("../Results/Project2_PowerCurves.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(5, 4.5, 3, 1))

plot(r_seq,
     sapply(r_seq, function(r) pwr.r.test(n=N, r=r, sig.level=alpha1)$power),
     type = "l", lwd = 2, col = "steelblue",
     xlab = "Correlation |r|", ylab = "Power",
     main = "A. Aim 1: Correlation Test", ylim = c(0, 1), las = 1)
lines(r_seq,
      sapply(r_seq, function(r) pwr.r.test(n=N, r=r, sig.level=alpha2)$power),
      lwd = 2, col = "darkorange", lty = 2)
abline(h = 0.80, lty = 3, col = "gray40")
points(obs_r, rep(0.02, 4), pch = 17, col = "red", cex = 1.2)
text(obs_r, rep(0.09, 4), names(obs_r), cex = 0.65, col = "red", srt = 45, adj = 0)
legend("bottomright",
       legend = c(expression(alpha == 0.05), expression(alpha == 0.0125)),
       col = c("steelblue", "darkorange"), lty = c(1, 2), lwd = 2, bty = "n")
mtext(sprintf("N = %d, two-sided", N), side = 3, line = 0.3, cex = 0.8)

plot(f2_seq,
     sapply(f2_seq, function(f2) pwr.f2.test(u=u, v=v, f2=f2, sig.level=alpha1)$power),
     type = "l", lwd = 2, col = "steelblue",
     xlab = expression(paste("Effect size ", f^2)), ylab = "Power",
     main = "B. Aim 2: Interaction F-test", ylim = c(0, 1), las = 1)
lines(f2_seq,
      sapply(f2_seq, function(f2) pwr.f2.test(u=u, v=v, f2=f2, sig.level=alpha2)$power),
      lwd = 2, col = "darkorange", lty = 2)
abline(h = 0.80, lty = 3, col = "gray40")
abline(v = c(0.02, 0.15), lty = 3, col = "gray70")
text(c(0.025, 0.155), c(0.15, 0.15), c("small", "medium"), cex = 0.7, col = "gray50")
legend("bottomright",
       legend = c(expression(alpha == 0.05), expression(alpha == 0.0125)),
       col = c("steelblue", "darkorange"), lty = c(1, 2), lwd = 2, bty = "n")
mtext(sprintf("N = %d, u = %d, v = %d", N, u, v), side = 3, line = 0.3, cex = 0.8)

dev.off()