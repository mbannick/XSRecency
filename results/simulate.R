rm(list=ls())
library(XSRecency)
library(data.table)
library(xtable)
library(ggplot2)
library(ggbeeswarm)
library(latex2exp)
source("./results/sim-helpers.R")

set.seed(100)

# number of simulations
n_sims <- 5000

# number screened
n <- 5000

# time being considered
big_T <- 1

# prevalence
p <- 0.29

# baseline incidence
inc <- 0.032

for(type in c("constant", "linear", "exponential")){
  for(phi in c(1, 2, 3)){

    if(type == "constant"){
      inc.function <- c.incidence
      infection.function <- c.infections
      rho <- NA
    } else if(type == "linear"){
      inc.function <- l.incidence
      infection.function <- l.infections
      rho <- 0.0028
    } else {
      inc.function <- e.incidence
      infection.function <- e.infections
      rho <- 0.07
    }

    if(phi == 1){
      phi.func <- phi.character.1
      frr <- 0
      mdri <- 142
    } else if(phi == 2){
      phi.func <- phi.character.1
      frr <- 0.015
      mdri <- 142
    } else {
      phi.func <- phi.character.2
      frr <- 0.015
      mdri <- 142
    }

    sim <- simulate(n_sims=n_sims, n=n,
                    inc.function=inc.function,
                    infection.function=infection.function,
                    baseline_incidence=inc, prevalence=p, rho=rho,
                    phi.func=phi.func, frr=frr, mdri=mdri, big_T=big_T)
    assign(paste(type, phi, sep="_"), sim)
    summ <- summarize.simulation(sim)
    assign(paste(type, phi, "summary", sep="_"), summ)
    rm(sim)
    rm(summ)
  }
}

data <- data.table()

for(phi in c(1, 2, 3)){
  for(type in c("constant", "linear")){
    sim <- get(paste(type, phi, "summary", sep="_"))

    if(type == "constant") tname <- "Constant"
    if(type == "linear") tname <- "Non-constant"

    if(phi == 1) pname <- "Zero"
    if(phi == 2) pname <- "Constant"
    if(phi == 3) pname <- "Non-constant"

    row <- c(
      tname,
      pname,
      sprintf('%.4f', sim$snap_true$bias),
      sprintf('%.3f', sim$snap_true$cover),

      sprintf('%.4f', sim$snap_est$bias),
      sprintf('%.3f', sim$snap_est$cover),

      sprintf('%.4f', sim$adj_true$bias),
      sprintf('%.3f', sim$adj_true$cover),

      sprintf('%.4f', sim$adj_est$bias),
      sprintf('%.3f', sim$adj_est$cover)
    )
    data <- rbind(data, t(row))
  }
}

# TABLE

addtorow <- list()
addtorow$pos <- list(0, 0, 0, 0, 0, 0)
addtorow$command <- c("\\multicolumn{2}{c}{Setting} & \\multicolumn{4}{c}{Snapshot Estimator (1)} & \\multicolumn{4}{c}{Kassanjee Estimator (2)} \\\\\n",
                      "\\hline \\\\\n",
                      "Incidence & FRR & \\multicolumn{2}{c}{$\\mu$} & \\multicolumn{2}{c}{$\\hat{\\mu}$} & \\multicolumn{2}{c}{$\\Omega_T$, $\\beta_T$} & \\multicolumn{2}{c}{$\\hat{\\Omega_T}$, $\\hat{\\beta_T}$} \\\\\n",
                      "\\hline \\\\\n",
                      "& & Bias & Cov & Bias & Cov & Bias & Cov & Bias & Cov \\\\\n",
                      "\\hline \\\\\n")
tab <- xtable(data, align=rep("c", 11), digits=2, caption="Simulation results.")
print(tab, add.to.row = addtorow,
      include.colnames = FALSE, include.rownames = FALSE)

# FIGURE A -- BIAS

all.fig.data <- data.table()

for(phi in c(1, 2, 3)){
  for(type in c("constant", "linear")){

    if(type == "constant") tname <- "A. Constant Incidence"
    if(type == "linear") tname <- "B. Non-constant Incidence"

    if(phi == 1) pname <- "1. Zero FRR"
    if(phi == 2) pname <- "2. Constant FRR"
    if(phi == 3) pname <- "3. Non-constant FRR"

    sim <- get(paste(type, phi, sep="_"))
    fig.data <- do.call(cbind, sim) %>% data.table
    names(fig.data) <- names(sim)
    fig.data[, phi := pname]
    fig.data[, type := tname]
    fig.data[, num := .I]

    all.fig.data <- rbind(all.fig.data, fig.data)

  }
}

sub <- all.fig.data[, .(type, phi, snap_true_est, snap_est_est,
                        adj_true_est, adj_est_est)]
setnames(sub,
         c("snap_true_est", "snap_est_est",
           "adj_true_est", "adj_est_est"),
         c("Snapshot,\nTrue", "Snapshot,\nEstimated",
           "Kassanjee,\nTrue", "Kassanjee,\nEstimated"))
sub <- melt(sub, id.vars=c("type", "phi"))

p <- ggplot(data=sub, mapping=aes(x=variable, y=value, fill=variable)) +
  geom_boxplot(outlier.alpha=0.5, lwd=0.2, outlier.size=0.5) +
  geom_hline(yintercept=mean(all.fig.data$truth), linetype='dashed', lwd=0.25) +
  facet_wrap(phi ~ type, nrow=3, ncol=2) +
  theme_minimal() +
  scale_fill_brewer(palette="Set2") +
  scale_color_brewer(palette="Set2") +
  ylab("Estimate") + xlab("Estimator") + labs(fill="Estimator") +
  theme(legend.position="bottom") +
  # ggtitle(TeX("Bias of Estimators for $\\lambda_0$")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(panel.grid.minor = element_line(size = 0.25),
        panel.grid.major = element_line(size = 0.25)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

pdf(file="/Users/marlena/OneDrive/Documents/2020_2021/RA/simulation-plots-20200107/bias-3.pdf",
    height=9, width=7)
p
dev.off()

# FIGURE B -- COVERAGE

cov.data <- data.table()

for(phi in c(1, 2, 3)){
  for(type in c("constant", "linear")){

    if(type == "constant") tname <- "A. Constant Inc."
    if(type == "linear") tname <- "B. Non-constant Inc."

    if(phi == 1) pname <- "1. Zero FRR"
    if(phi == 2) pname <- "2. Constant FRR"
    if(phi == 3) pname <- "3. Non-constant FRR"

    sim <- get(paste(type, phi, "summary", sep="_"))

    snap_true <- sim$snap_true$cover
    snap_est <- sim$snap_est$cover
    adj_true <- sim$adj_true$cover
    adj_est <- sim$adj_est$cover

    sv <- data.table(tname=tname, pname=pname,
                     snap_true=snap_true, snap_est=snap_est,
                     adj_true=adj_true, adj_est=adj_est)

    cov.data <- rbind(cov.data, sv)

  }
}

setnames(cov.data,
         c("snap_true", "snap_est",
           "adj_true", "adj_est"),
         c("Snapshot,\nTrue", "Snapshot,\nEstimated",
           "Kassanjee,\nTrue", "Kassanjee,\nEstimated"))
cov.data <- data.table::melt(cov.data, id.vars=c("tname", "pname"))
cov.data[, name := paste(tname, pname, sep="\n")]

p2 <- ggplot(data=cov.data, mapping=aes(x=name, y=value, color=variable, group=variable)) +
  geom_beeswarm(priority='random', cex=3.5, size=3, groupOnX=T) +
  # geom_line(linetype='dashed') +
  # facet_wrap(~tname, nrow=2) +
  geom_hline(yintercept=0.95, linetype='dashed', lwd=0.25) +
  theme_minimal() +
  scale_color_brewer(palette="Set2") +
  ylab("Estimate") + xlab("") + labs(color="Estimator") +
  theme(legend.position="bottom") +
  theme(panel.grid.minor = element_line(size = 0.25),
        panel.grid.major = element_line(size = 0.25))

pdf(file="/Users/marlena/OneDrive/Documents/2020_2021/RA/simulation-plots-20200107/coverage-summary-2.pdf",
    height=7, width=10)
p2
dev.off()

# FIGURE C -- ASSAY ESTIMATES

melt.est <- data.table::melt(all.fig.data, id.vars=c("type", "phi", "truth", "num"),
                             measure.vars=patterns("est$", cols=names(all.fig.data)),
                             variable.name="estimator", value.name="estimate")
melt.est[, estimator := lapply(.SD, function(x) gsub("_est$", "", x)), .SDcols="estimator"]
melt.var <- data.table::melt(all.fig.data, id.vars=c("type", "phi", "truth", "num"),
                             measure.vars=patterns("var$", cols=names(all.fig.data)),
                             variable.name="estimator", value.name="variance")
melt.var[, estimator := lapply(.SD, function(x) gsub("_var$", "", x)), .SDcols="estimator"]
melted <- merge(melt.est, melt.var, by=c("type", "phi", "truth", "num", "estimator"))

melted[, lower := estimate - qnorm(0.975) * variance ** 0.5]
melted[, upper := estimate + qnorm(0.975) * variance ** 0.5]
melted[, cover := (lower < truth) & (upper > truth)]
melted.b <- melted[, lapply(.SD, mean), .SDcols="estimate", by=c("type", "phi", "truth", "estimator")]

melted.small <- melted[num <= 200]

test <- melted.small[type == "A. Constant Incidence" & phi == "1. Zero FRR"]
test.b <- melted.b[type == "A. Constant Incidence" & phi == "1. Zero FRR"]
truth <- unique(test$truth)
ggplot(data=test, mapping=aes(x=num, y=estimate, color=cover)) +
  geom_pointrange(aes(ymin=lower, ymax=upper, color=cover), size=0.1) +
  geom_hline(yintercept=truth) +
  geom_hline(data=test.b, aes(yintercept=estimate), linetype='dashed') +
  theme_minimal() +
  theme(legend.position="bottom") +
  theme(panel.grid.minor = element_line(size = 0.25),
        panel.grid.major = element_line(size = 0.25)) +
  scale_color_brewer(palette="Set1") +
  labs("Coverage") + xlab("Simulation") + ylab("Estimate") +
  facet_wrap(~ estimator, nrow=4)

# MAKE THIS INTO A GRID WITH GEOM_GRID
# ALSO PUT THE ESTIMATORS IN THE RIGHT ORDER
