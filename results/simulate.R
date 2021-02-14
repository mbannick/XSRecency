rm(list=ls())
library(XSRecency)
library(data.table)
library(xtable)
library(ggplot2)
library(ggbeeswarm)
library(latex2exp)
library(magrittr)
library(gridExtra)
setwd("~/repos/XSRecency/")
source("./results/sim-helpers.R")

date <- format(Sys.time(), "%d-%m-%y-%h")
out_dir <- paste0("/Users/marlena/OneDrive/Documents/2020_2021/RA/simulation-plots-", date, "/")
dir.create(out_dir)

re.run <- TRUE

set.seed(100)

# number of simulations
n_sims <- 10

# number screened
n <- 5000

# prevalence
p <- 0.29

# baseline incidence
inc <- 0.032

# Which assay settings to use (1 or 2)
# setting <- 1

# out_dir <- paste0(out_dir, "setting_", setting, "-")
out_dir <- paste0(out_dir, "nsims_", n_sims)
dir.create(out_dir)

if(re.run){
  for(setting in c(1, 2)){
    for(type in c("constant", "linear", "exponential")){
      cat("Working on ", type, "\n")

      for(phi in c(1, 2, 3)){
        cat("Working on ", phi, "\n")
        set.seed(100)

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

        if(phi != 3) next

        if(setting == 1){
          window <- 71 # 45 # 80 # 45 # 71
          shadow <- 80 # 250 # 25 # 250 # 237
        } else {
          window <- 248
          shadow <- 306
        }
        params <- get.gamma.params(window=window/356.25, shadow=shadow/365.25)
        phit <- function(t, ...) 1-pgamma(t, shape = params[1], rate = params[2])

        if(setting == 1) ttime <- 2
        if(setting == 1) tval <- phit(2)

        if(setting == 2) tval <- 0.02
        if(setting == 2) ttime <- uniroot(function(t) phit(t) - tval, interval=c(0, 12))$root

        phit.const <- function(t, ...) phit(t)*(t <= ttime) + tval*(t > ttime)
        phit.const.dnorm <- function(t, ...) phit.const(t) + dnorm(t-7, mean=0, sd=1) / 8

        if(phi == 1){
          phi.func <- phit
          frr <- 0
        } else if(phi == 2){
          phi.func <- phit.const
          frr <- 0.015
        } else {
          phi.func <- phit.const.dnorm
          frr <- 0.015
        }

        sim <- simulate(n_sims=n_sims, n=n,
                        inc.function=inc.function,
                        infection.function=infection.function,
                        baseline_incidence=inc, prevalence=p, rho=rho,
                        phi.func=phi.func, frr=frr, window=window, shadow=shadow,
                        big_T=2)
        assign(paste(type, phi, setting, sep="_"), sim)
        summ <- summarize.simulation(sim)
        assign(paste(type, phi, setting, "summary", sep="_"), summ)
        rm(sim)
        rm(summ)
      }
    }
  }

  # for(setting in c(1, 2)){
  #   for(phi in c(1, 2, 3)){
  #     for(type in c("constant", "linear", "exponential")){
  #       name <- paste(type, phi, setting, "summary", sep="_")
  #       sim <- get(name)
  #       save(sim, file=paste0(name, ".Rdata"))
  #
  #       name <- paste(type, phi, setting, sep="_")
  #       sim <- get(name)
  #       save(sim, file=paste0(name, ".Rdata"))
  #     }
  #   }
  # }
}
# for(setting in c(1, 2)){
#   for(phi in c(1, 2, 3)){
#     for(type in c("constant", "linear", "exponential")){
#       name <- paste(type, phi, setting, "summary", sep="_")
#       load(file=paste0(name, ".Rdata"))
#
#       name <- paste(type, phi, setting, sep="_")
#       load(file=paste0(name, ".Rdata"))
#     }
#   }
# }

data <- data.table()

for(setting in c(1, 2)){
  for(phi in c(1, 2, 3)){
    for(type in c("constant", "linear", "exponential")){
      sim <- get(paste(type, phi, setting, "summary", sep="_"))

      if(type == "constant") tname <- "Constant"
      if(type == "linear") tname <- "Linear"
      if(type == "exponential") tname <- "Exponential"

      if(phi == 1) pname <- "(1) Zero"
      if(phi == 2) pname <- "(2) Constant"
      if(phi == 3) pname <- "(3) Non-constant"

      if(setting == 1) sname <- "Brookmeyer et al. 2013"
      if(setting == 2) sname <- "Laeyendecker et al. 2018"

      row <- c(
        # sname,
        tname,
        pname,
        sprintf('%.2f', sim$snap_bias * 100),
        sprintf('%.2f', sim$snap_true$bias * 100),
        sprintf('%.2f', sim$snap_true$se * 100),
        sprintf('%.2f', sim$snap_true$see * 100),
        sprintf('%05.2f', sim$snap_true$cover * 100),

        sprintf('%.2f', sim$snap_est$bias * 100),
        sprintf('%.2f', sim$snap_est$se * 100),
        sprintf('%.2f', sim$snap_est$see * 100),
        sprintf('%05.2f', sim$snap_est$cover * 100),

        sprintf('%.2f', sim$adj_bias * 100),
        sprintf('%.2f', sim$adj_true$bias * 100),
        sprintf('%.2f', sim$adj_true$se * 100),
        sprintf('%.2f', sim$adj_true$see * 100),
        sprintf('%05.2f', sim$adj_true$cover * 100),

        sprintf('%.2f', sim$adj_est$bias * 100),
        sprintf('%.2f', sim$adj_est$se * 100),
        sprintf('%.2f', sim$adj_est$see * 100),
        sprintf('%05.2f', sim$adj_est$cover * 100)
      )
      data <- rbind(data, t(row))
    }
  }
}

# TABLE

addtorow <- list()
addtorow$pos <- list(0, 0, 0, 0, 0, 0, 0, 3, 6, 9, 9, 9, 12, 15, 18)
addtorow$command <- c("\\multicolumn{2}{c}{Setting} & \\multicolumn{9}{c}{Snapshot Estimator (1)} & \\multicolumn{9}{c}{Kassanjee Estimator (2)} \\\\\n",
                      "\\hline ",
                      "Incidence & FRR & \\multicolumn{5}{c}{$\\mu$} & \\multicolumn{4}{c}{$\\hat{\\mu}$} & \\multicolumn{5}{c}{$\\Omega_{T^*}$, $\\beta_{T^*}$} & \\multicolumn{4}{c}{$\\hat{\\Omega}_{T^*}$, $\\hat{\\beta}_{T^*}$} \\\\\n",
                      "\\hline ",
                      "& & E[Bias] & Bias & SE & SEE & Cov & Bias & SE & SEE & Cov & E[Bias] & Bias & SE & SEE & Cov & Bias & SE & SEE & Cov \\\\\n",
                      "\\hline ",
                      "\\multicolumn{19}{c}{Brookmeyer et al. 2013} \\\\\n",
                      "\\hline ",
                      "\\hline ",
                      "\\hline ",
                      "\\multicolumn{19}{c}{Laeyendecker et al. 2018} \\\\\n",
                      "\\hline ",
                      "\\hline ",
                      "\\hline ",
                      "\\hline ")
tab <- xtable(data, align=rep("c", 21), digits=2, caption="Simulation results.")
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
  ylab("Estimate") + xlab("Estimator") + labs(fill="Estimator") +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(panel.grid.minor = element_line(size = 0.25),
        panel.grid.major = element_line(size = 0.25)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

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

# FIGURE C -- DETAILED ESTIMATOR

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

melted.small <- melted[num <= 100]

plots <- list()
i <- 1

for(phi in c(1, 2, 3)){
  for(type in c("constant", "linear")){

    if(type == "constant") tname <- "A. Constant Incidence"
    if(type == "linear") tname <- "B. Non-constant Incidence"

    if(phi == 1) pname <- "1. Zero FRR"
    if(phi == 2) pname <- "2. Constant FRR"
    if(phi == 3) pname <- "3. Non-constant FRR"

    lab <- i %in% c(5, 6)

    sub.melted <- copy(melted[type == tname & phi == pname])
    sub.data <- copy(melted.small[type == tname & phi == pname])
    b.data <- copy(melted.b[type == tname & phi == pname])
    truth <- unique(sub.melted$truth)

    pp <- ggplot(data=sub.data, mapping=aes(x=num, y=estimate, color=cover)) +
      geom_pointrange(aes(ymin=lower, ymax=upper, color=cover), size=0.3, fatten=0.5) +
      geom_hline(yintercept=truth) +
      geom_hline(data=b.data, aes(yintercept=estimate), linetype='dashed',
                 lwd=0.5) +
      theme_minimal() +
      theme(legend.position="none") +
      theme(panel.grid.minor = element_line(size = 0.25),
            panel.grid.major = element_line(size = 0.25)) +
      scale_color_brewer(palette="Set1") +
      labs("Coverage") + xlab("Simulation") + ylab("Estimate") +
      facet_wrap(~ factor(estimator,
                          levels=c("snap_true", "snap_est", "adj_true", "adj_est")),
                 nrow=4, strip.position="right") +
      ggtitle(paste0(tname, "\n", pname))
    if(lab){
      pp <- pp + theme(axis.title.x=element_blank())
    }

    plots[[i]] <- pp
    i <- i + 1

  }
}

# FIGURE D -- RECENCY ASSAY ESTIMATES
cat("Working on recency assay simulations.")

phi_sims <- list()

for(phi in c(1, 2, 3)){

  if(setting == 1){
    window <- 71
    shadow <- 237
  } else {
    window <- 248
    shadow <- 306
  }

  if(phi == 1){
    phi.func <- phi.character.1
    frr <- 0
  } else if(phi == 2){
    phi.func <- phi.character.1
    frr <- 0.015
  } else {
    phi.func <- phi.character.2
    frr <- 0.015
  }

  if(phi == 1) pname <- "1. Zero FRR"
  if(phi == 2) pname <- "2. Constant FRR"
  if(phi == 3) pname <- "3. Non-constant FRR"

  simulations <- replicate(n_sims, assay.properties.sim(phi.func, window, frr, shadow), simplify="matrix") %>% t
  simulations <- data.table(simulations)

  simulations[, case := pname]
  simulations[, constant := phi %in% c(1, 2)]
  simulations[, frr := frr]
  simulations[, window := window]
  mdri <- true.mdri(phi.func, window=window, frr=frr, shadow=shadow)
  simulations[, mdri := mdri]
  simulations[, num := .I]
  phi_sims[[phi]] <- simulations
}

df <- rbindlist(phi_sims)

columns <- colnames(df)
cols <- lapply(columns, function(x) unlist(df[[x]]))
df <- do.call(cbind, cols) %>% data.table
names(df) <- columns

mu <- melt(df, id.vars=c("case", "constant", "frr", "window", "mdri", "num"), measure.vars="mu_est")
omega <- melt(df, id.vars=c("case", "constant", "frr", "window", "mdri", "num"), measure.vars="omega_est")
beta <- melt(df, id.vars=c("case", "constant", "frr", "window", "mdri", "num"), measure.vars="beta_est")

df <- rbindlist(list(mu, omega, beta))
df[, true := ifelse(variable == "beta_est", as.numeric(frr), as.numeric(window)/365.25)]
df[, value := as.numeric(value)]
df[variable == "mu_est", variable := "Mean Window Period"]
df[variable == "omega_est", variable := "MDRI"]
df[variable == "beta_est", variable := "FRR"]

df.truth <- df[, lapply(.SD, mean), .SDcols="true", by=c("case", "variable")]

rp <- ggplot(data=df) + geom_boxplot(aes(y=value)) +
  geom_hline(data=df.truth, aes(yintercept=true), color='red') +
  facet_grid(variable ~ case, scales="free_y") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Output all of the plots

pdf(file=paste0(out_dir, "/bias-3.pdf"),
    height=9, width=7)
p
dev.off()

pdf(file=paste0(out_dir, "/coverage-summary-2.pdf"),
    height=7, width=10)
p2
dev.off()

pdf(file=paste0(out_dir, "/estimator-detailed.pdf"),
    height=11, width=9)
grid.arrange(grobs=plots, nrow=3, ncol=2)
dev.off()

pdf(file=paste0(out_dir, "/recency.pdf"),
    height=6, width=6)
rp
dev.off()
