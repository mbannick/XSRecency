# ----------------------------------------------------------------
# TABLES AND FIGURES FOR CROSS-SECTIONAL RECENCY ASSAY PERFORMANCE
# FEB 2021
# ----------------------------------------------------------------
# ----------------------------------------------------------------

rm(list=ls())

library(data.table)
library(xtable)
library(ggplot2)
library(ggh4x)
library(RColorBrewer)

# READ IN VERSIONED RESULTS ---------------------------------

version <- "~/Documents/FileZilla/xs-recent/15-02-21-Feb/"
detail <- fread(paste0(version, "detail.csv"))
summ <- fread(paste0(version , "summary.csv"))

# TABLE RESULTS ---------------------------------------------

summ[, tname := ifelse(itype == "constant", "Constant", ifelse(itype == "linear", "Linear", "Exponential"))]
detail[, tname := ifelse(itype == "constant", "Constant", ifelse(itype == "linear", "Linear", "Exponential"))]

summ[is.na(phi_frr) & is.na(phi_tfrr) & is.na(phi_norm_mu), pname := "Zero"]
summ[(!is.na(phi_frr) | !is.na(phi_tfrr)) & is.na(phi_norm_mu), pname := "Constant"]
summ[(!is.na(phi_frr) | !is.na(phi_tfrr)) & !is.na(phi_norm_mu), pname := "Non-constant"]
detail[is.na(phi_frr) & is.na(phi_tfrr) & is.na(phi_norm_mu), pname := "Zero"]
detail[(!is.na(phi_frr) | !is.na(phi_tfrr)) & is.na(phi_norm_mu), pname := "Constant"]
detail[(!is.na(phi_frr) | !is.na(phi_tfrr)) & !is.na(phi_norm_mu), pname := "Non-constant"]

summ[window == 71, sname := "Brookmeyer et al. 2013"]
summ[window == 248, sname := "Laeyendecker et al. 2018"]
detail[window == 71, sname := "Brookmeyer et al. 2013"]
detail[window == 248, sname := "Laeyendecker et al. 2018"]

settings <- summ$sname %>% unique
phis <- summ$pname %>% unique
epis <- summ$tname %>% unique

data <- data.table()

for(setting in settings){
  for(phi in phis){
    for(epi in epis){

      const <- 100

      this.row <- summ[(sname == setting) & (tname == epi) & (pname == phi)]

      snap.exp <- this.row[estimator == "snap_true", expected_bias]
      adj.exp <- this.row[estimator == "adj_true", expected_bias]

      row <- c(
        epi,
        phi,
        sprintf('%.2f', this.row[estimator == "snap_true", expected_bias] * 100),
        sprintf('%.2f', this.row[estimator == "snap_true", bias] * 100),
        sprintf('%.2f', this.row[estimator == "snap_true", se] * 100),
        sprintf('%.2f', this.row[estimator == "snap_true", see] * 100),
        sprintf('%05.2f', this.row[estimator == "snap_true", cover] * 100),

        sprintf('%.2f', this.row[estimator == "snap_est", bias] * 100),
        sprintf('%.2f', this.row[estimator == "snap_est", se] * 100),
        sprintf('%.2f', this.row[estimator == "snap_est", see] * 100),
        sprintf('%05.2f', this.row[estimator == "snap_est", cover] * 100),

        sprintf('%.2f', this.row[estimator == "adj_true", expected_bias] * 100),
        sprintf('%.2f', this.row[estimator == "adj_true", bias] * 100),
        sprintf('%.2f', this.row[estimator == "adj_true", se] * 100),
        sprintf('%.2f', this.row[estimator == "adj_true", see] * 100),
        sprintf('%05.2f', this.row[estimator == "adj_true", cover] * 100),

        sprintf('%.2f', this.row[estimator == "adj_est", bias] * 100),
        sprintf('%.2f', this.row[estimator == "adj_est", se] * 100),
        sprintf('%.2f', this.row[estimator == "adj_est", see] * 100),
        sprintf('%05.2f', this.row[estimator == "adj_est", cover] * 100)
      )
      data <- rbind(data, t(row))
    }
  }
}

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

# FIGURE OF RESULTS ----------------------------------------------------

detail[estimator == "snap_true", name := "Snapshot Truth"]
detail[estimator == "adj_true", name := "Adjusted Truth"]
detail[estimator == "snap_est", name := "Snapshot Estimated"]
detail[estimator == "adj_est", name := "Adjusted Estimated"]

detail[, name := factor(name,
                        levels=c("Snapshot Truth",
                                 "Snapshot Estimated",
                                 "Adjusted Truth",
                                 "Adjusted Estimated"))]
detail[, pname := factor(pname,
                         levels=c("Zero", "Constant", "Non-constant"))]
detail[, tname := factor(tname,
                         levels=c("Constant", "Linear", "Exponential"))]

detail[, cover_pct := lapply(.SD, mean), .SDcols="cover",
       by=c("tname", "pname", "sname", "name")]

detail[cover_pct < 0.98 & cover_pct >= 0.94, cover_group := "94-98"]
detail[cover_pct < 0.94 & cover_pct >= 0.90, cover_group := "90-94"]
detail[cover_pct < 0.90 & cover_pct >= 0.86, cover_group := "86-90"]
detail[cover_pct < 0.86 & cover_pct >= 0.82, cover_group := "82-86"]
detail[cover_pct < 0.82 & cover_pct >= 0.78, cover_group := "78-82"]
detail[cover_pct < 0.78 & cover_pct >= 0.74, cover_group := "74-78"]

pdf("boxplots.pdf", height=8, width=14)
ggplot() + geom_boxplot(data=detail, aes(x=tname, y=estimate,
                                         fill=cover_group, color=cover_group)) +
  facet_nested(sname ~name + pname, scales="free_y") +
  geom_hline(yintercept=unique(detail$truth),
             color='black', linetype='dashed') +
  scale_colour_manual(values=rev(brewer.pal(7,"BuPu"))) +
  scale_fill_manual(values=rev(brewer.pal(7, "BuPu"))) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(color="Coverage", fill="Coverage") +
  xlab("Incidence") +
  ylab("Estimate")
dev.off()

# PHI FIGURE ----------------------------------

pdf('RAplot.pdf',height=6,width=10)
par(mfrow=c(1,2))

for(i in 1:2){

  if(i == 1){
    mdri <- 71 # 45 # 71
    shadow <- 80 # 250 # 237
  } else {
    mdri <- 248
    shadow <- 306
  }

  t = seq(0, 12, 0.01)
  params <- get.gamma.params(window=mdri/356.25, shadow=shadow/365.25)

  phit <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])

  if(i == 1) ttime <- 2
  if(i == 1) tval <- phit(2)

  if(i == 2) tval <- 0.02
  if(i == 2) ttime <- uniroot(function(t) phit(t) - tval, interval=c(0, 12))$root

  phit.const <- function(t) phit(t)*(t <= ttime) + tval*(t > ttime)
  phit.const.dnorm <- function(t) phit.const(t) + dnorm(t-7, mean=0, sd=1) / 8

  plot(t, phit(t), col='black', type='l', ylab=expression(phi(t)), xlab=expression(t))
  abline(h=tval, lty='dashed')

  ###### MDRI at 2 year 142/365, FRR 1.5%: 1.5% + Gamma distribution with mean 142/365.25-1.5%
  lines(t, phit.const(t), col='red')

  ##### Non-constant FRR -- peak at 7 years
  lines(t, phit.const.dnorm(t), col='blue')

  if(i == 2){
    pgon <- c(seq(2, ttime, 0.01))
    pgon.phi <- 1-pgamma(pgon, shape = params[1], rate = params[2])
    polygon(x=c(pgon, rev(pgon)), y=c(rep(tval, length(pgon)), rev(pgon.phi)), col='lightgrey')
  }

  legend('topright',c('(1)','(2)','(3)'),lty=rep(1, 4),col=c("black", "red", "blue"), cex=0.75)
  # abline(h=0, lty='dashed')
  abline(v=2, lty='dashed')
}

dev.off()

# INCIDENCE FIGURE ----------------------------

Year <- seq(2005,2018,1)
Prevalence <- c(24.6,22.9,24.7,26.1,27.2,34,29.4,32.4,
                30.1,31.4,31.0,30.6,26.6,17.2)/100
plot(Year,Prevalence)

Incidence <- c(0,3.8,3.1,3.9,5.4,5.1,5.3,5.1,
               4.6,4.0,3.7,3.3,3.0,4.0)
plot(Year,Incidence)

# Year 2011 - 2018
Year1 = Year[Year>=2011]
Incidence1 = Incidence[Year>=2011]

res_linear = lm(Incidence1~Year1)
# rho = 0.28 (*10^-2) lambda(2018) = 3.2

res_exp = lm(log(Incidence1)~Year1)
# rho = 0.07, lambda(2018) = 3.2

Prevalence1 = Prevalence[Year>=2011]
mean(Prevalence1) # prevalence 0.29

pdf("incidence-plot.pdf",height=8,width=8)
par(mfrow=c(1,1))
plot(Year1,Incidence1/100,xlab='Year',ylab='Incidence')
lines(Year1,rep(3.2 /100,length(Year1)),lty=2)
lines(Year1, 0.032 + 0.0028 * (2018-Year1))
lines(Year1,(3.2 + 0.28*(2018-Year1))/100,lty=2,col='red')
lines(Year1,3.2 *exp(0.07*(2018-Year1))/100,lty=2,col='blue')
legend('topright',c('Constant','Linear','Exponential'),lty=rep(2,3),col=c(1,2,4))
dev.off()
