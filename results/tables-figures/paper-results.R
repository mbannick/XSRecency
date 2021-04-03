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
library(magrittr)
source("~/repos/XSRecency/R/phi-functions.R")
source("~/repos/XSRecency/R/data-generator.R")

# READ IN VERSIONED RESULTS ---------------------------------

# ORIGINAL
# version <- "~/Documents/FileZilla/xs-recent/15-02-21-Feb/"

# GAMMA PHI FIX
# version <- "~/Documents/FileZilla/xs-recent/25-03-21-21-3/"

# DUONG + GAMMA PHI FIX (MAIN TABLE)
version <- "~/Documents/FileZilla/xs-recent/25-03-21-21-2/"
# These are identical versions
# version <- "~/Documents/FileZilla/xs-recent/29-03-21-18-5000/"

# DUONG + GAMMA PHI FIX + 10 YEAR DUONG
# version <- "~/Documents/FileZilla/xs-recent/29-03-21-16-2/"

# DUONG + GAMMA PHI FIX + N 2000
# version <- "~/Documents/FileZilla/xs-recent/29-03-21-16-2000/"

# DUONG + GAMMA PHI FIX + N 10000
# version <- "~/Documents/FileZilla/xs-recent/29-03-21-16-10000/"

# DUONG + GAMMA PHI FIX + N 5000 + SENSITIVITY EPI
# version <- "~/Documents/FileZilla/xs-recent/02-04-21-09/"

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

summ[window == 101, sname := "1 A-C"]
summ[window == 248, sname := "2 A-C"]
detail[window == 101, sname := "1 A-C"]
detail[window == 248, sname := "2 A-C"]

summ[pname == "Zero" & sname == "1 A-C", assay := "1A"]
summ[pname == "Constant" & sname == "1 A-C", assay := "1B"]
summ[pname == "Non-constant" & sname == "1 A-C", assay := "1C"]
summ[pname == "Zero" & sname == "2 A-C", assay := "2A"]
summ[pname == "Constant" & sname == "2 A-C", assay := "2B"]
summ[pname == "Non-constant" & sname == "2 A-C", assay := "2C"]

detail[pname == "Zero" & sname == "1 A-C", assay := "1A"]
detail[pname == "Constant" & sname == "1 A-C", assay := "1B"]
detail[pname == "Non-constant" & sname == "1 A-C", assay := "1C"]
detail[pname == "Zero" & sname == "2 A-C", assay := "2A"]
detail[pname == "Constant" & sname == "2 A-C", assay := "2B"]
detail[pname == "Non-constant" & sname == "2 A-C", assay := "2C"]

settings <- summ$sname %>% unique
phis <- summ$pname %>% unique
epis <- summ$tname %>% unique

data <- data.table()

for(setting in settings){
  print(setting)
  for(phi in phis){
    print(phi)
    for(epi in epis){
      print(epi)

      const <- 100

      this.row <- summ[(sname == setting) & (tname == epi) & (pname == phi)]
      assay <- this.row[, assay][1]

      snap.exp <- this.row[estimator == "snap_true", expected_bias]
      adj.exp <- this.row[estimator == "adj_true", expected_bias]

      row <- c(
        epi,
        assay,

        "$\\times$",
        # sprintf('%.2f', this.row[estimator == "snap_true", bias] * 100),
        # sprintf('%.2f', this.row[estimator == "snap_true", se] * 100),
        # sprintf('%.2f', this.row[estimator == "snap_true", see] * 100),
        # sprintf('%05.2f', this.row[estimator == "snap_true", cover] * 100),

        sprintf('%.2f', this.row[estimator == "snap_est", bias] * 100),
        sprintf('%.2f', this.row[estimator == "snap_est", se] * 100),
        sprintf('%.2f', this.row[estimator == "snap_est", see] * 100),
        sprintf('%05.2f', this.row[estimator == "snap_est", cover] * 100),

        "$\\times$",
        # sprintf('%.2f', this.row[estimator == "adj_true", bias] * 100),
        # sprintf('%.2f', this.row[estimator == "adj_true", se] * 100),
        # sprintf('%.2f', this.row[estimator == "adj_true", see] * 100),
        # sprintf('%05.2f', this.row[estimator == "adj_true", cover] * 100),

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
addtorow$pos <- list(0, 0, 0, 0, 0, 0, 0, 0, 3, 6, 9, 9, 9, 12, 15)
addtorow$command <- c("\\multicolumn{2}{c}{Setting} & \\multicolumn{5}{c}{Snapshot Estimator (1)} & \\multicolumn{5}{c}{Kassanjee Estimator (2)} \\\\\n",
                      "\\hline ",
                      "Incidence & Assay & \\multicolumn{5}{c}{$\\hat{\\mu}$} & \\multicolumn{5}{c}{$\\hat{\\Omega}_{T^*}$, $\\hat{\\beta}_{T^*}$} \\\\\n",
                      "\\hline ",
                      "& & Asm. & Bias & SE & SEE & Cov & Asm. & Bias & SE & SEE & Cov \\\\\n",
                      "\\hline ",
                      "\\hline ",
                      "\\multicolumn{11}{c}{Recency Assay 1A-C} \\\\\n",
                      "\\hline ",
                      "\\hline ",
                      "\\hline ",
                      "\\multicolumn{11}{c}{Recency Assay 2A-C} \\\\\n",
                      "\\hline ",
                      "\\hline ",
                      "\\hline ")
tab <- xtable(data, align=rep("c", 13), digits=2, caption="Simulation results.")
print(tab, add.to.row = addtorow,
      include.colnames = FALSE, include.rownames = FALSE,
      sanitize.text.function=identity)

# FIGURE OF RESULTS ----------------------------------------------------

detail <- detail[estimator %in% c("snap_est", "adj_est")]

detail[estimator == "snap_est", name := "Snapshot"]
detail[estimator == "adj_est", name := "Adjusted"]

detail[, name := factor(name,
                        levels=c("Snapshot",
                                 "Adjusted"))]
detail[, tname := factor(tname,
                         levels=c("Constant", "Linear", "Exponential"))]

detail[, cover_pct := lapply(.SD, mean), .SDcols="cover",
       by=c("tname", "pname", "sname", "name")]

detail[cover_pct < 0.99 & cover_pct >= 0.97, cover_group := "97-99"]
detail[cover_pct < 0.97 & cover_pct >= 0.93, cover_group := "93-97"]
detail[cover_pct < 0.93 & cover_pct >= 0.90, cover_group := "90-93"]
detail[cover_pct < 0.90, cover_group := "<90"]

detail[pname == "Zero" & sname == "1 A-C", assay := "A"]
detail[pname == "Constant" & sname == "1 A-C", assay := "B"]
detail[pname == "Non-constant" & sname == "1 A-C", assay := "C"]
detail[pname == "Zero" & sname == "2 A-C", assay := "A"]
detail[pname == "Constant" & sname == "2 A-C", assay := "B"]
detail[pname == "Non-constant" & sname == "2 A-C", assay := "C"]

pdf("boxplots.pdf", height=6, width=10)
ggplot() + geom_boxplot(data=detail, aes(x=name, y=estimate,
                                         fill=cover_group)) +
  facet_nested(sname ~ assay + tname, scales="free_y") +
  geom_hline(yintercept=unique(detail$truth),
             color='black', linetype='dashed') +
  scale_fill_manual(values=(brewer.pal(7, "BuPu"))[3:7]) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(color="Coverage", fill="Coverage") +
  xlab("Incidence") +
  ylab("Estimate")
dev.off()

pdf("boxplots-2.pdf", height=6, width=10)
ggplot() + geom_boxplot(data=detail, aes(x=name, y=estimate,
                                         color=name)) +
  facet_nested(sname ~ assay + tname, scales="free_y") +
  geom_hline(yintercept=unique(detail$truth),
             color='black', linetype='dashed') +
  scale_colour_manual(values=c("#ffa75e", "#5ea4ff")) +
  scale_fill_manual(values=c("#ffa75e", "#5ea4ff")) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom") +
  labs(color="Estimator") +
  ylab("Estimate") + xlab("")
dev.off()

# PHI FIGURE ----------------------------------

pdf('RAplot.pdf',height=6,width=10)
par(mfrow=c(1,2))

for(i in 1:2){

  if(i == 1){
    mdri <- 101 # 45 # 71
    shadow <- 194 # 250 # 237
  } else {
    mdri <- 248
    shadow <- 306
  }

  t = seq(0, 12, 1e-3)
  params <- get.gamma.params(window=mdri/365.25, shadow=shadow/365.25)

  phit <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])

  if(i == 1) ttime <- 2
  if(i == 1) tval <- phit(2)
  if(i == 1) name <- "Assay 1"

  if(i == 2) tval <- 0.02
  if(i == 2) ttime <- uniroot(function(t) phit(t) - tval, interval=c(0, 12))$root
  if(i == 2) name <- "Assay 2"

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

  legend('topright', paste0(name, c("A", "B", "C")),
         lty=rep(1, 4),col=c("black", "red", "blue"), cex=0.75)
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

set.seed(10)
e <- runif(10000)
constant.inf <- infections.con(e, t=0, p=0.29, lambda_0=0.032, rho=0.0)
linear.inf <- infections.lin(e, t=0, p=0.29, lambda_0=0.032, rho=0.0028)
expon.inf <- infections.exp(e, t=0, p=0.29, lambda_0=0.032, rho=0.07)

pdf("incidence-plot.pdf",height=5,width=8)
layout(matrix(c(1, 2, 1, 3, 1, 4), ncol=2, byrow=TRUE), c(2, 1), c(1, 1, 1))
plot(Year1,Incidence1/100,xlab='Year',ylab='Incidence')
lines(Year1,rep(3.2 /100,length(Year1)),lty=2)
lines(Year1, 0.032 + 0.0028 * (2018-Year1))
lines(Year1,(3.2 + 0.28*(2018-Year1))/100,lty=2,col='red')
lines(Year1,3.2 *exp(0.07*(2018-Year1))/100,lty=2,col='blue')
legend('topright',c('Constant','Linear','Exponential'),lty=rep(2,3),col=c(1,2,4))
hist(-constant.inf, freq=TRUE, xlab="Years Infected",
     main="Constant Incidence", xlim=c(0, 13), breaks=200, border="grey", col='grey')
hist(-linear.inf, freq=TRUE, xlab="Years Infected",
     main="Linear Incidence", xlim=c(0, 13), breaks=200, border="#ff5454", col='#ff5454')
hist(-expon.inf, freq=TRUE, xlab="Years Infected",
     main="Exponential Incidence", xlim=c(0, 13), breaks=200, border='#6675ff', col='#6675ff')
dev.off()

pdf("incidence-plot.pdf",height=5,width=8)
layout(matrix(c(1, 2, 1, 3, 1, 4), ncol=2, byrow=TRUE), c(2, 1), c(1, 1, 1))
plot(Year1,Incidence1/100,xlab='Year',ylab='Incidence')
lines(Year1,rep(3.2 /100,length(Year1)),lty=2)
lines(Year1, 0.032 + 0.0028 * (2018-Year1))
lines(Year1,(3.2 + 0.28*(2018-Year1))/100,lty=2,col='red')
lines(Year1,3.2 *exp(0.07*(2018-Year1))/100,lty=2,col='blue')
legend('topright',c('Constant','Linear','Exponential'),lty=rep(2,3),col=c(1,2,4))
hist(-constant.inf, freq=TRUE, xlab="Years Infected",
     main="Constant Incidence", xlim=c(0, 13), breaks=200, border="grey", col='grey')
hist(-linear.inf, freq=TRUE, xlab="Years Infected",
     main="Linear Incidence", xlim=c(0, 13), breaks=200, border="#ff5454", col='#ff5454')
hist(-expon.inf, freq=TRUE, xlab="Years Infected",
     main="Exponential Incidence", xlim=c(0, 13), breaks=200, border='#6675ff', col='#6675ff')
dev.off()
