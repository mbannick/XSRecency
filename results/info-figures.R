Year <- seq(2005,2018,1)
Prevalence <- c(24.6,22.9,24.7,26.1,27.2,34,29.4,32.4,30.1,31.4,31.0,30.6,26.6,17.2)/100
plot(Year,Prevalence)

Incidence <- c(0,3.8,3.1,3.9,5.4,5.1,5.3,5.1,4.6,4.0,3.7,3.3,3.0,4.0)
plot(Year,Incidence)

### Year 2011 - 2018
Year1 = Year[Year>=2011]
Incidence1 = Incidence[Year>=2011]
res_linear = lm(Incidence1~Year1)
# rho = 0.28 (*10^-2) lambda(2018) = 3.2
Prevalence1 = Prevalence[Year>=2011]
mean(Prevalence1) # prevalence 0.29

res_exp = lm(log(Incidence1)~Year1)
# rho = 0.07, lambda(2018) = 3.2
# pdf('/Users/rita-gaofei/Desktop/RA_project_Longitudinal/Incidence_plot.pdf',height=8,width=8)
par(mfrow=c(1,1))
plot(Year1,Incidence1/100,xlab='Year',ylab='Incidence')
lines(Year1,rep(3.2 /100,length(Year1)),lty=2)
lines(Year1, 0.032 + 0.0028 * (2018-Year1))
lines(Year1,(3.2 + 0.28*(2018-Year1))/100,lty=2,col='red')
lines(Year1,3.2 *exp(0.07*(2018-Year1))/100,lty=2,col='blue')
legend('topright',c('Constant','Linear','Exponential'),lty=rep(2,3),col=c(1,2,4))
# dev.off()

# NEW ANALYSIS WITHOUT YEAR 2018
Year <- seq(2005,2017,1)
Prevalence <- c(24.6,22.9,24.7,26.1,27.2,34,29.4,32.4,30.1,31.4,31.0,30.6,26.6)/100
plot(Year,Prevalence)

Incidence <- c(0,3.8,3.1,3.9,5.4,5.1,5.3,5.1,4.6,4.0,3.7,3.3,3.0)
plot(Year,Incidence)

### Year 2011 - 2018
Year1 = Year[Year>=2011]
Incidence1 = Incidence[Year>=2011]
res_linear = lm(Incidence1~Year1)
# rho = 0.4 (*10^-2) lambda(2018) = 2.9
Prevalence1 = Prevalence[Year>=2011]
mean(Prevalence1) # prevalence 0.30

res_exp = lm(log(Incidence1)~Year1)
# rho = 0.10, lambda(2018) = 3.0

## GO WITH lambda = 3.0
# pdf('/Users/rita-gaofei/Desktop/RA_project_Longitudinal/Incidence_plot.pdf',height=8,width=8)
pdf('~/OneDrive/Documents/2020_2021/RA/new-inc.pdf',height=6,width=10)
par(mfrow=c(1,1))
plot(Year1,Incidence1/100,xlab='Year',ylab='Incidence',
     main=c("linear rho = 0.004, exp rho = 0.10, lambda0 = 0.003"))
lines(Year1,rep(3.0 /100,length(Year1)),lty=2)
lines(Year1, 0.03 + 0.004 * (2017-Year1))
lines(Year1,(3.0 + 0.4*(2017-Year1))/100,lty=2,col='red')
lines(Year1,3.0 *exp(0.1*(2017-Year1))/100,lty=2,col='blue')
legend('topright',c('Constant','Linear','Exponential'),lty=rep(2,3),col=c(1,2,4))
dev.off()



