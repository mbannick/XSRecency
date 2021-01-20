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
lines(Year1,(3.2 + 0.28*(2018-Year1))/100,lty=2,col='red')
lines(Year1,3.2 *exp(0.07*(2018-Year1))/100,lty=2,col='blue')
legend('topright',c('Constant','Linear','Exponential'),lty=rep(2,3),col=c(1,2,4))
# dev.off()

library(latex2exp)

pdf('~/OneDrive/Documents/2020_2021/RA/RAplot.pdf',height=8,width=8)
par(mfrow=c(1,1))

mdri <- 142
shadow <- 150

params <- get.gamma.params(window=mdri/356.25, shadow=shadow/365.25)

###### MDRI at 2 year 142/365, FRR 1.5%: 1.5% + Gamma distribution with mean 142/365.25-1.5%
alpha = 1; x = (142/365.25-0.015)/(1-0.015); beta=alpha/x
t = seq(0,12,0.01)
phit = (1-pgamma(t, shape = params[1], rate = params[2]))*(1-0.015)+0.015
plot(t,phit,type='l',ylab=expression(phi(t)), col='red')

###### window period 142/365, FRR 0%: Gamma distribution with mean 142/(365.25*2)
alpha = 1; x = 142/365.25; beta=alpha/x
phit = 1-pgamma(t, shape = params[1], rate = params[2])
lines(t,phit,col='black')

##### Non-constant FRR -- peak at 7 years
phit = (1-pgamma(t, shape=params[1], rate=params[2]))*(1 - 0.015) + 0.015 + dnorm(t-7, mean=0, sd=1) / 8
lines(t, phit, col='blue')

legend('topright',c('(1)','(2)','(3)'),lty=rep(1,4),col=c("black", "red", "blue"))
# abline(h=0, lty='dashed')
abline(v=2, lty='dashed')
dev.off()

