# INCIDENCE FIGURE ----------------------------

source("R/data-generator.R")

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
