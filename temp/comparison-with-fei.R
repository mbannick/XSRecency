
set.seed(2346)

# PARAMETERS --------
bigT <- 1
tau <- 10
lambda <- 0.03
prevalence <- 0.2
n_sims <- 1000
N <- 1000
a <- 0.5
b <- 1.5
q <- 0.4
rho <- 0

dat <- generate.raw.data(n_sims=n_sims, n=N, prevalence=prevalence)

sim <- simulate.recent(sim_data=dat, infection.function=infections.con,
                       baseline_incidence=lambda, prevalence=prevalence, rho=0,
                       phi.func=truephifunc_MDRI,
                       summarize=TRUE,
                       ptest.dist=function(n) runif(n, a, b),
                       ptest.prob=q, bigT=bigT)

# Get assay parameters simulation based on external data simulation
# assay <- assay.properties.nsim(100, phi.func=truephifunc_MDRI, bigT=1, tau=8)
# Calculate true assay parameters
true_frr <- true.frr(phi.func=truephifunc_MDRI, bigT=bigT, tau=tau)
true_mdri <- true.window.mdri(phi.func=truephifunc_MDRI, maxT=bigT)
true_window <- true.window.mdri(phi.func=truephifunc_MDRI, maxT=tau)

snapest <- snapshot.estimate(n_r=sim$n_r, n_n=sim$n_n, mu=true_window)
est <- adjusted.estimate(n_r=sim$n_r, n_n=sim$n_n, n_p=sim$n_p,
                         omega=true_mdri, beta=true_frr, big_T=bigT)
est.pt <- adjusted.estimate.pt(n_r_pt=sim$n_r_pt, n_n=sim$n_n, n_p=sim$n_p,
                     omega=true_mdri, beta=true_frr, big_T=bigT,
                     num_beta=sim$num_beta,
                     den_omega=sim$den_omega,
                     den_beta=sim$den_beta)

# MY CODE
est_sd <- sd(est)
est.pt_sd <- sd(est.pt)
my_code <- c(q, a, b, mean(est), mean(est.pt), est_sd, est.pt_sd, (est.pt_sd/est_sd)**2)

# FEI'S CODE

epidemic_para = epidemic_para_func(p=prevalence,lambda=lambda,rho=rho)

res1 = simu_const_truephi_TestTime(seed=2346, N, epidemic_para, nrep=n_sims, q, a, b)
res = apply(res1,2,sd,na.rm=T)
fei_code = c(q, a, b, apply(res1,2,mean,na.rm=T), res, (res[2]/res[1])^2)

res <- rbind(my_code, fei_code)
colnames(res) = c('q','a','b','Est','Est_new','sd','sd_new','RE')
print(res)

Timefun = caltimefun_const(lambda,prevalence)
res = cal_MDRI(0,10, truephifunc_MDRI,Time=1)
MDRI = res$MDRI; FRR = res$FRR
res = NULL

##### Generate Recency Data #####
set.seed(10)
data = generate_recency_TestTime1(N, prevalence, truephifunc_MDRI,Timefun,q,a,b)

Npos = data$Npos; Nneg = data$Nneg; Nrec = data$Nrec
Rec = data$Rec; Test = data$Test; TestTime = data$TestTime; Delta = data$Delta
time <- data$time

inc_res = calInc_MDRI(Nrec, Nneg, N, MDRI, MDRI_sd = 0, FRR, FRR_sd=0, Time=bigT, CI = F)
inc_res2 = calInc_TestTime(Nneg, Rec, Test, TestTime, Delta, truephifunc_MDRI,betaT=FRR, Tstar=bigT)
res = c(inc_res$est, inc_res2$est)

inc_res_me <- adjusted.estimate(n_r=Nrec, n_n=Nneg, n_p=Npos,
                         omega=MDRI, beta=FRR, big_T=bigT)

newTestTime <- TestTime
newTestTime[Test == 0] <- 0
Delta2 <- time >= newTestTime

ri_star <- (Delta2 == 0) & (newTestTime <= bigT)
ri_tild <- 1 - ((newTestTime > bigT) & (Delta2 == 1))
ri_new <- (Rec | ri_star) & ri_tild
Rec_new_me <- sum(ri_new)

# So the problem with mine is that I wasn't including people who didn't have a test.

Tstar <- 1
Rec_new = Rec - Test*Rec*(TestTime>Tstar)*Delta + Test*(1-Rec)*(TestTime<=Tstar)*(1-Delta)
sum(Rec_new)

# recent_ti <- lapply(ptest_times, function(x) -x <= bigT)
recent_ti <- newTestTime <= bigT
int_phi_ti <- sapply(newTestTime, function(x) x - integrate(f=truephifunc_MDRI, lower=0, upper=x)$value)
recent_int_ti <- recent_ti * int_phi_ti
int_phi_ti_ti <- (1 - recent_ti) * (TestTime)

num_beta <- sum(recent_ti)
# This is \sum I(T_i \leq T^*) * int_0^{T_i} (1 - \phi(u)) du
den_omega <- sum(recent_int_ti)
# This is I(T_i > T^*) * T_i
den_beta <- sum(int_phi_ti_ti)

est.pt <- adjusted.estimate.pt(n_r_pt=Rec_new_me,
                               n_n=Nneg,
                               n_p=Npos,
                               omega=true_mdri,
                               beta=true_frr,
                               big_T=bigT,
                               num_beta=num_beta,
                               den_omega=den_omega,
                               den_beta=den_beta)

numerator <- Rec_new_me - FRR * num_beta
denominator <- Nneg * (MDRI - FRR * bigT + (den_omega + FRR * den_beta)/Npos)

val <- numerator / denominator
mean(est.pt)
