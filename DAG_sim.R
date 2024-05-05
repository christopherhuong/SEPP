

library(simDAG)

# initialize empty DAG
dag <- empty_dag()

# parent nodes
dag <- dag + 
  node("age", type="rnorm", mean=50, sd=4) +
  node("sex", type="rbernoulli", p=0.5)

# define relations between parent and child nodes
dag <- dag +
  node("bmi", type="gaussian", parents=c("sex","age"),
      betas=c(1.1, 0.4), intercept=12, error=2) + 
  # bmi = 12 + sex*1.1 + age*0.4 + e~N(0,2)
  node("death", type="binomial", parents=c("age", "bmi"),
       betas=c(0.1, 0.3), intercept=-15)
  # death = -15 + age*0.1 + bmi*0.3

plot(dag)

# simulate data
sim_dat <- sim_from_dag(dag=dag, n_sim = 1e4)

head(sim_dat, 5)

# fit glm()

mod_bmi <- glm(bmi ~ age + sex, data=sim_dat, family="gaussian")
mod_bmi |> summary()
mod_death <- glm(death ~ age + bmi, data=sim_dat, family="binomial")
mod_death |> summary()

# estimate dag from data
dag1 <- empty_dag() +
  node("age", type="rnorm") +
  node("sex", type="rbernoulli") +
  node("bmi", type="gaussian", parents=c("sex", "age")) +
  node("death", type="binomial", parents=c("age", "bmi"))

est_dag <- dag_from_data(dag=dag1, data=sim_dat)
sim_dat1 <- sim_from_dag(est_dag$dag, n_sim=1e4)





# time-varying covariates -------------------------------------------------

dag <- empty_dag() +
  
  node("age", type="rnorm", mean=50, sd=4) + #parent = rnorm
  
  node("sex", type="rbernoulli", p=0.5) + #parent = rbernoulli
  
  node("bmi_t1", type="gaussian", betas=c(1.1, 0.4), parents=c("age", "sex"),
       intercept=12, error=2) +
  
  node("death_t1", type="binomial", parents=c("age", "sex", "bmi_t1"),
       betas=c(0.1, 0.3, 0.1), intercept=-15) +
  
  node("bmi_t2", type="gaussian", parents="bmi_t1", betas=c(1.1),
       intercept=0, error=2) +
  
  node("death_t2", type="binomial", betas=c(0.1, 0.3), 
       parents=c("age", "bmi_t2"), intercept=-15)
       

plot(dag)

sim_dat <- sim_from_dag(dag=dag, n_sim = 1e4)


head(sim_dat)


# discrime-time simulation ------------------------------------------------













# dagR --------------------------------------------------------------------
require("dagR")

# exposure = 0
# outcome = -1
# covariates = 1:total

confound_dag <- dag.init(
  covs = c(1), # specify that single covariate (common cause) is known
  symbols = c('X',  # symbol of primary exposure
              'C',  # common cause
              'O'), # outcome
  arcs = c(0, -1, # arc from primary exposure to outcome
           1, 0, # covariate to exposure
           1, -1), # covariate to outcome
  x.name = "Exposure",
  y.name = "Outcome",
  cov.names = c("Common cause"),
  noxy = T
  
)

confound_dag
summary(confound_dag)

confound_fig <- dag.draw(confound_dag,
                         numbering=F,
                         legend=T,
                         noxy=2)


confound_sim <- dag.sim(
  confound_dag,
  n=100,
 binary = c(1, # x = binary
            0, # covariate = continuous
            0),  # outcome = continuous
 mu = c(0.5, #prevalence p of primary exposure
        120, #mu of C,
        10), #mu of O)
 stdev = c(3, #sd of random noise when generating X
           20, #sd Gaussian random noise of C
           1), #sd of rGaussian andom noise of O
 b = c(5, #direct effect of X->0
       log(2), #direct effect *exp(b)* of C->X
       1.5), #direct effect b of C-> O
       naming=2
       
           
 
)

head(confound_sim)



















