`
library(tidyverse)
library(terra)
library(greta)
library(bayesplot)


# Read in the data and filter 5 countries with more data

standard_WHO_susc_test<- read_csv("data/1_standard-WHO-susc-test_complex-subgroup.csv")

standard_WHO_susc_test1<-standard_WHO_susc_test |> 
  select("Country", "Site name" ,"Start year" ,"End year","Insecticide tested",
         "Insecticide class","No. mosquitoes tested",
         "Percent mortality") |> 
  janitor::clean_names() |> 
  filter(start_year==end_year) |> 
  group_by(country,site_name) |> 
  summarise(test_conducted=n()) |> 
  arrange(desc(test_conducted)) |> 
  filter(test_conducted>21) |> 
  ungroup()

# Create a vector site and country names
site_vector<-standard_WHO_susc_test1$site_name
country_vector <-standard_WHO_susc_test1$country

# 2010 is the minimum year in the data

baseline_year=2010

# Baseline insecticide class (Pyrthroid)


# Set baseline year
standard_WHO_susc_test2 <-standard_WHO_susc_test |> 
  select("Country", "Site name" ,"Start year" ,"year"="End year","Insecticide tested",
         "Insecticide class","No. mosquitoes tested",
         "Percent mortality",Longitude,Latitude) |> 
  janitor::clean_names() |> 
  filter(site_name %in% site_vector) |> 
  select(-start_year) |> 
  filter(
    # drop any from before the baseline
    year >= baseline_year) |> 
  
  readr::type_convert() |> 
  mutate(year=as.numeric(year),
    # create an index to the simulation year (in 1-indexed integers)
    year_id = year - baseline_year + 1) |> 
  na.omit() |> 
  distinct(country,site_name,year,insecticide_tested,
  insecticide_class,no_mosquitoes_tested,
  percent_mortality,year_id) |> 
  filter(country=="Burkina Faso",
         insecticide_tested=="deltamethrin",site_name=="VK7")

#To check which sites have more tests
standard_WHO_susc_test2|> 
  group_by(country,site_name) |>
  count(site_name) |> 
  arrange(desc(n)) 

# Read in the ITN access data
itn_access <-read_csv("data/National Unit-data.csv") |> 
  select("country"=Name,Year,Value) |> 
  janitor:: clean_names() |> 
  distinct(country,year,value) |> 
  filter (country=="Burkina Faso")

## Combine net use by country and  per year

# Calculate number that died

standard_WHO_susc_test2.1 <-standard_WHO_susc_test2 |> 
  inner_join(itn_access, by=c("year","country")) |>
  filter(!is.na(value)) |> 
  mutate(no_mosquitoes_tested=as.numeric(no_mosquitoes_tested),
         died=round(percent_mortality/100*no_mosquitoes_tested),0) 


# create indices to categorical vectors
classes <- unique(standard_WHO_susc_test2.1$insecticide_class)
types <- unique(standard_WHO_susc_test2.1$insecticide_tested)
countries <- unique(standard_WHO_susc_test2.1$country)
site <-unique(standard_WHO_susc_test2.1$site_name)
years <- baseline_year - 1 + sort(unique(standard_WHO_susc_test2.1$year_id))


#Create id for categorical variables

standard_WHO_susc_test3<- standard_WHO_susc_test2.1|> 
  mutate(
    
    country_id = match(country, countries),
    class_id = match(insecticide_class, classes),
    type_id = match(insecticide_tested, types),
    site_id=match(site_name, site)
  ) |> 
  na.omit() |> 
  mutate(unique_id=row_number()) |> 
  select(unique_id,insecticide_class,insecticide_tested,site_name,year, no_mosquitoes_tested,
         percent_mortality,died,net_access=value) |> 
  filter(!percent_mortality %in% c(1,3,18))

# Dispersion parameter

rho <- function(suc) {
  # fraction susceptible
  q <-standard_WHO_susc_test3$percent_mortality
  # fraction resistant
  p <- 1 - q
  # fraction remaining susceptible
  q / (q + p )
}


betabinomial_p_rho <- function(N, p, rho) {
  
  # model the observation (betabinomial) sd as a multiplier on the binomial sd,
  # accounting for additional error due to nonindependent sampling of individuals
  # from the population. This is based on the INLA parameterisation
  
  # solve for a and b:
  #   p = a / (a + b)
  #   rho = 1 / (a + b + 1)
  a <- p * (1 / rho- 1)
  b <- a * (1 - p) / p
  
  # define betabinomial according to the greta interface
  beta_binomial(size = N, alpha = a, beta = b)
  
}

# Disperson to 1


rho_year <- normal(0, 0.025,
                   truncation = c(0, 1),
                   dim = 1)



# Initial susceptibility

p_1 <- normal(
  mean = 0.8,
  sd = 0.3,
  truncation = c(0, 1)
)
susc_matrix <- as.matrix(standard_WHO_susc_test3[, 6:9])
# 
# 
x <- as_data(susc_matrix[, "net_access"])
# 
y <- susc_matrix[, "died"]
y1 <- as_data(y[1])
y2 <- as_data(y[2])
y3 <- as_data(y[3])
y4 <- as_data(y[4])

# 
N <- as_data(susc_matrix[,"no_mosquitoes_tested"])


# convert from susceptibility (0-1) to the unconstrained scale for linear modelling

eta_1 <- log(p_1 / (1 - p_1))



# variables & priors

# intercept1 <- normal(0, 1)
# 
# intercept2 <- normal(0, 1)

coef <- normal(0, 1)

# sd <- normal(0, 1, truncation = c(0,Inf))

# mu1<-init_susc_relative

eta_2 <- eta_1 + coef * x[1]
eta_3 <- eta_2 + coef * x[2]
eta_4 <- eta_3 + coef * x[3]


p_2 <- ilogit(eta_2)
p_3 <- ilogit(eta_3)
p_4 <- ilogit(eta_4)

distribution(y1) <- binomial(N[1], p_1)
distribution(y2) <- binomial(N[2], p_2)
distribution(y3) <- binomial(N[3], p_3)
distribution(y4) <- binomial(N[4], p_4)

# Likelihood
# 
# distribution(y) <- betabinomial_p_rho(N, mu2,rho_year)

 sims <- calculate(sd, nsim = 1000)
hist(sims$sd, breaks = 100)
 abline(v = susc_matrix[1,])


#Write eight models for each year


 m1 <- model(p_1, coef)

 # MCMC draws
 
draws <- mcmc(
  model = m1,
  n_samples = 5000,
  thin = 5,
  warmup = 5000,
  chains = 5
)

# Trace Plot

mcmc_trace(
  x = draws)


# TO DO

#1 Autoregressive Model
#2 Evolution Mode;




















#Posterior simulations


post_sims <- calculate(p_1,coef,
                       nsim = 100,
                       values = draws)

param_sim_mat <- cbind(post_sims$p_1,
                       post_sims$coef)



plot(susc_matrix[3,] ~ susc_matrix[4,] , data = susc_matrix)
for (i in 1:nrow(param_sim_mat)) {
  abline(a = param_sim_mat[i, 1],
         b = param_sim_mat[i, 2])
}






 