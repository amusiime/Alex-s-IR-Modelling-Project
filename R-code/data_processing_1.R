
library(tidyverse)
library(terra)
library(greta)

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


#j-insecticide
#l-location
#t-time

# Linear regression
# A simple, one-variable Bayesian linear regression model using the attitude data
# 
# # variables & priors
# int <- normal(0, 10)
# coef <- normal(0, 10)
# sd <- cauchy(0, 3, truncation = c(0, Inf))
# 
# # linear predictor
# mu <- int + coef * attitude$complaints
# 
# # observation model
# distribution(attitude$rating) <- normal(mu, sd)


# Initial susceptibility

initial_susc_wt <- normal(
  mean = 80,
  sd = 30,
  truncation = c(0, 100)
)


# convert from relative (0-1) to the constrained scale for Initial susceptibility
init_susc_relative <- ilogit(initial_susc_wt)


susc_matrix <- as.matrix(standard_WHO_susc_test3[, 6:9])


# variables & priors
#intercept <- normal(0, 10,truncation = c(0, 100))
coef <- normal(0, 10,dim = ncol(susc_matrix))
sd <- cauchy(0, 3, truncation = c(0, 100))

#linear predictor
mu <- init_susc_relative + coef * standard_WHO_susc_test3$net_access

# observation model
distribution(standard_WHO_susc_test3$died) <- normal(mu, sd)


#Write eight models for each year


# Dispersion parameter

rho <- function(suc) {
  # fraction susceptible
  q <-suc
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
# Likelihood

distribution(standard_WHO_susc_test3$died) <- betabinomial_p_rho(standard_WHO_susc_test3$no_mosquitoes_tested,
                                                                 p=standard_WHO_susc_test3$percent_mortality,
                                                                 rho =rho_year)

m <- model(coef,
           sd,
           mu)


draws <- mcmc(
  model = m,
  n_samples = 1000,
  thin = 5,
  warmup = 1000,
  chains = 5
)





library(bayesplot)
mcmc_trace(
  x = draws)

# new data, from 1995, with hierarchical initial state
# user    system   elapsed 
# 23331.754  9113.193  5618.095 

save.image(file = "temporary/fitted_model.RData")

# check convergence
coda::gelman.diag(draws,
                  autoburnin = FALSE,
                  multivariate = FALSE)

 standard_WHO_susc_test3|> 
  filter(country=="Senegal" & year=="2013") |> view()
 
 
 
  #HR2 deletions
  # Self medication fevers
 # Behavior 
 # 
  
standard_WHO_susc_test3 |> 
  count(site_name)
