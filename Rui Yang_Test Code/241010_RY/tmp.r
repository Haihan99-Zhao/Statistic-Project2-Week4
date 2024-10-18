# setwd(here::here())
rm(list = ls())

start_time <- Sys.time()

## Read the data

engcov <- read.table("engcov.txt", header = TRUE, sep = " ")

my_data <- engcov[1:150, ]

julian_vecter <- my_data$julian
nhs_vecter <- my_data$nhs

t <- julian_vecter
deaths <- nhs_vecter

  
estimate_pois <- rep(0, 310)

round <- 0
for (i in t){
  round <- round + 1
  if (deaths[round] != 0){
    estimate_pois[i] <- rpois(1, lambda = deaths[round])
  }
  
}

print(estimate_pois)
  




















end_time <- Sys.time()
execution_time <- end_time - start_time

print(execution_time)




