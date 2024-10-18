setwd(here::here())
rm(list = ls())

start_time <- Sys.time()

## Read the data

engcov <- read.table("engcov.txt", header = TRUE, sep = " ")

my_data <- engcov[1:150, ]


## Algorithm step 1.

duration_80 <- 1:80
duration_80_density <- dlnorm(duration_80, meanlog = 3.152, sdlog = 0.451, log = FALSE)

duration_80_density_normalise <- duration_80_density / sum(duration_80_density)

## Algorithm step 2.

julian_vecter <- my_data$julian
nhs_vecter <- my_data$nhs

victim_number <- sum(nhs_vecter)    # victim_number = 29442

death_duration <- rep(0, times = victim_number)

death_duration <- sample(duration_80, size = victim_number, replace = TRUE, prob = duration_80_density_normalise)

julian_vecter <- julian_vecter   # To avoid negative days after change

set.seed(1)

cal_P <- function(estimate_data, real_data){

  tmp1 <- (real_data - estimate_data)^2
  tmp2_vec <- rep(1, times = length(estimate_data))
  tmp2 <- pmax(tmp2_vec, estimate_data)
  tmp3 <- tmp1 / tmp2
  p <- sum(tmp3)
  
  return(p)
}

deconv <- function(t, deaths, n.rep = 100, bs = FALSE, t0 = NULL){
  
  length_of_year <- 310
  
  deaths_t <- rep(0, times = length_of_year)
  round <- 1
  for (i in t){
    deaths_t[i] <- deaths[round]
    round <- round + 1
  }

  if (is.null(t0)) {

    duration_80_density <- dlnorm(1:80, meanlog = 3.152, sdlog = 0.451, log = FALSE)
    duration_80_density_normalise <- duration_80_density / sum(duration_80_density)

    death_duration <- sample(1:80, size = 29442, replace = TRUE, prob = duration_80_density_normalise)

    t0 <- rep(0, times = 29442)
    
    round <- 1
    index_of_i <- 0
    for (i in deaths_t){
      index_of_i <- index_of_i + 1
      if (i != 0){
        for (j in 1:i){
          t0[round] = index_of_i - death_duration[round]
          round <- round + 1
        }
      }
    }
  }
  
  t0 <- t0
  
  inft <- matrix(rep(0, length_of_year * 29442), nrow = length_of_year, ncol = 29442)
  P <- rep(0, n.rep)

  for (i in 1:n.rep){
    
    if (i <= 50){
      change_step <- c(-8, -4, -2, -1, 1, 2, 4, 8)
    }else if (i <= 75){
      change_step <- c(-4, -2, -1, 1, 2, 4)
    }else{
      change_step <- c(-2, -1, 1, 2)
    }
    
    last_day_of_death <- max(t)
    
    death_duration <- sample(duration_80, size = 29442, replace = TRUE, prob = duration_80_density_normalise)
    estimate_death_date <- t0 + death_duration
    
    estimate_death_date[estimate_death_date > last_day_of_death] <- last_day_of_death   # No one can death before get infect
    
    estimate_death <- tabulate(estimate_death_date, nbins = length_of_year)
    p0 <- cal_P(estimate_death, deaths_t)

    shuffle <- sample(1:length(t0), length(t0), replace = FALSE)
    step <- sample(change_step, length(t0), replace = TRUE)
    
    round <- 1
    for (j in shuffle){
      
      # estimate_death_tmp <- estimate_death
      # 
      # old_date <- t0[j] + death_duration[j]
      # 
      # if (old_date > last_day_of_death){
      #   old_date <- last_day_of_death
      # }
      # 
      # new_date <- old_date + step[round]
      # 
      # if (new_date > last_day_of_death){
      #   
      #   new_date <- last_day_of_death
      #   
      #   estimate_death_tmp[old_date] <- estimate_death_tmp[old_date] - 1
      #   estimate_death_tmp[new_date] <- estimate_death_tmp[new_date] + 1
      #   
      # } else{
      #   
      #   estimate_death_tmp[old_date] <- estimate_death_tmp[old_date] - 1
      #   estimate_death_tmp[new_date] <- estimate_death_tmp[new_date] + 1
      #   
      # }
      
      estimate_death_tmp <- estimate_death

      old_date <- t0[j] + death_duration[j]
      new_date <- old_date + step[round]
      
      # print(old_date)
      # print(new_date)
      # print(length(estimate_death_tmp))
      
      estimate_death_tmp[old_date] <- estimate_death_tmp[old_date] - 1
      estimate_death_tmp[new_date] <- estimate_death_tmp[new_date] + 1
      
      pi <- cal_P(estimate_death_tmp, deaths_t)

      if (pi <= p0){
        p0 <- pi
        t0[j] <- t0[j] + step[round]
        estimate_death <- estimate_death_tmp
      }
      round <- round + 1
    }
    
    P[i] <- p0
    inft[i,] <- t0
    
    if ( i == 1){
      print(p0)
    }
    
    if ( i %% 10 == 0){
      print(p0)
    }
    
  }
  
  
  t1 <- t0
  # return(t1)
  return(list(t0 = t0, inft = inft, P = P))
}


t <- julian_vecter
deaths <- nhs_vecter

output <- deconv(t, deaths, n.rep = 120, bs=FALSE, t0=NULL)

# output <- deconv(t, deaths, n.rep = 200, bs=FALSE, t0=output)
# aaa <- tabulate(output, nbins = 310)
# aaa

t0 <- output$t0
inft <- output$inft
P <- output$P

# new1 <- output + death_duration
# 
# print(new1)


end_time <- Sys.time()
execution_time <- end_time - start_time

print(execution_time)




