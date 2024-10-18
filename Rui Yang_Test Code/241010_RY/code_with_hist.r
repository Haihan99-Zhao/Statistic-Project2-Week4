setwd(here::here())
rm(list = ls())

start_time <- Sys.time()



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
  
  duration_80_density <- dlnorm(1:80, meanlog = 3.152, sdlog = 0.451, log = FALSE)
  duration_80_density_normalise <- duration_80_density / sum(duration_80_density)
  
  if (is.null(t0)) {
    
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
  
  inft <- matrix(rep(0, length_of_year * 29442), nrow = length_of_year, ncol = 29442)
  P <- rep(0, n.rep)
  
  for (i in 1:n.rep){
    
    last_day_of_death <- max(t)
    
    death_duration <- sample(1:80, size = 29442, replace = TRUE, prob = duration_80_density_normalise)
    estimate_death_date <- t0 + death_duration
    
    estimate_death_date[estimate_death_date > last_day_of_death] <- last_day_of_death   # No one can death before get infect
    
    estimate_death <- tabulate(estimate_death_date, nbins = length_of_year)
    
    mean_real_data <- sum(deaths) / length(t)
    
    if (i <= 50){
      change_step <- c(-8, -4, -2, -1, 1, 2, 4, 8)
    }else if (i <= 75){
      change_step <- c(-4, -2, -1, 1, 2, 4)
    }else{
      change_step <- c(-2, -1, 1, 2)
    }
    
    if (bs){
      tmp <- rpois(29442, lambda = mean_real_data)
      deaths_t <- tabulate(tmp, nbins = length_of_year)
    }
    
    p0 <- cal_P(estimate_death, deaths_t)
    
    shuffle <- sample(1:length(t0), length(t0), replace = FALSE)
    step <- sample(change_step, length(t0), replace = TRUE)
    
    round <- 1
    for (j in shuffle){
      
      estimate_death_tmp <- estimate_death
      
      old_date <- t0[j] + death_duration[j]
      new_date <- old_date + step[round]
      
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
    
    
    if (bs == FALSE){
      ## plot to observe convergence situation
      
      layout(matrix(c(1, 2), 1, 2))
      
      # Fig 1.
      # fig1_y_range <- c(1, 1600)
      # plot(1:310, estimate_death, main = "Death: Estimate (blue) & Real (green) | Estimate incidence (orange)", type = "l", col = "blue", 
      #      xlab = "date", ylab = "Number of Victims", ylim = fig1_y_range)
      # lines(1:310, deaths_t,type = "l", col = "darkgreen")
      # Estimate_incidence <- tabulate(t0, nbins = length_of_year)
      # lines(1:310, Estimate_incidence,type = "l", col = "darkorange")
      
      hist(estimate_death, main = "Overlapping Histograms", xlab = "Values", col = rgb(0, 0, 1, 0.5), xlim = c(0, 310), ylim = c(0, 1600), breaks = 300)
      hist(deaths_t, col = rgb(1, 0, 0, 0.5), add = TRUE)
      
      # Fig 2.
      fig2_x_range <- c(1, n.rep)
      plot(1:i, P[1:i], main = "P-Value", type = "l", col = "darkred",
           xlab = "date", ylab = "Value of P", xlim = fig2_x_range)
    } else {
      fig1_y_range <- c(1, 1600)
      plot(1:310, estimate_death, main = "Death: Estimate (blue) & Simulate (green) | Estimate incidence (orange)", type = "l", col = "blue", 
           xlab = "date", ylab = "Number of Victims", ylim = fig1_y_range)
      lines(1:310, deaths_t,type = "l", col = "darkgreen")
      Estimate_incidence <- tabulate(t0, nbins = length_of_year)
      lines(1:310, Estimate_incidence,type = "l", col = "darkorange")
      
      # Fig 2.
      fig2_x_range <- c(1, n.rep)
      plot(1:i, P[1:i], main = "P-Value", type = "l", col = "darkred",
           xlab = "date", ylab = "Value of P", xlim = fig2_x_range)
    }
    
    
  }
  return(list(t0 = t0, inft = inft, P = P))
}



####### Main ########

set.seed(1)

## Read the data

engcov <- read.table("engcov.txt", header = TRUE, sep = " ")

my_data <- engcov[1:150, ]

julian_vecter <- my_data$julian
nhs_vecter <- my_data$nhs

t <- julian_vecter
deaths <- nhs_vecter

output <- deconv(t, deaths, n.rep = 120, bs=FALSE, t0=NULL)

t0_0 <- output$t0
inft_0 <- output$inft
P_0 <- output$P

output1 <- deconv(t, deaths, n.rep = 1, bs=TRUE, t0=t0_0)

t0_1 <- output1$t0
inft_1 <- output1$inft
P_1 <- output1$P




end_time <- Sys.time()
execution_time <- end_time - start_time
print(execution_time)




