rm(list = ls())
start_time <- Sys.time()


cal_P <- function(real_data, estimate_data){
  tmp1 <- (real_data - estimate_data)^2
  tmp2_vec <- rep(1, times = length(estimate_data))
  tmp2 <- pmax(tmp2_vec, estimate_data)
  tmp3 <- tmp1 / tmp2
  p <- sum(tmp3)
  
  return(p)
}


deconv <- function(t,deaths,n.rep=100,bs=FALSE,t0=NULL){
  length_max <- 310
  inft <- matrix(0, nrow = length_max, ncol = n.rep)
  P <- rep(0, n.rep)

  #deaths_t means deaths occurred at day t
  #Following codes replace the function "tabulate"
  deaths_t <- rep(0, times = length_max)
  round <- 1
  for (i in t){
    deaths_t[i] <- deaths[round]
    round <- round + 1
  }
  
  #Step 1 
  #Find the vector of probabilities and normalise.
  ill_dlnorm <- dlnorm(1 : 80, meanlog = 3.152, sdlog = 0.451) 
  prob_ill <- ill_dlnorm / sum(ill_dlnorm)
  #No problem
  
  if (is.null(t0)) {
    ill_duration <- rep(0,29442)
    t0 <- rep(0, 29442)
    
    
    #      round = 1
    #      n_index <- 1
    #      for (i in t){
    #       if(i <= 80){
    #          for (n in 1:deaths[n_index]){
    #            ill_duration[round] <- sample(1:i, size= 1, prob= prob_ill[1:i])
    #            round = round+1
    #         }
    #         n_index <- 1 + n_index
    #       }
    #       else{
    #         for (n in 1:deaths[n_index]){
    #           ill_duration[round] <- sample(1:80, size= 1, prob= prob_ill)
    #           round = round+1
    #         }
    #         n_index <- 1 + n_index
    #       }
    #     }
    #print(ill_duration)
    ill_duration = sample(1:80, size = length(t0), prob = prob_ill, replace = TRUE)
    
    round <- 1
    index_of_i <- 0
    for (i in deaths_t){
      index_of_i <- index_of_i + 1
      if (i != 0){
        for (j in 1:i){
          t0[round] = index_of_i - ill_duration[round]
          round <- round + 1
        }
      }
    }
  }
  
  #!!!!!!!!!!!!!!!!!!!!!Step 3.
  #Repeat following steps n.rep times
  
 if(n.rep > 0) 
  for (j in 1:n.rep) {
    new_ill_duration <- rep(0 , 29442)
    abnormal_infection_date1 <- which(t0 > max(t))
    t0[abnormal_infection_date1] <- max(t)
    abnormal_infection_date2 <- which(t0 < 0)
    t0[abnormal_infection_date2] <- 1
    # round = 0
    # n_index <- 0
    # for (i in t){
    #   n_index <- 1 + n_index
    #   if(i <= 80){
    #     for (n in 1:deaths[n_index]){
    #       round = round+1
    #       new_ill_duration[round] <- sample(1:i, size= 1, prob= prob_ill[1:i])
    #     }
    #   }
    #   else{
    #     for (n in 1:deaths[n_index]){
    #       round = round+1
    #       new_ill_duration[round] <- sample(1:80, size= 1, prob= prob_ill)
    #     }
    #   }
    # }
    
    
    new_ill_duration <-  sample(1:80, size = length(t0), prob = prob_ill, replace = TRUE)
    
    estimate_death_date <- t0 + new_ill_duration
    
    #abnormal_death_date <- which(estimate_death_date >max(t))
    #estimate_death_date[abnormal_death_date] <- max(t)
    
    d_simulation <- tabulate(estimate_death_date, nbins = length_max)
    
    #Count the number of deaths in each days
    d_real <- deaths_t
    P_value <- cal_P(d_real, d_simulation)
    #create new t0 vector
    
    new_t0 <- c()
    if (j <= 50){
      random_steps <- c(-8, -4, -2, -1, 1, 2, 4, 8)
    }else if (j <= 75){
      random_steps <- c(-4, -2, -1, 1, 2, 4)
    }else{
      random_steps <- c(-2, -1, 1, 2)
    }
    
    steps_list <- sample(random_steps, size = length(t0), replace = TRUE)
    shuffle_location<-sample(1:length(t0), size = length(t0), replace = FALSE )
    
    round <- 1
    
    #迭代每个p，
    for (i in shuffle_location){
      
      new_d_simulation <- d_simulation
      
      old_esdeath_date <- estimate_death_date[i]
      
      new_t0 <- t0[i] + steps_list[round]
      new_esdeath_date <- new_t0 + new_ill_duration[i]
      
      new_d_simulation[old_esdeath_date] <- new_d_simulation[old_esdeath_date] - 1
      new_d_simulation[new_esdeath_date] <- new_d_simulation[new_esdeath_date] + 1
      
      new_P_value <- cal_P(d_real,new_d_simulation)
      
      if (new_P_value <= P_value){
        P_value <- new_P_value
        t0[i] <- new_t0
        d_simulation <- new_d_simulation
      }
      round <- round + 1
    }
    
    P[j] <- P_value
    inft[1:length_max,j] <- tabulate(t0, nbins = length_max)
    
    
    #plot the graphs
    layout(matrix(c(1, 2), 1, 2))
    
    fig1_y_range <- c(1, 1600)
    
    plot(1:length_max, d_real, main = "Basic situation", type = "l", col = "blue", 
         xlab = "date", ylab = "Number of Victims", ylim = fig1_y_range)
    deaths_t <- deaths_t
    lines(1:length_max, d_simulation,type = "l", col = "darkgreen")
    Estimate_incidence <- tabulate(t0, nbins = length_max)
    lines(1:length_max, Estimate_incidence,type = "l", col = "darkorange")
    
    legend("topright",                                      # Define legend position
           legend = c("Estimate Death", "Real Death", "Estimate Incidence"),  # Define legend text
           col = c("blue", "darkgreen", "darkorange"),      # Define color of legend position
           lty = 1,
           cex = 0.8)                                         # Define types of different lines
    
    fig2_x_range <- c(1, n.rep)
    plot(1:i, P[1:i], main = "P-Value", type = "l", col = "darkred",
         xlab = "date", ylab = "Value of P", xlim = fig2_x_range)
  }
  
  if(bs == TRUE){
    for(i in 1:length(deaths)){
      boots_P <- rep(0, 150)
      
      expected_value <- deaths[i]
      poisson_values <- rpois (n = 29442, lambda = expected_value)
      abnormal_poisson_loca1 <- which(poisson_values < 1)
      abnormal_poisson_loca2 <- which(poisson_values > max(t))
      poisson_values[abnormal_poisson_loca1] <- 1
      poisson_values[abnormal_poisson_loca2] <- 211
      boots_d_real <- tabulate(poisson_values)
      
      boots_P_value <- cal_P(boots_d_real, d_simulation)
      
      random_steps <- c(-4, -2, -1, 1, 2, 4)
      boots_step_list <- sample(random_steps, size = length(t0), replace = TRUE)
      boots_shuffle_location <- sample(1:length(t0), size = length(t0), replace = FALSE)
      
      #boots_ill_duration <- sample(1:80, size = length(t0), replace = TRUE,prob = prob_ill)
      boots_ill_duration <- rep(0, size = 29442)    
       round = 1
           n_index <- 1
           for (i in t){
            if(i <= 80){
               for (n in 1:deaths[n_index]){
                 boots_ill_duration[round] <- sample(1:i, size= 1, prob= prob_ill[1:i])
                 round = round+1
              }
              n_index <- 1 + n_index
            }
            else{
              for (n in 1:deaths[n_index]){
                boots_ill_duration[round] <- sample(1:80, size= 1, prob= prob_ill)
                round = round+1
              }
              n_index <- 1 + n_index
            }
           }
      print(boots_ill_duration)
      
      round <- 1
      for(k in boots_shuffle_location){
        boots_d_simulation <- d_simulation
        
        old_esdeath_date <- t0[k] + boots_ill_duration[k] 
        new_t0 <- t0[k] + boots_step_list[round]
        new_esdeath_date <- new_t0 + boots_ill_duration[k]
        
        boots_d_simulation[old_esdeath_date] <- boots_d_simulation[old_esdeath_date] - 1
        boots_d_simulation[new_esdeath_date] <- boots_d_simulation[new_esdeath_date] + 1
        print(boots_d_simulation)
        new_P_value <- cal_P(boots_d_real,boots_d_simulation)

        if(new_P_value <= boots_P_value){
          t0[k] <- new_t0
          d_simulation <- boots_d_simulation
          boots_P_value <- new_P_value
        }
        round = round + 1
      }
      boots_P[i] <- boots_P_value
      
    }
  }
  
  
  
  
  return(list(P = P, inft = inft, t0 = t0))
}

setwd("C:\\Program Files\\RStudio")
conv_data <- read.table("engcov.txt", header = TRUE, sep="")

set.seed(15)
t_location <- which(conv_data[1:150,2] > 0)
t <- conv_data[t_location,3]
deaths <- conv_data[t_location,2]

# output <- deconv(t, deaths, bs=FALSE, t0=NULL)
# t0_0 <- output$t0
# inft_0 <- output$inft
# P_0 <- output$P

output1 <- deconv(t, deaths, n.rep = 100, bs=TRUE, t0=NULL)

t0_1 <- output1$t0
inft_1 <- output1$inft
P_1 <- output1$P

end_time <- Sys.time()
execution_time <- end_time - start_time
print(execution_time)