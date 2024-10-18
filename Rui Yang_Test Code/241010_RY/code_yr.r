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
  
  inft_out <- matrix(rep(0, length_of_year * n.rep), nrow = length_of_year, ncol = n.rep)
  P_out <- rep(0, n.rep)

  for (i in 1:n.rep){
    
    last_day_of_death <- max(t)

    death_duration <- sample(1:80, size = 29442, replace = TRUE, prob = duration_80_density_normalise)
    estimate_death_date <- t0 + death_duration
    
    # indices_outrip <- which(estimate_death_date > last_day_of_death)
    # estimate_death_date[indices_outrip] <- last_day_of_death   # No one can death before get infect
    # gap <- last_day_of_death * rep(1, times = length(indices_outrip)) - (t0[indices_outrip] + death_duration[death_duration])
    # t0[indices_outrip] <- t0[indices_outrip] - gap

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

      change_step <- c(-2, -1, 1, 2)

      estimate_pois <- rep(0, length_of_year)

      round <- 0

      for (k in t){
        round <- round + 1
        if (deaths[round] != 0){
          estimate_pois[k] <- rpois(1, lambda = deaths[round])
        }

      }

      deaths_t <- estimate_pois

    }
    

    p0 <- cal_P(estimate_death, deaths_t)
    # print(p0)

    shuffle <- sample(1:length(t0), length(t0), replace = FALSE)
    step <- sample(change_step, length(t0), replace = TRUE)
    
    round <- 1
    

    t0_tmp <- t0

    
    for (j in shuffle){
      
      estimate_death_tmp <- estimate_death
      
      old_date <- t0_tmp[j] + death_duration[j]
      
      if (old_date >= length_of_year){
        old_date <- length_of_year

        t0_tmp[j] <- length_of_year - death_duration[j]

        
      }
      
      estimate_death_tmp[old_date] <- estimate_death_tmp[old_date] - 1
      
      new_date <- old_date + step[round]
      
      if (new_date >= length_of_year){
        next
      }
      
      if (new_date <= 1){
        next
      }
      
      estimate_death_tmp[new_date] <- estimate_death_tmp[new_date] + 1

      pi <- cal_P(estimate_death_tmp, deaths_t)

      if (pi < p0){
        p0 <- pi
        

        t0_tmp[j] <- t0_tmp[j] + step[round]

        
        estimate_death <- estimate_death_tmp
      }
      round <- round + 1
    }
    
    if(bs != TRUE){
      t0 <- t0_tmp
    }
    
    
    # print(p0)
    
    P_out[i] <- p0
    inft_out[,i] <- tabulate(t0_tmp, nbins = length_of_year)
    
    
    if (bs){
      
      layout(matrix(c(1, 2), 1, 2))

      fig1_y_range <- c(1, 1800)
      plot(1:310, estimate_death, main = "[bs]: Quantify model Uncertainty.", type = "l", col = "blue",
           xlab = "date", ylab = "Number of Victims", ylim = fig1_y_range)
      lines(1:310, deaths_t,type = "l", col = "darkgreen")
      Estimate_incidence <- tabulate(t0, nbins = length_of_year)
      lines(1:310, Estimate_incidence,type = "l", col = "darkorange")

      legend("topright",                                      # 图例位置
             legend = c("Estimate Death", "Simulate Death", "Convergent Incidence"),  # 图例文本
             col = c("blue", "darkgreen", "darkorange"),      # 颜色
             lty = 1)                                         # 线型

      # Fig 2.
      fig2_x_range <- c(1, n.rep)
      fig2_y_range <- c(0, 1000)
      plot(1:i, P_out[1:i], main = "[bs]: P-Value", type = "l", col = "darkred",
           xlab = "date", ylab = "Value of P", xlim = fig2_x_range, ylim = fig2_y_range)
      
    } else {
      
      ## plot to observe convergence situation
      
      layout(matrix(c(1, 2), 1, 2))
      
      # Fig 1.
      fig1_y_range <- c(1, 1800)
      plot(1:310, estimate_death, main = "[Training]: Convergence Situation", type = "l", col = "blue", 
           xlab = "date", ylab = "Number of Victims", ylim = fig1_y_range)
      lines(1:310, deaths_t,type = "l", col = "darkgreen")
      Estimate_incidence <- tabulate(t0, nbins = length_of_year)
      lines(1:310, Estimate_incidence,type = "l", col = "darkorange")
      
      legend("topright",                                      # 图例位置
             legend = c("Estimate Death", "Real Death", "Estimate Incidence"),  # 图例文本
             col = c("blue", "darkgreen", "darkorange"),      # 颜色
             lty = 1)                                         # 线型
      
      # Fig 2.
      fig2_x_range <- c(1, n.rep)
      plot(1:i, P_out[1:i], main = "[Training]: P-Value", type = "l", col = "darkred",
           xlab = "date", ylab = "Value of P", xlim = fig2_x_range)

    }

    
  }
  return(list(t0 = t0, inft = inft_out, P = P_out))
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

n.rep <- 100

output <- deconv(t, deaths, n.rep = 100, bs=FALSE, t0=NULL)

t0_0 <- output$t0
inft_0 <- output$inft
P_0 <- output$P

output1 <- deconv(t, deaths, n.rep = 100, bs=TRUE, t0=t0_0)

t0_1 <- output1$t0
inft_1 <- output1$inft
P_1 <- output1$P



layout(matrix(c(1, 1), 1, 1))

color_0 <- adjustcolor("darkorange", alpha.f = 0.01)

plot(inft_1[,1], type = "l",lty = 2, col = color_0, ylim = range(inft_1), 
     xlab = "date", ylab = "Number of Voctims", main = "Model Result and Its Uncertainty")
for(i in 2:100) {
  lines(inft_1[,i], lty = 2, col = color_0)
}

color_1 <- adjustcolor("darkorange", alpha.f = 1)
lines(inft_0[, ncol(inft_0)], lty = 1, col = color_1)


length_of_year <- 310

deaths_t <- rep(0, times = length_of_year)
round <- 1
for (i in t){
  deaths_t[i] <- deaths[round]
  round <- round + 1
}

color_2 <- adjustcolor("darkgreen", alpha.f = 1)

lines(1:310, deaths_t,type = "l", lty = 1, col = color_2)

color_3 <- adjustcolor("red", alpha.f = 0.7)

abline(v = 84, col = color_3, lwd = 2, lty = 1)


legend("topright",
       legend = c("Convergent Incidence", "Real Death Data", "Uncertainty interval", "UK Lockdown Date"),
       col = c(color_1, color_2, adjustcolor("darkorange", alpha.f = 0.3), color_3),
       lty = 1)



end_time <- Sys.time()
execution_time <- end_time - start_time
print(execution_time)





