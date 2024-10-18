rm(list = ls())
start_time <- Sys.time()#


cal_P <- function(real_data, estimate_data){
  # P = sum_{i} ((d_i - d_i^s)^2 / max(1, d_i^s))
  # The purpose of writing this function is to make the main function look cleaner.

  # The aim of this function is to calculate the modified Pearson statistic
  # for two input values. The final of the function outputs the p-value,
  # which is used to quantify the degree of fit.
  
  tmpt1 <- (real_data - estimate_data)^2 # The square operation on the numerator.
  
  tmpt2 <- pmax(1, estimate_data) # The "max" operation on the denominator.
  # It is noticeable that using the "pmax" function can return the larger value 
  # from each pair of elements in two vectors. By avoiding the use of a "for" loop, 
  # you can significantly speed up the execution.

  tmpt3 <- tmpt1 / tmpt2

  p <- sum(tmpt3) # Sum up each element in the vector.
  
  return(p) # The return value is the function is p.
}




deconv <- function(t, deaths, n.rep=100, bs = FALSE, t0 = NULL){
  # The [deconv] gets [t: the date occurs death] and [deaths: the number of victims in "t" day]
  # as input. If we input [t0: the estimated incidence trajectory] as variable, then we use it
  # in following steps. Elif we do not input it, then function will generate one.
  # The parameter [bs] serve as an indicator that whether we will do the bootstrap to evaluate
  # the uncertainty of the model. The [n.rep] is the iteration times of the function.
  
  # The [deconv] function will have 3 output: 
  # (1) [P]: record the P_value we calculated in each iteration. 
  # (2) [inft]: record the [t0] (which may changed) after each iteration.  
  # (3) [t0]: final t0 after the last iteration.

  # Step 0:
  # *name a series of vectors and fill them with 0.*
  # *generate vector"deaths_t" using a "for" loop*

  length_max <- 310
  inft <- matrix(0, nrow = length_max, ncol = n.rep)
  # The function "matrix" generates a matrix with "length_max" rows and "n.rep" columns.
  # [inft] is a matrix, each column records the estimate number of victims in 310 days 
  # after each times of iteration.

  P <- rep(0, times = n.rep) # store each p-value for each iteration.
  deaths_t <- rep(0, times = length_max) # store number of real death data from day 1-130. 
  # The function "rep" generates a new vector by filling it with 0 repeated times and
  # parameter "times" define the length of new vector

  round <- 1 # temporary variable, indicate the index of [deaths_t] start from 1.
  for (i in t){
    # The "deaths_t[i]" represents the number of deaths on i-th day.
    # eg. "deaths_t[1]" = 1 means one person has died at day1. 

    deaths_t[i] <- deaths[round]
    # It is important to note that the real death data starting from 62 and ending at 211, so
    # we only assign value to "deaths_t" from 62-th to 211-th elements.
    # The death counts for the remaining days would be 0.

    round <- round + 1
    # Since the data in "t" and "deaths" correspond one-to-one, the i-th data taken
    # from "t" corresponds to the i-th data in "deaths".

    # The purpose of adding 1 each time is to when "i" moves to the next time point,
    # the number of deaths should move to the next data.
  }
  

  
  # Step 1*:
  # Using the known distribution and given parameters("meanlog", "sdlog"),
  # generate the probabilities for the illness duration and normalize them.

  ill_dlnorm <- dlnorm(1 : 80, meanlog = 3.152, sdlog = 0.451) 
  # Function "dlnorm" calculates the PDF of the log-normal distribution with two paraments.
  # This code generates the probabilities for the range of 1 to 80 under the log-normal
  # distribution with the specified meanlog and sdlog. Store in variable "ill_dlnorm".

  prob_ill <- ill_dlnorm / sum(ill_dlnorm)
  # Normalize the obtained probabilities so that their sum equals 1.


  # Step 2**:
  # If we haven't assigned "t0", it can be generated using the following approach.
  # Subtract the generated infection-to-death period from the actual death date to get the 
  # initial infection time.

  if (is.null(t0)){ 
    ill_duration <- rep(0,29442) # Store the time for infection-to-death in "ill_duration".
    t0 <- rep(0, 29442) # Store the initial time for infection in t0.


    ill_duration = sample(1:80, size = length(t0), prob = prob_ill, replace = TRUE)
    # Generate infection duration by sampling the data between 1 and 80 under specific prob.
    # sample(a, size = ?, pro = ?, replace = TRUE/FALSE) ----->
    # The "sample" function gets elements from vector "a" based on "prob" probabilities
    # to generates a vector of "size" size. Parameter "replace" considering whether to sample 
    # with or without replacement.

    death_date_list <- rep(t, deaths) # Store death date of all victims arranged by date.
    #rep(vec,vec) means repeatly input t[i] with deaths[i] times into "death_date_list".

    t0 <- death_date_list - ill_duration
  }
  
  # Step 3***:
  # Using the obtained "t0" ( determined by a specific step or given), estimate 
  # the death time, and fit it with the actual death time.The goodness of fit is represented
  # by the Pearson statistic(using "cal_P" function), and this process is repeated 100 times.
  # !!!If bs == TRUE, we need to replace the true death date in following codes
  boots_P <- rep(0, 150)

  new_ill_duration <- rep(0 , 29442) # Store the time for new infection-to-death.
  for (j in 1:n.rep) {
    # The outer loop runs 100 times, updating "t0" by comparing the Pearson statistic values.
    # In each iteration, we first need to shift the current "t0" by a step to the left and right,
    # and then re-estimate the infection-to-date duration to obtain the estimated death date.
    # Finally, use the "tabulate" function to determine the number of deaths per day, and compare
    # the resulting p-value. Based on this, update both "t0" and the p-value

    # Step 3.1*:
    # When we use bootstrapping, which means "bs == TRUE" as function input. We need to plot the 
    # initial time(t0) for each iteration.
    # 
    if(bs){
      poisson_values <- rep(0,length_max)
      
      for(i in t){
        expected_value <- deaths_t[i]
        poisson_values[i] <- rpois (n = 1, lambda = expected_value)
      } 
    }
        
    new_ill_duration <-  sample(1:80, size = length(t0), prob = prob_ill, replace = TRUE)
    estimate_death_date <- t0 + new_ill_duration 
    # The estimated "t0" plus the randomly generated infection time gives the estimated death date.

    d_simulation <- tabulate(estimate_death_date, nbins = length_max)
    # The number of daily death victims is stored in d_simulation.   
    # The "tabulate" function counts the times of each unique number in the input vector and stores
    # the counts in the corresponding elements and "nbins" limits the maximum range of counts.
    # eg. 5 appears three times in the vector, the 3rd element in the return value of "tabulate"
    #     will be 5. The occurences of 311 will not be returned by "tabulate".

    d_real <- deaths_t
    P_value <- cal_P(d_real, d_simulation)
    # Use the external function "cal_P" to calculate the Pearson statistic.

    if (j <= 50){
      random_steps <- c(-8, -4, -2, -1, 1, 2, 4, 8)
    }
    else if (j <= 75){
      random_steps <- c(-4, -2, -1, 1, 2, 4)
    }
    else{
      random_steps <- c(-2, -1, 1, 2)
    }
    # Define different step-vectors based on the number of iterations. The reason is that as the
    # number of iterations increases, the fit of the expected death time improves, and the required
    # step size for adjustments will decrease accordingly.

    steps_list <- sample(random_steps, size = length(t0), replace = TRUE)
    # Store the random step sizes in "step_list". It's noticeable that with different iteration
    # counts, the vector sampled by "sample" function will also vary.
    
    round <- 1 # temporary variable, indicate the index of [deaths_t] start from 1.
    for (i in sample(1:length(t0), size = length(t0), replace = FALSE)){
          
      tmpt_d_simulation <- d_simulation
      tmpt_t0 <- t0[i] + steps_list[round]
      tmpt_old_esdeath_date <- estimate_death_date[i]
          

      tmpt_new_esdeath_date <- tmpt_t0 + new_ill_duration[i]
          
      if(tmpt_old_esdeath_date > 310){
        tmpt_old_esdeath_date <- 310
        t0[i] <- t0[i] - (tmpt_old_esdeath_date - 310)
      }
          
      if (tmpt_new_esdeath_date < 1){
        next
      }
      if (tmpt_new_esdeath_date > 310){
        next
      }
          
          

          
      tmpt_d_simulation[tmpt_old_esdeath_date] <- tmpt_d_simulation[tmpt_old_esdeath_date] - 1
      tmpt_d_simulation[tmpt_new_esdeath_date] <- tmpt_d_simulation[tmpt_new_esdeath_date] + 1
          
          
          
      tmpt_P_value <- cal_P(d_real, tmpt_d_simulation)
          
      if (tmpt_P_value < P_value){
        P_value <- tmpt_P_value
        t0[i] <- tmpt_t0
        d_simulation <- tmpt_d_simulation
      }
      round <- round + 1
    }
        
    P[j] <- P_value
    inft[1:length_max, j] <- tabulate(t0, nbins = length_max)
        
        
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
             xlab = "date", ylab = "Value of P", xlim = fig2_x_range, ylim = c(1,500))
  }
  return(list(P = P, inft = inft, t0 = t0, boots <- boots_P))
     }



setwd("C:\\Program Files\\RStudio")
conv_data <- read.table("engcov.txt", header = TRUE, sep="")



set.seed(8535)

conv_data <- conv_data[1:150,]

t <- conv_data$julian
deaths <- conv_data$nhs

output <- deconv(t, deaths, n.rep = 100, bs = FALSE, t0 = NULL)
t0_0 <- output$t0
inft_0 <- output$inft
P_0 <- output$P

output1 <- deconv(t, deaths, n.rep =  100, bs = TRUE, t0 = t0_0)
boots <- output1$boots
t0_1 <- output1$t0
inft_1 <- output1$inft
P_1 <- output1$P

end_time <- Sys.time()
execution_time <- end_time - start_time
print(execution_time)
