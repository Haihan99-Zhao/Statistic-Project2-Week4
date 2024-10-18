# Our Group Member are：
# Name: Rui Yang    | Student ID: S2747080
# Name: Haihan Zhao | Student ID: S2668314
# Name: Di Wu       | Student ID: S2636080

# The workflow of our group is: 
# (1) Conduct multiple group discussions to clarify the assignment's work content
#     and specific details. ->
# (2) Each member independently completes all parts of the assignment (to ensure
#     that every member gets to practice with R programming). ->
# (3) Multiple discussions on the code details, and finalize the submission version.

# The GitHub URL of our repository is:
# 


start_time <- Sys.time() # Set the starting point for timing.

#setwd("C:\\Program Files\\RStudio") # Set up working environment.
setwd(here::here())



# The purpose of writing this [cal_P] is to make the main function look cleaner.

# The aim of this function is to calculate the modified Pearson statistic
# for two input values. The final of the function outputs the p-value,
# which is used to quantify the degree of fit.

cal_P <- function(real_data, estimate_data){
  # P = sum_{i} ((d_i - d_i^s)^2 / max(1, d_i^s))
  
  tmpt1 <- (real_data - estimate_data)^2 # The square operation on the numerator.
  
  tmpt2 <- pmax(1, estimate_data) # The "max" operation on the denominator.
  # It is noticeable that using the "pmax" function can return the larger value 
  # from each pair of elements in two vectors. By avoiding the use of a "for" loop, 
  # you can significantly speed up the execution.

  tmpt3 <- tmpt1 / tmpt2

  p <- sum(tmpt3) # Sum up each element in the vector.
  
  return(p) # The return value is the function is p.
}


# The [deconv] gets [t]: "the date occurs death" and [deaths]: "the number of 
# victims in "t[i]" day" as input. 

# If we input [t0]: "the estimated incidence trajectory" as variable, then we use
# it in following steps. Else if we do not input it, then function will generate one.

# The parameter [bs] serve as an indicator that whether we will do the bootstrap 
# to evaluate the uncertainty of the model.

# The [n.rep] is the iteration times of the function.

# The [deconv] function will have 3 output: 
# (1) [P]: record the P_value we calculated in each iteration. 
# (2) [inft]: record the "incidence-date" (which may changed) after each iteration.  
# (3) [t0]: final incidence date of each individual after the last iteration.

deconv <- function(t, deaths, n.rep=100, bs = FALSE, t0 = NULL){

  # Step 0:
  # *difine a series of vectors and fill them with 0.*
  # *generate vector"deaths_t" using a "for" loop*

  length_max <- 310
  inft_output <- matrix(0, nrow = length_max, ncol = n.rep)
  # The function "matrix" generates a matrix with "length_max" rows and "n.rep" columns.
  # [inft_output] is a matrix, each column records the estimate number of victims  
  # in 310 days after each times of iteration.

  P_output <- rep(0, times = n.rep) # store each p-value for each iteration.
  deaths_t <- rep(0, times = length_max) # store number of real death data from day 1-310. 
  # The function "rep" generates a new vector by filling it with 0 repeated times and
  # parameter "times" define the length of new vector
  
  deaths_t[t]<- deaths
  # The "deaths_t[t]" represents the number of deaths on i-th day(i for all
  # elements in vector t)
  

  
  # Step 1*:
  # *Using the known distribution and given parameters("meanlog", "sdlog"),
  # generate the probabilities for the illness duration and normalize them.*

  ill_dlnorm <- dlnorm(1 : 80, meanlog = 3.152, sdlog = 0.451) 
  # Function "dlnorm" calculates the PDF of the log-normal distribution with two paraments.
  # This code generates the probabilities for the range of 1 to 80 under the log-normal
  # distribution with the specified meanlog and sdlog. Store in variable "ill_dlnorm".

  prob_ill <- ill_dlnorm / sum(ill_dlnorm)
  # Normalize the obtained probabilities so that their sum equals 1.


  # Step 2**:
  # *If we haven't assigned "t0", it can be generated using the following approach:
  # Subtract the generated infection-to-death period from the actual death date to get the 
  # initial infection time.*

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
  # *Using the obtained "t0" (determined by a specific step or given), estimate 
  # the death date, and fit it with the actual death date. The goodness of fit is
  # represented by the Pearson statistic(using "cal_P" function), and this process
  # is repeated 100 times.*

  # **If using bootstrapping, at this step, we will use the actual daily death victims
  # data as the parameter for the Poisson simulation. A single number will be  
  # generated in the place of real daily number of victims, which will then enter
  # the next loop.

  new_ill_duration <- rep(0 , 29442) # Store the time for new infection-to-death.
  for (j in 1:n.rep) {
    # The outer loop runs 100 times, updating "t0" by comparing the Pearson statistic values.
    # In each iteration, we first need to shift the current "t0" by a step to the left and right,
    # and then re-estimate the infection-to-date duration to obtain the estimated death date.
    # Finally, use the "tabulate" function to determine the number of deaths per day, and compare
    # the resulting p-value. Based on this, update both "t0" and the p-value
        
    new_ill_duration <-  sample(1:80, size = length(t0), prob = prob_ill, replace = TRUE)
    estimate_death_date <- t0 + new_ill_duration 
    # The estimated "t0" plus the randomly generated infection time gives the first 
    # prediction of death date.

    d_simulation <- tabulate(estimate_death_date, nbins = length_max)
    # The number of daily death victims is stored in d_simulation.   
    # The "tabulate" function counts the times of each unique number in the input
    # vector and stores the counts in the corresponding elements and "nbins" limits
    # the maximum range of counts.
    
    # eg. 5 appears 3 times in the vector, the 5th element in the return value of 
    # "tabulate" will be 3. The occurrence above 311 will not be returned by "tabulate",
    # because the "nbins" parameter confine the length of output.

    if (j <= 50){
      random_steps <- c(-8, -4, -2, -1, 1, 2, 4, 8)
    }
    else if (j <= 75){
      random_steps <- c(-4, -2, -1, 1, 2, 4)
    }
    else{
      random_steps <- c(-2, -1, 1, 2)
    }
    # Define different step-vectors based on the number of iterations. 
    # The reason is that as the number of iterations increases, we are getting
    # closer and closer to convergent state, and we need to make our movement
    # smaller to achieve accurate adjustment.

    d_real <- deaths_t # "d_real" means the real number of victims for each day. 
    # The length is 310. "deaths_t" stores number of real death data from day 1-130.
    #  If "bs == FALSE", we use the real data to calculate the Pearson statistic.
    
    

    # Step 3.1*:
    # When we use bootstrapping, which means "bs == TRUE" as function input, we 
    # need to replace the actual death data. The method is to use the actual death 
    #  count for each day as the expected value for the Poisson simulation.

    # if(bs) == if(bs == TRUE)
    if(bs){
      poisson_values <- rep(0,length_max)
      round<- 1 # temporary variable, indicate the index of [deaths] start from 1.

      for(i in t){  
        expected_value <- deaths[round] # real data as expected value
        poisson_values[i] <- rpois (n = 1, lambda = expected_value)
        # "t" and "deaths" correspond to recording the time and the number of deaths,
        #  respectively.
        
        # "for(i in t)": "i" will automatically start from the first element, and
        #  it corresponds to the first element of "deaths" as well.
        round <- round + 1 # implement sequential input.
      }
      d_real <- poisson_values # replace the real death date.
    }
    
    
    P_value <- cal_P(d_real, d_simulation)
    # Use the defined function "cal_P" to calculate the Pearson statistic.
    
    steps_list <- sample(random_steps, size = length(t0), replace = TRUE)
    # Store the random step sizes in "step_list". It's noticeable that with different
    # iteration counts, the vector sampled by "sample" function will also vary.

    original_t0 <- t0
    # The following code will change "t0", so store the "original t0" beforehand.
    round <- 1 # Temporary variable, indicate the index of [steps_list] start from 1.
    shuffle <- sample(1:length(t0), size = length(t0), replace = FALSE)
    # To avoid sequential dependence, process "t0" in a random order. 
    # "shuffle" refers to a random sequence of index.



    # Step 3.2*:
    # In this part, we need to recalculate the estimated death date with different
    # steps. Then, we update the vector of deaths per day. Finally, by comparing
    # the results, we update the "t0", the vector of deaths per day and P_value 
    # that produce the smaller Pearson statistic.
    for (i in shuffle){
          
      tmpt_d_simulation <- d_simulation 
      # Temporary variable, decide whether to keep it, based on Pearson statistic.
      
      tmpt_t0 <- t0[i] + steps_list[round] # perform the step adjustment.
      tmpt_old_esdeath_date <- estimate_death_date[i]
      tmpt_new_esdeath_date <- tmpt_t0 + new_ill_duration[i]

      # "tmpt_d_simulation" is the vector stored deaths per day, After moving step,
      #  the victim's death date is updated from old date to a new one. So the count
      #  in old date's position should be decreased by one, and new date's position
      #  should be increased by one.
      tmpt_d_simulation[tmpt_old_esdeath_date] <- tmpt_d_simulation[tmpt_old_esdeath_date] - 1
    
      # In some cases, the change in the Pearson statistic can be directly assessed.
      if (tmpt_new_esdeath_date < 1){
        next
        # If "new death date < 1", which affect outcome larger in Pearson statistic.
        # This is beacause we don't have real death data before day 1.  
        # And to ensure no errors occur. Discard in advance and jump into next loop.
      }
      if (tmpt_new_esdeath_date > 310){
        next
        # If "new death date > 310", which affect outcome larger in Pearson statistic.
        # This is beacause we don't have real death data after day 310.  
        # And to ensure no errors occur. Discard in advance and jump into next loop.
      }
      
      tmpt_d_simulation[tmpt_new_esdeath_date] <- tmpt_d_simulation[tmpt_new_esdeath_date] + 1
          
      tmpt_P_value <- cal_P(d_real, tmpt_d_simulation) # the Pearson statistic after adjusting
          
      if (tmpt_P_value < P_value){
        P_value <- tmpt_P_value
        t0[i] <- tmpt_t0
        d_simulation <- tmpt_d_simulation
        # Update the fit with a better Pearson statistic, initial infected time and deaths per day.
      }
      round <- round + 1 # Going next index of "steps_list"
    }
    
        
    P_output[j] <- P_value # After fully updating "t0" each time, store the Pearson statistic.  
    inft_output[1:length_max, j] <- tabulate(t0, nbins = length_max)
    # After fully updating "t0" each time, store the current estimated incidence



    # Step 4****:plot the graphs
    # *It is worth noting that this plotting occurs within the loop,, so each time 
    #  a plot is generated, the changes can be clearly observed.*
    if (bs){
      layout(matrix(c(1, 2), 1, 2)) # Generate a 1x2 layout for displaying plots.

      # Fig1 showing estimated deaths and new infection each day;simulated deaths
      #  after bootstrapping

      fig1_y_range <- c(1, 1800) # Set up a limits vector
      plot(1:length_max, d_simulation, main = "[bs]: Quantify model Uncertainty.n", type = "l",  
           col = "blue", xlab = "date", ylab = "Number of Victims", ylim = fig1_y_range)
      # "plot" function will draw a graph on the first place in created layout:
      # The 1st vector input defines the range of the x-axis.
      # The 2nd vector input defines the data points for each points on the x-axis.
      # In this grap, we draw a line showing current simulated deaths per day.
      # The 3rd "main" defines the title of this graph.
      # The 4th "type" defines lines types, and "l" means *Line chart.
      # The 5th "col" defines the color of line.
      # The "xlab" and "ylab" label the two axes.
      # The "ylim" means the range displayed on the y-axis

      lines(1:length_max, d_real,type = "l", col = "darkgreen")
      # Draw a new line on the existing graph. This line is showing the real deaths per day.

      estimate_incidence <- tabulate(t0, nbins = length_max)
      lines(1:length_max, estimate_incidence,type = "l", col = "darkorange")
      # This line is showing the estimated incidence per day.
      
      legend("topright",legend = c("Estimate Death", "simulate Death", "Convergent Incidence with subtle change "),  
              col = c("blue", "darkgreen", "darkorange"), lty = 1, cex = 0.8)                                 
      # "legend" function adds a legend to the existing  graph.
      # The 1st string input defines legend position. "col" defines color. 
      # "lty = 1" defines the type of lines showing is solid line.
      # "cex = 0.8" defines reducing the legend size by 80%.

      # Fig2 showing P-value after bootstrapping, only slight fluctuations during iteration.
      fig2_x_range <- c(1, n.rep)
      fig2_y_range <- c(0, 1000)
      plot(1:j, P_output[1:j], main = "[bs]: P-Value", type = "l", col = "darkred",
           xlab = "date", ylab = "Value of P", xlim = fig2_x_range, ylim = fig2_y_range)
      }
    else{     
      layout(matrix(c(1, 2), 1, 2)) # Create a layout to plot.
      
      # Fig1
      fig1_y_range <- c(1, 1800) 
      plot(1:length_max, d_simulation, main = "[Training]: Convergence Situation", type = "l",  
           col = "blue", xlab = "date", ylab = "Number of Victims", ylim = fig1_y_range)
      # This line is showing estimated death per day.

      lines(1:length_max, d_real,type = "l", col = "darkgreen")
      # This line is showing the real number of deaths per day.

      
      estimate_incidence <- tabulate(t0, nbins = length_max)
      lines(1:length_max, estimate_incidence,type = "l", col = "darkorange")
      # This line is showing the estimated incidence per day.

      legend("topright", legend = c("Estimate Death", "Real Death", "Estimate Incidence"),                                      
              col = c("blue", "darkgreen", "darkorange"), lty = 1, cex = 0.8)
          
    
      # Fig2 shows the changes of Pearson statistic and focus on the convergence state.
      fig2_x_range <- c(1, n.rep)
      plot(1:j, P_output[1:j], main = "[Training]: P-Value", type = "l", col = "darkred",
          xlab = "date", ylab = "Value of P", xlim = fig2_x_range)
      }

    if(bs == TRUE) {t0 <- original_t0}
    # When we use bootstrapping, the "t0" of the inner loop would be updated, 
    # but it doesn't go the next loop. In the other words, in the loop of
    # n.rep = 100, the initial "t0" is always the "t0" in the convergence state.
    

  }
  return(list(P = P_output, inft = inft_output, t0 = t0, d_simulation = d_simulation))
  # "return" function returns a range of variable names in the list, 
  #  which are stored in the defined function as a list form.
}


#********************** MAIN ****************************


# Step 5*****:
# *Read the data in a table, and utilize hte function.

conv_data <- read.table("engcov.txt", header = TRUE, sep="")
# "read.table" reads the file in the form of a table
# "header =TRUE" means the first column in file is names.
# "sep = "" " means read data using "sep"

set.seed(70) # Set up a random seed
conv_data <- conv_data[1:150,]

t <- conv_data$julian # Read the "julian" vector from the "conv_data" matrix.
deaths <- conv_data$nhs # Read the "deaths" vector from the "conv_data" matrix.

output <- deconv(t, deaths, n.rep = 100, bs = FALSE, t0 = NULL)

t0_0 <- output$t0 # Store the result of "t0" returned from the function into "t0_0" variable.
inft_0 <- output$inft # Store the matrix of inft outcome
P_0 <- output$P # store the vector of Pearson statistic.

output1 <- deconv(t, deaths, n.rep =  100, bs = TRUE, t0 = t0_0)
# Use bootstrapping and set the t0 at converged state.

t0_1 <- output1$t0
inft_1 <- output1$inft
P_1 <- output1$P



# Step 6: 
# *Visualize our result*

layout(1)
# Create a drawing board with size of 1. (We just need one figure)

# [inft] record the changeable case of convergent Incidence date.
# It records 100 times of re-sample, so we plot them 100 times and 
# use the density of lines to represent the model's uncertainty.

# Since they are very "dense", we Choose "black" with a transparency of 0.1
# to show them. 
color_0 <- adjustcolor("black", alpha.f = 0.1)
# "Plot" one column of data and followed by "lines" with other data can make
#  all the lines shown in one single figure.

plot(inft_1[,1], type = "l",lty = 2, col = color_0, ylim = range(inft_1), 
     xlab = "date", ylab = "Number of Voctims", main = "Model Result and Its Uncertainty")
for(i in 2:100) {
  lines(inft_1[,i], lty = 2, col = color_0)
}

# We choose an orange line to represent our convergent incidence, which record 
# in the last column of [inft_0] (which is the one of the output after we
# run the [deconv] for the first time)

color_1 <- adjustcolor("darkorange", alpha.f = 1)
lines(inft_0[, ncol(inft_0)], lty = 1, col = color_1)

# The same with the operation in [deconv], we combined the
# [julian] and [nhs] into one vector [deaths_t]. 
# deaths_t[n] = m, means there are m victims in the n-th day . 

length_max <- 310

deaths_t <- rep(0, times = length_max)
round <- 1
for (i in t){
  deaths_t[i] <- deaths[round]
  round <- round + 1
}

# the [deaths_t] represents the real death data,
# We choose the darkgreen line to show it.

color_2 <- adjustcolor("darkgreen", alpha.f = 1)
lines(1:310, deaths_t,type = "l", lwd = 2, lty = 2, col = color_2)

# day 84 is the first day of UK lockdown, we draw a Vertical
# line with red color to show it. Combine all lines into evaluate
# will show some policy implications.

color_3 <- adjustcolor("red", alpha.f = 0.7)
abline(v = 84, col = color_3, lwd = 2, lty = 1)

# As what we have done in [deconv], we plot the death-trajectory estimated by 
# convergent model (with blue line). Combine it with green line (real 
# death-trajectory) to show the fit of our model.

color_4 <- adjustcolor("blue", alpha.f = 0.6)

# The convergent estimated death stored in [d_simulation] label in first output.
d_simulation_0 <- output$d_simulation

lines(1:310, d_simulation_0, type = "l", lty = 2, col = color_4)



legend("topright",
       legend = c("Convergent Incidence", "Real Death Trajectory", "Uncertainty interval",
                  "UK Lockdown Date", "Estimated Death Trajectory "),
       col = c(color_1, color_2, adjustcolor("black", alpha.f = 0.3), color_3, color_4),
       lty=c(1, 2, 1, 1, 2))





# Obtain the time now, subtract it with start time
# to show the Program run time

end_time <- Sys.time()
execution_time <- end_time - start_time
print(execution_time)
