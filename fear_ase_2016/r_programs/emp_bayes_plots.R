setwd("/home/agerken/mclab/cegs_ase_paper/r_data")

emp_bayes_M= read.csv('emp_bayes_to_compare_M_20150218.csv', head=TRUE)
attach(emp_bayes_M)

#r101 only
r101_M <- emp_bayes_M[ which(line=='r101'),]
r101_M <- within(r101_M, line_over_total <- Line_sum / Total_count)
r101_M <- within(r101_M, tester_over_total <- Tester_sum / Total_count)
attach(r101_M)

#plot summary stuff
plot(q5_mean_theta, Line_over_tester)
plot(q5_mean_theta, Tester_over_line)
plot(q5_mean_theta, Line_Average)
plot(q5_mean_theta, Tester_Average)
plot(q5_mean_theta, Line_sum)
plot(q5_mean_theta, Tester_sum)
plot(q5_mean_theta, Total_count)

#test directionality of results
plot(q5_mean_theta, line_over_total)
plot(q5_mean_theta, tester_over_total)

#correlations of q4 and q5 thetas
plot(q4_mean_theta, q5_mean_theta)

#density plots
q4 <- density(q4_mean_theta) # returns the density data
plot(q4) # plots the results 
q5 = density(q5_mean_theta)
plot(q5)
q6 = density(q6_mean_theta)
plot(q6)

#set up multiple density plot function
plot.multi.dens <- function(s)
{
  junk.x = NULL
  junk.y = NULL
  for(i in 1:length(s))
  {
    junk.x = c(junk.x, density(s[[i]])$x)
    junk.y = c(junk.y, density(s[[i]])$y)
  }
  xr <- range(junk.x)
  yr <- range(junk.y)
  plot(density(s[[1]]), xlim = xr, ylim = yr, main = "")
  for(i in 1:length(s))
  {
    lines(density(s[[i]]), xlim = xr, ylim = yr, col = i, lty=i)
  }
}

#plot multiple density plots
x <- q4_mean_theta
y <- q5_mean_theta
z <- q6_mean_theta
plot.multi.dens(list(x,y,z))
legend("topleft", legend=c("q4", "q5", "q6"), lty=c(1,2,3), col=c('black', 'red', 'green'), bty='n', cex=.75)
abline(v=0.5)


