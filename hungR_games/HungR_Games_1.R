library(data.table)

# Using a Gaussian generator, build a function to return a simulation of a discrete Brownian motion
# for a given vector t.

leconte_mb <-function(t){
  return(cumsum(sqrt(t-shift(t, fill=0))*rnorm(length(t))))
}