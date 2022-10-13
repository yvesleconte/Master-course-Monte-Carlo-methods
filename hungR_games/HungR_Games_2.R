
# Build a function that takes a real matrix x of size (n,p) and a constant m as arguments to return
# a vector v in R^p such that for all k s.t. 1<k<p v_k = f(x_k) where x_k is the column k of x
# (f defined as in the annex)

leconte_f_bis <- function(y,m){
  return(m*exp((1/length(y))*sum(cumsum(y)%*%(1/(length(y):1)^2)))*sin(prod(sqrt((1+pi*y)/(1+2*pi))))^2*prod(y*exp(-y*y/2)))
}

leconte_f <- function(x,m){
  return(apply(t(x),1,leconte_f_bis,m))
}
