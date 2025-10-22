q_scad <- function(lambda,par){
  b <- 3.7
  term1 <- 1
  term2 <- ((b*lambda-par)*as.numeric((b*lambda-par)>0))/((b-1)*lambda)
  out <- lambda*(term1*as.numeric(par <= lambda)+term2*as.numeric(par > lambda))
  return(out)
}
