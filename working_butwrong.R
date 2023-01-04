x1 = 0.4365980
p1 = 0.0329664
x2 = 17.9241816
p2 = 1 - p1

a = 0.05884
b = 4.298
c = 21.8

w1 = p1
w2 = p2

x = c(x1, x2, p1)

cost_function = function(x){
  x1 = x[1]
  x2 = x[2]
  p1 = x[3]
  p2 = 1 - p1
  a = 0.05884
  b = 4.298
  c = 21.8
  value = objfunction(x1, x2, 
                      p1, p2,
                      a, b, c)
  return(value)
}

objfunction = function(x1, x2, 
                       w1, w2,
                       a, b, c,
                       iter){
  A = information_m(x1, x2, 
                    w1, w2, 
                    a, b, c)
  cvec = c(-c/(a^2), c/(b^2), (1/a) - (1/b))
  # value = -t(cvec) %*% (A + (10^(-6) * diag(3)))
  # value = pracma::mrdivide(-(cvec), (A + (10^(-6) * diag(3))))
  value =  (-(cvec)) %*%  solve(A + (10^(-6) * diag(3))) 
  
  con1 = x[1] > x[2]
  con2 = p[1] < 0
  return(value - 1e100 * (con1 + con2))
}

information_m = function(x1, x2, 
                         p1, p2,
                         a, b, c){
  x = c(x1, x2)
  p = c(p1, p2)
  
  # sapply(1:2, function(x))
  test = c()
  for(i in 1:2){
    FFr = sqrt(p[i]) * c(x[i] * c * exp(-a*x[i]), 
            -x[i] * c * exp(-b * x[i]),
            exp(-b  * x[i]) - exp(-a * x[i]))
    test = rbind(test, FFr)
  }
  A = t(test) %*% test
  return(A)
}














