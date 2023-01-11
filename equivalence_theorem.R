#Equivalence theorem

x1 = 0.4365980
p1 = 0.0329664
x2 = 17.9241816
p2 = 1 - p1

x1 = 0.179288
p1 = 0.60616
x2 = 3.56583
p2 = 1 - p1

a = 0.05884
b = 4.298
c = 21.8

w1 = p1
w2 = p2

lb = 0
ub = 20
maxf = aaa$optimumValue
maxdiff = -1000
N = 1000
mb = lb + (ub - lb) / 4
h1 = (ub - mb) / N
h2 = (mb - lb) / (2 * N)
x = seq(from = 0, to = 20, length.out = 1000)
differ = rep(10000, 1000)
objfunction_c1 = function(x,
                          a, b, c,
                          IA){
  fx = c(x*c*exp(-a*x), -x * c * exp(-b*x), exp(-b*x)-exp(-a*x))
  e = a-b
  f = log(a/b)
  cvector = c((e/a-f)/(e^2), (f - e/b)/(e^2), 0)
  # value =  (-(cvec)) %*%  solve(A + (10^(-6) * diag(3)))  %*% cvec
  value = (t(fx) %*% IA %*% cvector)^2 - t(cvector) %*% IA %*% cvector
  return(value)
}

for(i in 1: (N + 2*N)){
  if(i <= (2*N)){
    x[i] = lb + h2*i
  } else{
    x[i] = mb + h1 * (i - 2*N)
  }
  
  A = information_m(x1, x2, 
                    p1, p2, 
                    a, b, c)
  IA = solve(A + (10^(-6) * diag(3)))
  differ[i] = objfunction_c1(x[i], 
                             a,b,c,
                             IA)
  if(differ[i] > maxdiff){
    maxdiff = differ[i]
  }
}

maxdiff = abs(maxdiff)
eff = 1.0 - maxdiff / maxf
eff


plot(x, differ)




