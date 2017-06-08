#Fundamental Simulation Function
##Generate p-value function
compute_pvalue = function(Yk, sigk, n, K){
  nk = rep(n, K) # assume evenly distributed number of subgroups nk
  # initialize the parameters
  tao = NULL # causal estimand
  Yob = NULL
  sig = NULL
  V = NULL # neyman estimators of variance
  t = NULL # t-statiscs
  # p-value matrix
  p = matrix(0, ncol = K, nrow = n)
  # gen YK obs 
  for (j in 1:n){
    for (i in 1:K){ 
      # derive the sample mean and sample variance 
      Yob[i] = rnorm(1, Yk[i], sigk[i] / sqrt(n))
      sig[i] = sigk[i] * rchisq(1, nk[i] - 1) / (nk[i] - 1)
      Y0 = Yob[1] # non-treat mean
      tao[i] = Yob[i] - Y0 
      V[i] = sig[i]/nk[i] + sig[1]/nk[1]
      # calculate the t-statistics
      t[i] = (tao[i] - 0) / sqrt(V[i])
      # get the two-sided p-value
      p[j,i] = (1 - pt(abs(t[i]), nk[i] - 2)) * 2
    }
  }
  return(p)
}


##Generate p-combine function
compute_pcombine <- function(Yk, sigk, n, K){
  p_comb_val <- NULL
  pmatrix <- compute_pvalue(Yk, sigk, n, K)
  for (k in 1:n){
  p_comb_val[k] = min(pmatrix[k,-1]) #minimum
  }
  return (p_comb_val)
}


##Generate power function
compute_power <- function(Yk, sigk, n, K, p_cv, pcombine_func = compute_pcombine){
  
  power = 0
  pcomb = pcombine_func(Yk, sigk, n, K)
  rejec = 0
  rec = 0
  
  for (i in 1:n){
    if(pcomb[i] <= p_cv){
      rejec = rejec + 1
    }else{
      rec = rec + 1
    }
  }
  power = rejec / (rejec + rec)
}

##80% Power critical ratio calculating function
compute_ctratio = function(Y, n, K, r1_start = 0.01, r1_step =0.001, r1_ceil = 1, rk = 1){
# set the starting critical value
r1_critical = r1_start

power = ratio_power(Y, r1_critical, rk, n, K, compute_pcombine, compute_power)
  if (r1_start > r1_ceil){
    print("Error: critical value does not exist")
  }else{
    r1_start = r1_start + r1_step
    power_ahead = ratio_power(Y, r1_start, rk, n, K, compute_pcombine, compute_power)
    if (power < 0.8 | power >= power_ahead){
      return(compute_ctratio(Y, n, K, r1_start = r1_start, r1_step = r1_step, r1_ceil = r1_ceil, rk = rk))
    }else{
    return (r1_critical)  
    } 
  }  
}

##Binary critical ratio calculation
compute_binary_ctratio = function(Y, n, K, delta, r1_start = 0, r1_step =0.1, r1_ceil = 1, rk = 1){
  critical = (r1_start + r1_ceil)/2

  power = ratio_power(Y, critical, rk, n, K, compute_pcombine, compute_power)
  
  # check for the monotonicity
  critical_next = critical + r1_step
  
  power_ahead = ratio_power(Y, critical_next, rk, n, K, compute_pcombine, compute_power)
  diff = abs(power - 0.8)
  while (diff > delta) {
    print(diff)
    print(power)
    print(critical)
   if (power > power_ahead | power < 0.8 ){
     return(compute_binary_ctratio(Y, n, K, delta = delta, r1_start = critical, r1_step = r1_step, r1_ceil = r1_ceil, rk = rk))
   }else{
     if (power > 0.8){
     return(compute_binary_ctratio(Y, n, K, delta = delta, r1_start = r1_start, r1_step = r1_step, r1_ceil = critical, rk = rk))
   }else{
     print("Error: does not exists")
   }   
   }
  }  
  return(critical)
}

