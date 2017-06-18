#Fundamental Simulation Function
##Generate p-value function
compute_pvalue = function(Yk, sigk, N, K){
  n <- N/K
  nk = rep(n, K) # assume evenly distributed number of subgroups nk
  # initialize the parameters
  tao = NULL # causal estimand
  Yob = NULL
  sig = NULL
  V = NULL # neyman estimators of variance
  t = NULL # t-statiscs
  # p-value matrix
  p = NULL
  # gen YK obs 
    for (i in 1:K){ 
      # derive the sample mean and sample variance 
      Yob[i] = rnorm(1, Yk[i], sigk[i] * sqrt(K/N))
      sig[i] = sigk[i] * rchisq(1, nk[i] - 1) / (nk[i] - 1)
      Y0 = Yob[1] # non-treat mean
      tao[i] = Yob[i] - Y0 
      V[i] = sig[i]/nk[i] + sig[1]/nk[1]
      # calculate the t-statistics
      t[i] = (tao[i] - 0) / sqrt(V[i])
      # get the two-sided p-value
      p[i] = (1 - pt(abs(t[i]), nk[i] - 2)) * 2
    }
  return(p)
}


##Generate p-combine function
compute_pcombine <- function(Yk, sigk, N, K, ptype, n){
  p_comb_val <- NULL
  for (k in 1:n){
    pmatrix <- compute_pvalue(Yk, sigk, N, K)
    pmatrix <- pmatrix[-1]
    if (ptype == "min"){
      p_comb_val[k] = min(pmatrix)
    }else{
      if (ptype == "max"){
        p_comb_val[k] = max(pmatrix)  
      }else{
        if (ptype == "median"){
          p_comb_val[k] = median(pmatrix)
        }else{
          if (ptype == "mean"){
            p_comb_val[k] = mean(pmatrix)
          }else{
            print("Error: not choosing the existent type")
          }
        }
      }
    }
  }
  return (p_comb_val)
}


##Generate power function
compute_power <- function(Yk, sigk, N, K, p_cv, ptype, R = N/K, n){
  power = 0
  pcomb = compute_pcombine(Yk, sigk, N, K, ptype, n)
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
  return(power)
}

##80% Power critical ratio calculating function
compute_ctratio = function(Y, N, K, r1_start = 0.01, r1_step =0.001, r1_ceil = 1, rk = 1){
  # set the starting critical value
  r1_critical = r1_start
  
  power = ratio_power(Y, r1_critical, rk, N, K, compute_pcombine, compute_power)
  if (r1_start > r1_ceil){
    print("Error: critical value does not exist")
  }else{
    r1_start = r1_start + r1_step
    power_ahead = ratio_power(Y, r1_start, rk, N, K, compute_pcombine, compute_power)
    if (power < 0.8 | power >= power_ahead){
      return(compute_ctratio(Y, n, K, r1_start = r1_start, r1_step = r1_step, r1_ceil = r1_ceil, rk = rk))
    }else{
      return (r1_critical)  
    } 
  }  
}

##Binary critical ratio calculation
compute_binary_ctratio = function(Y, N, K, delta, r1_start = 0, r1_step = 0.15, r1_ceil = 1, rk = 1){
  critical = (r1_start + r1_ceil)/2
  
  power = ratio_power(Y, critical, rk, N, K, compute_pcombine, compute_power)
  
  # check for the monotonicity
  critical_next = critical + r1_step
  
  power_ahead = ratio_power(Y, critical_next, rk, N, K, compute_pcombine, compute_power)
  diff = abs(power - 0.8)
  while (diff > delta) {
    print(diff)
    print(power)
    print(critical)
    if (r1_start == r1_ceil ){
      ## avoid infinite loop case, reduce the floor r1_start
    critical = r1_start/2
    return(compute_binary_ctratio(Y, N, K, delta = delta, r1_start = critical, r1_step = r1_step, r1_ceil = r1_ceil, rk = rk))
    }else{
    if (power > power_ahead | power < 0.8 ){
      return(compute_binary_ctratio(Y, N, K, delta = delta, r1_start = critical, r1_step = r1_step, r1_ceil = r1_ceil, rk = rk))
    }else{
      if (power > 0.8){
        return(compute_binary_ctratio(Y, N, K, delta = delta, r1_start = r1_start, r1_step = r1_step, r1_ceil = critical, rk = rk))
      }else{
        print("Error: does not exists")
        }
      }   
    }
  }  
  return(critical)
}
