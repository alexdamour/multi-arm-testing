#Fundamental Simulation Function
##Generate p-value function
compute_pvalue = function(Yk, sigk, N, K, n = 1000){
  nk = N/K # assume evenly distributed number of subgroups nk
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
    for (i in 1: K){ 
      # derive the sample mean and sample variance 
      Yob[i] = rnorm(1, Yk[i], sigk[i] * sqrt(K/N))
      sig[i] = sigk[i] * rchisq(1, nk - 1) / (nk - 1)
      Y0 = Yob[1] # non-treat mean
      tao[i] = Yob[i] - Y0 
      V[i] = sig[i]/nk + sig[1]/nk
      # calculate the t-statistics
      t[i] = (tao[i] - 0) / sqrt(V[i])
      # get the two-sided p-value
      p[j,i] = (1 - pt(abs(t[i]), nk - 2)) * 2
    }
  }
  return(p)
}


##Generate p-combine function
compute_pcombine <- function(Yk, sigk, N, K, ptype = "min", n =1000){
  # p-combined value array
  p_comb_val <- NULL
  pmatrix <- compute_pvalue(Yk, sigk, N, K, n)
  n = nrow(pmatrix)
  for (k in 1:n){
    if (ptype == "min"){
      p_comb_val[k] = min(pmatrix[k,-1])
    }else{
      if (ptype == "max"){
        p_comb_val[k] = max(pmatrix[k,-1])  
      }else{
        if (ptype == "median"){
          p_comb_val[k] = median(pmatrix[k,-1])
        }else{
          if (ptype == "mean"){
            p_comb_val[k] = mean(pmatrix[k,-1])
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
compute_power <- function(Yk, sigk, N, K, p_cv, ptype = "min", R = 1000){
  # R denotes repeated experiment times for calculating power
  
  power = 0
  pcomb = compute_pcombine(Yk, sigk, N, K, ptype, n = R)
  
  # the success rate 
  rejec = 0
  rec = 0
  
  for (i in 1:R){
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

##Multi-dimension Y generate function
multi_y = function(epsilon, K){
Ymv_temp <- rmvnorm(1, rep(0, K-1), diag(K-1))
Ymv <- c(0, Ymv_temp)
Ynorm <- Ymv * epsilon / sqrt(sum(Ymv^2))
return(Ynorm)
}


##Binary level calculation 
compute_binary_level = function(Y, sigk, N, K, delta, p_cv, eps_start = 0, eps_ceil = 1, eps_rep = 100, power_rep = 100, ptype = "min" ){

  # the iteration value
  critical = (eps_start + eps_ceil)/2
  
  power_exp = NULL
  for (exp in 1:eps_rep){
    Y_e = multi_y(critical, K)
    power_exp[exp] = compute_power(Y_e, sigk, N, K, p_cv, ptype = ptype, R = power_rep)
  }
  power = mean(power_exp)
  
  # control for the difference
  diff = abs(power - 0.8)
  
  while (diff > delta) {
    print(diff)
    print(power)
    print(critical)
    
    if (power < 0.8 ){
      return(compute_binary_level(Y, sigk, N, K, delta = delta, p_cv = p_cv, eps_start = critical, eps_ceil = eps_ceil, ptype = ptype))
    }else if (power >= 0.8){
        return(compute_binary_level(Y, sigk, N, K, delta = delta, p_cv = p_cv,eps_start = eps_start, eps_ceil = critical, ptype = ptype))
        }
    }
  return(critical)
}



