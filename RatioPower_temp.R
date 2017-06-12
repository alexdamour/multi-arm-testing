# define a function that returns the power value given 
# ratio r1 = y1/sig1 and ratio rk = sigk/sig1
ratio_power = function(Y, r1, rk, N, K, compute_pcombine, compute_power, knum = 1, ktype = "homo", ptype = "min", power_rep = N){
  # the treatment level y1
  y1 = Y[2]
  sig1 = y1/r1
  # derive the null mean Y_null
  Y_null = rep(0, K)
  
  # derive the standard deviation sig array 
  sigk = NULL
  if (ktype == "homo"){
    sigk = rep(rk* y1/r1, K)
  }else{
    if (ktype == "rk-hetero"){
      sigk = rk * y1/r1
    }else{
      if(ktype == "sig-hetero"){
        sigk_rest = exp(rnorm(K-1, log(rk*sig1)))
        sigk = c(sigk_rest[1], 0, sigk_rest[2:K-1])     
      }else{
        print("Error: not choosing the existent ktype-type")
      }
    }
  }
  
  sigk[2] = y1/r1 
  
  # compute the corresponding p_power 
  power = NULL
  
  # get the experiment number to calculate power
  R = power_rep/K
  
  # derive the critical value p_cv for null hypothesis
  p_cv = quantile(compute_pcombine(Y_null, sigk, N, K, ptype), 0.05)
  # get the sparse effect 
  for (j in 1:knum){
    power[j] = compute_power(Y, sigk, N, K, p_cv, ptype, R)
  }
  avg_power = mean(power)
  
  return(avg_power) 
}

# minimal sigma case 
sig_stat_gen = function(mean, K, ratio_type = "min", sig_type = "uniform"){
  if (sig_type == "uniform"){
    rk = runif(K, max = 2 * mean)
  }else{
    if (sig_type == "log-normal"){
      rk = exp(rnorm(K, log(mean)))
    }else{
      print("Error: not choosing the existent sig-type")
    }
  }
  # get the rest of sigmas apart from sig1
  rk_rest = c(rk[1], rk[3:K])
  
  if (ratio_type == "min"){
    rk_item = min(rk_rest)
  }else{
    if (ratio_type == "max" ){
      rk_item = max(rk_rest)  
    }else{
      if (ratio_type == "median"){
        rk_item = median(rk_rest)
      }else{
        if (ratio_type == "mean"){
          rk_item = mean(rk_rest)
        }else{
          print("Error: not choosing the existent type")
        }
      }
    }
  }
  return( list(rk_item = rk_item,  rk = rk) )
}

power_ratio_test = function(Y, N, K, n_r1, n_rk, ktype = "homo", knum = 1, ratio_type, sig_type, ratio_rep = 50, step_r1 = 0.1, step_rk = 0.1){
  
  power_matrix= matrix(0, nrow = n_r1, ncol = n_rk)
  r1_item = NULL
  rk_item = NULL
  
  if (ktype == "homo" | ktype == "sig-hetero"){
    for (i in 1:n_r1){ # Row
      r1_i = i * step_r1
      r1_item[i] = r1_i
      for (j in 1:n_rk){ # Column
        rk_j = j * step_rk
        rk_item[j] = rk_j
        power_matrix[i,j] = ratio_power(Y, r1_i, rk_j, N, K, compute_pcombine, compute_power, ktype = ktype) 
      }
    }
  }else{
    if (ktype == "rk-hetero"){
      for (i in 1:n_r1){
        r1_i = i * step_r1
        r1_item[i] = r1_i
        ratio_exp = NULL
        for (j in 1:n_rk){ # Column
          mean_j = j * step_rk
          rk_item_exp = NULL
          power_exp = NULL
          for (e in 1:ratio_rep){ # repeated experiments
            statis_exp = sig_stat_gen(mean_j, K, ratio_type = ratio_type, sig_type = sig_type)
            # get the rk_item sample
            rk_item_exp[e] = statis_exp$rk_item 
            # get the corresponding power
            rk_exp = statis_exp$rk
            power_exp[e] = ratio_power(Y, r1_i, rk_exp, N, K, compute_pcombine, compute_power, ktype = ktype, knum = knum)
          }
          # get the mean of rk_item and power for repeated experiments
          rk_item[j] = mean(rk_item_exp)  
          power_matrix[i,j] = mean(power_exp)}
      }
    }else{
      if (ktype == "rk-hetero"){
        for (i in 1:n_r1){
          r1_i = i * step_r1
          r1_item[i] = r1_i
          ratio_exp = NULL
          for (j in 1:n_rk){ # Column
            mean_j = j * step_rk
            rk_item_exp = NULL
            power_exp = NULL
            for (e in 1:ratio_rep){ # repeated experiments
              statis_exp = sig_stat_gen(mean_j, K, ratio_type = ratio_type, sig_type = sig_type)
              # get the rk_item sample
              rk_item_exp[e] = statis_exp$rk_item 
              # get the corresponding power
              rk_exp = statis_exp$rk
              power_exp[e] = ratio_power(Y, r1_i, rk_exp, n, K, compute_pcombine, compute_power, ktype = ktype, knum = knum)
            }
            # get the mean of rk_item and power for repeated experiments
            rk_item[j] = mean(rk_item_exp)  
            power_matrix[i,j] = mean(power_exp)
          }
        }
      }else{
        print("Error: not choosing the appropriate ktype")
      }
    }
  }
  return( list(power_matrix = power_matrix, rk = rk_item, r1 = r1_item) )
}

