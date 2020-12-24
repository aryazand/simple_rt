
##----General Epidemic Modeling Functions

seir_model <- function(time, state, parameters) {

  with(as.list(c(state, parameters)), {  # tell R to unpack variable names from the state and parameters inputs

    # Calculating the total population size N (the sum of the number of people in each compartment)
    N <- S+E+I+R

    # Defining lambda as a function of beta and I:
    lambda <- beta * I/N

    # The differential equations
    dS <- (-lambda * S)
    dE <- (lambda * S)  - (delta * E)
    dI <- (delta * E) - (gamma * I)
    dR <- (gamma * I)

    return(list(c(dS, dE, dI, dR)))
  })

}

gi_distribution <- function(delta_val, gamma_val) {

  mean_latent_period = 1/delta_val
  mean_infectious_period = 1/gamma_val
  mean_gi = mean_latent_period + mean_infectious_period


  # Per park et al, 2018 - Appendix A:
  ### squared-coefficient-of-variance (SCV) = 1 - (2/(mean_gi^2*delta_val*gamma_val))
  ### shape = 1/SCV

  # coefficient-of-variance = sd/mean
  # sd of gamma distribution = sqrt(shape)*scale
  # mean of gamma distribution = shape*scale
  # thus CV = 1/sqrt(shape)

  SCV = 1 - (2/(mean_gi^2*delta_val*gamma_val))
  shape_gi = 1/SCV
  scale_gi = mean_gi/shape_gi
  sd_gi = sqrt(SCV)*mean_gi

  return(c(shape = shape_gi, scale = scale_gi, sd = sd_gi))
}

logistic_func = function(n1, n2, length.out, k) {
  L = n2 - n1
  x = map_dbl(seq(-4,4, length.out = length.out), function(x) L/(1 + exp(k*(-x))))
  x = x + n1
  return(x)
}

create_r0_vector <- function(r0, period_lengths, transition_lengths) {

  time_lengths <- vector(class(period_lengths), length(c(period_lengths, transition_lengths)))
  time_lengths[c(TRUE, FALSE)] <- period_lengths
  time_lengths[c(FALSE, TRUE)] <- transition_lengths

  r0 = rep(r0, each=2)

  r0_vector = pmap(list(head(r0, -1), tail(r0,-1), time_lengths), function(x,y,z) logistic_func(n1=x,n2=y,length.out=z, k=1))
  r0_vector = r0_vector %>% unlist()

  return(r0_vector)
}

model_pandemic <- function(initial_state_values, mean_gi, gamma_val, r0, period_lengths, transition_lengths) {

  time_length = sum(c(period_lengths, transition_lengths))
  delta_val = 1/(mean_gi - 1/gamma_val)

  # Generational Interval Distribution
  gi_gamma_param = gi_distribution(delta_val, gamma_val)

  # Model Parameters
  r0_vector = create_r0_vector(r0, period_lengths, transition_lengths)
  parameters <- c(beta = r0_vector*gamma_val,
                  delta = rep(delta_val, time_length),
                  gamma = rep(gamma_val, time_length))

  parameters <- matrix(parameters, ncol=3)
  colnames(parameters) = c("beta", "delta", "gamma")

  # Vector storing the sequence of timesteps to solve the model at
  times <- c(0,1)  # from 0 to 100 days in daily intervals

  # MODEL OUTPUT (solving the differential equations):
  # For-loop through each time in the model
  # Solving the differential equations using the ode integration algorithm
  output <- vector(mode = "list", length=time_length)

  state_values = initial_state_values
  for(i in 1:time_length) {
    time_step = as.data.frame(ode(y = state_values,
                                  times = times,
                                  func = seir_model,
                                  parms = parameters[i,]))

    state_values = unlist(time_step[2,2:5])

    if(time_step$S[2] < 1000 | time_step$I[2] < 1) {
      break
    }

    output[[i]] = time_step

  }

  output = bind_rows(output)
  output = output[seq(2, nrow(output), by=2),]
  output$time = seq(1, nrow(output))

  # calculate Rt
  output$r0 = r0_vector[1:nrow(output)]
  N = output$S + output$E + output$I + output$R

  output <- output %>%
    mutate(proportion_susceptible = S/N) %>%
    mutate(r_eff = r0 * proportion_susceptible) %>%
    mutate(total_infected = round(E + I + R)) %>%
    mutate(new_cases = c(NA, diff(total_infected)))


  output = as.tibble(output)
  attr(output, "metadata") = c(mean_gi = mean_gi, gamma_val = gamma_val)

  return(output)
}

model_stochastic_parameters <- function(initial_state_values, mean_gi_range, r0_range, num_periods, period_lengths_range, transition_lengths_range,
                                        min_simLength, min_newCases) {

  i = 0
  sim_data = data.frame()
  while(nrow(sim_data) < min_simLength | any(sim_data$new_cases < min_newCases, na.rm = T)) {
    # i = i + 1
    # print(paste("trial...",i))

    # Setup generation-interval parameters
    mean_gi = sample(mean_gi_range, 1)
    gamma_val = 1/sample(1:(mean_gi-1), 1)
    delta_val = 1/(mean_gi - 1/gamma_val)
    gi_distribution_params = gi_distribution(delta_val, gamma_val)

    # Setup r0
    r0 = sample(r0_range, num_periods)

    # Period and transition lengths
    sample.vec <- function(x, ...) x[sample(length(x), ...)]
    period_lengths = sample.vec(period_lengths_range, size = num_periods, replace = T)
    transition_lengths = sample.vec(transition_lengths_range, size = num_periods-1, replace = T)

    # Run Model
    sim_data = model_pandemic(initial_state_values = initial_state_values,
                              mean_gi = mean_gi,
                              gamma_val=gamma_val,
                              r0 = r0,
                              period_lengths = period_lengths,
                              transition_lengths = transition_lengths)
  }

  return(sim_data)
}

model_stochastic_parameters_multiwrapper <- function(simulation_parameters) {

  # Set seed
  myseed = 06161989
  set.seed(myseed)

  epidemic_simulations =
    with(simulation_parameters,
         replicate(Nsims,
                   model_stochastic_parameters(initial_state_values = initial_state_values,
                                               mean_gi_range = seq(mean_gi_range[1], mean_gi_range[2], mean_gi_range[3]),
                                               r0_range = seq(r0_range[1], r0_range[2], r0_range[3]),
                                               num_periods = num_periods,
                                               period_lengths = seq(period_lengths[1], period_lengths[2], period_lengths[3]),
                                               transition_lengths =  seq(transition_lengths[1], transition_lengths[2], transition_lengths[3]),
                                               min_simLength = min_simLength,
                                               min_newCases = min_newCases), simplify = F)
    )
}

##----Esimate Rt functions

estimates_rt.simple <- function(time, Is, mean_gi, tau) {

  df = tibble(time, Is)

  # Basic Ratio
  df = df %>%
    mutate(new_cases_sum = roll_sum(Is,  tau, align="center", fill = c(NA, NA, NA))) %>%
    mutate(new_cases_sum_gi_days_ago = lag(new_cases_sum, mean_gi)) %>%
    mutate(rt.simple_mean = (new_cases_sum+0.1)/(new_cases_sum_gi_days_ago+0.1))

  # Confidence Interval.using beta representation of Clopper-Pearson method
  # (instead of using poisson.test(), which uses binom.test indirectly)

  p.L <- function(x, n, alpha = 0.025) {
    p = qbeta(alpha, x, n - x + 1)
    p[x == 0] = 0
    return(p)
  }
  p.U <- function(x, n, alpha = 0.025) {
    p = qbeta(1 - alpha, x + 1, n - x)
    p[x == n] = 1
    return(p)
  }

  CI.L = p.L(x = df$new_cases_sum, n = df$new_cases_sum + df$new_cases_sum_gi_days_ago)
  CI.U = p.U(x = df$new_cases_sum, n = df$new_cases_sum + df$new_cases_sum_gi_days_ago)

  df$rt.simple_quantile_025 = CI.L/(1-CI.L)
  df$rt.simple_quantile_975 = CI.U/(1-CI.U)

  # Select columns to output
  df = df %>% dplyr::select(time, starts_with("rt."))

  return(df)
}

estimates_rt.cori <- function(time, Is, mean_gi, sd_gi, tau) {

  r_estimates = EpiEstim::estimate_R(Is,
                                     method = "parametric_si",
                                     config = make_config(list(mean_si = mean_gi,
                                                               std_si = sd_gi,
                                                               t_start = seq(2, length(Is) - tau),
                                                               t_end = seq(2+tau, length(Is))))
  )

  df = tibble(r_estimates$R) %>%
    select(t_start, t_end,
           rt.cori_mean = `Mean(R)`,
           rt.cori_quantile_025 = `Quantile.0.025(R)`,
           rt.cori_quantile_975 = `Quantile.0.975(R)`) %>%
    mutate(time = time[ceiling((t_start + t_end)/2)]) %>%
    select(-t_start, -t_end)

  return(df)
}

estimates_rt <- function(data, mean_gi, gamma_val, tau, error.mean = 0, error.sd = 0) {

  gamma_val = unname(gamma_val)
  delta_val = 1/(mean_gi - 1/gamma_val) %>% unname()

  # Calculate Generational Interval Distribution
  gi_gamma_param = gi_distribution(delta_val, gamma_val)

  # Calculate rt.simple_mean
  df.rt_simple <- estimates_rt.simple(time = data$time[-1], Is = data$new_cases[-1], mean_gi = mean_gi + error.mean, tau = tau)
  data = left_join(data, df.rt_simple, by = "time")


  # Calculate rt.cori
  df.rt_cori <- estimates_rt.cori(time = data$time[-1], Is = data$new_cases[-1], mean_gi = mean_gi + error.mean,
                                  sd_gi = gi_gamma_param["sd"] + error.sd, tau = tau)

  data = full_join(data, df.rt_cori, by="time")


  r_eff = data$r_eff
  data = data %>%
    mutate_at(vars(ends_with("_mean")), list(error = function(x) x - r_eff)) %>%
    mutate_at(vars(ends_with("_quantile_975")), list(error = function(x) x > r_eff)) %>%
    mutate_at(vars(ends_with("_quantile_025")), list(error = function(x) x < r_eff)) %>%
    pivot_longer(matches("Quantile_[[:digit:]]{3}_error"), names_to = c("Method", "Quantile"),
               names_pattern = "rt.([[:alpha:]]+)_(quantile_[[:digit:]]{3})", values_to = "value") %>%
    pivot_wider(names_from = "Quantile", values_from = "value") %>%
    mutate(CI_error = quantile_975 + quantile_025 == 2) %>%
    select(-starts_with("quantile_")) %>%
    pivot_wider(names_from = "Method", names_glue = "rt.{Method}_in_CI", values_from = "CI_error")


  return(data)
}

