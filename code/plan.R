plan <- drake_plan(
  # Load simulation parameters
  simulation_parameters = list(
    mean_gi_range = c(4,12,1),
    r0_range = c(0.5,1.8,0.01),
    num_periods = 3,
    period_lengths = c(50,50,1),
    transition_lengths = c(14,30,1),
    Nsims = 100,
    initial_state_values = c(S = 2e7-100,
                              E = 0,
                              I = 1000,
                              R = 0),
    min_simLength = 150,
    min_newCases = 20
  ),

  # Run Epidemic Simulations
  epidemic_simulations = model_stochastic_parameters_multiwrapper(simulation_parameters),

  # Estimate Rt from simulation with accurate parameters
  rt_estimates = map(epidemic_simulations, function(i) estimates_rt(i,
                      mean_gi = attributes(i)$metadata["mean_gi"],
                       gamma_val = attributes(i)$metadata["gamma_val"], tau = 7,
                       error.mean = 0, error.sd = 0)) %>% bind_rows(., .id = "simulation_run"),

  # Estimate Rt from simulation with inaccurate mean
  rt_estimates_meanerr = map(epidemic_simulations, function(i) estimates_rt(i,
                               mean_gi = attributes(i)$metadata["mean_gi"],
                               gamma_val = attributes(i)$metadata["gamma_val"], tau = 7,
                               error.mean = 2, error.sd = 0)) %>% bind_rows(., .id = "simulation_run"),

  # Estimate Rt from simulation with inaccurate sd
  rt_estimates_sderr = map(epidemic_simulations, function(i) estimates_rt(i,
                           mean_gi = attributes(i)$metadata["mean_gi"],
                           gamma_val = attributes(i)$metadata["gamma_val"], tau = 7,
                           error.mean = 0, error.sd = 3)) %>% bind_rows(., .id = "simulation_run"),

  # Real-world data
  COVID19_params = list(mean_gi = 4.83, sd_gi = 1.72),
  national_data_raw = download_covid19_data(),
  national_data = national_data_raw %>%
    group_by(country) %>% nest() %>%
    mutate(rt.simple = map(data, function(x) estimates_rt.simple(time = x$time, Is = x$cases_new, mean_gi = round(COVID19_params$mean_gi), tau = 7))) %>%
    mutate(rt.cori = map(data, function(x) estimates_rt.cori(time = x$time, Is = x$cases_new, mean_gi = COVID19_params$mean_gi, sd_gi = COVID19_params$sd_gi, tau = 7))) %>%
    mutate(data = map2(data, rt.simple, full_join, by = "time")) %>% select(-rt.simple) %>%
    mutate(data = map2(data, rt.cori, full_join, by = "time")) %>% select(-rt.cori) %>%
    unnest(data) %>%
    mutate(estimation_diff = abs(rt.cori_mean - rt.simple_mean)),
)
