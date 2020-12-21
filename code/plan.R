run_simulation <- drake_plan(
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
                              R = 0)
  ),

  # Run Epidemic Simulations
  epidemic_simulations = model_stochastic_parameters_multiwrapper(simulation_parameters),

  # Estimate Rt from simulation with accurate parameters
  rt_estimates_params_acc = map(epidemic_simulations, function(i) estimates_rt(i,
                                                    mean_gi = attributes(i)$metadata["mean_gi"],
                                                     gamma_val = attributes(i)$metadata["gamma_val"], tau = 7,
                                                     error.mean = 0, error.sd = 0)),

  # Estimate Rt from simulation with inaccurate parameters
  rt_estimates_params_inacc = map(epidemic_simulations, function(i) estimates_rt(i,
                                   mean_gi = attributes(i)$metadata["mean_gi"],
                                   gamma_val = attributes(i)$metadata["gamma_val"], tau = 7,
                                   error.mean = 0, error.sd = 2))
)
