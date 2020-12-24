source("./code/packages.R")
source("./code/functions.R")
source("./code/plan.R")
source("./code/create_figures.R")
source("./code/COVID19_data_analysis.R")

drake_config(
  plan,
  verbose = 2
)
