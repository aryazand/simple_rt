
download_covid19_data = function() {

  national_data = get_national_data(totals = F, source = "WHO")

  national_data = national_data %>%
    mutate(time = date) %>% select(-date) %>%
    filter(time > as.Date("2020-04-01") & time < as.Date("2020-12-01"))

  return(national_data)
}


