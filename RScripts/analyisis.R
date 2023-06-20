#run_report will output two csvs
run_report = function(root_path) {
  
  dalys_by_sim = function(data_) {
    data_ %>%
    group_by(sim.id) %>% 
    group_map(
      function(data,group) {
        data %>% 
          mutate(sim.id = group[1,1]) %>% 
          EMODAnalyzeR::calculate.DALY(life_expectancy = 66) %>%
          mutate(sim.id = group$sim.id) 
      } 
    ) %>% bind_rows
  }
    
  get_infections_and_py = function (data) {
    infections = data %>% mutate(Year = floor(Year)) %>% group_by(Year, sim.id) %>% summarize(Infections=sum(Newly.Infected * pop_scaling_factor))
    py_on_treatment = data %>% filter(HasIntervention.PrEP==1) %>% group_by(Year, sim.id) %>% summarize(Population=sum(Population * pop_scaling_factor)) %>% mutate(Year = floor(Year)) %>% group_by(Year, sim.id) %>% summarize(Population=mean(Population))
    inner_join(infections, py_on_treatment, by=c("Year","sim.id"))
  }
  # todo - offset year just like daly calculation
  infections_averted_over_time = function(baseline, intervention) {
    intervention = calculate.pop_scaling_factor(intervention, 2009, reference_population = 33868111)
    baseline_data = get_infections_and_py(baseline)
    intervention_data = get_infections_and_py(intervention)
    intervention_data$infections.averted = baseline_data$Infections - intervention_data$Infections
    results = intervention_data %>% group_by(Year) %>% summarize(infections.averted = median(infections.averted), py_on_treatment = median(Population))
    results %>% ggplot() + 
      geom_point(aes(x=Year, y=infections.averted))
    data.frame(infections.averted=sum(results$infections.averted), py_on_treatment=sum(results$py_on_treatment))
    
  }
  
  
  
  calc_max_price = function(baseline, intervention, cost_art = 237) {
    intervention = calculate.pop_scaling_factor(intervention, 2009, reference_population = 33868111)
    py = get_infections_and_py(intervention)
    dbs_intervention= dalys_by_sim(intervention)
    dbs_baseline = dalys_by_sim(baseline)
    dbs_intervention$dalys_averted = dbs_baseline$daly_future_discounted - dbs_intervention$daly_future_discounted
    dbs_intervention$on_art_averted = dbs_baseline$discount_factor * (dbs_baseline$on_art - dbs_intervention$on_art)
    py_and_dbs = inner_join(dbs_intervention %>% rename(Year = year), py, by=c("Year", "sim.id"))
    wastage = 0.05
    cost_delivery = (6 + 0.05) * 2 # cost of delivery + cost of syringe
    py_and_dbs %>% mutate(discounted_py = Population * discount_factor) %>%
      mutate(discounted_wasted_py = discounted_py / (1-wastage)) %>%
      summarize (dalys_per_py = sum(dalys_averted) / sum(discounted_wasted_py), art_per_py = sum(on_art_averted) / sum(discounted_wasted_py)) %>% 
      mutate(max_cost_without_art = 500*dalys_per_py) %>% 
      mutate(max_cost_with_art = max_cost_without_art + cost_art*art_per_py) %>%
      mutate(max_cost_with_art_minus_delivery = max_cost_with_art - cost_delivery) %>%
      mutate(max_cost_with_art_minus_delivery_10PctDemandGen = max_cost_with_art_minus_delivery * 0.9) %>% 
      mutate(final_max_cost_per_dose = max_cost_with_art_minus_delivery_10PctDemandGen / 2)
  }
  
  calc_coverage = function(intervention) {
    intervention %>% filter(Year==2040) %>% pivot_wider(id_cols="Year",names_from=c("HasIntervention.PrEP"), values_from=c("Population", "Infected"),values_fn=sum)
  }

  dirs = list.dirs(root_path,recursive = FALSE)
  baseline_dir = Filter(function(.) {grepl("baseline", .)}, dirs)
  experiment_dirs = Filter(function(.) {!grepl("baseline", .)}, dirs)
  
  baseline = read.simulation.results(paste0(baseline_dir,"/ReportHIVByAgeAndGender"), "baseline", stratify_columns = c("Year","Gender","Age","HasIntervention.PrEP"), 
                                               summarize_columns = c("Newly.Infected","Newly.Tested.Positive", "Newly.Tested.Negative", "Population", 
                                                                     "Infected", "On_ART", "Died", "Died_from_HIV",
                                                                     "Tested.Ever", "Diagnosed"), 
                                               min_age_inclusive=15, max_age_inclusive=64)
  baseline = calculate.pop_scaling_factor(baseline, 2009, reference_population = 33868111)
  
  inf_averted_fun = function(path) {
    experiment = 
      paste0(path,"/ReportHIVByAgeAndGender") %>%
      read.simulation.results("", stratify_columns = c("Year","Gender","Age","HasIntervention.PrEP"), 
                              summarize_columns = c("Newly.Infected","Newly.Tested.Positive", "Newly.Tested.Negative", "Population", 
                                                    "Infected", "On_ART", "Died", "Died_from_HIV",
                                                    "Tested.Ever", "Diagnosed"), 
                              min_age_inclusive=15, max_age_inclusive=64)
    
    infections_averted_over_time(baseline, experiment) %>% mutate(experiment = path)
  }
  
  EMODAnalyzeR::bigpurple.add_slurm_to_path()
  bigpurple_opts = list(partition = "a100_short", time = "12:00:00")
  inf_averted = slurmR::Slurm_lapply(as.list(experiment_dirs), 
                                     inf_averted_fun, 
                                     sbatch_opt=bigpurple_opts, 
                                     njobs = length(experiment_dirs), 
                                     export = c("dalys_by_sim", "infections_averted_over_time", "baseline","get_infections_and_py"))
  
  inf_averted %>% bind_rows %>% 
    mutate(experiment = stringr::str_split(experiment, "/") %>% map(~ .[[9]]) %>% unlist()) %>% 
    write.csv(file=paste0(root_path,"/infections_averted.csv"))
    #ggplot + geom_point(aes(x=py_on_treatment ,y=infections.averted)) + scale_x_continuous(expand = c(0, 0)) + 
    #scale_y_continuous(expand = c(0, 0), limits = c(0, 4.5e6)) + geom_text(aes(x=py_on_treatment ,y=infections.averted, label=experiment, hjust=0))
  
  max_price_fun <- function(path) {
    experiment = 
      paste0(path,"/ReportHIVByAgeAndGender") %>%
      read.simulation.results("", stratify_columns = c("Year","Gender","Age","HasIntervention.PrEP"), 
                              summarize_columns = c("Newly.Infected","Newly.Tested.Positive", "Newly.Tested.Negative", "Population", 
                                                    "Infected", "On_ART", "Died", "Died_from_HIV",
                                                    "Tested.Ever", "Diagnosed"), 
                              min_age_inclusive=15, max_age_inclusive=64)
    calc_max_price(baseline, experiment) %>% mutate(experiment = path)
  }
  
  
  max_prices = slurmR::Slurm_lapply(as.list(experiment_dirs), 
                                    max_price_fun, 
                                    sbatch_opt=bigpurple_opts, 
                                    njobs = length(experiment_dirs), 
                                    export = c("dalys_by_sim", "calc_max_price", "baseline","get_infections_and_py"))
  
  max_prices %>% bind_rows %>% 
      mutate(experiment = stringr::str_split(experiment, "/") %>% map(~ .[[9]]) %>% unlist()) %>% 
      write.csv(file=paste0(root_path,"/max_prices.csv"))
}

run_report( root_path = "/gpfs/data/bershteynlab/EMOD/kaftad01/202306_SA_NewEXE/output_2023_06_19_17_25")
run_report( root_path = "/gpfs/data/bershteynlab/EMOD/kaftad01/202306_Nyanza_gates_len/output_2023_06_18_12_22")
run_report( root_path = "/gpfs/data/bershteynlab/EMOD/kaftad01/202306_Zimbabwe_NewEXE/output_2023_06_19_20_08")
