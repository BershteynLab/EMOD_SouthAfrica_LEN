library(purrr)
library(rsample)

#run_report will output two csvs
run_report = function(root_path, cost_art, life_expectancy, pop_scale_param_inst) {

  emodplot.incidence = 
    function (data_baseline, intervention_df, date.start, date.end) 
    {
      data.incidence <- EMODAnalyzeR::calculate.incidence(data_baseline)
      y.lim.max <- min(max((data.incidence %>% filter(Year > date.start, Year < date.end))$incidence * 1.2), 1)
      intervention.incidence <- EMODAnalyzeR::calculate.incidence(intervention_df)
      all_data_incidence = rbind(data.incidence, intervention.incidence)
      p <- emodplot.by_gender(all_data_incidence, date.start, date.end, 
                              "incidence") + 
                              scale_y_continuous(labels = scales::percent_format(accuracy = .1), 
                                                 breaks = seq(0, y.lim.max, 0.005), limits = c(0, y.lim.max)) + 
                              ylab("HIV Incidence (%)")
      return(p)
    }
  
  dalys_by_sim = function(data_) {
    data_ %>%
    group_by(sim.id) %>% 
    group_map(
      function(data,group) {
        data %>% 
          mutate(sim.id = group[1,1]) %>% 
          EMODAnalyzeR::calculate.DALY(life_expectancy = life_expectancy) %>%
          mutate(sim.id = group$sim.id) 
      } 
    ) %>% bind_rows
  }
    
  get_infections_and_py = function (data) {
    infections = data %>% 
                  mutate(Year = floor(Year)) %>% 
                  group_by(Year, sim.id) %>% 
                  summarize(Infections=sum(Newly.Infected * pop_scaling_factor))
    py_on_treatment = data %>% filter(HasIntervention.PrEP==1) %>% 
      group_by(Year, sim.id) %>% 
      summarize(Population=sum(Population * pop_scaling_factor)) %>% 
      mutate(Year = floor(Year)) %>% 
      group_by(Year, sim.id) %>% 
      summarize(py_on_treatment=mean(Population))
    py_total = data %>% 
      group_by(Year, sim.id) %>% 
      summarize(Population=sum(Population * pop_scaling_factor)) %>% 
      mutate(Year = floor(Year)) %>% 
      group_by(Year, sim.id) %>% 
      summarize(py_total=mean(Population))
    inner_join(infections, py_on_treatment, by=c("Year","sim.id")) %>% inner_join(py_total, by=c("Year","sim.id"))
  }
  # todo - offset year just like daly calculation
  infections_averted_over_time = function(baseline, intervention) {
    baseline_data = get_infections_and_py(baseline)
    intervention_data = get_infections_and_py(intervention)
    intervention_data$infections.averted = baseline_data$Infections - intervention_data$Infections
    intervention_data$infections.baseline = baseline_data$Infections
    results = intervention_data %>% 
              group_by(Year) %>%
              summarize(infections.averted = median(infections.averted),
                        py_on_treatment = median(py_on_treatment),
                        py_total = median(py_total),
                        baseline_infections = median(infections.baseline))
    results_time_on_treatment = results %>% filter(Year >= 2025, Year < 2045)
    results %>% ggplot() + 
      geom_point(aes(x=Year, y=infections.averted))
    data.frame(infections.averted=sum(results$infections.averted), 
               py_on_treatment=sum(results$py_on_treatment),
               percent_coverage=sum(results_time_on_treatment$py_on_treatment)/sum(results_time_on_treatment$py_total),
               percent_infections_averted = sum(results_time_on_treatment$infections.averted)/sum(results_time_on_treatment$baseline_infections))
    
  }
  
  
  
  calc_max_price = function(baseline, intervention, cost_art) {
    intervention = calculate.pop_scaling_factor(intervention, 
                                                pop_scale_param_inst$year, 
                                                reference_population = pop_scale_param_inst$population,
                                                age_max_inclusive = pop_scale_param_inst$age_max_inc,
                                                age_min_inclusive = pop_scale_param_inst$age_min_inc)
    py = get_infections_and_py(intervention)
    dbs_intervention= dalys_by_sim(intervention)
    dbs_baseline = dalys_by_sim(baseline)
    dbs_intervention$dalys_averted = dbs_baseline$daly_future_discounted - dbs_intervention$daly_future_discounted
    dbs_intervention$on_art_averted = dbs_baseline$discount_factor * (dbs_baseline$on_art - dbs_intervention$on_art)
    py_and_dbs = inner_join(dbs_intervention %>% rename(Year = year), py, by=c("Year", "sim.id"))
    cost_delivery = (6 + 0.05 + 2.5) * 2 # cost of delivery + cost of syringe + COST OF test
    summary_fun <- function (.data) {
      .data %>% summarize(dalys_per_py = sum(dalys_averted) / sum(discounted_py), 
                          art_per_py = sum(on_art_averted) / sum(discounted_py) ) %>%
                mutate(max_cost_without_art = 500*dalys_per_py) %>% 
                mutate(max_cost_with_art = max_cost_without_art + cost_art*art_per_py) %>%
                mutate(demand_generation = .1 * max_cost_with_art) %>%
                mutate(max_cost_with_art_minus_dg = max_cost_with_art - demand_generation) %>%
                mutate(max_cost_with_art_minus_dg_delivery = max_cost_with_art_minus_dg - cost_delivery) %>%
                mutate(max_cost_per_dose = max_cost_with_art_minus_dg_delivery / 2) %>%
                mutate(wastage_per_dose = 0.05*max_cost_per_dose) %>%
                mutate(max_cost_per_dose_minus_wastage = max_cost_per_dose - wastage_per_dose)
    }
    bootstraps = 
      py_and_dbs %>% 
        mutate(discounted_py = py_on_treatment * discount_factor) %>% group_by(sim.id) %>%
        summarize (dalys_averted = sum(dalys_averted),
                   discounted_py = sum(discounted_py),
                   on_art_averted = sum(on_art_averted),
                   undiscounted_py = sum(py_on_treatment)) %>% 
        bootstraps( times = 500)
    resamples <- map_dfr(bootstraps$splits, function(.data) {.data %>% analysis() %>% summary_fun()} )
    summary <- rbind(      resamples %>% summarize_all(function(data) {quantile(data, 0.5)}) %>% mutate(aggregation = "Median"),
                           resamples %>% summarize_all(function(data) {quantile(data, 0.025)}) %>% mutate(aggregation = "Lower Bound (95% CI)"),
                           resamples %>% summarize_all(function(data) {quantile(data, 0.975)}) %>% mutate(aggregation = "Upper Bound (95% CI)"))
    summary        
  }
  
  calc_coverage = function(intervention) {
    intervention %>% filter(Year==2040) %>% pivot_wider(id_cols="Year",names_from=c("HasIntervention.PrEP"), values_from=c("Population", "Infected"),values_fn=sum)
  }

  dirs = list.dirs(root_path,recursive = FALSE)
  baseline_dir = Filter(function(.) {grepl("-baseline", .)}, dirs)
  experiment_dirs = Filter(function(.) {!grepl("-baseline", .)}, dirs)
  
  baseline = read.simulation.results(paste0(baseline_dir,"/ReportHIVByAgeAndGender"), "baseline", stratify_columns = c("Year","Gender","Age","HasIntervention.PrEP"), 
                                               summarize_columns = c("Newly.Infected","Newly.Tested.Positive", "Newly.Tested.Negative", "Population", 
                                                                     "Infected", "On_ART", "Died", "Died_from_HIV",
                                                                     "Tested.Ever", "Diagnosed"), 
                                               min_age_inclusive=15, max_age_inclusive=64)
  baseline = calculate.pop_scaling_factor(baseline, 
                                          pop_scale_param_inst$year, 
                                          reference_population = pop_scale_param_inst$population,
                                          age_max_inclusive = pop_scale_param_inst$age_max_inc,
                                          age_min_inclusive = pop_scale_param_inst$age_min_inc)
  
  
  inf_averted_fun = function(path) {
    experiment = 
      paste0(path,"/ReportHIVByAgeAndGender") %>%
      read.simulation.results("", stratify_columns = c("Year","Gender","Age","HasIntervention.PrEP"), 
                              summarize_columns = c("Newly.Infected","Newly.Tested.Positive", "Newly.Tested.Negative", "Population", 
                                                    "Infected", "On_ART", "Died", "Died_from_HIV",
                                                    "Tested.Ever", "Diagnosed"), 
                              min_age_inclusive=15, max_age_inclusive=64)
    experiment =    calculate.pop_scaling_factor(experiment, 
                                                pop_scale_param_inst$year, 
                                                reference_population = pop_scale_param_inst$population,
                                                age_max_inclusive = pop_scale_param_inst$age_max_inc,
                                                age_min_inclusive = pop_scale_param_inst$age_min_inc)
    plt = experiment %>% 
      filter(Gender==1, Age < 50) %>% 
      group_by(Year, HasIntervention.PrEP, sim.id) %>% 
      summarize(Population=sum(Population)) %>% 
      pivot_wider(id_cols=c("Year","sim.id"),names_from=c("HasIntervention.PrEP"), values_from=Population) %>% 
      mutate(PctCoverage=`1`/(`1`+`0`)) %>% 
      ggplot() + geom_point(aes(x=Year, y=PctCoverage))
    ggplot2::ggsave(paste0(path,"/coverage_women.png"),plot=plt)
    
    plt = experiment %>% 
      filter( Age < 50) %>% 
      group_by(Year, HasIntervention.PrEP, sim.id) %>% 
      summarize(Population=sum(Population)) %>% 
      pivot_wider(id_cols=c("Year","sim.id"),names_from=c("HasIntervention.PrEP"), values_from=Population) %>% 
      mutate(PctCoverage=`1`/(`1`+`0`)) %>% 
      ggplot() + geom_point(aes(x=Year, y=PctCoverage))
    ggplot2::ggsave(paste0(path,"/coverage_pop.png"),plot=plt)
    infections_averted_over_time(baseline %>% filter(Age<50), experiment %>% filter(Age<50)) %>% mutate(experiment = path)
  }
  
  EMODAnalyzeR::bigpurple.add_slurm_to_path()
  bigpurple_opts = list(partition = "a100_short", time = "12:00:00")
  inf_averted = slurmR::Slurm_lapply(as.list(experiment_dirs), 
                                     inf_averted_fun, 
                                     sbatch_opt=bigpurple_opts, 
                                     njobs = length(experiment_dirs), 
                                     export = c("dalys_by_sim", "infections_averted_over_time", "baseline","get_infections_and_py", "pop_scale_param_inst"))
  
  inf_averted %>% bind_rows %>% 
    mutate(experiment = stringr::str_split(experiment, "/") %>% map(~ .[[9]]) %>% unlist()) %>% 
    write.csv(file=paste0(root_path,"/infections_averted.csv"))
    #ggplot + geom_point(aes(x=py_on_treatment ,y=infections.averted)) + scale_x_continuous(expand = c(0, 0)) + 
    #scale_y_continuous(expand = c(0, 0), limits = c(0, 4.5e6)) + geom_text(aes(x=py_on_treatment ,y=infections.averted, label=experiment, hjust=0))
  
  max_price_fun <- function(path) {
    experiment = 
      paste0(path,"/ReportHIVByAgeAndGender") %>%
      read.simulation.results("Intervention", stratify_columns = c("Year","Gender","Age","HasIntervention.PrEP"), 
                              summarize_columns = c("Newly.Infected","Newly.Tested.Positive", "Newly.Tested.Negative", "Population", 
                                                    "Infected", "On_ART", "Died", "Died_from_HIV",
                                                    "Tested.Ever", "Diagnosed"), 
                              min_age_inclusive=15, max_age_inclusive=64)
    plt = emodplot.incidence(baseline, experiment, 2020,2059)
    ggplot2::ggsave(paste0(path,"/incidence.png"),plot=plt)
    calc_max_price(baseline, experiment, cost_art) %>% mutate(experiment = path)
  }
  
  
  max_prices = slurmR::Slurm_lapply(as.list(experiment_dirs), 
                                    max_price_fun, 
                                    sbatch_opt=bigpurple_opts, 
                                    njobs = length(experiment_dirs), 
                                    export = c("dalys_by_sim", "calc_max_price", "baseline","get_infections_and_py", "emodplot.incidence","pop_scale_param_inst","cost_art","life_expectancy"))
  
  max_prices %>% bind_rows %>% 
      mutate(experiment = stringr::str_split(experiment, "/") %>% map(~ .[[9]]) %>% unlist()) %>% 
      write.csv(file=paste0(root_path,"/max_prices.csv"))
}

run_report( root_path = "/gpfs/data/bershteynlab/EMOD/kaftad01/202306_SA_NewEXE/output_2023_06_21_14_20", 
            cost_art=187, 
            life_expectancy = 66, 
            pop_scale_param_inst = pop_scale_params(year=2009,population=33868111, age_min_inc=15, age_max_inc= 64))
run_report( root_path = "/gpfs/data/bershteynlab/EMOD/kaftad01/202306_Nyanza_gates_len/output_consolidated", 
            cost_art=193, 
            life_expectancy = 70,
            pop_scale_param_inst = pop_scale_params(year=2009, population=2891150, age_min_inc = 15, age_max_inc = 64))
run_report( root_path = "/gpfs/data/bershteynlab/EMOD/kaftad01/202306_Zimbabwe_NewEXE/output_2023_06_22_04_04", 
            cost_art=94, 
            life_expectancy = 66,
            pop_scale_param_inst = pop_scale_params(year=2020, population=9140126, age_min_inc = 15, age_max_inc=64))
