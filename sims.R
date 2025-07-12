library(dplyr)
library(tidyr)
library(ggplot2)
library(brms)
library(rlang)

#install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
library(cmdstanr)
#install_cmdstan()
# Read and source all .R files in the "R/" directory
r_files <- list.files("R/", pattern = "\\.R$", full.names = TRUE)

# Source each file
invisible(lapply(r_files, source))

#Define doses and their corresponding true probabilities 
doses <- c(.1, .13, .15, .17, .2)
set.seed(122)

#1st thresh MTD
probs     = c(0.13, 0.18, 0.22, 0.24, 0.29)
set.seed(122)

#2nd thresh MTD
probs     = c(0.1, 0.15, 0.18, 0.22, 0.27)
set.seed(229)

#3rd thresh MTD
probs     = c(0.1, 0.12, 0.15, 0.18, 0.23)
set.seed(122)

#4th thresh MTD
probs     = c(0.1, 0.115, 0.135, 0.15, 0.18)
set.seed(122)

#5th thresh MTD
probs     = c(0.1, 0.112, 0.128, 0.14, 0.15)
set.seed(12)


### If patient fails safety discharge criteria/threshold at hour 24 an d 48 over 3 visits, they fail to qualify for analysis
### 15-20% don't play nice with MTX; rule of 3 
### still collect data on their clearance for threshold adjusted discharge 
### if at least 6 visit 
### simulate with lost patient visits vs adding. How much does it hurt us osing patients
### 



#Define the true model 
dlt_model <- true_model(doses, 
                        probs, 
                        target_dlt_rate = 0.15,
                        plot = TRUE, 
                        sigma_re = .3, 
                        n_patients_spaghetti= 40 )


 # Simulate patient data
# 1 every 2 weeks per site in the study 
visits <- sim_data(n_patients=40,
                   mean_enrollment_gap = 14,
                   mean_visits = 12)



# Get prior paramaters based on doses and probs
prior_info = get_prior(
  doses = c(0.1, .13, .15, .17, .20),
  probs = c(0.1, .12, .15, .23, .30),
  true_mtd=0.15,
  true_dlt_rate = 0.15,
  force_through_mtd=T,
  int_prior = .5,
  dose_scale_factor = 10,
  slope_prior = .5,
  link="logit",
  plot=T
)

#IRB Buffer built in 
## at least 15 visits have to have occurred and if not push back 30 days
    
# Plot Prior Distribution Dose_response curve assumed (i.e. model skeleton)
plot_prior(prior_info,
           dlt_model)


## Run Single trial 
results <- run_trial(
  visits = visits,
  prior_info=prior_info,
  dlt_model=dlt_model,
  target_dlt_rate = 0.15,
  starting_dose = .1,
  decision_interval = 60, 
  min_data_per_update=1,
  safety_lead_in = 180, #lets assume 180 
  max_overdose_prob     = 0.50,
  overdose_eval_threshold = .16
)

#### look at power at 180 day estimates across simulations 

###### Plot MTD Estimates over Trial duration
plot_mtd_over_time(results)
plot_est_dr(results)


#############################################
############## Run in Parallel ##############
#############################################
library(parallel)

# Set number of simulations and cores
n_sims   <- 1000
n_cores  <- max(1, detectCores() - 2)
sims_per_core <- ceiling(n_sims / n_cores)

# Define all `probs` vectors and patient counts
probs_list <- list(
  c(0.13, 0.18, 0.21, 0.25, 0.30),
  c(0.11, 0.15, 0.18, 0.22, 0.27),
  c(0.10, 0.12, 0.15, 0.19, 0.24),
  c(0.10, 0.115, 0.130, 0.15, 0.19),
  c(0.10, 0.112, 0.128, 0.14, 0.15)
)

patient_counts <- c(40, 46)

output_dir <- "sim_outputs"
dir.create(output_dir, showWarnings = FALSE)

# Batching helper
split_sim_ids <- function(n, k) {
  split(seq_len(n), cut(seq_len(n), breaks = k, labels = FALSE))
}

# Outer scenario loop
for (probs_idx in seq_along(probs_list)) {
  for (n_patients in patient_counts) {
    message("Running scenario: probs=", probs_idx, " | n=", n_patients)

    probs <- probs_list[[probs_idx]]
    
    # Set up cluster
    cl <- makeCluster(n_cores)
    
    # Export necessary objects and packages
    clusterExport(cl, varlist = c("true_model", "sim_data", "get_prior", "run_trial", "probs", "n_patients"), envir = environment())
    clusterEvalQ(cl, library(brms))  # add others if needed
    
    # Define batched simulation runner
    run_batch <- function(sim_ids) {
      lapply(sim_ids, function(sim_id) {
        tryCatch({
          dlt_model <- true_model(
            doses     = c(0.1, 0.13, 0.15, 0.17, 0.20),
            probs     = probs,
            target_dlt_rate = 0.15,
            sigma_re  = 0.3,
            plot      = F,
            dose_scale_factor = 10,
            n_patients_spaghetti = n_patients
          )
          
          visits <- sim_data(n_patients = n_patients,
                             mean_enrollment_gap = 14,
                             mean_visits = 11.7)
          
          prior_info <- get_prior(
            doses = c(0.1, 0.13, 0.15, 0.17, 0.20),
            probs = c(0.08, 0.125, 0.15, 0.19, 0.24),
            true_mtd = 0.15,
            true_dlt_rate = 0.15,
            int_prior = .8,
            slope_prior = 0.8,
            link = "logit",
            dose_scale_factor = 10,
            plot = T,
            force_through_mtd =T
          )
          
          result <- run_trial(
            visits = visits,
            prior_info = prior_info,
            dlt_model = dlt_model,
            target_dlt_rate = 0.15,
            starting_dose = 0.15,
            decision_interval = 60,
            min_data_per_update = 20,
            safety_lead_in = 180,
            max_overdose_prob = 0.50,
            overdose_eval_threshold = 0.16
          )
          
          list(
            mtd_estimates       = result$mtd_estimates,
            final_dose_response = result$final_dose_response,
            true_mtd            = result$true_mtd,
            final_mtd           = result$final_mtd  
          )
        }, error = function(e) {
          message("    Sim ", sim_id, " failed: ", conditionMessage(e))
          return(NULL)
        })
      })
    }
    
    # Split work among workers
    sim_batches <- split_sim_ids(n_sims, n_cores)
    batch_results <- parLapply(cl, sim_batches, run_batch)
    stopCluster(cl)
    
    # Flatten the list
    sim_results <- do.call(c, batch_results)
    
    # Save results
    outfile <- file.path(output_dir, paste0("results_probs", probs_idx, "_n", n_patients, ".rds"))
    saveRDS(sim_results, outfile)
  }
}



library(dplyr)
library(tibble)
library(purrr)
library(readr)
library(stringr)

# Define doses used
dose_levels <- c(0.1, 0.13, 0.15, 0.17, 0.20)

# Path to saved simulation result files
output_dir <- "sim_outputs"
rds_files <- list.files(output_dir, pattern = "^results_probs.*\\.rds$", full.names = TRUE)

# Function to extract scenario info from filename
parse_filename <- function(filename) {
  # Example filename: "results_probs2_n40.rds"
  tibble(
    file = filename,
    probs_idx = as.integer(str_match(basename(filename), "probs(\\d+)")[, 2]),
    n_patients = as.integer(str_match(basename(filename), "_n(\\d+)")[, 2])
  )
}

# Load and summarize one file
summarize_file <- function(file_row) {
  sims <- readRDS(file_row$file)
  tibble(
    final_mtd = map_chr(sims, ~ if (!is.null(.x)) as.character(.x$final_mtd) else NA_character_),
    probs_idx = file_row$probs_idx,
    n_patients = file_row$n_patients
  ) %>%
    filter(!is.na(final_mtd))
}

# Combine all simulations
summary_all <- rds_files %>%
  map_df(~ summarize_file(parse_filename(.x)))

# Tabulate MTD selection rates
dose_summary <- summary_all %>%
  count(probs_idx, n_patients, final_mtd) %>%
  group_by(probs_idx, n_patients) %>%
  mutate(pct_selected = n / sum(n)) %>%
  ungroup() %>%
  arrange(probs_idx, n_patients, desc(pct_selected))

# Convert dose to numeric
dose_summary <- dose_summary %>%
  mutate(final_mtd = as.numeric(final_mtd))

# Print or return summary
dose_summary


























# Define true model
dlt_model <- true_model(
  doses     = c(0.1, 0.13, 0.15, 0.17, 0.20),
  # probs     = c(0.13, 0.175, 0.2, 0.24, 0.29),
  # probs     = c(0.1, 0.15, 0.175, 0.20, 0.25),
  #probs     = c(0.1, 0.12, 0.15, 0.19, 0.25),
  # probs     = c(0.1, 0.115, 0.135, 0.15, 0.18),
  probs     = c(0.1, 0.112, 0.128, 0.14, 0.15),
  target_dlt_rate = 0.15,
  sigma_re  = 0.5,
  plot      = T,
  n_patients_spaghetti = 40,
  dose_scale_factor = 10
)

# Simulate data

visits <- sim_data(n_patients=40,
                     mean_enrollment_gap = 14,
                     mean_visits = 11)


# Get prior
prior_info <- get_prior(
  doses = c(0.1, 0.13, 0.15, 0.17, 0.20),
  probs = c(0.1, 0.12, 0.15, 0.19, 0.25),
  true_mtd = 0.15,
  true_dlt_rate = 0.15,
  int_prior = 0.3,
  force_through_mtd = T,
  slope_prior = .08,
  link = "logit",
  dose_scale_factor = 10,
  plot = T
)

# Run trial
results = run_trial(
  visits = visits,
  prior_info = prior_info,
  dlt_model = dlt_model,
  target_dlt_rate = 0.15,
  starting_dose = 0.15,
  decision_interval = 60,
  min_data_per_update = 1,
  safety_lead_in = 180,
  max_overdose_prob = 0.50,
  overdose_eval_threshold = 0.16
)

plot_mtd_over_time(results)
results$final_dose_response
plot_est_dr(results)
