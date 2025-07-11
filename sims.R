library(dplyr)
library(tidyr)
library(ggplot2)
library(brms)
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
  slope_prior = 3,
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
  max_overdose_prob     = 0.80,
  overdose_eval_threshold = .16
)

#### look at power at 180 day estimates across simulations 

###### Plot MTD Estimates over Trial duration
plot_mtd_over_time(results)


#### Plot estimated 
plot_est_dr(results)


#############################################
############## Run in Parallel ##############
#############################################

library(parallel)

# Set number of simulations and cores
n_sims   <- 10
n_cores  <- detectCores() - 1

# Define all `probs` vectors and patient counts
probs_list <- list(
  c(0.13, 0.18, 0.21, 0.25, 0.30),
  c(0.11, 0.15, 0.175, 0.20, 0.25),
  c(0.10, 0.12, 0.15, 0.19, 0.25),
  c(0.10, 0.115, 0.135, 0.15, 0.18),
  c(0.10, 0.112, 0.128, 0.14, 0.15)
)
patient_counts <- c(40, 46)

# Output directory
output_dir <- "sim_outputs"
dir.create(output_dir, showWarnings = FALSE)

# Loop over all scenario combinations
for (probs_idx in seq_along(probs_list)) {
  for (n_patients in patient_counts) {
    
    message("Running scenario: probs=", probs_idx, " | n=", n_patients)
    
    # Pull fixed probs for this scenario
    probs <- probs_list[[probs_idx]]
    
    # Define wrapper for parallelized simulation
    run_single_trial <- function(sim_id) {
      tryCatch({
        message("  Sim ", sim_id)
        
        dlt_model <- true_model(
          doses     = c(0.1, 0.13, 0.15, 0.17, 0.20),
          probs     = probs,
          target_dlt_rate = 0.15,
          sigma_re  = 0.3,
          plot      = FALSE,
          dose_scale_factor = 10,
          n_patients_spaghetti = n_patients
        )
        
        visits <- sim_data(n_patients = n_patients,
                           mean_enrollment_gap = 14,
                           mean_visits = 11.7)
        
        prior_info <- get_prior(
          doses = c(0.1, 0.13, 0.15, 0.17, 0.20),
          probs = c(0.1, 0.12, 0.15, 0.19, 0.25),  # fixed prior
          true_mtd = 0.15,
          true_dlt_rate = 0.15,
          int_prior = 0.3,
          slope_prior = 0.08,
          link = "logit",
          dose_scale_factor = 10,
          plot = FALSE,
          force_through_mtd = TRUE
        )
        
        result <- run_trial(
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
        
        # Keep only essential output
        list(
          mtd_estimates       = result$mtd_estimates,
          final_dose_response = result$final_dose_response,
          true_mtd            = result$true_mtd
        )
        
      }, error = function(e) {
        message("    Sim ", sim_id, " failed: ", conditionMessage(e))
        return(NULL)
      })
    }
    
    # Run simulations in parallel
    sim_results <- mclapply(1:n_sims, run_single_trial, mc.cores = n_cores)
    
    # Save compiled results list for this scenario
    outfile <- file.path(output_dir, paste0("results_probs", probs_idx, "_n", n_patients, ".rds"))
    saveRDS(sim_results, outfile)
  }
}





# Extract final MTDs from each result
final_mtds <- sapply(results_list, function(res) {
  if (is.null(res)) return(NA)  # in case of try-catch failures
  tail(res$mtd_estimates$mtd, 1)
})


results_list[5]
# Tabulate distribution
mtd_distribution <- table(final_mtds, useNA = "ifany")

# Convert to data frame
mtd_df <- as.data.frame(mtd_distribution)
names(mtd_df) <- c("Dose", "Count")

# Add proportions
mtd_df$Proportion <- mtd_df$Count / sum(mtd_df$Count, na.rm = TRUE)

# View results
print(mtd_df)







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
