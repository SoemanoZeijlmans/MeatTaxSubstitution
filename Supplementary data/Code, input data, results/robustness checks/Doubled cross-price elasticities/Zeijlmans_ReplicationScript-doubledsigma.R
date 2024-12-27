# Replication script for tax-induced animal product substitution and multiple externalities
# School of Business and Economics, Vrije Universiteit Amsterdam
# (c) Soemano Zeijlmans, June 2024
# Contact: soem@pm.me

rm(list=ls()) #Empty environment.
cat("\014")
cat ("Replication script result for tax-induced animal product substitution and multiple externalities, MSc STREEM VU, Soemano Zeijlmans, June 2024 \n")

## Setup ----
# Note: Choose UTF-8 encoding when prompted.
# Note: specify tax type and level by scrolling to the bottom of this R script!
# Load necessary library
library(nleqslv)  # For solving non-linear equations. Information at CRAN: https://cran.r-project.org/web/packages/nleqslv/index.html
library(ggplot2)  # For making the threshold plot

# Disable scientific number notation
options(scipen = 999)

# Create output folders
if (!dir.exists("output")) {
  dir.create("output")
}
if (!dir.exists("output/scenarios")) {
  dir.create("output/scenarios")
}
if (!dir.exists("output/thresholds")) {
  dir.create("output/thresholds")
}
if (!dir.exists("intermediate-files")) {
  dir.create("intermediate-files")
}

# Import the data
data <- read.csv("Zeijlmans_InputData_doubled-sigma.csv", sep = ",")  
# Make sure to set the working directory to the folder of this file.
# If the R script and the dataset are stored in the same folder, and you are using RStudio:
# Menu -> Session -> Set Working Directory -> To Source File Location

# Initial guesses for prices from the CSV file. The algorithm needs initial estimates to work. 
# These initial guesses are the pre-tax estimates for improved calculation efficiency.
initial_guesses <- data$P_eq_0

## Calculate scaling parameters for demand and supply ----

# Calculate L_D and append to the data frame
# First, write L_D this as a function.
calculate_L_D <- function(P_eq_0, Q_eq_0, epsilon, sigma) {
  n <- length(P_eq_0)
  L_D <- numeric(n)
  names(L_D) <- names(P_eq_0)
  for (i in 1:n) {
    j <- setdiff(1:n, i)
    product_term_initial <- prod(P_eq_0[j]^sigma[i, j])
    L_D[i] <- Q_eq_0[i] / ((P_eq_0[i]^epsilon[i]) * product_term_initial)
  }
  return(L_D)
}

# Append L_D to df
data$L_D <- calculate_L_D(
  P_eq_0 = data$P_eq_0,
  Q_eq_0 = data$Q_eq_0,
  epsilon = data$epsilon,
  sigma = as.matrix(data[, grep("^sigma_", names(data))])
)

# Find value for L_S (no function needed, simpler calculation), and
# append L_S to df
data$L_S <- data$Q_eq_0/data$P_eq_0^data$eta



## Define functions ----

# Find equilibrium prices after tax
equilibrium_function <- function(P, data) {
  equations <- numeric(nrow(data))  # Initialize vector to hold equations
  for (i in 1:nrow(data)) {
    L_D <- data$L_D[i]  # Demand parameter for product i
    L_S <- data$L_S[i]  # Supply parameter for product i
    tau <- data$tau[i]  # Specific tau value for product i
    epsilon <- data$epsilon[i]  # Price elasticity of demand for product i
    eta <- data$eta[i]  # Price elasticity of supply for product i
    # Construct a product term for substitutes
    product_term <- 1
    for (j in 1:nrow(data)) {
      if (j != i) {
        sigma <- data[[paste0("sigma_", data$product[j])]][i]  # Cross-price elasticity
        product_term <- product_term * P[j]^sigma
      }
    }
    # Construct the equation for product i
    equations[i] <- L_D * (P[i])^epsilon * product_term - L_S * (P[i] - tau)^eta
  }
  return(equations)
}


# Define the demand function
calculate_Q_D <- function(P_D, L_D, epsilon, sigma) {
  n <- length(P_D)
  Q_D <- numeric(n)
  names(Q_D) <- names(P_D)
  for (i in 1:n) {
    product_term <- 1
    for (j in setdiff(1:n, i)) {
      product_term <- product_term * (P_D[j]^sigma[i, j])
    }
    Q_D[i] <- L_D[i] * (P_D[i]^epsilon[i]) * product_term
  }
  return(Q_D)
}

# Define a function for calculating total carbon emissions
calculate_E <- function(Q, CI) {
  n <- length(Q)
  E <- numeric(n)
  names(E) <- names(Q)
  for (i in 1:n) {
    E[i] <- CI[i]*Q[i]
  }
  return(E)
} #Calculates emissions for given consumption level

calculate_A <- function(Q, AI, type) {
  n <- length(Q)
  A <- numeric(n)
  names(A) <- names(Q)
  if (type == "DD") {
    AI = data$deaths_direct
  } else if (type == "TD") {AI = data$deaths_total
  } else if (type == "DT") {AI = data$time_direct
  } else if (type == "TT") {AI = data$time_total
  }
  for (i in 1:n) {
    A[i] <- AI[i]*Q[i]
  }
  return(A)
} #Calculates animal impacts for given consumption level

calculate_TR <- function(Q, tau) {
  n <- length(Q)
  TR <- numeric(n)
  names(TR) <- names(Q)
  for (i in 1:n) {
    TR[i] <- tau[i]*Q[i]
  }
  return(TR)
}

calculate_Theta <- function(deltaE, deltaA) {
  mapply(function(e, a) {
    if ((e > 0 && a > 0) || (e < 0 && a < 0) || e == 0 || a == 0) { # Check for edge cases
      return(NaN)
    } else {
      Theta <- abs(e) / abs(a) / 1000
      return(Theta)
    }
  }, deltaE, deltaA)
}
#deltaE likely negative, deltaA likely positive for carbon taxation.
# Input is for total, not per product group. No vector input but individual numbers.
# Answers the question: how important should you find animal welfare relative to 
#   emissions reductions for this policy to have net zero external costs?
# If both deltas are positive or negative, or one is zero, the result of this function is undefined.


## Define function for output table ----
scenario_result <- function(taxlevel, taxtype, data){
  n <- length(data$uid)
  
  # Calculate the tax level per product
  if (taxtype == "carbon") {
    tau = data$CI * taxlevel/1000 #Divide by 1000 because the carbon price is per tonne, not per kg.
  } else if (taxtype == "meat") {tau = rep(taxlevel, n)
  } else if (taxtype == "slaughter") {tau = taxlevel * data$deaths_direct
  } else if (taxtype == "no") {tau = 0
  } else stop ("ERROR: Undefined tax type")
  data$tau <- tau #Writes tax levels to the dataframe, because the equilibrium function refers to the dataframe for the tax level.
  initial_guesses <- data$P_eq_0 + tau #Update initial guesses for improved calculation efficiency
  
  # Calculate equilibrium prices (for producers)
  result <- nleqslv(initial_guesses, 
                    equilibrium_function, 
                    data = data,
                    method = "Newton")
  if (result$termcd == 1) {  # Check if the solution converged
    P_1_cons <- result$x  # Extract the equilibrium prices from result
    names(P_1_cons) <- data$product  # Assign product names to the prices
    cat("SUCCESS: Numerical solution for new equilibrium prices found using Newton's method.\n")
  } else { # When solution didn't converge
    stop("FAILURE: Error in post-tax equilibrium calculation. Failed to find equilibrium prices.")
    P_1_cons <- rep(NaN, n)
  }
  
  # Calculate new consumer prices
  P_1_prod = P_1_cons - tau
  delta_Pp_abs = P_1_prod - data$P_eq_0
  delta_Pp_rel = delta_Pp_abs / data$P_eq_0
  delta_P_abs = P_1_cons - data$P_eq_0
  delta_P_rel = delta_P_abs / data$P_eq_0
  
  # Calculate new quantity
  Q_D_1 <- calculate_Q_D(
    P_D = P_1_cons,
    L_D = data$L_D,
    epsilon = data$epsilon,
    sigma = as.matrix(data[, grep("^sigma_", names(data))])
  )
  delta_Q_abs = Q_D_1 - data$Q_eq_0
  delta_Q_rel = delta_Q_abs / data$Q_eq_0
  
  # Calculate new carbon emissions
  E_1 <- calculate_E(
    Q=Q_D_1,
    CI=data$CI
  )
  E_0 = calculate_E(data$Q_eq_0, data$CI)
  delta_E_abs = E_1 - E_0
  delta_E_rel = delta_E_abs/E_0
  
  # Calculate new animal welfare impacts: direct deaths
  AW_DD_1 <- calculate_A(
    Q = Q_D_1,
    AI = data$deaths_direct,
    type = "DD"
  )
  AW_DD_0 = calculate_A(
    Q = data$Q_eq_0,
    AI = data$deaths_direct,
    type = "DD"
  )
  delta_AW_DD_abs = AW_DD_1 - AW_DD_0
  delta_AW_DD_rel = delta_AW_DD_abs / AW_DD_0
  # Calculate new animal welfare impacts: total deaths
  AW_TD_1 <- calculate_A(
    Q = Q_D_1,
    AI = data$deaths_total,
    type = "TD"
  )
  AW_TD_0 = calculate_A(
    Q = data$Q_eq_0,
    AI = data$deaths_direct,
    type = "TD"
  )
  delta_AW_TD_abs = AW_TD_1 - AW_TD_0
  delta_AW_TD_rel = delta_AW_TD_abs / AW_TD_0
  # Calculate new animal welfare impacts: direct time
  AW_DT_1 <- calculate_A(
    Q = Q_D_1,
    AI = data$time_direct,
    type = "DT"
  )
  AW_DT_0 = calculate_A(
    Q = data$Q_eq_0,
    AI = data$deaths_direct,
    type = "DT"
  )
  delta_AW_DT_abs = AW_DT_1 - AW_DT_0
  delta_AW_DT_rel = delta_AW_DT_abs / AW_DT_0
  # Calculate new animal welfare impacts: total time
  AW_TT_1 <- calculate_A(
    Q = Q_D_1,
    AI = data$time_total,
    type = "TT"
  )
  AW_TT_0 = calculate_A(
    Q = data$Q_eq_0,
    AI = data$deaths_direct,
    type = "TT"
  )
  delta_AW_TT_abs = AW_TT_1 - AW_TT_0
  delta_AW_TT_rel = delta_AW_TT_abs / AW_TT_0
  
  # Calculate threshold values (wrong, takes abs number not vector! Always returns NaN)
  Theta_DD <- calculate_Theta( #Direct deaths
    deltaE = delta_E_abs,
    deltaA = delta_AW_DD_abs
  )
  Theta_TD <- calculate_Theta( #Total deaths
    deltaE = delta_E_abs,
    deltaA = delta_AW_TD_abs
  )
  Theta_DT <- calculate_Theta( #Direct time
    deltaE = delta_E_abs,
    deltaA = delta_AW_DT_abs
  )
  Theta_TT <- calculate_Theta( #Total time
    deltaE = delta_E_abs,
    deltaA = delta_AW_TT_abs
  )
  
  # Calculate tax revenue
  TR <- calculate_TR(
    Q=Q_D_1,
    tau=tau
  )
  
  # Make output data frame for individual products
  ScenarioResult <- data.frame(
    DepVar = c("Tax (USD/kg)", "Δ Consumer price (USD/kg)", "Δ Consumer price (%)","Δ Producer price (USD/kg)", "Δ Producer price (%)", "Δ Consumption (kg)", "Δ Consumption (%)", "Tax revenue (USD)",
               "Δ Emissions (kg CO2-eq.)", "Δ Emissions (%)",
               "Δ Direct deaths", "Δ Direct deaths (%)",
               "Δ Total deaths", "Δ Total deaths (%)",
               "Δ Direct time in agriculture (days)", "Δ Direct time in agriculture (%)",
               "Δ Total time in agriculture (days)", "Δ Total time in agriculture (%)",
               "SCAU threshold (direct deaths)", "SCAU threshold (total deaths)", "SCAU threshold (direct time)", "SCAU threshold (total time)"),
    beef = c(tau[1], delta_P_abs[1], delta_P_rel[1], delta_Pp_abs[1], delta_Pp_rel[1], delta_Q_abs[1], delta_Q_rel[1], TR[1], delta_E_abs[1], delta_E_rel[1], 
             delta_AW_DD_abs[1], delta_AW_DD_rel[1], delta_AW_TD_abs[1], delta_AW_TD_rel[1], 
             delta_AW_DT_abs[1], delta_AW_DT_rel[1], delta_AW_TT_abs[1], delta_AW_TT_rel[1],
             Theta_DD[1], Theta_TD[1], Theta_DT[1], Theta_TT[1]),
    chicken = c(tau[2], delta_P_abs[2], delta_P_rel[2], delta_Pp_abs[2], delta_Pp_rel[2], delta_Q_abs[2], delta_Q_rel[2], TR[2], delta_E_abs[2], delta_E_rel[2], 
                delta_AW_DD_abs[2], delta_AW_DD_rel[2], delta_AW_TD_abs[2], delta_AW_TD_rel[2], 
                delta_AW_DT_abs[2], delta_AW_DT_rel[2], delta_AW_TT_abs[2], delta_AW_TT_rel[2],
                Theta_DD[2], Theta_TD[2], Theta_DT[2], Theta_TT[2]),  
    pork = c(tau[3], delta_P_abs[3], delta_P_rel[3], delta_Pp_abs[3], delta_Pp_rel[3], delta_Q_abs[3], delta_Q_rel[3], TR[3], delta_E_abs[3], delta_E_rel[3], 
             delta_AW_DD_abs[3], delta_AW_DD_rel[3], delta_AW_TD_abs[3], delta_AW_TD_rel[3], 
             delta_AW_DT_abs[3], delta_AW_DT_rel[3], delta_AW_TT_abs[3], delta_AW_TT_rel[3],
             Theta_DD[3], Theta_TD[3], Theta_DT[3], Theta_TT[3]),
    farmedfish = c(tau[4], delta_P_abs[4], delta_P_rel[4], delta_Pp_abs[4], delta_Pp_rel[4], delta_Q_abs[4], delta_Q_rel[4], TR[4], delta_E_abs[4], delta_E_rel[4], 
                   delta_AW_DD_abs[4], delta_AW_DD_rel[4], delta_AW_TD_abs[4], delta_AW_TD_rel[4], 
                   delta_AW_DT_abs[4], delta_AW_DT_rel[4], delta_AW_TT_abs[4], delta_AW_TT_rel[4],
                   Theta_DD[4], Theta_TD[4], Theta_DT[4], Theta_TT[4])
  )
  rownames(ScenarioResult) <- ScenarioResult$DepVar
  ScenarioResult <- ScenarioResult[, -1]
  
  # Add a column for totals
  ScenarioResult$total <- NA
  # Sum up columns that can simply be summed up
  rows_to_sum <- c("Δ Consumption (kg)", "Tax revenue (USD)", "Δ Emissions (kg CO2-eq.)", "Δ Direct deaths", 
                   "Δ Total deaths", "Δ Direct time in agriculture (days)", "Δ Total time in agriculture (days)") #Which columns should be summed?
  for (row_name in rows_to_sum) {
    ScenarioResult[row_name, "total"] <- sum(ScenarioResult[row_name, c("beef", "chicken", "pork", "farmedfish")], na.rm = TRUE)
  }
  # Calculate columns that cannot simply be summed up
  delta_Q_rel_total <- sum(delta_Q_abs) / sum(data$Q_eq_0)
  ScenarioResult["Δ Consumption (%)", "total"] <- delta_Q_rel_total
  delta_E_rel_total <- sum(delta_E_abs) / sum(data$Q_eq_0 * data$CI)
  ScenarioResult["Δ Emissions (%)", "total"] <- delta_E_rel_total
  # Same, but for animal welfare indicators:
  delta_AW_DD_rel_total <- sum(delta_AW_DD_abs) / sum(data$Q_eq_0 * data$deaths_direct)
  ScenarioResult["Δ Direct deaths (%)", "total"] <- delta_AW_DD_rel_total
  delta_AW_TD_rel_total <- sum(delta_AW_TD_abs) / sum(data$Q_eq_0 * data$deaths_total)
  ScenarioResult["Δ Total deaths (%)", "total"] <- delta_AW_TD_rel_total
  delta_AW_DT_rel_total <- sum(delta_AW_DT_abs) / sum(data$Q_eq_0 * data$time_direct)
  ScenarioResult["Δ Direct time in agriculture (%)", "total"] <- delta_AW_DT_rel_total
  delta_AW_TT_rel_total <- sum(delta_AW_TT_abs) / sum(data$Q_eq_0 * data$time_total)
  ScenarioResult["Δ Total time in agriculture (%)", "total"] <- delta_AW_TT_rel_total
  # Same, but for the threshold values:
  #direct deaths
  Theta_DD_total <- calculate_Theta(
    deltaE = sum(delta_E_abs),
    deltaA = sum(delta_AW_DD_abs)
  )
  ScenarioResult["SCAU threshold (direct deaths)", "total"] <- Theta_DD_total
  #total deaths
  Theta_TD_total <- calculate_Theta(
    deltaE = sum(delta_E_abs),
    deltaA = sum(delta_AW_TD_abs)
  )
  ScenarioResult["SCAU threshold (total deaths)", "total"] <- Theta_TD_total
  #direct time in ag
  Theta_DT_total <- calculate_Theta(
    deltaE = sum(delta_E_abs),
    deltaA = sum(delta_AW_DT_abs)
  )
  ScenarioResult["SCAU threshold (direct time)", "total"] <- Theta_DT_total
  #total time in ag
  Theta_TT_total <- calculate_Theta(
    deltaE = sum(delta_E_abs),
    deltaA = sum(delta_AW_TT_abs)
  )
  ScenarioResult["SCAU threshold (total time)", "total"] <- Theta_TT_total
  
  #Return output table for this scenario
  return (ScenarioResult)
  
}

## Calculate and export specified scenarios ----

#Set tax scenarios to be exported
carbontaxes = c(10, 25, 50)
meattaxes = c(1, 2.5, 5)
slaughtertaxes = c(1, 2.5, 5)
taxtypes = c("carbon", "meat", "slaughter")

cat("\n----    PART 1: MAKING TABLES FOR SPECIFIED TAX SCENARIOS.   ----\n\n")

for (taxtype in taxtypes) {
  if (taxtype == "carbon") {
    taxlevels = carbontaxes
  } else if (taxtype == "meat") {
    taxlevels = meattaxes
  } else if (taxtype == "slaughter") {
    taxlevels = slaughtertaxes
  }
  cat ("CALCULATING", toupper(taxtype), "TAX RESULTS.\n\n")
  for (taxlevel in taxlevels) {
    cat("Calculating", taxtype, "tax at level", taxlevel, "\n")
    ScenarioName <- paste("Scenario -", taxtype, "tax of", taxlevel, "USD")
    FileName <- paste("output/scenarios/", ScenarioName, ".csv", sep="")
    Scenario <- scenario_result(taxlevel = taxlevel, taxtype = taxtype, data=data)
    write.csv(Scenario, file = FileName, row.names=TRUE)
    cat("File exported:", FileName, "\n")
  }
  cat("\n")
}


cat("Tax scenario tables exported.\n")

## Find threshold values ----
cat("\n----    PART 2: FINDING THRESHOLD VALUES FOR SCAU    ----\n")

TRandExt_list_carbon <- list()
TRandExt_list_slaughter <- list()

TRandExt_function <- function(taxlevel, taxtype, data){
  n <- length(data$uid)
  # Calculate the tax level per product
  if (taxtype == "carbon") {
    tau = data$CI * taxlevel/1000 #Divide by 1000 because the carbon price is per tonne, not per kg.
  } else if (taxtype == "meat") {tau = rep(taxlevel, n)
  } else if (taxtype == "slaughter") {tau = taxlevel * data$deaths_direct
  } else if (taxtype == "no") {tau = 0
  } else stop ("ERROR: Undefined tax type")
  data$tau <- tau #Writes tax levels to the dataframe, because the equilibrium function refers to the dataframe for the tax level.
  initial_guesses <- data$P_eq_0 + tau #Update initial guesses for improved calculation efficiency
  
  # Calculate equilibrium prices (for producers)
  result <- nleqslv(initial_guesses, 
                    equilibrium_function, 
                    data = data,
                    method = "Newton")
  if (result$termcd == 1) {  # Check if the solution converged
    P_1_cons <- result$x  # Extract the equilibrium prices from result
    names(P_1_cons) <- data$product  # Assign product names to the prices
    #cat("SUCCESS: Numerical solution for new equilibrium prices found using Newton's method.\n")
  } else { # When solution didn't converge
    #cat("FAILURE: Error in post-tax equilibrium calculation. Failed to find equilibrium prices. Scenario: $", taxlevel, " USD ", taxtype, " tax.")
    P_1_cons <- result$x
    names(P_1_cons) <- data$product
  }
  
  # Calculate new consumer prices
  P_1_prod = P_1_cons - tau
  delta_P_abs = P_1_cons - data$P_eq_0
  delta_P_rel = delta_P_abs / data$P_eq_0
  
  # Calculate new quantity
  Q_D_1 <- calculate_Q_D(
    P_D = P_1_cons,
    L_D = data$L_D,
    epsilon = data$epsilon,
    sigma = as.matrix(data[, grep("^sigma_", names(data))])
  )
  delta_Q_abs = Q_D_1 - data$Q_eq_0
  delta_Q_rel = delta_Q_abs / data$Q_eq_0
  
  # Calculate new carbon emissions
  E_1 <- calculate_E(
    Q=Q_D_1,
    CI=data$CI
  )
  E_0 = calculate_E(data$Q_eq_0, data$CI)
  delta_E_abs = E_1 - E_0
  delta_E_rel = delta_E_abs/E_0
  
  # Calculate new animal welfare impacts: direct deaths
  AW_DD_1 <- calculate_A(
    Q = Q_D_1,
    AI = data$deaths_direct,
    type = "DD"
  )
  AW_DD_0 = calculate_A(
    Q = data$Q_eq_0,
    AI = data$deaths_direct,
    type = "DD"
  )
  delta_AW_DD_abs = AW_DD_1 - AW_DD_0
  delta_AW_DD_rel = delta_AW_DD_abs / AW_DD_0
  
  # Calculate tax revenue
  TR <- calculate_TR(
    Q=Q_D_1,
    tau=tau
  )
  
  TRtotal = sum(TR)
  delta_E_abs_total = sum(delta_E_abs)
  delta_AW_DD_abs_total = sum(delta_AW_DD_abs)

  #Return output vector #changed
  if (result$termcd == 1) {
    return (c(taxlevel, TRtotal, delta_E_abs_total, delta_AW_DD_abs_total))
  } else return (c(taxlevel, NaN, NaN, NaN))
}

cat("\nCreating table for tax revenue of carbon taxation... (This may take a while.)\n")
#carbon tax table
carbonsteps <- c(seq(0.01, 15, by=0.01), seq(15.01, 30, 0.01), seq(30.1, 150, 0.1) ) #set iterations
i<-0
for (taxlevel in carbonsteps) {
  TRandExt_vector <- TRandExt_function(taxlevel, "carbon", data=data)
  TRandExt_list_carbon[[length(TRandExt_list_carbon)+1]] <- TRandExt_vector
  i <- i+1
  cat("\r", round(i/length(carbonsteps)*100), "% completed.", sep="")
}
df_TRandExt_carbon <- do.call(rbind, TRandExt_list_carbon)
colnames(df_TRandExt_carbon) <- c("CarbonTaxLevel", "CarbonTaxRevenue", "CarbonTaxDeltaCarbon", "CarbonTaxDeltaDD")
df_TRandExt_carbon <- as.data.frame(df_TRandExt_carbon)
write.csv(df_TRandExt_carbon, file = "intermediate-files/TaxRevenuesCarbonTax.csv")
cat("\nCarbon tax revenue table exported.\n")

cat("\nCreating table for tax revenue of slaughter taxation... (This may take a while.)\n")
#slaughter tax table
slaughtersteps <- c(seq(0.001, 1, by=0.001), seq(1.001, 3, 0.001), seq(3.01, 30, 0.01) )
i <- 0
for (taxlevel in slaughtersteps) {
  TRandExt_vector <- TRandExt_function(taxlevel, "slaughter", data=data)
  TRandExt_list_slaughter[[length(TRandExt_list_slaughter)+1]] <- TRandExt_vector
  i <- i+1
  cat("\r", round(i/length(slaughtersteps)*100), "% completed.", sep="")
}
df_TRandExt_slaughter <- do.call(rbind, TRandExt_list_slaughter)
colnames(df_TRandExt_slaughter) <- c("SlaughterTaxLevel", "SlaughterTaxRevenue", "SlaughterTaxDeltaCarbon", "SlaughterTaxDeltaDD")
df_TRandExt_slaughter <- as.data.frame(df_TRandExt_slaughter)
write.csv(df_TRandExt_slaughter, file = "intermediate-files/TaxRevenuesSlaughterTax.csv")
cat("\nSlaughter tax revenue table exported.\n")

#Initialise new DataFrame for matched values.
TaxRevenueDF <- data.frame(TR = seq(1, 250, 10)) #This one is for plotting because TR=1 is the first sensible result.
TaxRevenueDFround <- data.frame(TR = seq(0, 250, 10)) #This one is for the table because round values are better than 11, 21, 31, etc.

# Function to find the closest value
find_closest <- function(target, vec) {
  vec[which.min(abs(vec - target))]
}

cat("\nMatching carbon tax revenues with slaughter tax revenues...\n")

# Create new columns for the closest matches for the figure
# Add columns for the closest matches from climate data
TaxRevenueDF$ClosestCarbonTaxRevenue <- sapply(TaxRevenueDF$TR, function(x) find_closest(x, df_TRandExt_carbon$CarbonTaxRevenue))
TaxRevenueDF$CarbonTaxLevel <- sapply(TaxRevenueDF$ClosestCarbonTaxRevenue, function(x) df_TRandExt_carbon$CarbonTaxLevel[which(df_TRandExt_carbon$CarbonTaxRevenue == x)])
TaxRevenueDF$CarbonTaxDeltaCarbon <- sapply(TaxRevenueDF$ClosestCarbonTaxRevenue, function(x) df_TRandExt_carbon$CarbonTaxDeltaCarbon[which(df_TRandExt_carbon$CarbonTaxRevenue == x)])
TaxRevenueDF$CarbonTaxDeltaDD <- sapply(TaxRevenueDF$ClosestCarbonTaxRevenue, function(x) df_TRandExt_carbon$CarbonTaxDeltaDD[which(df_TRandExt_carbon$CarbonTaxRevenue == x)])

# Add columns for the closest matches from slaughter data
TaxRevenueDF$ClosestSlaughterTaxRevenue <- sapply(TaxRevenueDF$TR, function(x) find_closest(x, df_TRandExt_slaughter$SlaughterTaxRevenue))
TaxRevenueDF$SlaughterTaxLevel <- sapply(TaxRevenueDF$ClosestSlaughterTaxRevenue, function(x) df_TRandExt_slaughter$SlaughterTaxLevel[which(df_TRandExt_slaughter$SlaughterTaxRevenue == x)])
TaxRevenueDF$SlaughterTaxDeltaCarbon <- sapply(TaxRevenueDF$ClosestSlaughterTaxRevenue, function(x) df_TRandExt_slaughter$SlaughterTaxDeltaCarbon[which(df_TRandExt_slaughter$SlaughterTaxRevenue == x)])
TaxRevenueDF$SlaughterTaxDeltaDD <- sapply(TaxRevenueDF$ClosestSlaughterTaxRevenue, function(x) df_TRandExt_slaughter$SlaughterTaxDeltaDD[which(df_TRandExt_slaughter$SlaughterTaxRevenue == x)])

# Create new columns for the closest matches for the table
# Add columns for the closest matches from climate data
TaxRevenueDFround$ClosestCarbonTaxRevenue <- sapply(TaxRevenueDFround$TR, function(x) find_closest(x, df_TRandExt_carbon$CarbonTaxRevenue))
TaxRevenueDFround$CarbonTaxLevel <- sapply(TaxRevenueDFround$ClosestCarbonTaxRevenue, function(x) df_TRandExt_carbon$CarbonTaxLevel[which(df_TRandExt_carbon$CarbonTaxRevenue == x)])
TaxRevenueDFround$CarbonTaxDeltaCarbon <- sapply(TaxRevenueDFround$ClosestCarbonTaxRevenue, function(x) df_TRandExt_carbon$CarbonTaxDeltaCarbon[which(df_TRandExt_carbon$CarbonTaxRevenue == x)])
TaxRevenueDFround$CarbonTaxDeltaDD <- sapply(TaxRevenueDFround$ClosestCarbonTaxRevenue, function(x) df_TRandExt_carbon$CarbonTaxDeltaDD[which(df_TRandExt_carbon$CarbonTaxRevenue == x)])

# Add columns for the closest matches from slaughter data
TaxRevenueDFround$ClosestSlaughterTaxRevenue <- sapply(TaxRevenueDFround$TR, function(x) find_closest(x, df_TRandExt_slaughter$SlaughterTaxRevenue))
TaxRevenueDFround$SlaughterTaxLevel <- sapply(TaxRevenueDFround$ClosestSlaughterTaxRevenue, function(x) df_TRandExt_slaughter$SlaughterTaxLevel[which(df_TRandExt_slaughter$SlaughterTaxRevenue == x)])
TaxRevenueDFround$SlaughterTaxDeltaCarbon <- sapply(TaxRevenueDFround$ClosestSlaughterTaxRevenue, function(x) df_TRandExt_slaughter$SlaughterTaxDeltaCarbon[which(df_TRandExt_slaughter$SlaughterTaxRevenue == x)])
TaxRevenueDFround$SlaughterTaxDeltaDD <- sapply(TaxRevenueDFround$ClosestSlaughterTaxRevenue, function(x) df_TRandExt_slaughter$SlaughterTaxDeltaDD[which(df_TRandExt_slaughter$SlaughterTaxRevenue == x)])



cat("Calculating threshold values...\n")
TaxRevenueDF$SCAUThreshold <- (TaxRevenueDF$CarbonTaxDeltaCarbon - TaxRevenueDF$SlaughterTaxDeltaCarbon) / (TaxRevenueDF$SlaughterTaxDeltaDD - TaxRevenueDF$CarbonTaxDeltaDD) / 1000
TaxRevenueDFround$SCAUThreshold <- (TaxRevenueDFround$CarbonTaxDeltaCarbon - TaxRevenueDFround$SlaughterTaxDeltaCarbon) / (TaxRevenueDFround$SlaughterTaxDeltaDD - TaxRevenueDFround$CarbonTaxDeltaDD) / 1000
  #note: Divide by 1000, because the carbon emissions column is in kg, and the SCC is in tonnes.

#Export intermediate file
write.csv(TaxRevenueDF, file = "intermediate-files/MatchedTaxRevenues.csv")

cat("Making and exporting simplified table...\n")
# Make a simple shorter table
TaxRevenueDFSimple <- TaxRevenueDFround[TaxRevenueDFround$TR %% 20 == 0, ]
TaxRevenueDFSimple <- subset(TaxRevenueDFSimple, select = -c(ClosestCarbonTaxRevenue, ClosestSlaughterTaxRevenue, CarbonTaxDeltaCarbon, CarbonTaxDeltaDD, SlaughterTaxDeltaCarbon, SlaughterTaxDeltaDD))

# Export simple table as CSV
write.csv(TaxRevenueDFSimple, file = "output/thresholds/SimplifiedThresholdTable.csv")

# Make a plot
cat("Exporting figure with threshold values...\n")
png("output/thresholds/ThresholdValuesSCAU.png", width = 1800, height = 1800, res=300) #Specify the file

windowsFonts(A = windowsFont("Garamond"))
plot(TaxRevenueDF$TR, TaxRevenueDF$SCAUThreshold, type = "b", pch = 19, col = "black", 
     xlab = "Tax revenue (TR, USD/cap./annum)", 
     ylab = expression(SCAU~Threshold~(theta[SCAU]*", relative to SCC")), 
     main = "Threshold values for the\nSocial Cost of Animal Use (SCAU)",
     family = "A",
     font = 1
     )
text(x = max(TaxRevenueDF$TR) - 75, 
     y = max(TaxRevenueDF$SCAUThreshold) - 0.02, 
     labels = "Revenue-equivalent slaughter\ntax reduces external costs\nmore than a carbon tax.", 
     pos = 3, 
     cex = 0.8,
     family = "A",
     font = 1
     )
text(x = min(TaxRevenueDF$TR) + 60, 
     y = min(TaxRevenueDF$SCAUThreshold) + 0.015, 
     labels = "Revenue-equivalent carbon\ntax reduces external costs\nmore than a slaughter tax.", 
     pos = 3, 
     cex = 0.8,
     family = "A",
     font = 1
     )
dev.off()

## Footer ----
cat("\nEnd of script.\n")
# End of script.