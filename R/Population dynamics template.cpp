


   // Data variable declarations:
   DATA_FACTOR(cw);                             // Vector of observed carapace widths (mm).
   DATA_FACTOR(year);                           // Vector of osurvey year.
   DATA_FACTOR(maturity);                       // Vector of observed crab maturities.
   DATA_FACTOR(shell);                          // Vector of observed crab shell conditions.
    
   // 

   (1-p_moult)                                  // Skip moulting probability.
   p_moult * (1-p_maturity)                     // Skip to immaturity.
   p_moult * p_maturity                         // Skip to maturity.
    
   // Growth model parameters:
   PARAMETER(alpha_immature);                   // Hiatt intercept parameter for immature to immature growth.
   PARAMETER_VECTOR(beta_immature);             // Hiatt slope parameter for immature to immature growth.  
   PARAMETER(transition_immature);              // SPLM transition parameter for immature to immature growth. 
   PARAMETER(window_immature);                  // SPLM window parameter for immature to immature growth. 
   PARAMETER_VECTOR(log_sigma_growth);          // Parameters to model Variability of
   
   PARAMETER(alpha_mature);                     // Hiatt intercept parameter for immature to mature growth.
   PARAMETER(beta_mature);                      // Hiatt slope parameter for immature to mature growth.
 
   
   PARAMETER_VECTOR(logit_p_moult);
   PARAMETER_VECTOR(logit_p_maturity);
   
   
   // Mean growth for immatures:
   mu_immature = alpha_immature + beta_immature[0] * cws[i] +  window_immature * (beta_immature[1] - beta_immature[0]) * log(1 + exp((cws[i] - transition_immature) / window_immature));
   
   // Error around immature growth:
   sigma_immature = exp(sigma_growth[0] + sigma_growth[1] * mu_immature);  
   
   
   p_moult                                      // Probability of moulting
   p_maturity                                   // Probability of moulting to maturity.  
   
     