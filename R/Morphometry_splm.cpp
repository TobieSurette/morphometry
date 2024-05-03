// Male chela morphometry analysis, Gaussian regression mixture: 

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() (){
   // Data variable declarations:
   DATA_VECTOR(y);                              // Vector of observed chela heights (mm).
   DATA_VECTOR(cw);                             // Vector of observed carapace widths (mm).
   DATA_FACTOR(station);                        // Vector of station IDs.
   DATA_VECTOR(cw_imm);                         // Vector of observed carapace widths (mm) for small immature crab.
   
   // Parameter variable declarations:
   PARAMETER(log_alpha_immature);               // Log-linear allometric immature intercept.
   PARAMETER_VECTOR(log_beta_immature);         // Log-linear allometric immature slope.
   PARAMETER(log_alpha_mature);                 // Log-linear allometric mature intercept.   
   PARAMETER(log_beta_mature);                  // Log-linear allometric mature slope.
   PARAMETER(log_precision_immature);            // Precision parameter for immature SPL model.
   Type precision_immature = exp(log_precision_immature);
   PARAMETER(transition_immature);              // Transition point for immature SPL model.

   PARAMETER_VECTOR(station_effect_immature);   // Vector of station random effects for immatures.
   PARAMETER(log_sigma_station_immature);       // Log-scale station error parameter for immatures.
   PARAMETER_VECTOR(station_effect_mature);     // Vector of station random effects for matures.
   PARAMETER(log_sigma_station_mature);         // Log-scale station error parameter for matures.
      
   // Chela error parameters:       
   PARAMETER(log_sigma);                        // Log-scale common error parameter.
   Type sigma = exp(log_sigma);                 // Chela height error parameter.
   PARAMETER(log_sigma_outlier);                // Log-scale addiotnal variability of data outliers.
   Type sigma_outlier = exp(log_sigma_outlier); // Additional variability of data outliers.

   // Outlier proportions parameters:
   PARAMETER(alpha_outlier);                    // Logit-scale mean outlier rate.
     
   // Maturity proportions model:
   PARAMETER(eta_alpha);                        // Logit-linear intercept for mixture proportions.
   PARAMETER_VECTOR(eta_beta);                  // Logit-linear slopes for mixture proportions.
   PARAMETER_VECTOR(eta_transition);            // Logit-linear transition sizes for mixture proportions.
   PARAMETER(log_eta_precision);                // Logit-linear transition precision for mixture proportions.
   Type eta_precision = exp(log_eta_precision); // Logit-linear transition precision for mixture proportions.
      
   // Calculated variables:
   Type res = 0;                                   // Negative log-likelihood accumulator.
   int n_obs = y.size();                           // Number of length-frequency categories. 
   int n_station = station_effect_immature.size(); // Number of sampling stations.
   int n_imm = cw_imm.size();                      // Number of sampling stations.
   
   // Prior over station effects:
   /// Type sigma_station = exp(log_sigma_station); // Station effect error parameter.
   //for (int j = 0; j < n_station; j++){
   //   res -= dnorm(station_effect[j], Type(0.0), sigma_station, true);
   //}
   
   // Prior over station effects:
   Type sigma_station_immature = exp(log_sigma_station_immature); // Station effect error parameter for immatures.
   Type sigma_station_mature = exp(log_sigma_station_mature);     // Station effect error parameter for matures.
   for (int j = 0; j < n_station; j++){
      res -= dnorm(station_effect_immature[j], Type(0.0), sigma_station_immature, true);
      res -= dnorm(station_effect_mature[j], Type(0.0), sigma_station_mature, true);
   }
   
   // Likelihood functions:
   Type mu_immature = 0;
   Type mu_mature = 0;
   Type logit_p_mature = 0;
   Type p_mature = 0;
   Type p_outlier = 0; 
   for (int i = 0; i < n_obs; i++){  
      // Allometric relation for immatures SPL(1):
      mu_immature = log_alpha_immature + log_beta_immature[0] * log(cw[i]) +
                    precision_immature * (log_beta_immature[1] - log_beta_immature[0]) * log(1 + exp((log(cw[i]) - transition_immature) / precision_immature)) + 
                    station_effect_immature[station[i]-1]; 
      
      // Allometric relation for matures:
      mu_mature = log_alpha_mature + log_beta_mature * log(cw[i]) + station_effect_mature[station[i]-1];
      
      // Outlier component proportions:
      p_outlier = 1 / (1 + exp(-alpha_outlier));      
      
      // Proportion of mature crab:
      logit_p_mature = eta_alpha + eta_beta[0] * cw[i] + 
                       eta_precision * (eta_beta[1] - eta_beta[0]) * log(1 + exp((cw[i] - eta_transition[0]) / eta_precision)) + 
                       eta_precision * (eta_beta[2] - eta_beta[1]) * log(1 + exp((cw[i] - eta_transition[1]) / eta_precision));
      p_mature = Type(1) / (Type(1) + exp(-logit_p_mature));
      
      // Gaussian mixture density:                                            
      res -= log((1-p_mature) * (p_outlier * dnorm(log(y[i]), mu_immature, sigma + sigma_outlier, false) + 
                           (1 - p_outlier) * dnorm(log(y[i]), mu_immature, sigma, false)) +                
                   (p_mature) * (p_outlier * dnorm(log(y[i]), mu_mature, sigma + sigma_outlier, false) + 
                           (1 - p_outlier) * dnorm(log(y[i]), mu_mature, sigma, false)));    
   }   
  
   // Log-likelohood contribution of small immature crab:
   for (int i = 0; i < n_imm; i++){ 
      // Proportion of mature crab:
      logit_p_mature = eta_alpha + eta_beta[0] * cw_imm[i] + 
                       eta_precision * (eta_beta[1] - eta_beta[0]) * log(1 + exp((cw_imm[i] - eta_transition[0]) / eta_precision)) + 
                       eta_precision * (eta_beta[2] - eta_beta[1]) * log(1 + exp((cw_imm[i] - eta_transition[1]) / eta_precision));
      p_mature = Type(1) / (Type(1) + exp(-logit_p_mature));
   
      // Bernouilli log-density:
      res -= log(1-p_mature); 
   }
     
   return res;
}  

  