// Male chela morphometry analysis, Gaussian regression mixture: 

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() (){
   // Data variable declarations:
   DATA_VECTOR(y);                          // Vector of observed chela heights (mm).
   DATA_VECTOR(cw);                         // Vector of observed carapace widths (mm).
   DATA_FACTOR(station);                    // Vector of station IDs.
   
   // Parameter variable declarations:
   PARAMETER(log_alpha_immature);           // Log-linear allometric immature intercept.
   PARAMETER(log_alpha_mature);             // Log-linear allometric mature intercept.
   PARAMETER(log_beta_immature);            // Log-linear allometric immature slope.
   PARAMETER(log_beta_mature);              // Log-linear allometric mature slope.
   PARAMETER(log_sigma_immature);           // Log-scale immature error parameter.
   PARAMETER(log_sigma_mature);             // Log-scale mature error parameter.
   PARAMETER_VECTOR(eta);                   // Logit-linear slope and intercept parameters for mixture proportions.
   
   PARAMETER(log_sigma_kurtotic);           // Extra error parameter to model outliers.
   PARAMETER(logit_p_kurtotic);             // Logit-linear proportion of kurtotic data.
   
   Type sigma_kurtotic = exp(log_sigma_kurtotic);      // Extra error parameter to model outliers.
   Type p_kurtotic = 1 / (1 + exp(-logit_p_kurtotic)); // Proportion of kurtotic outliers.
   
   // Station effects:
   PARAMETER_VECTOR(station_immature_effect);                     // Vector of station random effects.
   PARAMETER(log_sigma_station_immature);                         // Log-scale station error parameter.
   Type sigma_station_immature = exp(log_sigma_station_immature); // Station effect error parameter.
   PARAMETER_VECTOR(station_mature_effect);                       // Vector of station random effects.
   PARAMETER(log_sigma_station_mature);                           // Log-scale station error parameter.
   Type sigma_station_mature = exp(log_sigma_station_mature);     // Station effect error parameter.
      
   // Calculated variables:
   Type res = 0;                            // Negative log-likelihood accumulator.
   int n_obs = y.size();                    // Number of length-frequency categories. 
   Type sigma_immature = exp(log_sigma_immature); // Immature error parameter.
   Type sigma_mature = exp(log_sigma_mature);     // Mature error parameter.
   int n_station = station_immature_effect.size();   // Number of sampling stations.
   
   // Prior over station effects:
   for (int j = 0; j < n_station; j++){
      res -= dnorm(station_immature_effect[j], Type(0.0), sigma_station_immature, true);
   }
   for (int j = 0; j < n_station; j++){
      res -= dnorm(station_mature_effect[j], Type(0.0), sigma_station_mature, true);
   }
      
   // Likelihood functions:
   Type mu_immature = 0;
   Type mu_mature = 0;
   Type p = 0;
   for (int i = 0; i < n_obs; i++){  
      // Proportion of mature crab:
      p = Type(1) / (Type(1) + exp(-eta[0] -eta[1] * cw[i]));  
      
      // Allometric immature relation:
      mu_immature = log_alpha_immature + log_beta_immature * log(cw[i]) + station_immature_effect[station[i]-1];
      
      // Allometric mature relation:
      mu_mature = log_alpha_mature + log_beta_mature * log(cw[i]) + station_mature_effect[station[i]-1];
      
      // Gaussian kurtotic mixture density:
      res -= log((1-p) * ((1-p_kurtotic) * dnorm(log(y[i]), mu_immature, sigma_immature, false) +                 
                            (p_kurtotic) * dnorm(log(y[i]), mu_immature, sigma_immature + sigma_kurtotic, false)) +
                   (p) * ((1-p_kurtotic) * dnorm(log(y[i]), mu_mature, sigma_mature, false) +
                            (p_kurtotic) * dnorm(log(y[i]), mu_mature, sigma_mature + sigma_kurtotic, false)));
   }   
  
   return res;
}  

  