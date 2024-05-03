// Observer male chela morphometry analysis: 

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() (){
   // Data variable declarations:
   DATA_VECTOR(y);        // Vector of observed chela heights (mm).
   DATA_VECTOR(cw);       // Vector of observed carapace widths (mm).
   DATA_FACTOR(observer); // Vector of observer IDs.
   DATA_FACTOR(grid);     // Vector of fishing grid IDs.
   DATA_FACTOR(week);     // Vector of fishing week IDs.
   DATA_FACTOR(trip);     // Vector of observer trip IDs.
      
   // Chela morphometry parameters:
   PARAMETER(log_alpha_immature);               // Log-linear allometric immature intercept.
   PARAMETER(log_beta_immature);                // Log-linear allometric immature slope.
   PARAMETER(log_alpha_mature);                 // Log-linear allometric mature intercept.   
   PARAMETER(log_beta_mature);                  // Log-linear allometric mature slope.
   PARAMETER(mu_observer);                      // Mean observer bias.
   PARAMETER_VECTOR(observer_effect);           // Chela height observer effects.
   PARAMETER(log_sigma_observer);               // Log-scale error parameter on observer chela effects.  
   PARAMETER_VECTOR(trip_effect);               // Observer chela deviation fishing trip random effect.
   PARAMETER(log_sigma_trip);                   // Log-scale error parameter on fishing trip effects.
      
   // Chela error parameters:       
   PARAMETER(log_sigma);                        // Log-scale common error parameter.
   PARAMETER_VECTOR(observer_error_effect);     // Observer error random effect. 
   PARAMETER(log_sigma_observer_error);         // Observer error random effect error parameter.
   PARAMETER(log_sigma_outlier);                // Log-scale addiotnal variability of data outliers.
   Type sigma_outlier = exp(log_sigma_outlier); // Additional variability of data outliers.
   
   // Outlier proportions parameters:
   PARAMETER(alpha_outlier);                    // Logit-scale mean outlier rate.
   PARAMETER_VECTOR(observer_outlier_effect);   // Observer outlier proportion effects.
   PARAMETER(log_sigma_observer_outlier);       // Log-scale error parameter on observer outlier proportion effects.
   
   // Logit-linear parameters for mature mixture proportions:
   PARAMETER(alpha_p);                          // Logit-linear intercept.
   PARAMETER(beta_p);                           // Logit-linear slope.           
   PARAMETER_VECTOR(observer_effect_p);         // Chela height observer effects on mature proportions.
   PARAMETER(log_sigma_observer_p);             // Log-scale error parameter on observer mature proportions.
   PARAMETER_VECTOR(observer_effect_beta_p);    // Chela height observer effects on mature proportions slope.
   PARAMETER(log_sigma_observer_beta_p);        // Log-scale error parameter on observer mature proportions slope.
   PARAMETER_VECTOR(grid_effect_p);             // Fishing grid effect.
   PARAMETER(log_sigma_grid_p);                 // Log-scale error parameter on fishing grid effects.
   PARAMETER_VECTOR(week_effect_p);             // Fishing week effect.
   PARAMETER(log_sigma_week_p);                 // Log-scale error parameter on fishing week effects.
 
   // Calculated variables:
   Type res = 0;                                // Negative log-likelihood accumulator.
   int n_obs = y.size();                        // Number of observations. 
   int n_observer = observer_effect.size();     // Number of observers.
   int n_grid = grid_effect_p.size();           // Number of fishing grids.
   int n_week = week_effect_p.size();           // Number of fishing weeks.
   int n_trip = trip_effect.size();             // Number of observer trips.
      
   // Priors over various observer effects:
   Type sigma_observer = exp(log_sigma_observer);                 // Observer chela bias error.
   Type sigma_observer_p = exp(log_sigma_observer_p);             // Observer maturity proportions error.
   Type sigma_observer_beta_p = exp(log_sigma_observer_beta_p);   // Observer maturity proportions error slope.
   Type sigma_observer_error = exp(log_sigma_observer_error);     // Observer chela variability error.
   Type sigma_observer_outlier = exp(log_sigma_observer_outlier); // Observer outlier proportion error.
   for (int j = 0; j < n_observer; j++){
      res -= dnorm(observer_effect[j], Type(0.0), sigma_observer, true);                 // Chela bias effects.
      res -= dnorm(observer_effect_p[j], Type(0.0), sigma_observer_p, true);             // Maturity proportions.
      res -= dnorm(observer_effect_beta_p[j], Type(0.0), sigma_observer_beta_p, true);   // Maturity proportions slope.
      res -= dnorm(observer_error_effect[j], Type(0.0), sigma_observer_error, true);     // Chela error effects.
      res -= dnorm(observer_outlier_effect[j], Type(0.0), sigma_observer_outlier, true); // Outlier proportion effects.
   }
               
   // Priors over fishing grid effects:
   Type sigma_grid_p = exp(log_sigma_grid_p);    // Error parameter on fishing grid effects.
   for (int j = 0; j < n_grid; j++){
      res -= dnorm(grid_effect_p[j], Type(0.0), sigma_grid_p, true);
   }
   
   // Priors over fishing week effects:
   Type sigma_week_p = exp(log_sigma_week_p);    // Error parameter on fishing week effects.
   for (int j = 0; j < n_week; j++){
      res -= dnorm(week_effect_p[j], Type(0.0), sigma_week_p, true);
   }
   
   // Priors over fishing trip effects:
   Type sigma_trip = exp(log_sigma_trip);    // Error parameter on fishing trip effects.
   for (int j = 0; j < n_trip; j++){
      res -= dnorm(trip_effect[j], Type(0.0), sigma_trip, true);
   }
        
   // Likelihood functions:
   Type logit_p_mature = 0;
   Type p_mature = 0;
   Type mu_immature = 0;
   Type mu_mature = 0;
   Type sigma = 0;
   Type logit_p_outlier = 0;
   Type p_outlier = 0;
   Type y_mod = 0;
   for (int i = 0; i < n_obs; i++){        
      // Maturity proportion model:
      logit_p_mature = alpha_p + 
                       (beta_p + observer_effect_beta_p[observer[i]-1]) * cw[i] + 
                       grid_effect_p[grid[i]-1] + 
                       week_effect_p[week[i]-1] + 
                       observer_effect_p[observer[i]-1];
      p_mature = Type(1) / (Type(1) + exp(-logit_p_mature));
      
      // Allometric relations by maturity:
      mu_immature = log_alpha_immature + log_beta_immature * log(cw[i]) ;
      mu_mature = log_alpha_mature + log_beta_mature * log(cw[i]);
                                                                 
      // Chela observation error:                                            
      sigma = exp(log_sigma + observer_error_effect[observer[i]-1]);
      
      // Outlier component proportions:
      logit_p_outlier = alpha_outlier + observer_outlier_effect[observer[i]-1];
      p_outlier = 1 / (1 + exp(-logit_p_outlier));

      // Offset chela height:
      y_mod = y[i] - mu_observer - observer_effect[observer[i]-1] - trip_effect[trip[i]-1];
      
      // Gaussian mixture density:                                            
      res -= log((1-p_mature) * (p_outlier * dnorm(log(y_mod), mu_immature, sigma + sigma_outlier, false) + 
                       (Type(1)-p_outlier) * dnorm(log(y_mod), mu_immature, sigma, false)) +                
                   (p_mature) * (p_outlier * dnorm(log(y_mod), mu_mature, sigma + sigma_outlier, false) + 
                       (Type(1)-p_outlier) * dnorm(log(y_mod), mu_mature, sigma, false)));                  
   }   
  
   return res;
}  

  