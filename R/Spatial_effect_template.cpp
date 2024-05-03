   DATA_VECTOR(grid_x);
   DATA_VECTOR(grid_y);
   
   // Grid centroid distance matrix:
   matrix<Type> D_grid(n_grid, n_grid); 
   for (i = 0; i < n_grid; i++){
      dist_grid(i,i) = Type(0);
      for (j = 0; j < i; j++){
         D_grid(i,j) = sqrt((grid_x[i] - grid_x[j]) *  (grid_x[i] - grid_x[j]) + (grid_y[i] - grid_y[j]) *  (grid_y[i] - grid_y[j]));
         D_grid(j,i) = D_grid(i,j);
      }
   }
   
   // Grid spatial effect:
   matrix<Type> Sigma_grid(n_grid, n_grid); 
   for (i = 0; i < n_grid; i++){
      Sigma_grid_p(i,i) = sigma_grid + sigma_grid_scale;
      for (j = 0; j < i; j++){
         Sigma_grid_p(i,j) = sigma_grid_scale * exp(-D_grid(i,j) / range_grid); // Exponential covariance function.               
         Sigma_grid_p(j,i) = Sigma_grid_p(i,j);
      }
   }

   using namespace density;
   MVNORM_t<Type> neg_log_density_grid_p(Sigma_grid_p);
   res += neg_log_density_grid(grid_spatial_effect_p);

    