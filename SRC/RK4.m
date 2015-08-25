function [State_out, t_out] = RK4(State_init, t_init, t_fin, dt, RHS, out_period)
   
   % Coefficients for time integration
   coeff_1_2 = 1/2 * dt;
   coeff_1_3 = 1/3 * dt;
   coeff_1_6 = 1/6 * dt;
   coeff_1   =  1  * dt;
   
   % Change dt sign but not in the coeffs above
   dt = sign(t_fin - t_init) * dt;

   % Init output variables
   n_step = floor( (t_fin - t_init) / dt );
   l_out = out_period > 0;
   if( l_out )
      n_out = floor( n_step / out_period ) + 1;
      if mod(n_step, out_period) > 0
         n_out = n_out + 1;
      end
   else
      n_out = 2;
   end
   State_out = NaN( [size(State_init), n_out] );
   t_out     = NaN( [1               , n_out] );
   State_out(:,:,1) = State_init;
   t_out    (1    ) = t_init    ;
   k_out = 1;
   
   % monitor progress
   fprintf( '%6.2f %%',  0 );
   
   % Init time loop
   State_cur = State_init;
   t_cur = t_init;
   RHS_tmp = NaN(size(State_init));
   
   for k_step = 1:n_step
      
      % Step 1
      RHS_tmp = RHS(t_cur, State_cur);
      State_next = State_cur + coeff_1_6 * RHS_tmp;
      % Step 2
      State_tmp = State_cur + coeff_1_2 * RHS_tmp;
      RHS_tmp = RHS(t_cur + dt/2, State_tmp);
      State_next = State_next    + coeff_1_3 * RHS_tmp;
      % Step 3
      State_tmp = State_cur + coeff_1_2 * RHS_tmp;
      RHS_tmp = RHS(t_cur + dt/2, State_tmp);
      State_next = State_next    + coeff_1_3 * RHS_tmp;
      % Step 4
      State_tmp = State_cur + coeff_1   * RHS_tmp;
      RHS_tmp = RHS(t_cur + dt, State_tmp);
      State_cur = State_next + coeff_1_6 * RHS_tmp;
      
      t_cur = t_cur + dt;
      
      % Output
      if( l_out  & ( mod(k_step, out_period) == 0 ) )
         k_out = k_out + 1;
         State_out(:,:,k_out) = State_cur;
         t_out(k_out) = t_cur;
      end
      
      % monitor progress
      fprintf( '\b\b\b\b\b\b\b\b%6.2f %%',  k_step / n_step * 100 );
     
   end
   
   % Output last step if necessary
   if ( k_out < n_out )
      State_out(:,:,n_out) = State_cur;
      t_out(n_out) = t_cur;
   end
   
end