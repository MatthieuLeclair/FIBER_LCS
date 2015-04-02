function [Fx Fy] = generate_RO_flow(param, file, l_load, l_save)

   % Setup dimensions
   a = param.a; %radius in m
   h = param.h;  %tank width in m - was 9.9cm
   c = param.c; %initial position of the rod/rotor on the c axis (y=0)

   if( ~l_load )
      
      fprintf('\nComputing rotor oscillator analytical function on a grid\n')
      
      x_range = param.x_range;
      y_range = param.y_range;
      dx      = param.dx     ;
      dy      = param.dy     ;

      x_grid = x_range(1):dx:x_range(2);
      y_grid = y_range(1):dy:y_range(2);

      [X,Y] = ndgrid(x_grid, y_grid);

      %limits of integration
      kmax=450;
      kmin=10.^-8;%0.1;%

      term1dX = - 0.5 * (exp(pi*Y) + 2*cos((pi*(c + X))/2).*exp((pi*Y)/2) + 1)           ...
                .* ( ( pi*exp((pi*Y)/2) .* sin((pi*(c - X))/2) )                      ...
                     ./ ( exp(pi*Y) + 2*cos((pi*(c + X))/2) .* exp((pi*Y)/2) + 1 )    ...
                     - ( pi*sin((pi*(c + X))/2) .* exp((pi*Y)/2)                      ...
                         .*(exp(pi*Y) - 2*exp((pi*Y)/2).*cos((pi*(c - X))/2) + 1) )   ...
                     ./ ( exp(pi*Y) + 2*cos((pi*(c + X))/2).*exp((pi*Y)/2) + 1 ).^2 )    ...
                ./ ( exp(pi*Y) - 2*exp((pi*Y)/2).*cos((pi*(c - X))/2) + 1 );

      term1dY = 0.5 * ( exp(pi*Y) + 2*cos((pi*(c + X))/2).*exp((pi*Y)/2) + 1 )           ...
                .* ( ( pi*exp(pi*Y) - pi*exp((pi*Y)/2).*cos((pi*(c - X))/2) )         ...
                     ./ ( exp(pi*Y) + 2*cos((pi*(c + X))/2).*exp((pi*Y)/2) + 1 )      ...
                     - ( ( pi*exp(pi*Y) + pi*cos((pi*(c + X))/2).*exp((pi*Y)/2) )     ...
                         .* ( exp(pi*Y) - 2*exp((pi*Y)/2).*cos((pi*(c - X))/2) + 1) ) ...
                     ./ ( exp(pi*Y) + 2*cos((pi*(c + X))/2).*exp((pi*Y)/2) + 1 ).^2 )    ...
                ./ ( exp(pi*Y) - 2*exp((pi*Y)/2).*cos((pi*(c - X))/2) + 1 );


      % term 2
      Qfx = zeros(length(kmin),length(kmax));
      Qfy = zeros(length(kmin),length(kmax));

      ni = numel(x_grid);
      nj = numel(y_grid);  
      for i=1:ni
         % monitor progress
         fprintf( '\b\b\b\b\b\b\b\b%6.2f %%', i / ni * 100 );
         for j = 1:nj
            
            x = x_grid(i);
            y = y_grid(j);
            
            fdiffx = (@(k)cos(k*y).*...
                      ((sinh(k*c).*(2*cosh(k*x) - 2*k.*cosh(k*x).*coth(k) + ...
                                    2*k*x.*sinh(k*x)))./(2*k - sinh(2*k)) - ...
                       (cosh(k*c).*(2*sinh(k*x) - 2*k.*sinh(k*x).*tanh(k) + ...
                                    2*k*x.*cosh(k*x)))./(2*k + sinh(2*k))));

            fdiffy = (@(k)k.*sin(k*y).*...
                      ((2*cosh(c*k).*(x*sinh(k*x) - cosh(k*x).*tanh(k)))./...
                       (2*k + sinh(2*k)) - ... 
                       (2*sinh(c*k).*(x*cosh(k*x) - sinh(k*x).*coth(k)))./...
                       (2*k - sinh(2*k))));
            
            Qfx(i,j)=quadgk(fdiffx,kmin,kmax,'WayPoints',0);
            Qfy(i,j)=quadgk(fdiffy,kmin,kmax,'WayPoints',0);
            
         end
      end
      
      psidx=term1dX + Qfx;
      [m, n] = find(isinf(psidx));
      psidx(m,n) = 0;
      clear m n
      [m, n] = find(isnan(psidx));
      psidx(m,n) = 0;
      clear m n

      psidy=term1dY+Qfy;
      [m, n] = find(isinf(psidy));
      psidy(m,n) = 0;
      clear m n
      [m, n] = find(isnan(psidy));
      psidy(m,n) = 0;
      clear m n

      Fx = griddedInterpolant(X, Y, psidx, 'cubic'); 
      Fy = griddedInterpolant(X, Y, psidy, 'cubic');
      
      fprintf('\n')
      
   else
      
      save_obj = matfile(file, 'Writable', false);
      Fx = save_obj.Fx;
      Fy = save_obj.Fy;
      
   end

   if( l_save )
      save_obj = matfile(file, 'Writable', true);
      save_obj.Fx = Fx;
      save_obj.Fy = Fy;
   end

end