function [Ux0 Uy0 varargout] = generate_RO_flow(param, file, l_load, l_save, l_grad)
 
   if( l_load )
      
      % Load interpolants from file or ...
      % ==================================
      
      save_obj = matfile(file, 'Writable', false);
      Ux0 = save_obj.Ux0;
      Uy0 = save_obj.Uy0;
      if( l_grad )
         varargout{1} = save_obj.dUx0_dx;
         varargout{2} = save_obj.dUx0_dy;
         varargout{3} = save_obj.dUy0_dx;
         varargout{4} = save_obj.dUy0_dy;
      end
      
   else

      % ... Compute interpolants
      % ========================
      
      if( l_grad )
         fprintf('\nComputing rotor oscillator analytical function and gradients on a grid\n')
      else
         fprintf('\nComputing rotor oscillator analytical function on a grid\n')
      end
      
      % RO flow parameters
      % ------------------
      a = param.a; %radius in m
      h = param.h; %tank width in m - was 9.9cm
      c = param.c; %initial position of the rod/rotor on the c axis (y=0)
      
      % Define grid
      % -----------
      x_range = param.x_range;
      y_range = param.y_range;
      dx      = param.dx     ;
      dy      = param.dy     ;
      dcl     = param.dcl    ;

      x_grid = x_range(1):dx:x_range(2);
      y_grid = y_range(1):dy:y_range(2);

      [XG,YG] = ndgrid(x_grid, y_grid);
      sz_grid = size(XG);
      n_points = numel(XG);
      XG = reshape(XG, n_points, 1);
      YG = reshape(YG, n_points, 1);
      
      if( l_grad )
         X = [XG ; XG+dcl ; XG     ; XG-dcl ; XG    ];
         Y = [YG ; YG     ; YG+dcl ; YG     ; YG-dcl];
         n_points_tot = 5*n_points;
      else
         X = XG;
         Y = YG;
         n_points_tot = n_points;
      end

      % Compute term 1
      % --------------
      psidx = - 0.5 * (exp(pi*Y) + 2*cos((pi*(c + X))/2).*exp((pi*Y)/2) + 1)        ...
              .* ( ( pi*exp((pi*Y)/2) .* sin((pi*(c - X))/2) )                      ...
                   ./ ( exp(pi*Y) + 2*cos((pi*(c + X))/2) .* exp((pi*Y)/2) + 1 )    ...
                   - ( pi*sin((pi*(c + X))/2) .* exp((pi*Y)/2)                      ...
                       .*(exp(pi*Y) - 2*exp((pi*Y)/2).*cos((pi*(c - X))/2) + 1) )   ...
                   ./ ( exp(pi*Y) + 2*cos((pi*(c + X))/2).*exp((pi*Y)/2) + 1 ).^2 ) ...
              ./ ( exp(pi*Y) - 2*exp((pi*Y)/2).*cos((pi*(c - X))/2) + 1 );

      psidy = 0.5 * ( exp(pi*Y) + 2*cos((pi*(c + X))/2).*exp((pi*Y)/2) + 1 )        ...
              .* ( ( pi*exp(pi*Y) - pi*exp((pi*Y)/2).*cos((pi*(c - X))/2) )         ...
                   ./ ( exp(pi*Y) + 2*cos((pi*(c + X))/2).*exp((pi*Y)/2) + 1 )      ...
                   - ( ( pi*exp(pi*Y) + pi*cos((pi*(c + X))/2).*exp((pi*Y)/2) )     ...
                       .* ( exp(pi*Y) - 2*exp((pi*Y)/2).*cos((pi*(c - X))/2) + 1) ) ...
                   ./ ( exp(pi*Y) + 2*cos((pi*(c + X))/2).*exp((pi*Y)/2) + 1 ).^2 ) ...
              ./ ( exp(pi*Y) - 2*exp((pi*Y)/2).*cos((pi*(c - X))/2) + 1 );


      % Compute term 2
      % --------------
      Qfx = zeros(n_points_tot,1);
      Qfy = zeros(n_points_tot,1);
      %limits of integration
      kmax=450;
      kmin=10.^-8;%0.1;%
      for k = 1:n_points_tot
         fprintf( '\b\b\b\b\b\b\b\b%6.2f %%', k / n_points_tot * 100 );
         x = X(k);
         y = Y(k);
         fdiffx = @(r) cos(r*y).*...
                  ((sinh(r*c).*(2*cosh(r*x) - 2*r.*cosh(r*x).*coth(r) + ...
                                2*r*x.*sinh(r*x)))./(2*r - sinh(2*r)) - ...
                   (cosh(r*c).*(2*sinh(r*x) - 2*r.*sinh(r*x).*tanh(r) + ...
                                2*r*x.*cosh(r*x)))./(2*r + sinh(2*r)));
         fdiffy =  @(r)r.*sin(r*y).*...
                   ((2*cosh(c*r).*(x*sinh(r*x) - cosh(r*x).*tanh(r)))./...
                    (2*r + sinh(2*r)) - ... 
                    (2*sinh(c*r).*(x*cosh(r*x) - sinh(r*x).*coth(r)))./...
                    (2*r - sinh(2*r)));
         
         psidx(k) = psidx(k) + quadgk(fdiffx,kmin,kmax,'WayPoints',0);
         psidy(k) = psidy(k) + quadgk(fdiffy,kmin,kmax,'WayPoints',0);
      end
      
      % psidx(isinf(psidx)|isnan(psidx)) = 0;
      % psidy(isinf(psidy)|isnan(psidy)) = 0;

      % Create interpolants
      % -------------------
      XG = reshape(XG, sz_grid);
      YG = reshape(YG, sz_grid);
      
      Ux0 = griddedInterpolant(XG, YG, reshape(-psidy(1:n_points),sz_grid), 'cubic');
      Uy0 = griddedInterpolant(XG, YG, reshape( psidx(1:n_points),sz_grid), 'cubic');
      
      if( l_grad )
         % dUx0_dx
         field_tmp = reshape( psidy(3*n_points+1:4*n_points) - psidy(1*n_points+1:2*n_points) , sz_grid ) / dcl;
         varargout{1} = griddedInterpolant(XG, YG, reshape(field_tmp,sz_grid), 'cubic');
         % dUx0_dy
         field_tmp = reshape( psidy(4*n_points+1:5*n_points) - psidy(2*n_points+1:3*n_points) , sz_grid ) / dcl;
         varargout{2} = griddedInterpolant(XG, YG, reshape(field_tmp,sz_grid), 'cubic');
         % dUy0_dx
         field_tmp = reshape( psidx(1*n_points+1:2*n_points) - psidx(3*n_points+1:4*n_points) , sz_grid ) / dcl;
         varargout{3} = griddedInterpolant(XG, YG, reshape(field_tmp,sz_grid), 'cubic');
         % dUy0_dy
         field_tmp = reshape( psidx(2*n_points+1:3*n_points) - psidx(4*n_points+1:5*n_points) , sz_grid ) / dcl;
         varargout{4} = griddedInterpolant(XG, YG, reshape(field_tmp,sz_grid), 'cubic');
         
         clear field_tmp
      end
      
      fprintf('\n')
      
      % Write interpolants to file if required
      % --------------------------------------
      if( l_save )
         if(exist(file) == 2) delete(file); end
         save_obj = matfile(file, 'Writable', true);
         save_obj.Ux0 = Ux0;
         save_obj.Uy0 = Uy0;
         if( l_grad )
            save_obj.dUx0_dx = varargout{1};
            save_obj.dUx0_dy = varargout{2};
            save_obj.dUy0_dx = varargout{3};
            save_obj.dUy0_dy = varargout{4};
         end
      end
      
   end

   

end