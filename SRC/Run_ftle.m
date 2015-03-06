clear
close all

% ======================================================================================
%                                       USER INIT
% ======================================================================================

% Flow_field
% ----------
flow_field = 'double_gyre';
% double gyre flow parameters
DG_parameters.A = 0.1        ;  % Amplitude of the stream function
DG_parameters.e = 0.25       ;  % Magnitude of the transversal oscillation
DG_parameters.omega = 2*pi/10;  % frequency of gyre oscillation
% - ML - Only 'double_gyre' so far, add 'rotor_oscillator' later

% Flow_field_handle = @double_gyre_field;   % Write similar function to use other velocity fields

% Particle type
% -------------
particle_type = 'passive';   % Choose 'fiber' or 'passive tracer'
l_Cauchy_Green = true  ;   % Compute CG tensor if true otherwise only advect particles (no cluster)

% Parameters for fiber particles
% ------------------------------
fiber_parameters.eps = 2    ;     % ellipse aspect ratio: epsilon =  b/a
fiber_parameters.St  = 1.e-3;     % particle Stokes number
fiber_parameters.Re  = 1.e-3;     % Reynolds number

% RHS_handle = @RHS_fiber;   % RHS for fiber shaped particles
% RHS_handle = @RHS_passive_test;   % RHS por pasive particles (could be faster with only 2 eq. instead of 5)

% Result file
% -----------
result_file = 'DG_Fiber_at_local_vel';

% Time info
% ---------
t_init = 0.4;        % Initial time
t_fin = t_init - 20; % Final time
out_period = -1;     % Intermediate state output every out_period time steps
                     % Put negative value for no intermediate output => only initial and final states

DT = 1.e-1;          % RK4 time step

% Initial conditions
% ------------------
x_min = 0 ; x_max = 2 ;   % x domain boundaries
y_min = 0 ; y_max = 1 ;   % y domain boundaries
dx = 0.01 ; dy = 0.01 ;   % Resolution

init_at_rest = false  ;   % if false, start at local fluid velocity

cluster_size = 1.e-5  ;   % Size of the cluster for computing flow map gradient

% Figure info
% -----------
final_figure_plot = true;
format = 'png';
visibility = 'on';
l_export = false;
resolution = 680;
ttl = 'Fiber start at local velocity';


% ======================================================================================
%                                      PROGRAM INIT
% ======================================================================================


% Program init checks
% ===================
% Check flow field and particle type
% ----------------------------------
if( ~ismember(flow_field, {'double_gyre', 'rotor_oscillator'} ) )
   fprintf('ERROR: Choose between the folloing available flow fields:\n   %s or %s\n', ...
           '''double_gyre''', '''rotor_oscillator''')
   return
end
if( ~ismember(particle_type, {'fiber', 'passive'} ) )
   fprintf('ERROR: Choose between the folloing available particle types:\n   %s or %s\n', ...
           '''fiber''', '''passive''')
   return
end
% Forward or backward integration
% -------------------------------
delta = sign(t_fin - t_init);

% Initial conditions
% ==================
% Positions
% ---------
% Initial particle positions
X_0 = x_min+dx:dx:x_max-dx;
Y_0 = y_min+dy:dy:y_max-dy;
[X_0, Y_0] = meshgrid(X_0, Y_0);
sz_ini = size(X_0);   % Store size for later use

n_part = numel(X_0);
sz = [n_part, 1];
X_0 = reshape(X_0, sz);
Y_0 = reshape(Y_0, sz);

% Initial cluster points positions
if( l_Cauchy_Green )
   X_cl1 = X_0 + cluster_size/2 ; Y_cl1 = Y_0                  ;
   X_cl2 = X_0                  ; Y_cl2 = Y_0 + cluster_size/2 ;
   X_cl3 = X_0 - cluster_size/2 ; Y_cl3 = Y_0                  ;
   X_cl4 = X_0                  ; Y_cl4 = Y_0 - cluster_size/2 ;
   
   X_0 = [X_0; X_cl1; X_cl2; X_cl3; X_cl4];
   Y_0 = [Y_0; Y_cl1; Y_cl2; Y_cl3; Y_cl4];
   clear X_cl1 X_cl2 X_cl3 X_cl4 Y_cl1 Y_cl2 Y_cl3 Y_cl4
   sz = size(X_0);
end

% Additional fiber state variables
% --------------------------------
if( strcmp(particle_type, 'fiber') )
   % Velocity
   % --------
   if init_at_rest
      X_dot_0 = zeros(sz);
      Y_dot_0 = zeros(sz);
   else   % init with local velocity
      switch flow_field
        case 'double_gyre'
          [X_dot_0, Y_dot_0]  = double_gyre_flow(t_init, [X_0 , Y_0], delta, DG_parameters);
          X_dot_0 = 0.9 * X_dot_0;
          Y_dot_0 = 0.9 * Y_dot_0;
        case 'rotor_oscillator'
          fprintf('ERROR: rotor_oscillator not coded so far\n')
          return
      end
   end

   % Angle
   % -----
   % Theta_0 = zeros(sz);
   Theta_0 = rand(sz)*2*pi;
   % switch flow_field
   %   case 'double_gyre'
   %     [~, ~, ~, ~, du, ~] = double_gyre_flow_and_grad( t, [X_0 , Y_0], delta, , DG_parameters )
   %   case 'rotor_oscillator'
   %     fprintf('ERROR: rotor_oscillator not coded so far\n')
   %     return
   % end
   % Theta_0 = 1/2 * atan2(du(:,4)+du(:,3), du(:,1)-du(:,2));

   % State vector
   % ------------
   State_init = [X_0 Y_0 X_dot_0 Y_dot_0 Theta_0];
   clear X_0 Y_0 X_dot_0 Y_dot_0 Theta_0
else
   State_init = [X_0 Y_0];
   clear X_0 Y_0
end


% ======================================================================================
%                               INTEGRATE TRAJECTORIES
% ======================================================================================

% Determine flow field and RHS functions
% --------------------------------------
switch flow_field
  case 'double_gyre'
    switch particle_type
      case 'passive'
        flow_field_handle = @(t, XY) double_gyre_flow( t, XY, delta, DG_parameters );
        RHS = @(t, State) RHS_passive( t, State, flow_field_handle );
      case 'fiber'
        flow_field_handle = @(t, XY) double_gyre_flow_and_grad( t, XY, delta, DG_parameters );
        RHS = @(t, State) RHS_fiber( t, State, flow_field_handle, fiber_parameters );
    end
  otherwise
    fprintf('ERROR: Only double gyre coded so far\n')
    return
end

% Integrate particle trajectories
% -------------------------------
[State_out, t_out] = RK4(State_init, t_init, t_fin, DT, RHS, out_period);


% ======================================================================================
%                             CAUCHY-GREEN TENSOR INVARIANTS
% ======================================================================================

% Recover State variables
n_out = size(State_out, 3);
X = reshape( State_out(1:n_part,1,:), [sz_ini n_out] );
Y = reshape( State_out(1:n_part,2,:), [sz_ini n_out] );
if( strcmp(particle_type, 'fiber') )
   X_dot = reshape( State_out(1:n_part,3,:), [sz_ini n_out] );
   Y_dot = reshape( State_out(1:n_part,4,:), [sz_ini n_out] );
   Theta = reshape( State_out(1:n_part,5,:), [sz_ini n_out] );
end

% Compute the Cauchy Green invariants
% -----------------------------------
if( l_Cauchy_Green )
   X_cl = squeeze( State_out(n_part+1:end,1,:) );
   Y_cl = squeeze( State_out(n_part+1:end,2,:) );
   clear State_out
   lda1 = NaN( [sz_ini n_out] );
   lda2 = NaN( [sz_ini n_out] );
   xi1x = NaN( [sz_ini n_out] );
   xi1y = NaN( [sz_ini n_out] );
   xi2x = NaN( [sz_ini n_out] );
   xi2y = NaN( [sz_ini n_out] );
   for k_out = 1:n_out
      % Flow map gradient matrix components
      M11 = ( X_cl(0*n_part+1:1*n_part,k_out) - X_cl(2*n_part+1:3*n_part,k_out) ) / cluster_size ;
      M21 = ( Y_cl(0*n_part+1:1*n_part,k_out) - Y_cl(2*n_part+1:3*n_part,k_out) ) / cluster_size ;
      M12 = ( X_cl(1*n_part+1:2*n_part,k_out) - X_cl(3*n_part+1:4*n_part,k_out) ) / cluster_size ;
      M22 = ( Y_cl(1*n_part+1:2*n_part,k_out) - Y_cl(3*n_part+1:4*n_part,k_out) ) / cluster_size ;
      
      % Eigen values and eigen spaces of the CG tensor
      [l1, l2, x1x, x1y, x2x, x2y] = CauchyGreen(M11, M12, M21, M22);
      
      % Store variables at each output step
      X   (:,:,k_out) = reshape( X_cl(1:n_part,k_out), sz_ini );
      Y   (:,:,k_out) = reshape( Y_cl(1:n_part,k_out), sz_ini );
      lda1(:,:,k_out) = reshape( l1                   , sz_ini );
      lda2(:,:,k_out) = reshape( l2                   , sz_ini );
      xi1x(:,:,k_out) = reshape( x1x                  , sz_ini );
      xi1y(:,:,k_out) = reshape( x1y                  , sz_ini );
      xi2x(:,:,k_out) = reshape( x2x                  , sz_ini );
      xi2y(:,:,k_out) = reshape( x2y                  , sz_ini );
   end
else
   clear State_out
end

% ======================================================================================
%                                      SAVE RESULTS
% ======================================================================================

result_obj = matfile(['../OUTPUT/' result_file], 'Writable', true);
result_obj.X     = X     ;
result_obj.Y     = Y     ;
result_obj.time  = t_out ;
if( strcmp(particle_type, 'fiber') )
   result_obj.X_dot            = X_dot            ;
   result_obj.Y_dot            = Y_dot            ;
   result_obj.Theta            = Theta            ;
   result_obj.fiber_parameters = fiber_parameters ;
end
if( strcmp(flow_field, 'double_gyre') )
   result_obj.DG_parameters = DG_parameters;
end
if( l_Cauchy_Green )
   result_obj.lda1  = lda1  ;
   result_obj.lda2  = lda2  ;
   result_obj.xi1x  = xi1x  ;
   result_obj.xi1y  = xi1y  ;
   result_obj.xi2x  = xi2x  ;
   result_obj.xi2y  = xi2y  ;
end
% ======================================================================================
%                                       PLOT
% ======================================================================================

if( final_figure_plot & l_Cauchy_Green )
   
   fig = plot_ftle( X(1,:,1), Y(:,1,1), log(lda2(:,:,end))/(2*abs(t_fin-t_init)), ttl, visibility );
   
   if( l_export )
      full_file_name = ['../FIGURES/' result_file '.' format];
      switch format
        case 'pdf'
          print(fig,'-dpdf', '-painters', full_file_name);
        case 'eps'
          print(fig,'-depsc2', '-painters', full_file_name);
        case 'jpg'
          print(fig,'-djpeg', ['-r' int2str(resolution)], full_file_name);
        case 'png'
          print(fig,'-dpng' , ['-r' int2str(resolution)], full_file_name);
        case 'fig'
          saveas( fig, full_file_name);
      end
   end
   
end

