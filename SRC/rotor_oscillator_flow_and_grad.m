function [Ux, Uy, DUxDt, DUyDt, du, rotation] = rotor_oscillator_flow_and_grad(t, XY, delta, Ux0, Uy0, dUx0_dx, dUx0_dy, dUy0_dx, dUy0_dy, param)
   
   % Amplitude
   switch param.amp_type
     case 0   % constant amplitude
       amp = param.amp;
     case 1   % Sine wave
       amp = param.amp * sin( param.amp_omega*t + param.amp_phi );
     case 2   % Square wave   
       fprintf('ERROR: Rotor oscillator flow and fiber-type particles not compatible with square amplitude')
   end
   
   % Equivalent position for flow at rest
   X_pseudo = XY(:,1);
   Y_pseudo = XY(:,2) - param.eps * sin(param.omega*t);
   
   % Velocity and Velocity spatial gradients without amplitude and time direction yet
   Ux      = Ux0    ( X_pseudo , Y_pseudo );
   Uy      = Uy0    ( X_pseudo , Y_pseudo );
   du(:,1) = dUx0_dx( X_pseudo , Y_pseudo );
   du(:,2) = dUy0_dy( X_pseudo , Y_pseudo );
   du(:,3) = dUx0_dy( X_pseudo , Y_pseudo );
   du(:,4) = dUy0_dx( X_pseudo , Y_pseudo );
   
   % Velocity time derivative
   dY_pseudo_dt = - param.omega * param.eps * cos(param.omega*t);
   dUx_dt = amp * du(:,3) * dY_pseudo_dt;
   dUy_dt = amp * du(:,2) * dY_pseudo_dt;
   if( param.amp_type == 1 )
       damp_dt = param.amp_omega * param.amp * cos( param.amp_omega*t + param.amp_phi );
       dUx_dt = dUx_dt + damp_dt * Ux; 
       dUy_dt = dUy_dt + damp_dt * Uy; 
   end
   
   % Multiply Velocity and Velocity spatial gradients by amplitude and time direction
   Ux = delta * amp * Ux;
   Uy = delta * amp * Uy;
   du = delta * amp * du;
   
   % Compute "full" velocity time derivative
   DUxDt = dUx_dt + Ux.*du(:,1) + Uy.*du(:,3);
   DUyDt = dUy_dt + Ux.*du(:,4) + Uy.*du(:,2);
   
   % Compute rotation
   rotation = 1/2 * ( du(:,4) - du(:,3) );

end