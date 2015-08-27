function [Ux, Uy] = rotor_oscillator_flow(t, XY, delta, Ux0, Uy0, param)
   
   % Amplitude
   switch param.amp_type
     case 0   % constant amplitude
       amp = param.amp;
     case 1   % Sine wave
       amp = param.amp * sin( param.amp_omega*t + param.amp_phi );
     case 2   % Square wave   
       amp = param.amp * sign( sin( param.amp_omega*t + param.amp_phi ) );
   end
   
   % Equivalent position for flow at rest
   X_pseudo = XY(:,1);
   Y_pseudo = XY(:,2) - param.eps * sin(param.omega*t);
   
   % Velocity
   Ux = delta * amp * Ux0( X_pseudo , Y_pseudo );
   Uy = delta * amp * Uy0( X_pseudo , Y_pseudo );
   
end