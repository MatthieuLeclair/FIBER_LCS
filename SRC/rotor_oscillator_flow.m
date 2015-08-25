function [Ux, Uy] = rotor_oscillator_flow(t, XY, delta, Ux0, Uy0, param)
   
   Ux = delta * Ux0( XY(:,1) , XY(:,2) - param.eps * sin(param.omega*t) );
   Uy = delta * Uy0( XY(:,1) , XY(:,2) - param.eps * sin(param.omega*t) );
   
   switch param.amp_type
     case 0   % constant amplitude
       amp = param.amp;
     case 1   % Sine wave
       amp = param.amp * sin( param.amp_omega*t + param.amp_phi );
     case 2   % Square wave   
       amp = param.amp * sign( sin( param.amp_omega*t + param.amp_phi ) );
   end
   
   Ux = amp * Ux;
   Uy = amp * Uy;
   
end