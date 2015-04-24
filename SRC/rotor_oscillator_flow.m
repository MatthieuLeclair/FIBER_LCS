function [Ux, Uy] = rotor_oscillator_flow(t, XY, delta, Fx, Fy, param)
   
   Ux = - delta * Fy( XY(:,1) , XY(:,2) - param.eps * sin(param.omega*t) );
   Uy =   delta * Fx( XY(:,1) , XY(:,2) - param.eps * sin(param.omega*t) );
   
   switch param.amp_type
     case 0   % constant amplitude
       amp = param.amp;
     case 1   % Sine wave
       amp = param.amp * sin( param.omswitch*t + param.amp_phi );
     case 2   % Square wave   
       amp = param.amp * sign( sin( param.omswitch*t + param.amp_phi ) );
   end
   
   Ux = amp * Ux;
   Uy = amp * Uy;
   
end