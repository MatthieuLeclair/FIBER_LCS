function [Ux, Uy] = rotor_oscillator_flow(t, XY, delta, Fx, Fy, param)
   
   Ux = - delta * Fy( XY(:,1) , XY(:,2) - param.eps * sin(param.omega*t) );
   Uy =   delta * Fx( XY(:,1) , XY(:,2) - param.eps * sin(param.omega*t) );
   
end