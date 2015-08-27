function [Ux, Uy] = double_gyre_flow( t, XY, delta, param )

% XY is the field with the x and y coordinates of points where u and v need
% to be determined

% Parameters
A     = param.A    ;    % Stream function amplitude
e     = param.e    ;    % Amplitud of the transversal oscillation
omega = param.omega;    % frequency of gyre oscillation

C = -A*pi^2*e*omega;    % Constant acceleration

a_t = e*sin(omega*t);
b_t = 1 - 2*e*sin(omega*t);

f = a_t.*XY(:,1).^2 + b_t.*XY(:,1);
dfdx = 2*a_t.*XY(:,1) + b_t;

% Double gyre flow stream function
% psi = A*sin(pi*f).*sin(pi*x(2));

% fluid velocity: u
Ux = delta * (-pi*A*sin(pi*f).*cos(pi*XY(:,2)));
Uy = delta * (pi*A*cos(pi*f).*sin(pi*XY(:,2)).*dfdx);

end