function [Ux, Uy, DUx_Dt, DUy_Dt, dUx_dx, dUx_dy, dUy_dx, dUy_dy, rotation] = double_gyre_flow_and_grad( t, XY, delta, param )

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
dfdx2 = 2*a_t;


% Double gyre flow stream function (from M.B.)
% psi = A*sin(pi*f).*sin(pi*x(2));

% fluid velocity: u
Ux = delta * (-pi*A*sin(pi*f).*cos(pi*XY(:,2)));
Uy = delta * (pi*A*cos(pi*f).*sin(pi*XY(:,2)).*dfdx);

%fluid gradients: du_x / dx, du_y / dy, du_x /dy and du_y /dx
dUx_dx = delta * (-pi^2*A*cos(pi*f).*dfdx.*cos(pi*XY(:,2)));
dUy_dy = - dUx_dx;
dUx_dy = delta * (pi^2*A*sin(pi*f).*sin(pi*XY(:,2)));
dUy_dx = delta * (pi*A*sin(pi*XY(:,2)).*( -pi*sin(pi*f).*dfdx.^2 + cos(pi*f).*dfdx2 ));

rotation = 1/2*(dUy_dx - dUx_dy);

% fluid accelaration: Du/Dt = du/dt + (u.nabla)u
dUx_dt = C*cos(pi*XY(:,2))*cos(omega*t).*cos(pi*f).*(XY(:,1).^2 - 2*XY(:,1));
dUy_dt = C*sin(pi*XY(:,2))*cos(omega*t).*( (XY(:,1).^2 - 2*XY(:,1)).*sin(pi*f).*(1 + 2*e*sin(omega*t).*(XY(:,1)-1)) - 2/pi.*cos(pi*f).*(XY(:,1) - 1) );

DUx_Dt = dUx_dt + Ux.*dUx_dx + Uy.*dUx_dy;
DUy_Dt = dUy_dt + Ux.*dUy_dx + Uy.*dUy_dy;

end