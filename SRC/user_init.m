% Flow_field
% ----------
flow_field = 'double_gyre';
% flow_field = 'rotor_oscillator';

% double gyre flow parameters
DG_parameters.A = 0.1         ;  % Amplitude of the stream function
DG_parameters.e = 0.25        ;  % Magnitude of the transversal oscillation
DG_parameters.omega = 2*pi/10 ;  % frequency of gyre oscillation

% rotor oscillator generator parameters
RO_generator_parameters.x_range = [-1 1]     ;   % add description here
RO_generator_parameters.y_range = [-2.5 2.5] ;   % add description here
RO_generator_parameters.dx      = 0.01       ;   % add description here
RO_generator_parameters.dy      = 0.01       ;   % add description here
RO_generator_parameters.a       = 0.003      ;   % add description here
RO_generator_parameters.h       = 0.01       ;   % add description here
RO_generator_parameters.c       = 0.54       ;   % add description here

% rotor oscillator flow parameters
RO_parameters.eps       = 1       ;   % amplitude of oscillation
RO_parameters.omega     = 2*pi    ;   % pulsation of oscillation
RO_parameters.amp_type  = 0       ;   % type of amplitude evolution: 0 => constant ; 1 => sinusoidal ; 2 => square
RO_parameters.amp       = 1       ;   % amplitude of the stream function
RO_parameters.amp_omega = 2*pi/10 ;   % amplitude pulsation
RO_parameters.amp_phi   = 0*pi    ;   % amplitude initial phase

% rotor oscillator generator options
RO_generator_file   = '../INPUT/RO_interpolant_01' ;  % file where interpolant is read or written 
RO_generator_save   = false  ;  % save or not RO interpolant
RO_generator_load   = true   ;  % load or not RO interpolant

% Particle type
% -------------
particle_type  = 'passive';   % Choose 'fiber' or 'passive'
% Parameters for fiber particles
fiber_parameters.eps = 8    ;     % ellipse aspect ratio: epsilon =  b/a
fiber_parameters.St  = 3.e-2;     % particle Stokes number
fiber_parameters.Re  = 1.e-3;     % Reynolds number

% Cauchy Green Tensor
% -------------------
l_Cauchy_Green = true;   % Compute CG tensor if true otherwise only advect particles (no cluster)

% Result file
% -----------
result_file = 'test';

% Time info
% ---------
t_init = 0;          % Initial time
t_fin = 20;          % Final time
out_period = -1;     % Intermediate state output every out_period time steps
                     %    Put negative value for no intermediate output => only initial and final states
DT = 1.e-1;          % RK4 time step

% Initial conditions
% ------------------
% Domain
x_min = 0  ; x_max = 2 ;   % x domain boundaries
y_min = 0  ; y_max = 1 ;   % y domain boundaries
dx = 0.01  ; dy = 0.01 ;   % Resolution
% Velocity
init_at_rest = true  ;   % if false, start at local fluid velocity
init_coeff   = 1     ;   % start at local fluid velocity x init_coeff
% Cluster size
cluster_size = 1.e-5  ;   % Size of the cluster for computing flow map gradient

% Figure info
% -----------
l_plot_ftle = true;
l_export    = false;
format      = 'png';
visibility  = 'on';
resolution  = 340;
ttl         = 'test figure';