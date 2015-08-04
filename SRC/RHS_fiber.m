function F = RHS_fiber(t, State, flow_field, param)
   
   % State is the state vector for all particles:
   % | X         |     | State(:,1) |
   % | Y         |     | State(:,2) |
   % | X_dot     |  =  | State(:,3) |
   % | Y_dot     |     | State(:,4) |
   % | Theta     |     | State(:,5) |
   
   % Decompose State vector (for code readability)
   X     = State(:,1);
   Y     = State(:,2);
   X_dot = State(:,3);
   Y_dot = State(:,4);
   Theta = State(:,5);
      
   % Coefficients
   
   L1 = (param.eps^2 - 1) / ( 2*((2*param.eps^2 - 1)*log(param.eps+sqrt(param.eps^2 - 1))/sqrt(param.eps^2 - 1))  - param.eps) ;
   L2 = (param.eps^2 - 1) / (   ((2*param.eps^2 - 3)*log(param.eps+sqrt(param.eps^2 - 1))/sqrt(param.eps^2 - 1))  + param.eps) ;
   
   pre_fact_1 = (8*param.Re)/(3*param.St);
   pre_fact_2 = 1/2*(1-param.eps^2)/(1+param.eps^2);
   
   % Flow properties
   [Ux, Uy, DUxDt, DUyDt, du, rotation] = flow_field(t, [X,Y]);
   
   % Right Hand Side
   COS = cos(Theta);
   SIN = sin(Theta);
   COS_2   = COS .* COS;
   SIN_2   = SIN .* SIN;
   COS_SIN = (L1-L2) * COS .* SIN;
   clear COS SIN
   
   F(:,1) = X_dot;
   F(:,2) = Y_dot;
   F(:,3) = DUxDt - pre_fact_1 * ( (L1*COS_2 + L2*SIN_2) .* (X_dot- Ux)  + COS_SIN .* (Y_dot- Uy) );
   F(:,4) = DUyDt - pre_fact_1 * ( (L1*SIN_2 + L2*COS_2) .* (Y_dot- Uy)  + COS_SIN .* (X_dot- Ux) );
   F(:,5) = rotation + pre_fact_2 * ( sin(2*Theta).*(du(:,1)-du(:,2)) - cos(2*Theta).*(du(:,4)+du(:,3)) );
      
end