function F = RHS_passive(t, State, flow_field)
   
   % State is the state vector for all particles:
   % | X |      | State(:,1) |
   % |   |  =   |            |
   % | Y |      | State(:,2) |
   
   % Flow properties
   [Ux, Uy] = flow_field( t, [State(:,1), State(:,2) ] );
   
   % Right Hand Side
   F(:,1) = Ux;
   F(:,2) = Uy;

end