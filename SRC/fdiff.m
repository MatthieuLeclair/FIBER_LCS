function f = fdiff(r,c,x,y)

   rx = r*x;
   ry = r*y;
   r2 = 2*r;
   
   cosh_rx = cosh(rx)  ; sinh_rx = sinh(rx)  ;
   cosh_ry = cosh(ry)  ; sinh_ry = sinh(ry)  ;
   
   tanh_r  = tanh(r);
   coth_r  = coth(r);
   sinh_2r = sinh(r2);
   
   d1 = cosh(r*c) ./ (sinh_2r + r2);
   d2 = sinh(r*c) ./ (sinh_2r - r2);
   
   a = 2*cos(ry)   .* ( d2.*( cosh_rx.*(r.*coth_r-1) + rx.*sinh_rx) - d1.*( sinh_rx.*(r.*tanh_r+1) + rx.*cosh_rx) );
   b = r2.*sin(ry) .* ( d2.*( x*cosh_rx - sinh_rx.*coth_r         ) + d1.*( x*sinh_rx - cosh_rx.*tanh_r         ) );
   
   f = complex(a,b);
   
end