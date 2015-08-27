function f = fdiff(r,c,x,y)

   rx = r*x;
   ry = r*y;
   r2 = 2*r;
   
   cosh_rx = cosh(rx)  ; sinh_rx = sinh(rx)  ;
   cosh_ry = cosh(ry)  ; sinh_ry = sinh(ry)  ;
   
   cosh_rc = cosh(r*c) ; sinh_rc = sinh(r*c) ;
   
   tanh_r  = tanh(r) ;
   coth_r  = coth(r) ;
   
   sinh_2r = sinh(r2);
   d1 = sinh_2r + r2;
   d2 = sinh_2r - r2;
   
   a = -2*cos(ry) .* ( ( cosh_rc.*( sinh_rx.*(r.*tanh_r+1) + rx.*cosh_rx) ) ./ d1 + ...
                       ( sinh_rc.*( cosh_rx.*(r.*coth_r-1) + rx.*sinh_rx) ) ./ d2       );
   
   b = r2.*sin(ry).* ( ( cosh_rc.*( x*sinh_rx - cosh_rx.*tanh_r         ) ) ./ d1 + ... 
                       ( sinh_rc.*( x*cosh_rx - sinh_rx.*coth_r         ) ) ./ d2       );
   
   f = complex(a,b);
   
end