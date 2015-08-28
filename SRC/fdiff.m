function f_comp = fdiff(r,c,x,y)

   rx = r*x;
   ry = r*y;
   r2 = 2*r;
   
   cosh_rx = cosh(rx);
   sinh_rx = sinh(rx);
   
   tanh_r  = tanh(r);
   coth_r  = coth(r);
   sinh_2r = sinh(r2);
   
   d1 = sinh(r*c) ./ (sinh_2r - r2);
   d2 = cosh(r*c) ./ (sinh_2r + r2);
   
   a =    cos(ry) .* ( d1.*( (r.*coth_r-1).*cosh_rx + rx.*sinh_rx) - d2.*( (r.*tanh_r+1).*sinh_rx + rx.*cosh_rx) );
   b = r.*sin(ry) .* ( d1.*( x*cosh_rx - coth_r.*sinh_rx         ) + d2.*( x*sinh_rx - tanh_r.*cosh_rx         ) );
   
   f_comp = complex(a,b);
   
end