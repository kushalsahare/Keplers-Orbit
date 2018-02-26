function K = Kepler_func(x)

r3 = (x(1)^2 + x(2)^2)^1.5;
K = [ x(3) ;
      x(4);
      -x(1)/r3 ;
      -x(2)/r3] ;
end