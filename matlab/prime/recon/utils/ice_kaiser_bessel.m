function value=ice_kaiser_bessel(x, W, beta)

arg = 1. - 4.*x*x/W/W;

if ( arg < 0.0)
    value=0.0;
else
    arg = beta * sqrt(arg);
    bess = bessi0( arg );
    value=bess/W;
end;
return;


function value=bessi0( x )
ax=abs(x);
if ( ax < 3.75 )
    y=x/3.75;
    y=y*y;
    ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
    
else
    y=3.75/ax;
    ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1...
        +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2...
        +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1...
        +y*0.392377e-2))))))));
    
end;
value=ans;
