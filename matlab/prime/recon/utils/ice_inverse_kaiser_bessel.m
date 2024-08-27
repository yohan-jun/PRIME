function value=ice_inverse_kaiser_bessel( x, W, beta)

arg = pi*pi * W*W * x*x - beta*beta;
if ( arg >= 0.0 )
    arg = sqrt( arg );
    if ( arg < 1.0e-4 )
        dsinc = 1.0 ;
    else
        dsinc = sin( arg ) / arg;
    end;
    invKB = dsinc * beta / sinh(beta);
    
else
    arg =arg.* -1.;
    arg = sqrt( arg );
    if ( arg < 1.0e-4 )
        dsinc = 1.0;
    else
        dsinc = sinh( arg ) / arg;
    end;
    invKB =sqrt(-1).*( dsinc * beta / sinh(beta));  %// make sure result = 1 when x = 0
end;

value=invKB;

return;