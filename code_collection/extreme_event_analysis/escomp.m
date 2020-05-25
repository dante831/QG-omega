function esat = escomp(temp)
    
    %{
    hlv = 2.500e6;
    rvgas = 461.50;
    T0 = 273.16;
    e0 = 610.78;

    esat = e0 * exp( -hlv/rvgas*(1.0./temp - 1.0/T0));
    %}

    %formula in Paul's code
    Tc = temp - 273.15;
    esat = 611.20*exp(17.67*Tc./(Tc+243.5));
end

