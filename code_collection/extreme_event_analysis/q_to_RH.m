function RH = q_to_RH(q, T, p)

    rdgas = 287.04;
    rvgas = 461.50;
    
    nu = rvgas / rdgas;

    pv = p .* q .* nu ./ (1 - q .* (1 - nu));% vapor pressure
    
    % saturation vapor pressure
    Tc = T - 273.15;
    pv_s = 611.20*exp(17.67*Tc./(Tc+243.5)); % Bolton 1980 equation 10

    %pv_s = comp(T); % saturation vapor pressure
    
    RH = pv ./ pv_s;

end
