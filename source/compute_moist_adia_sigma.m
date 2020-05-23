function [dtheta_dp_ma, theta] = compute_moist_adia_sigma(temp, p)

    % This code is slightly modified from https://pog.mit.edu/src/eff_stat_stab.m

    Rd       = 287.04;  % dry air gas constant, in J/kg/K
    Rv       = 461.50;  % water vapor gas constant, in J/kg/K
    cpd      = 1005.7;  % specific heat, dry air, in J/kg/K
    cpv      = 1870;    % specific heat, water vapor, in J/kg/K
    p0       = 1e5;     % reference pressure, Pa
    g        = 9.80665; % gravitational acceleration, m/s^2
    kappa    = Rd/cpd;
    gc_ratio = Rd/Rv;

    % saturation vapor pressure [Pa]
    Tc = temp - 273.15;
    es = 611.20*exp(17.67*Tc./(Tc+243.5));

    % latent heat of condensation [J/kg]
    Tc = temp - 273.15;
    L = (2.501 - 0.00237 * Tc) * 1e6; % Bolton 1980 equation 2

    % saturation mixing ratio and specific humidity 
    rs = gc_ratio*es./(p-es);
    qs = rs./(1+rs);

    % potential temperature
    exponent = kappa.*(1+rs./gc_ratio)./(1+rs.*cpv./cpd);
    theta    = temp.*(p0./p).^exponent;

    % density
    temp_virtual = temp.*(1.0+rs/gc_ratio)./(1.0+rs);
    rho          = p/Rd./temp_virtual;

    % moist adiabatic lapse rate
    malr = g/cpd*(1+rs)./(1+cpv/cpd.*rs) ...
            .*(1+L.*rs./Rd./temp)...
            ./(1+L.^2.*rs.*(1+rs/gc_ratio)./(Rv.*temp.^2.*(cpd+rs*cpv)));

    % derivative of potential temperature wrt pressure along a moist adiabat
    % (neglects small contribution from vertical variations of exponent)
    dtemp_dp_ma  = malr/g./rho;
    dtheta_dp_ma = dtemp_dp_ma.*theta./temp-exponent.*theta./p;

    % include the T/theta factor in the definition of dtheta_dp_ma
    dtheta_dp_ma = temp ./ theta .* dtheta_dp_ma;

end


