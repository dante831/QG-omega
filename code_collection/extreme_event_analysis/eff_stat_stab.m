function [dtheta_dp_ma, dtheta_dp_eff, theta] = eff_stat_stab(p, temp, lambda)

% calculates the effective static stability derived in O'Gorman, JAS, 2011, pages 75-90 according to equation 8 
%
% inputs are 1d vertical profiles of pressure and temperature and the asymmetry parameter lambda (default value 0.6):
% p       pressure (Pa)
% temp    temperature (K)
% lambda  asymmetry parameter defined by equation 5 of O'Gorman, JAS, 2011


 if nargin<3
  lambda = 0.6; % default value of lambda (appropriate for midlatitudes)
 end

 if ~isequal(size(temp), size(p)) | min(size(p))~=1
  error('temp and p must be one-dimensional and of same size')
 end

 % constants
 Rd       = 287.04;      % gas constant for dry air [J/kg/K]
 Rv       = 461.5;       % gas constant for water vapor [J/kg/K]
 cpd      = 1005.7;      % specific heat dry air [J/kg/K]
 cpv      = 1870;        % specific heat water vapor [J/kg/K]
 g        = 9.80665;     % gravitational acceleration [m/s^2]
 p0       = 1e5;         % reference pressure [Pa]
 kappa    = Rd/cpd;
 gc_ratio = Rd/Rv;

 % saturation vapor pressure [Pa]
 Tc = temp-273.15;
 es = 611.20*exp(17.67*Tc./(Tc+243.5)); % Bolton 1980 equation 10

 % latent heat of condensation [J/kg]
 L = (2.501-0.00237*Tc)*1e6; % Bolton 1980 equation 2

 % saturation mixing ratio
 rs = gc_ratio*es./(p-es);

 % saturation specific humidity 
 qs = rs./(1+rs);

 % potential temperature
 exponent = kappa.*(1+rs./gc_ratio)./(1+rs.*cpv./cpd);
 theta    = temp.*(p0./p).^exponent;

 % derivative of potential temperature with respect to pressure
 dtheta_dp = gradient(theta, p);
 %dtheta_dp = d_dp(theta, p); % which one is more accurate? 
 
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

 % effective static stability following equation 8 of O'Gorman, JAS, 2011
 dtheta_dp_eff = dtheta_dp-lambda.*dtheta_dp_ma;

end
