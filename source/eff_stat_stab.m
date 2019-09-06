function theta = potential_temp(p, temp, lambda)

% p       pressure (Pa)
% temp    temperature (K)

 if ~isequal(size(temp), size(p)) | min(size(p))~=1
  error('temp and p must be one-dimensional and of same size')
 end

 % constants
 Rd       = 287.04;      % gas constant for dry air [J/kg/K]
 Rv       = 461.5;       % gas constant for water vapor [J/kg/K]
 cpd      = 1005.7;      % specific heat dry air [J/kg/K]
 cpv      = 1870;        % specific heat water vapor [J/kg/K]
 p0       = 1e5;         % reference pressure [Pa]
 kappa    = Rd/cpd;
 gc_ratio = Rd/Rv;

 % saturation vapor pressure [Pa]
 Tc = temp-273.15;
 es = 611.20*exp(17.67*Tc./(Tc+243.5)); % Bolton 1980 equation 10

 % saturation mixing ratio
 rs = gc_ratio*es./(p-es);

 % potential temperature
 exponent = kappa.*(1+rs./gc_ratio)./(1+rs.*cpv./cpd);
 theta    = temp.*(p0./p).^exponent;

end
