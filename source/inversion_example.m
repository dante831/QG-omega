
% add functions from source/
addpath('./source/')

% define useful constants
Omega = 7.2921e-5;  % Earth's angular speed in rad/s 
R = 6371000.0;      % Earth's average radius from https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html, in meters
Ra = 287.04;        % specific gas constant, in J/(kg*K)
cp = 1005.0;        % heat capacity at 300K, in J/(kg*K)
kappa = Ra/cp;
dt = 3600.0 * 6.0;  % time difference, in seconds



% shrink_box parameters
box_min = 15; % values for CESM-LE data
box_max = 29;
box_max_x = box_max;
box_max_y = box_max;
computation_y = box_max_y;
computation_x = box_max_x;
p_min = 10000; % upper limit
p_max = 55000; % lower limit, events should at least have full non-NaN values at this level
window_x = 1;       % smoothing window. This translates to a (2*window_x+1)-by-(2*window_y+1) smoothing filter
window_y = 1;



% read data of an example event
nc_filename = './data/example_event.nc';
[T, ua, va, ug, vg, omega, omega_b, ...
 lon, lat, level, event_timespan] = read_ncfile(nc_filename);

p_start_0 = 1; % lowest level
p_start = find(level == p_max); % lower boundary start point
p_end = find(level == p_min); % highest level
level_indices = p_start_0 : p_end;

% add shrink_box, if range as large as box_max has NaN value, then shrink until box_min
[x_ind, y_ind, computation_x, computation_y, p_start, NaN_tag] = ...
        shrink_box(omega, box_min, p_max, computation_x, computation_y, p_start, p_end, level(level_indices));

% subselect all fields and indices according to the result from shrink_box
level_indices = p_start : p_end;
level = level(level_indices);
T       = T      (y_ind, x_ind, level_indices, :);
ug      = ug     (y_ind, x_ind, level_indices, :);
vg      = vg     (y_ind, x_ind, level_indices, :);
ua      = ua     (y_ind, x_ind, level_indices, :);
va      = va     (y_ind, x_ind, level_indices, :);
omega   = omega  (y_ind, x_ind, level_indices, :);
omega_b = omega_b(y_ind, x_ind, level_indices, :);
lat     = lat(y_ind);
lon     = lon(x_ind);
%event_latspan = event_latspan(y_ind);
%event_lonspan = event_lonspan(x_ind);

% define some other useful arrays
phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);

% smooth T field
T_smoothed = T;
for k = 1 : length(level)
    for t = 1 : length(event_timespan)
        T_smoothed(:, :, k, t) = smooth2a(T(:, :, k, t), window_x, window_y);
        % renormalize so that the smoothed field has the same variance as the previous field
        T_smoothed(:, :, k, t) = smooth_normalize(T(:, :, k, t), T_smoothed(:, :, k, t));
    end
end
T = T_smoothed;
clear('T_smoothed');

% define Corriolis parameters
f0 = 2 * Omega * sin(lat((1+end)/2) / 180 * 3.1415926);
beta = 2 * Omega / R * cos(lat((1+end)/2) / 180 * 3.1415926);

% define place-holders for the inversion process
[A, B, C, rhs, J, Qx, Qy, ...
 dug_dlambda, dvg_dlambda, dug_dphi, dvg_dphi, ...
 dT_dlambda, dT_dphi, sigma_accu, zeta] = ...
            deal(zeros(length(lat), length(lon), length(level), length(event_timespan)));
[sigma, dtheta_dp_ma, T_avg] = deal(zeros(length(level), length(event_timespan)));
omega_QG = zeros(length(lat), length(lon), length(level), length(event_timespan));

for t = 1 : length(event_timespan)

    T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));

    % effective static stability: sigma_eff
    temp_omega = squeeze(omega(:, :, :, t));
    omega_prime = temp_omega - repmat(mean(mean(temp_omega, 1), 2), ...
            [size(temp_omega(:, :, 1)), 1]);
    omega_up = temp_omega; omega_up(temp_omega > 0) = 0;
    omega_up_prime = omega_up - repmat(mean(mean(omega_up, 1), 2), ...
            [size(omega_up(:, :, 1)), 1]);
    % equation (5) in O'Gorman, 2011
    lambda_eff = squeeze(mean(mean(omega_up .* omega_up_prime, 1), 2) ./ ...
                     mean(mean(omega_up.^2, 1), 2));
    [~, ~, theta_avg] = eff_stat_stab(level, T_avg(:, t), lambda_eff);

    sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg) .* gradient(theta_avg, level);

    Level = repmat(reshape(level, 1, 1, length(level)), length(lat), length(lon), 1);
    theta = T(:, :, :, t) .* (1e5 ./ Level) .^ kappa;
    sigma_accu(:, :, :, t) = - Ra * T(:, :, :, t) ./ (Level .* theta) .* d_dp(theta, level);

end

%% Q vector form of the advective forcing
GEO_T = true;
dlambda    = mod(lon(2) - lon(1), 360) / 180.0 * 3.1415926; % x angular increment
dphi       = mod(lat(2) - lat(1), 360) / 180.0 * 3.1415926; % y angular increment
[A, B] = Q_vector(level, ug, vg, T, Phi, event_timespan, dphi, dlambda, f0, beta, GEO_T);

%% traditional form of the advective forcing
% get relative vorticity
for t = 1 : length(event_timespan)
    for k = 1 : length(level)
        zeta(:, :, k, t) = curl(phi, ug(:, :, k, t), vg(:, :, k, t), dphi, dlambda);
    end
end
% term A
[A_traditional, ~, ~, ~] = traditional_A(phi, level, ua, va, zeta, event_timespan, dphi, dlambda, f0, beta);
% term B
B_traditional = traditional_B(phi, level, ug, vg, T, event_timespan, dphi, dlambda, f0, GEO_T);

% term C (latent heating)
[C1, C2, J1, J2] = diabatic_heating(phi, T, level, ua, va, sigma_accu, omega, ...
                            event_timespan, dphi, dlambda, dt);
J = J1 + J2;
C = C1 + C2;

% right hand side
rhs = A + B + C;

% remove negative values in sigma_accu
[sigma_accu, tag_sigma] = sigma_remove_negative(sigma_accu, sigma);
if ~tag_sigma; disp('Warning: more than 10% points has negative stability!'); end

% smooth sigma_accu field
sigma_accu_smoothed = sigma_accu;
for k = 1 : length(level)
    for t = 1 : length(event_timespan)
        sigma_accu_smoothed(:, :, k, t) = smooth2a(sigma_accu(:, :, k, t), window_x, window_y);
        sigma_accu_smoothed(:, :, k, t) = smooth_normalize(sigma_accu(:, :, k, t), sigma_accu_smoothed(:, :, k, t));
    end
end
sigma_accu = sigma_accu_smoothed;
clear('sigma_accu_smoothed');

% do the inversion
for t = 1 : length(event_timespan)
    omega_QG(:, :, :, t) = SIP_inversion(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
                                        sigma_accu(:, :, :, t), dphi, dlambda, rhs(:, :, :, t), omega_b(:, :, :, t), false);
end

% load precipitation fields centered around this grid point
load('./data/event_precip.mat')
quantile_file_h = './data/precip_99.9th_quantile_h.mat';
temp = load(quantile_file_h, 'Q');
historical_quantile = temp.Q(177, 44); clear('temp')

% plot the inverted event, the same as figure 1 of Li and O'Gorman, 2019
plot_level = 50000;
plot_path = 'plots/';
filename = [plot_path, 'example_event'];
threeD_event_plot(0, historical_quantile, ...
                    lat, lon, level, omega_QG, event_precip(x_ind, y_ind, :), plot_level, event_timespan, ...
                    lon((1+end)/2), lat((1+end)/2), filename)


