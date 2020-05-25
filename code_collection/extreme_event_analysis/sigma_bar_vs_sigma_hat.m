
% composite variables vs averaged variables for a single ensemble member

% specify whether spatially verying sigma is used or not
ensemble_path = '~/CESM_LENS/exp/001_2D_varsig_Q_findcenter_correct_precip/';
if exist('OCEAN_DESERT') && OCEAN_DESERT
    ensemble_path = '~/CESM_LENS/exp/001_2D_varsig_Q_ocean_desert/';
end
run([ensemble_path, 'info.m']);
exit = false;
SMOOTH = false;
smooth_window_x = 1;
smooth_window_y = 1;
if SIGMA_2 || VORTICITY
    omega_range = [-1.0, 1.0];
else
    omega_range = [-1.0, 1.0];
end
change_range = [-0.5, 0.5];

if ~exist('sigma_tag')
    sigma_tag = false;
end

plot_level = 50000;
plot_level_name = [num2str(plot_level/100.), 'hPa'];
plevels = [100000, 97500, 95000, 92500, 90000, 87500, 85000, 82500, 80000, 77500, 75000, 70000, 65000, ...
            60000, 55000, 50000, 45000, 40000, 35000, 30000, 25000, 22500, 20000, 17500, 15000]';
plot_ind = plevels == plot_level;
Ra = 287.04;
R = 6371000.0;
cpd= 1005.7;
kappa = Ra / cpd;
figsize = [10 10 400 250];
figsize_zonal = [10 10 600 240];

% add functions from ../../source/
addpath('../../source/');

folder = [ensemble_path, 'event_output/'];
historical_file = 'precip_events_historical.mat';
rcp85_file = 'precip_events_rcp85.mat';

if strfind(folder, 'CESM')
    input_path = '/net/aimsir/archive1/ziweili/CESM_LENS/output/';
    box_max_x = 29;
    box_max_y = 29;
    temp_filename = ['/archive1/ziweili/CESM_LENS/data/', ...
                        'b.e11.B20TRC5CNBDRD.f09_g16.001.cam.h2.PRECT.1990010100Z-2005123118Z.nc'];
    omega_avg_h_name = ['../../output/omega_avg_', num2str(n_ensemble, '%.3d'), '_h.nc'];
    omega_avg_r_name = ['../../output/omega_avg_', num2str(n_ensemble, '%.3d'), '_r.nc'];
elseif strfind(folder, 'GFDL')
    %input_path = '/net/chickpea/volume1/scratch/ziweili/test1_GFDL/output_data/';
    input_path = '/net/aimsir/archive1/ziweili/test1_GFDL/data/';
    box_max_x = 19;
    box_max_y = 19;
    temp_filename = '/disk7/ziweili/test1_GFDL/input_data/orog_fx_GFDL-CM3_historical_r0i0p0.nc';
    omega_avg_h_name = '/net/aimsir/archive1/ziweili/test1_GFDL/data/omega_avg/omega_avg_historical.nc';
    omega_avg_r_name = '/net/aimsir/archive1/ziweili/test1_GFDL/data/omega_avg/omega_avg_rcp85.nc';
end
lat_series = ncread(temp_filename, 'lat');
lon_series = ncread(temp_filename, 'lon');

dlambda = (lon_series(2) - lon_series(1)) / 180.0 * pi;
dphi = (lat_series(2) - lat_series(1)) / 180.0 * pi;

Omega = 7.2921e-5;
f0 = 2 * Omega * sin(lat_series / 180 * 3.1415926);
if strfind(folder, 'CESM')
    string_1_1 = num2str(n_ensemble, '%.3d');
    string_1_2 = string_1_1;
elseif strfind(folder, 'GFDL')
    string_1_1 = 'historical';
    string_1_2 = 'rcp85';
end

if ~exist('events_sta_historical')
    events_sta_historical = load([folder, historical_file], '-mat', 'events_sta');
    num_event_historical = load([folder, historical_file], '-mat', 'num_event');
    num_event_discarded_h = load([folder, historical_file], '-mat', 'num_event_discarded');
    num_event_discarded_misalign_h = load([folder, historical_file], '-mat', 'num_event_discarded_misalign');
end
num_event_historical = load([folder, historical_file], '-mat', 'num_event');
if ~exist('events_sta_rcp85')
    events_sta_rcp85 = load([folder, rcp85_file], '-mat', 'events_sta');
    num_event_rcp85 = load([folder, rcp85_file], '-mat', 'num_event');
    num_event_discarded_r = load([folder, rcp85_file], '-mat', 'num_event_discarded');
    num_event_discarded_misalign_r = load([folder, rcp85_file], '-mat', 'num_event_discarded_misalign');
end
num_event_rcp85 = load([folder, rcp85_file], '-mat', 'num_event');
clear('events_sta');

[events_sta_historical, num_event_historical_2] = ...
    remove_by_omega(events_sta_historical, num_event_historical);
num_event_discarded_h.num_event_discarded = num_event_discarded_h.num_event_discarded + ...
        num_event_historical.num_event - num_event_historical_2.num_event;
num_event_historical = num_event_historical_2;
[events_sta_rcp85, num_event_rcp85_2] = ...
    remove_by_omega(events_sta_rcp85, num_event_rcp85);
num_event_discarded_r.num_event_discarded = num_event_discarded_r.num_event_discarded + ...
        num_event_rcp85.num_event - num_event_rcp85_2.num_event;
num_event_rcp85 = num_event_rcp85_2;
clear('num_event_historical_2', 'num_event_rcp85_2');

if ~exist('VORTICITY') || ...
   (exist('VORTICITY') && ~VORTICITY)
    threshold = 0;
    %[events_sta_historical, num_event_historical] = ...
    %    remove_by_Adv(events_sta_historical, num_event_historical, plot_level, threshold, string_1_1, input_path, SIGMA_2);
    %[events_sta_rcp85, num_event_rcp85] = ...
    %    remove_by_Adv(events_sta_rcp85, num_event_rcp85, plot_level, threshold, string_1_2, input_path, SIGMA_2);
end

latmax = 90; % degrees
latmin = -90;
latmax_plot = 60;
latmin_plot = -60;
lat = zeros(size(events_sta_historical.events_sta, 1), 1);
for j = 1 : length(lat)
    i = 1;
    while isempty(events_sta_historical.events_sta{j, i}) && ...
        i < size(events_sta_historical.events_sta, 2)
        i = i + 1;
    end
    if i ~= size(events_sta_historical.events_sta, 2)
        lat(j) = events_sta_historical.events_sta{j, i}(1).event_lat;
    else
        disp(['no data on j == ', num2str(j)])
        lat(j) = 2 * lat(j - 1) - lat(j - 2);
    end
end
if length(events_sta_historical.events_sta{j, i}(1).event_timespan) == 4
    timespan = 4;
elseif length(events_sta_historical.events_sta{j, i}(1).event_timespan) == 3
    timespan = 1;
end
sigma_eff_tag = ~isempty(events_sta_historical.events_sta{j, i}(1).sigma_eff);
if sigma_eff_tag
    [sigma_eff_h, sigma_eff_r] = deal(nan(length(lat_series), 1));
end

lonmax = 360; % degrees
lonmin = 0;
[~, lon_indices] = latlonindices(lat_series, lon_series, latmin, latmax, lonmin, lonmax);
[~, ~, lat_indices] = intersect(lat, lat_series);
lon = lon_series(lon_indices);
F0 = repmat(f0(lat_indices), [1, length(lon_indices), length(plevels)]);

[k2_h, k2_star_h, l2_h, m2_h, dtheta_dp_ma_h, dtheta_dp_ma_avg_h, ...
 sigma_h, sigma_star_h, omega_sigma_h, omega_k2_sigma_h, ...
 J_h, J_l2_h, d2_dp2_omega_QG_grid_h, omega_QG_h, omega_h, ...
 Adv_h, C_h, A1_h, A2_h, A3_h, B_h, num_event_historical, ...
 k2_m_h, T_adv_h, precip_h] = ...
    comp_composite_k2_l2_m2_grid_v1(timespan, events_sta_historical, num_event_historical, ...
    lat_series, lon_series, box_max_x, box_max_y, sigma_tag, plevels, string_1_1, input_path, SIGMA_2, VORTICITY, TRENBERTH);
[k2_r, k2_star_r, l2_r, m2_r, dtheta_dp_ma_r, dtheta_dp_ma_avg_r, ...
 sigma_r, sigma_star_r, omega_sigma_r, omega_k2_sigma_r, ...
 J_r, J_l2_r, d2_dp2_omega_QG_grid_r, omega_QG_r, omega_r, ...
 Adv_r, C_r, A1_r, A2_r, A3_r, B_r, num_event_rcp85, ...
 k2_m_r, T_adv_r, precip_r] = ...
    comp_composite_k2_l2_m2_grid_v1(timespan, events_sta_rcp85, num_event_rcp85, ...
    lat_series, lon_series, box_max_x, box_max_y, sigma_tag, plevels, string_1_2, input_path, SIGMA_2, VORTICITY, TRENBERTH);

[k2_nc_h, l2_nc_h, m2_nc_h, sigma_nc_h] = ... % '_nc_' means non-composite
        k2_l2_m2_grid(events_sta_historical, lat_series, lon_series, box_max_x, box_max_y, ...
        sigma_tag, plevels, input_path, SIGMA_2);

[k2_nc_r, l2_nc_r, m2_nc_r, sigma_nc_r] = ...
        k2_l2_m2_grid(events_sta_rcp85, lat_series, lon_series, box_max_x, box_max_y, ...
        sigma_tag, plevels, input_path, SIGMA_2);

plot_path_ensemble = ['plots_0.1/bar_vs_hat/'];
colors = get(gca,'colororder');
ind_lat = abs(lat) >= 30;
if exist('OCEAN_DESERT') && OCEAN_DESERT
    ind_lat = ones(size(lat));
end

%% composite k^2 vs averaged k^2
fig = figure('pos', [10, 10, 300, 270]);
temp_1 = k2_h(ind_lat, :, plot_ind);
temp_2 = k2_nc_h(ind_lat, :, plot_ind);
s1 = scatter(temp_2(:), temp_1(:), 1, 'o', 'filled', 'MarkerFaceColor', colors(1, :));
hold on
temp_1 = k2_r(ind_lat, :, plot_ind);
temp_2 = k2_nc_r(ind_lat, :, plot_ind);
s2 = scatter(temp_2(:), temp_1(:), 1.5, 'o', 'filled', 'MarkerFaceColor', colors(2, :));
plot([-10, 10], [-10, 10], 'k--');
legend([s1, s2], {'Historical', 'RCP8.5'}, 'location', 'southeast', 'interpreter', 'latex')
legend boxoff
xlim([0.0, 2e-10])
ylim([0.0, 2e-10])
xlabel('$\hat{k^2}$', 'interpreter', 'latex')
ylabel('$\overline{k^2}$', 'interpreter', 'latex')
set(gca, 'Position', [0.17, 0.15, 0.75, 0.75])
set(gca, 'TickLabelInterpreter','latex')
set(gca, 'TickDir', 'out');
figname = [plot_path_ensemble, 'k_hat_vs_k_bar_scatter_', plot_level_name];
print(fig, figname, '-depsc', '-r0', '-painters')
clf;

%% composite sigma vs averaged sigma
fig = figure('pos', [10, 10, 300, 270]);
temp_1 = sigma_h(ind_lat, :, plot_ind);
temp_2 = sigma_nc_h(ind_lat, :, plot_ind);
s1 = scatter(temp_2(:), temp_1(:), 1, 'o', 'filled', 'MarkerFaceColor', colors(1, :));
hold on
temp_1 = sigma_r(ind_lat, :, plot_ind);
temp_2 = sigma_nc_r(ind_lat, :, plot_ind);
s2 = scatter(temp_2(:), temp_1(:), 1.5, 'o', 'filled', 'MarkerFaceColor', colors(2, :));
plot([-10, 10], [-10, 10], 'k--');
legend([s1, s2], {'Historical', 'RCP8.5'}, 'location', 'southeast', 'interpreter', 'latex')
legend boxoff
xlim([0.0, 6e-6])
ylim([0.0, 6e-6])
xlabel('$\hat{\sigma}$', 'interpreter', 'latex')
ylabel('$\overline{\sigma}$', 'interpreter', 'latex')
set(gca, 'Position', [0.17, 0.15, 0.75, 0.75])
set(gca, 'TickLabelInterpreter','latex')
set(gca, 'TickDir', 'out');
figname = [plot_path_ensemble, 'sigma_hat_vs_sigma_bar_scatter_', plot_level_name];
print(fig, figname, '-depsc', '-r0', '-painters')
clf;

%% composite l2 vs averaged l2
fig = figure('pos', [10, 10, 300, 270]);
temp_1 = l2_h(ind_lat, :, plot_ind);
temp_2 = l2_nc_h(ind_lat, :, plot_ind);
s1 = scatter(temp_2(:), temp_1(:), 1, 'o', 'filled', 'MarkerFaceColor', colors(1, :));
hold on
temp_1 = l2_r(ind_lat, :, plot_ind);
temp_2 = l2_nc_r(ind_lat, :, plot_ind);
s2 = scatter(temp_2(:), temp_1(:), 1.5, 'o', 'filled', 'MarkerFaceColor', colors(2, :));
plot([-10, 10], [-10, 10], 'k--');
legend([s1, s2], {'Historical', 'RCP8.5'}, 'location', 'southeast', 'interpreter', 'latex')
legend boxoff
xlim([0.0, 2e-10])
ylim([0.0, 2e-10])
xlabel('$\hat{l^2}$', 'interpreter', 'latex')
ylabel('$\overline{l^2}$', 'interpreter', 'latex')
set(gca, 'Position', [0.17, 0.15, 0.75, 0.75])
set(gca, 'TickLabelInterpreter','latex')
set(gca, 'TickDir', 'out');
figname = [plot_path_ensemble, 'l2_hat_vs_l2_bar_scatter_', plot_level_name];
print(fig, figname, '-depsc', '-r0', '-painters')
clf;


%% composite m2 vs averaged m2
fig = figure('pos', [10, 10, 300, 270]);
temp_1 = m2_h(ind_lat, :, plot_ind);
temp_2 = m2_nc_h(ind_lat, :, plot_ind);
s1 = scatter(temp_2(:), temp_1(:), 1, 'o', 'filled', 'MarkerFaceColor', colors(1, :));
hold on
temp_1 = m2_r(ind_lat, :, plot_ind);
temp_2 = m2_nc_r(ind_lat, :, plot_ind);
s2 = scatter(temp_2(:), temp_1(:), 1.5, 'o', 'filled', 'MarkerFaceColor', colors(2, :));
plot([-10, 10], [-10, 10], 'k--');
legend([s1, s2], {'Historical', 'RCP8.5'}, 'location', 'southeast', 'interpreter', 'latex')
legend boxoff
xlim([-1e-8, 4e-8])
ylim([-1e-8, 4e-8])
xlabel('$\hat{m^2}$', 'interpreter', 'latex')
ylabel('$\overline{m^2}$', 'interpreter', 'latex')
set(gca, 'Position', [0.17, 0.15, 0.75, 0.75])
set(gca, 'TickLabelInterpreter','latex')
set(gca, 'TickDir', 'out');
figname = [plot_path_ensemble, 'm2_hat_vs_m2_bar_scatter_', plot_level_name];
print(fig, figname, '-depsc', '-r0', '-painters')
clf;



