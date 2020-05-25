% plot omega change at a given level between the two climate scenarios

% specify whether spatially verying sigma is used or not
info;
pwd_str = pwd;
exit = false;
SMOOTH = false;
REMOVE_OMEGA = true;
if ~isempty(strfind(pwd_str, 'ocean_desert')) || ~isempty(strfind(pwd_str, 'Adv_center'))
    GT_ZERO = false; % if this is true, we neglect events with omega_QG > 0 at 500hPa
else
    GT_ZERO = true;
end
GT_ZERO = true;
if REMOVE_OMEGA
    matfilename = 'plot_diagnostic.mat';
else
    matfilename = 'plot_diagnostic_no_omega_removal.mat';
end

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

if exist('OCEAN_DESERT') && OCEAN_DESERT
    plot_level = 87500;
else
    OCEAN_DESERT = false;
    plot_level = 50000;
end

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

folder = 'event_output/';
plot_path = 'event_plots/';
historical_file = 'precip_events_historical.mat';
rcp85_file = 'precip_events_rcp85.mat';
%historical_file = 'precip_events_historical_with_dtheta_dp.mat';
%rcp85_file = 'precip_events_rcp85_with_dtheta_dp.mat';

if strfind(pwd_str, 'CESM')
    input_path = '/net/aimsir/archive1/ziweili/CESM_LENS/output/';
    box_max_x = 29;
    box_max_y = 29;
    %temp_filename = ['/net/aimsir/archive1/ziweili/CESM_LENS/data/', ...
    %                    'b.e11.B20TRC5CNBDRD.f09_g16.001.cam.h2.PRECT.1990010100Z-2005123118Z.nc'];
    temp_filename = ['/archive1/ziweili/CESM_LENS/data/', ...
                        'b.e11.B20TRC5CNBDRD.f09_g16.001.cam.h2.PRECT.1990010100Z-2005123118Z.nc'];
    omega_avg_h_name = ['../../output/omega_avg_', num2str(n_ensemble, '%.3d'), '_h.nc'];
    omega_avg_r_name = ['../../output/omega_avg_', num2str(n_ensemble, '%.3d'), '_r.nc'];
elseif strfind(pwd_str, 'GFDL')
    %input_path = '/net/chickpea/volume1/scratch/ziweili/test1_GFDL/output_data/';
    input_path = '/net/aimsir/archive1/ziweili/test1_GFDL/data/';
    box_max_x = 19;
    box_max_y = 19;
    temp_filename = '/disk7/ziweili/test1_GFDL/input_data/orog_fx_GFDL-CM3_historical_r0i0p0.nc';
    %omega_avg_h_name = ['/net/chickpea/volume1/scratch/ziweili/test1_GFDL/omega/omega_avg_historical.nc'];
    %omega_avg_r_name = ['/net/chickpea/volume1/scratch/ziweili/test1_GFDL/omega/omega_avg_rcp85.nc'];
    omega_avg_h_name = '/net/aimsir/archive1/ziweili/test1_GFDL/data/omega_avg/omega_avg_historical.nc';
    omega_avg_r_name = '/net/aimsir/archive1/ziweili/test1_GFDL/data/omega_avg/omega_avg_rcp85.nc';
end
lat_series = ncread(temp_filename, 'lat');
lon_series = ncread(temp_filename, 'lon');

dlambda = (lon_series(2) - lon_series(1)) / 180.0 * pi;
dphi = (lat_series(2) - lat_series(1)) / 180.0 * pi;

Omega = 7.2921e-5;
f0 = 2 * Omega * sin(lat_series / 180 * 3.1415926);
if strfind(pwd_str, 'CESM')
    string_1_1 = num2str(n_ensemble, '%.3d');
    string_1_2 = string_1_1;
elseif strfind(pwd_str, 'GFDL')
    string_1_1 = 'historical';
    string_1_2 = 'rcp85';
end

if ~exist('events_sta_historical')
    events_sta_historical = load([folder, historical_file], '-mat', 'events_sta');
    num_event_historical = load([folder, historical_file], '-mat', 'num_event');
    num_event_discarded_h = load([folder, historical_file], '-mat', 'num_event_discarded');
    num_event_discarded_misalign_h = load([folder, historical_file], '-mat', 'num_event_discarded_misalign');
    %events_sta = calculate_dtheta_dp_ma_avg(events_sta, string_1_1, input_path);
    %matfilename = 'precip_events_historical_with_dtheta_dp.mat';
    %events_sta_historical = events_sta;
    %events_sta = events_sta.events_sta;
    %num_event = num_event.num_event;
    %num_event_historical = num_event;
    %save([folder, matfilename], 'num_event', 'plevels', 'lon_series', 'lat_series', 'events_sta', '-v7.3');
end
if ~exist('events_sta_rcp85')
    events_sta_rcp85 = load([folder, rcp85_file], '-mat', 'events_sta');
    num_event_rcp85 = load([folder, rcp85_file], '-mat', 'num_event');
    num_event_discarded_r = load([folder, rcp85_file], '-mat', 'num_event_discarded');
    num_event_discarded_misalign_r = load([folder, rcp85_file], '-mat', 'num_event_discarded_misalign');
    %events_sta = calculate_dtheta_dp_ma_avg(events_sta, string_1_2, input_path);
    %matfilename = 'precip_events_rcp85_with_dtheta_dp.mat';
    %events_sta_rcp85 = events_sta;
    %num_event_rcp85 = num_event;
    %events_sta = events_sta.events_sta;
    %num_event = num_event.num_event;
    %save([folder, matfilename], 'num_event', 'plevels', 'lon_series', 'lat_series', 'events_sta', '-v7.3');
end
clear('events_sta');

if REMOVE_OMEGA
    [events_sta_historical, num_event_historical_2] = ...
        remove_by_omega(events_sta_historical, num_event_historical, GT_ZERO, plot_level);
    num_event_discarded_h.num_event_discarded = num_event_discarded_h.num_event_discarded + ...
            num_event_historical.num_event - num_event_historical_2.num_event;
    num_event_historical = num_event_historical_2;
    [events_sta_rcp85, num_event_rcp85_2] = ...
        remove_by_omega(events_sta_rcp85, num_event_rcp85, GT_ZERO, plot_level);
    num_event_discarded_r.num_event_discarded = num_event_discarded_r.num_event_discarded + ...
            num_event_rcp85.num_event - num_event_rcp85_2.num_event;
    num_event_rcp85 = num_event_rcp85_2;
    clear('num_event_historical_2', 'num_event_rcp85_2');
end


if ~exist('VORTICITY') || ...
   (exist('VORTICITY') && ~VORTICITY)
    threshold = 0;
    %[events_sta_historical, num_event_historical] = ...
    %    remove_by_Adv(events_sta_historical, num_event_historical, plot_level, threshold, string_1_1, input_path, SIGMA_2);
    %[events_sta_rcp85, num_event_rcp85] = ...
    %    remove_by_Adv(events_sta_rcp85, num_event_rcp85, plot_level, threshold, string_1_2, input_path, SIGMA_2);
end

latmax_plot = 60;
latmin_plot = -60;

%{
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
        %lat(j) = 2 * lat(j - 1) - lat(j - 2);
        lat(j) = NaN;
    end
end
%}

first_event_index = min(find(~cellfun(@isempty, {events_sta_historical.events_sta{:}})));

if length(events_sta_historical.events_sta{first_event_index}(1).event_timespan) == 4
    timespan = 4;
elseif length(events_sta_historical.events_sta{first_event_index}(1).event_timespan) == 3
    timespan = 1;
end
sigma_eff_tag = ~isempty(events_sta_historical.events_sta{first_event_index}(1).sigma_eff);
if sigma_eff_tag
    [sigma_eff_h, sigma_eff_r] = deal(nan(length(lat_series), 1));
end

[latmax, latmin, lonmax, lonmin] = get_lat_lon_limits(pwd_str);

[lat_indices, lon_indices] = latlonindices(lat_series, lon_series, latmin, latmax, lonmin, lonmax);
%[~, ~, lat_indices] = intersect(lat, lat_series);
lon = lon_series(lon_indices);
lat = lat_series(lat_indices);
F0 = repmat(f0(lat_indices), [1, length(lon_indices), length(plevels)]);

[k2_h, k2_star_h, l2_h, m2_h, dtheta_dp_ma_h, dtheta_dp_ma_avg_h, ...
 sigma_h, sigma_star_h, omega_sigma_h, omega_k2_sigma_h, ...
 J_h, J_l2_h, d2_dp2_omega_QG_grid_h, omega_QG_h, omega_h, ...
 Adv_h, C_h, A1_h, A2_h, A3_h, B_h, num_event_historical, ...
 k2_m_h, T_adv_h, precip_h, omega_500hPa_std_h, omega_QG_500hPa_std_h] = ...
    comp_composite_k2_l2_m2_grid_v1(timespan, events_sta_historical, num_event_historical, ...
    lat_series, lon_series, box_max_x, box_max_y, sigma_tag, plevels, string_1_1, input_path, SIGMA_2, VORTICITY, TRENBERTH);
[k2_r, k2_star_r, l2_r, m2_r, dtheta_dp_ma_r, dtheta_dp_ma_avg_r, ...
 sigma_r, sigma_star_r, omega_sigma_r, omega_k2_sigma_r, ...
 J_r, J_l2_r, d2_dp2_omega_QG_grid_r, omega_QG_r, omega_r, ...
 Adv_r, C_r, A1_r, A2_r, A3_r, B_r, num_event_rcp85, ...
 k2_m_r, T_adv_r, precip_r, omega_500hPa_std_r, omega_QG_500hPa_std_r] = ...
    comp_composite_k2_l2_m2_grid_v1(timespan, events_sta_rcp85, num_event_rcp85, ...
    lat_series, lon_series, box_max_x, box_max_y, sigma_tag, plevels, string_1_2, input_path, SIGMA_2, VORTICITY, TRENBERTH);

% get non-composite variables
[k2_nc_h, l2_nc_h, m2_nc_h, sigma_nc_h] = ... % '_nc_' means non-composite
        k2_l2_m2_grid(events_sta_historical, lat_series, lon_series, box_max_x, box_max_y, ...
        sigma_tag, plevels, input_path, SIGMA_2);

[k2_nc_r, l2_nc_r, m2_nc_r, sigma_nc_r] = ...
        k2_l2_m2_grid(events_sta_rcp85, lat_series, lon_series, box_max_x, box_max_y, ...
        sigma_tag, plevels, input_path, SIGMA_2);

omega_avg_h = ncread(omega_avg_h_name, 'omega_avg');
omega_avg_r = ncread(omega_avg_r_name, 'omega_avg');
plevels_2 = ncread(omega_avg_h_name, 'level');

% get boundary contributions
if ~isempty(events_sta_historical.events_sta{first_event_index}(1).omega_QG_full)
    BOUNDARY = true;
    [omega_QG_interior_h, omega_QG_b_h, ...
     omega_QG_lateral_b_h, omega_QG_lateral_clim_b_h] = ...
        composite_boundary_contributions(timespan, events_sta_historical, num_event_historical, ...
        lat_series, lon_series, box_max_x, box_max_y, plevels, string_1_2, input_path);
    [omega_QG_interior_r, omega_QG_b_r, ...
     omega_QG_lateral_b_r, omega_QG_lateral_clim_b_r] = ...
        composite_boundary_contributions(timespan, events_sta_rcp85, num_event_rcp85, ...
        lat_series, lon_series, box_max_x, box_max_y, plevels, string_1_2, input_path);
else
    BOUNDARY = false;
end

%{
% boundary contributions (no need of this section, may safely delete)
if BOUNDARY
    omega_QG_interior_h = omega_QG_interior_h(:, :, plevels == plot_level);
    omega_QG_interior_r = omega_QG_interior_r(:, :, plevels == plot_level);
    omega_QG_b_h = omega_QG_b_h(:, :, plevels == plot_level);
    omega_QG_b_r = omega_QG_b_r(:, :, plevels == plot_level);
    omega_QG_lateral_b_h = omega_QG_lateral_b_h(:, :, plevels == plot_level);
    omega_QG_lateral_b_r = omega_QG_lateral_b_r(:, :, plevels == plot_level);
    omega_QG_lateral_clim_b_h = omega_QG_lateral_clim_b_h(:, :, plevels == plot_level);
    omega_QG_lateral_clim_b_r = omega_QG_lateral_clim_b_r(:, :, plevels == plot_level);
end
%}

[N_lat, N_lon] = size(events_sta_historical.events_sta);
[Le_h, Le_r, LR_h, LR_r] = deal(zeros(N_lat, N_lon));
%[sigma_center_h, sigma_center_r] = deal(zeros(N_lat, N_lon));

for j = 1 : N_lat
    disp(['calculating Rossby deformation radius LR, j = ', num2str(j)])
    for i = 1 : N_lon
        % historical
        n_event = length(events_sta_historical.events_sta{j, i}(:));
        if n_event == 0
            [LR_h(j, i)] = deal(NaN);
            %[LR_h(j, i), sigma_center_h(j, i)] = deal(NaN);
        else
            [temp_sigma_center, temp_C] = deal(zeros(n_event, 1));
            %{
            for n = 1 : n_event
                level = events_sta_historical.events_sta{j, i}(n).event_level;
                ind = level == plot_level;
                %temp_sigma_center(n) = mean(events_sta_historical.events_sta{j, i}(n).sigma_accu(ind, :));
                %temp_C(n)            = mean(events_sta_historical.events_sta{j, i}(n).C         (ind, :));
            end
            %}
            LR_h(j, i)           = sqrt(sigma_h(j, i, plot_ind)) * (max(plevels) - 20000) / abs(f0(lat_indices(j)));
            if ~isreal(LR_h(j, i))
                LR_h(j, i) = NaN;
            end
            %sigma_center_h(j, i) = mean(temp_sigma_center);
        end
        % rcp85
        n_event = length(events_sta_rcp85.events_sta{j, i}(:));
        if n_event == 0
            [LR_r(j, i)] = deal(NaN);
            %[LR_r(j, i), sigma_center_r(j, i)] = deal(NaN);
        else
            [temp_sigma_center, temp_C] = deal(zeros(n_event, 1));
            %{
            for n = 1 : n_event
                level = events_sta_rcp85.events_sta{j, i}(n).event_level;
                ind = level == plot_level;
                %temp_sigma_center(n) = mean(events_sta_rcp85.events_sta{j, i}(n).sigma_accu(ind, :));
                %temp_C(n)            = mean(events_sta_rcp85.events_sta{j, i}(n).C         (ind, :));
            end
            %}
            LR_r(j, i)           = sqrt(sigma_r(j, i, plot_ind)) * (max(plevels) - 200) / f0(lat_indices(j));
            if ~isreal(LR_r(j, i))
                LR_r(j, i) = NaN;
            end
            %sigma_center_r(j, i) = mean(temp_sigma_center);
        end
    end
end

[omega_500hPa_composite_grid_h, omega_500hPa_composite_grid_r] = deal(cell(size(events_sta_historical.events_sta)));

for j = 1 : N_lat
    event_latind = find(lat_series == lat(j));
    phi = lat_series(event_latind - (box_max_y - 1) / 2 : event_latind + (box_max_y - 1) / 2) / 180.0 * 3.1415926;
    for i = 1 : N_lon
        
        % historical
        n_event = length(events_sta_historical.events_sta{j, i}(:));
        if n_event == 0
            [Le_h(j, i)] = deal(NaN);
        else
            [omega_QG_500hPa_composite, omega_500hPa_composite, J_500hPa_composite, event_omega_avg] = ...
                    deal(zeros(box_max_y, box_max_x));
            counter = zeros(box_max_y, box_max_x);
            for n = 1 : n_event
                dim = size(events_sta_historical.events_sta{j, i}(n).omega_500hPa_full);
                y_ind = box_max_y/2 + (-dim(1)/2 + 1 :  dim(1)/2);
                x_ind = box_max_x/2 + (-dim(2)/2 + 1 :  dim(2)/2);
                omega_500hPa_composite   (y_ind, x_ind) = omega_500hPa_composite   (y_ind, x_ind) + ...
                    squeeze(mean(events_sta_historical.events_sta{j, i}(n).omega_500hPa_full, 4));
                omega_QG_500hPa_composite(y_ind, x_ind) = omega_QG_500hPa_composite(y_ind, x_ind) + ...
                    events_sta_historical.events_sta{j, i}(n).omega_QG_500hPa_full;
                J_500hPa_composite       (y_ind, x_ind) = J_500hPa_composite       (y_ind, x_ind) + ...
                    events_sta_historical.events_sta{j, i}(n).J_500hPa;
                counter(y_ind, x_ind) = counter(y_ind, x_ind) + 1;
            end
            omega_500hPa_composite      = omega_500hPa_composite    ./ counter;
            omega_QG_500hPa_composite   = omega_QG_500hPa_composite ./ counter;
            J_500hPa_composite          = J_500hPa_composite        ./ counter;
            omega_500hPa_composite_grid_h{j, i} = omega_500hPa_composite;
        
            n_lon = size(omega_QG_500hPa_composite, 2);
            n_lat = size(omega_QG_500hPa_composite, 1);
            temp_lon_indices = mod((lon_indices(i) - (n_lon - 1) / 2 : lon_indices(i) + (n_lon - 1) / 2) - 1, length(lon)) + 1;
            temp_lat_indices = (lat_indices(j) - (n_lat - 1) / 2 : lat_indices(j) + (n_lat - 1) / 2);
            % get the e-folding distance at 500hPa, consistent with Tandon et al, 2018
            event_omega_avg = omega_avg_h(temp_lon_indices, temp_lat_indices, plevels_2 == 50000)'; 
            % use omega_QG, since the event is centered on omega_QG field, not omega field
            %Le_h(j, i) = e_fold_v1(omega_QG_500hPa_composite, event_omega_avg, dphi, dlambda, lat(j)/180*pi, R); % e-folding distance in meters
            % prefer to use omega_QG_500hPa_composite to get Le because the event is centered on it
            Le_h(j, i) = e_fold_v1(omega_500hPa_composite, event_omega_avg, dphi, dlambda, lat(j)/180*pi, R);
        end

        % rcp85
        n_event = length(events_sta_rcp85.events_sta{j, i}(:));
        if n_event == 0
            [Le_h(j, i)] = deal(NaN);
        else
            [omega_QG_500hPa_composite, omega_500hPa_composite, J_500hPa_composite] = deal(zeros(box_max_y, box_max_x));
            counter = zeros(box_max_y, box_max_x);
            for n = 1 : n_event
                dim = size(events_sta_rcp85.events_sta{j, i}(n).omega_500hPa_full);
                y_ind = box_max_y/2 + (-dim(1)/2 + 1 :  dim(1)/2);
                x_ind = box_max_x/2 + (-dim(2)/2 + 1 :  dim(2)/2);
                omega_500hPa_composite   (y_ind, x_ind) = omega_500hPa_composite   (y_ind, x_ind) + ...
                    squeeze(mean(events_sta_rcp85.events_sta{j, i}(n).omega_500hPa_full, 4));
                omega_QG_500hPa_composite(y_ind, x_ind) = omega_QG_500hPa_composite(y_ind, x_ind) + ...
                    events_sta_rcp85.events_sta{j, i}(n).omega_QG_500hPa_full;
                J_500hPa_composite       (y_ind, x_ind) = J_500hPa_composite       (y_ind, x_ind) + ...
                    events_sta_rcp85.events_sta{j, i}(n).J_500hPa;
                counter(y_ind, x_ind) = counter(y_ind, x_ind) + 1;
            end
            omega_500hPa_composite      = omega_500hPa_composite         ./ counter;
            omega_QG_500hPa_composite   = omega_QG_500hPa_composite ./ counter;
            J_500hPa_composite          = J_500hPa_composite        ./ counter;
            omega_500hPa_composite_grid_r{j, i} = omega_500hPa_composite;
        
            n_lon = size(omega_QG_500hPa_composite, 2);
            n_lat = size(omega_QG_500hPa_composite, 1);
            temp_lon_indices = mod((lon_indices(i) - (n_lon - 1) / 2 : lon_indices(i) + (n_lon - 1) / 2) - 1, length(lon)) + 1;
            temp_lat_indices = (lat_indices(j) - (n_lat - 1) / 2 : lat_indices(j) + (n_lat - 1) / 2);
            % get the e-folding distance at 500hPa, consistent with Tandon et al, 2018
            event_omega_avg = omega_avg_r(temp_lon_indices, temp_lat_indices, plevels_2 == 50000)';
            % use omega_QG, since the event is centered on omega_QG field, not omega field
            %Le_r(j, i) = e_fold_v1(omega_QG_500hPa_composite, event_omega_avg, dphi, dlambda, lat(j)/180*pi, R); % e-folding distance in meters
            % prefer to use omega_QG_500hPa_composite to get Le because the event is centered on it
            Le_r(j, i) = e_fold_v1(omega_500hPa_composite, event_omega_avg, dphi, dlambda, lat(j)/180*pi, R);
        end
    end
end

%% calculate how close to saturation of model 005

if strfind(pwd_str, '005')
    [RH_h, RH_r] = deal(cell(size(events_sta_historical.events_sta)));
    for j = 1 : N_lat    
        event_latind = find(lat_series == lat(j));
        phi = lat_series(event_latind - (box_max_y - 1) / 2 : event_latind + (box_max_y - 1) / 2) / 180.0 * 3.1415926;            
        for i = 1 : N_lon
            % historical
            n_event = length(events_sta_historical.events_sta{j, i}(:));
            if n_event == 0
                RH_h{j, i} = [];
            else
                RH_h{j, i} = zeros(n_event, 1);
                for n = 1 : n_event
                    level = events_sta_historical.events_sta{j, i}(n).event_level;
                    ind = level == plot_level;
                    q1 = events_sta_historical.events_sta{j, i}(n).q(ind);
                    T1 = events_sta_historical.events_sta{j, i}(n).T(ind);
                    RH_h{j, i}(n) = q_to_RH(q1, T1, plot_level);
                end
            end

            % rcp85
            n_event = length(events_sta_rcp85.events_sta{j, i}(:));
            if n_event == 0
                RH_r{j, i} = [];
            else
                RH_r{j, i} = zeros(n_event, 1);
                for n = 1 : n_event
                    level = events_sta_rcp85.events_sta{j, i}(n).event_level;
                    ind = level == plot_level;
                    q1 = events_sta_rcp85.events_sta{j, i}(n).q(ind);
                    T1 = events_sta_rcp85.events_sta{j, i}(n).T(ind);
                    RH_r{j, i}(n) = q_to_RH(q1, T1, plot_level);
                end 
            end 
        end
    end

    [RH_h_mean, RH_r_mean, sat_h, sat_r, num_h, num_r] = deal(zeros(size(RH_h)));
    for j = 1 : length(RH_h(:))
        RH_h_mean(j) = nanmean(RH_h{j}(:));
        RH_r_mean(j) = nanmean(RH_r{j}(:));
        sat_h(j) = sum(RH_h{j}(:) > 0.90);
        sat_r(j) = sum(RH_r{j}(:) > 0.90);
    end
    % percentage
    sum(sat_h(:)) / sum(num_event_historical.num_event(:))
    sum(sat_r(:)) / sum(num_event_rcp85.num_event(:))
end






set(0,'DefaultFigureVisible','off');
figsize = [10, 10, 1100, 350];
[X, Y] = meshgrid(lon, lat);

temp = num_event_historical.num_event;
plot_title = 'historical number of events';
plot_filename = 'num_event_historical';
climit = [0, 50];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = num_event_rcp85.num_event;
plot_title = 'rcp85 number of events';
plot_filename = 'num_event_rcp85';
climit = [0, 50];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);


temp = Le_h;
plot_title = 'historical $Le$, (meters)';
plot_filename = 'Le_h';
climit = [0, 6.0e5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = Le_r;
plot_title = 'rcp85 $Le$, (meters)';
plot_filename = 'Le_r';
climit = [0, 6.0e5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);


temp = LR_h;
plot_title = 'historical Rossby deformation radius $L_R$ (meters)';
plot_filename = 'LR_h';
climit = [0, 5.0e6];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = (Le_r.^2 - Le_h.^2)./Le_h.^2;
plot_title = '$\Delta Le^2/Le^2$';
plot_filename = 'Le_change';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = omega_h(:, :, plot_ind);
plot_title = '$\omega$';
plot_filename = 'omega_h';
climit = omega_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = omega_r(:, :, plot_ind);
plot_title = '$\omega$';
plot_filename = 'omega_r';
climit = omega_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = m2_h(:, :, plot_ind);
plot_title = 'historical $m^2$';
plot_filename = 'm2_h';
climit = [-1e-8, 1e-8];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = m2_r(:, :, plot_ind);
plot_title = 'rcp85 $m^2$';
plot_filename = 'm2_r';
climit = [-1e-8, 1e-8];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);


temp = omega_r(:, :, plot_ind) - omega_h(:, :, plot_ind);
plot_title = '$\Delta\omega$ (Pa/s)';
plot_filename = 'omega_change';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = omega_QG_r(:, :, plot_ind) - omega_QG_h(:, :, plot_ind);
plot_title = '$\Delta\omega_{QG}$ (Pa/s)';
plot_filename = 'omega_QG_change';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = - (k2_r(:, :, plot_ind) - k2_h(:, :, plot_ind)) ./ k2_h(:, :, plot_ind);
plot_title = '$- \Delta k^2/k^2$';
plot_filename = 'k2_change';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = - (dtheta_dp_ma_r(:, :, plot_ind) - dtheta_dp_ma_h(:, :, plot_ind)) ./ dtheta_dp_ma_h(:, :, plot_ind);
plot_title = '$- \Delta \frac{T}{\theta}\frac{d\theta}{dp}|_{\theta^*}/\frac{T}{\theta}\frac{d\theta}{dp}|_{\theta^*}$';
plot_filename = 'dtheta_dp_ma_change';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = - (sigma_r(:, :, plot_ind) - sigma_h(:, :, plot_ind)) ./ sigma_h(:, :, plot_ind);
plot_title = '$- \Delta \sigma/\sigma$';
plot_filename = 'sigma_change';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = Adv_r(:, :, plot_ind) - Adv_h(:, :, plot_ind);
plot_title = '$\Delta RHS$';
plot_filename = 'Adv_change';
climit = [-1e-16, 1e-16];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = Adv_h(:, :, plot_ind);
plot_title = 'historical Adv';
plot_filename = 'Adv_h';
climit = [-1e-16, 1e-16];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = Adv_r(:, :, plot_ind);
plot_title = 'rcp85 Adv';
plot_filename = 'Adv_r';
climit = [-1e-16, 1e-16];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);



%% reconstruction of omega_QG
% first, test if the QG omega equation holds for the local grid analysis
% result: the equation holds well for historical and rcp85
error_h = - k2_h .* omega_QG_h .* sigma_h - m2_h .* F0.^2 .* omega_QG_h - Adv_h - C_h;
error_r = - k2_r .* omega_QG_r .* sigma_r - m2_r .* F0.^2 .* omega_QG_r - Adv_r - C_r;
figure('pos', figsize)
%surf(X, Y, error_h(:, :, plot_ind))
surf(X, Y, error_r(:, :, plot_ind))
view(2)
caxis([-1.0e-18, 1.0e-18])
colorbar
figure('pos', figsize)
%surf(X, Y, Adv_h(:, :, plot_ind))
surf(X, Y, Adv_r(:, :, plot_ind))
view(2)
caxis([-1.0e-18, 1.0e-18])
colorbar

% second, test if approximation omega*dtheta_dp|theta^* is a good approximation to J
% result: it's incredibly good!
temp = J_h(:, :, plot_ind);
plot_title = 'historical J';
plot_filename = 'J_h';
climit = [-1.0, 1.0];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = J_r(:, :, plot_ind);
plot_title = 'rcp85 J';
plot_filename = 'J_r';
climit = [-1.0, 1.0];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = omega_h(:, :, plot_ind) .* dtheta_dp_ma_h(:, :, plot_ind) * cpd;
plot_title = 'historical $c_p\frac{d\theta}{dp}|_{\theta^*}\omega$';
plot_filename = 'J_estimation_h';
climit = [0, 1.0];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = omega_r(:, :, plot_ind) .* dtheta_dp_ma_r(:, :, plot_ind) * cpd;
plot_title = 'rcp85 $c_p\frac{d\theta}{dp}|_{\theta^*}\omega$';
plot_filename = 'J_estimation_r';
climit = [0, 1.0];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

figure('pos', [10, 10, 500, 450])
%ind_lat_scatter = abs(lat) >= 30;
ind_lat_scatter = abs(lat) >= 30 & abs(lat) <= 50;
J_h_temp = J_h(ind_lat_scatter, :, plot_ind);
J_r_temp = J_r(ind_lat_scatter, :, plot_ind);
omega_dtheta_dp_h_temp = omega_h(ind_lat_scatter, :, plot_ind) .* dtheta_dp_ma_h(ind_lat_scatter, :, plot_ind) * cpd;
omega_dtheta_dp_r_temp = omega_r(ind_lat_scatter, :, plot_ind) .* dtheta_dp_ma_r(ind_lat_scatter, :, plot_ind) * cpd;
s1 = scatter(omega_dtheta_dp_h_temp(:), J_h_temp(:), 1.5, 'o', 'filled');
hold on
s2 = scatter(omega_dtheta_dp_r_temp(:), J_r_temp(:), 1.5, 'o', 'filled');
xlim([-0.05, 1.2])
ylim([-0.05, 1.2])
legend([s1, s2], {'Historical', 'RCP8.5'}, 'location', 'southeast', 'interpreter', 'latex')
xlabel('$c_p\frac{d\theta}{dp}|_\theta^*\omega$', 'interpreter', 'latex')
ylabel('$J$', 'interpreter', 'latex')
title('J vs $c_p\frac{d\theta}{dp}|_\theta^*\omega$', 'interpreter', 'latex')
saveas(gca, [plot_path, 'J_estimation_scatter'], 'png')
clf;

% third, reconstruction of omega_QG
% result: it agrees with omega_QG perfectly for historical and rcp
Plevels = repmat(reshape(plevels', 1, 1, length(plevels)), size(Adv_h, 1), size(Adv_h, 2), 1);
if SIGMA_2
    omega_QG_h_rec =  - 1 ./ (k2_h .* sigma_h + m2_h .* F0.^2) .* Adv_h;
    omega_QG_r_rec =  - 1 ./ (k2_r .* sigma_r + m2_r .* F0.^2) .* Adv_r;
elseif VORTICITY
    omega_QG_h_rec =  - 1 ./ (m2_h .* F0.^2) .* Adv_h;
    omega_QG_r_rec =  - 1 ./ (m2_r .* F0.^2) .* Adv_r;
else
    omega_QG_h_rec =  - 1 ./ (k2_h .* sigma_h + m2_h .* F0.^2) .* (Adv_h + l2_h.*kappa./Plevels.*J_h);
    omega_QG_r_rec =  - 1 ./ (k2_r .* sigma_r + m2_r .* F0.^2) .* (Adv_r + l2_r.*kappa./Plevels.*J_r);
end

temp = omega_QG_h(:, :, plot_ind);
plot_title = 'historical $\omega_{QG}$';
plot_filename = 'omega_QG_h';
climit = omega_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = omega_QG_r(:, :, plot_ind);
plot_title = 'rcp85 $\omega_{QG}$';
plot_filename = 'omega_QG_r';
climit = omega_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

% comparison between omega and omega_QG in the two climates
figure('pos',figsize_zonal);
colors = get(gca,'colororder');
hold on;
grid on;
tags_ocean = ones([size(omega_h, 1), size(omega_h, 2)]);
level_ind = plevels == plot_level;
omega_h_mean_2    = nanmean(omega_h   (:, :, plot_ind).*tags_ocean, 2);
omega_r_mean_2    = nanmean(omega_r   (:, :, plot_ind).*tags_ocean, 2);
omega_QG_h_mean_2 = nanmean(omega_QG_h(:, :, plot_ind).*tags_ocean, 2);
omega_QG_r_mean_2 = nanmean(omega_QG_r(:, :, plot_ind).*tags_ocean, 2);

plot_x = lat_series(lat_indices);
p1 = plot(plot_x, - omega_h_mean_2);
p2 = plot(plot_x, - omega_QG_h_mean_2);
p3 = plot(plot_x, - omega_r_mean_2);
p4 = plot(plot_x, - omega_QG_r_mean_2);

lgd = legend( ...
[p1, p2, p3, p4], {...
'Historical $\omega$', 'Historical $\omega_{QG}$', ...
'RCP8.5 $\omega$'    , 'RCP8.5 $\omega_{QG}$'}, ...
'location', 'best', ...
'interpreter', 'latex');

title('Comparison of zonal averaged \omega');
axis([-70 70 -0.2 1.0]);
xlabel('Latitude');
ylabel('-\omega (Pa/s)')
hold off;

saveas(gca, [plot_path, 'zonal_omega_comparison_', num2str(plevels(level_ind)/100), 'hPa'], 'png');


temp = omega_QG_h_rec(:, :, plot_ind);
plot_title = 'reconstructed historical $\omega_{QG}$';
plot_filename = 'omega_QG_h_rec';
climit = omega_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = omega_QG_r_rec(:, :, plot_ind);
plot_title = 'reconstructed rcp85 $\omega_{QG}$';
plot_filename = 'omega_QG_r_rec';
climit = omega_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

% fourth, check l2 and k2
temp = k2_h(:, :, plot_ind);
plot_title = 'historical $k^2$';
plot_filename = 'k2_h';
climit = [-1.5e-10, 1.5e-10];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = k2_r(:, :, plot_ind);
plot_title = 'rcp85 $k^2$';
plot_filename = 'k2_r';
climit = [-1.5e-10, 1.5e-10];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = l2_h(:, :, plot_ind);
plot_title = 'historical $l^2$';
plot_filename = 'l2_h';
climit = [-1.5e-10, 1.5e-10];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = l2_r(:, :, plot_ind);
plot_title = 'rcp85 $l^2$';
plot_filename = 'l2_r';
climit = [-1.5e-10, 1.5e-10];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);


%% fifth, dry decomposition of omega_QG_h_rec

% reference:
% omega_QG_h_rec =  - 1 ./ (k2_h .* sigma_h + m2_h .* F0.^2) .* (Adv_h + l2_h.*kappa./Plevels.*J_h);

d_sigma = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* (sigma_r - sigma_h) .* k2_h .* omega_QG_h_rec;
d_k2    = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* sigma_h .* (k2_r - k2_h) .* omega_QG_h_rec;
d_l2    = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* kappa ./ Plevels .* (l2_r - l2_h) .* J_h; 
d_J     = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* kappa ./ Plevels .* l2_h .* (J_r - J_h);
% for DJF in sigma, l2 is not well-behaved
%d_l2    = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* kappa ./ Plevels .* (k2_r - k2_h) .* J_h;
%d_J     = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* kappa ./ Plevels .* k2_h .* (J_r - J_h);
if SIGMA_2
    d_J(:) = 0;
    d_l2(:) = 0;
end
if VORTICITY
    d_m2  = - 1 ./ (F0.^2.*m2_h) .* F0.^2 .* (m2_r - m2_h) .* omega_QG_h_rec;
    d_Adv = - 1 ./ (F0.^2.*m2_h) .* (Adv_r - Adv_h);
else
    d_m2    = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* F0.^2 .* (m2_r - m2_h) .* omega_QG_h_rec;
    d_Adv   = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* (Adv_r - Adv_h);
end
d_dtheta_dp_ma_omega = - 1./(sigma_h .* k2_h + F0.^2.*m2_h) .* ...
        kappa .* cpd ./ Plevels .* l2_h .* (omega_r .* dtheta_dp_ma_r - omega_h .* dtheta_dp_ma_h);

d_rec   = d_sigma + d_k2 + d_l2 + d_J + d_m2 + d_Adv;
if VORTICITY
    d_rec = d_m2 + d_Adv;
end
d_omega_QG = (omega_QG_r - omega_QG_h);
d_k2_l2 = d_k2 + d_l2;

% zonal mean plots

colors = get(gca,'colororder');
plot_x = lat_series(lat_indices);
figure('pos', figsize_zonal);
hold on;
plot(-90:90, zeros(size(-90:90)), 'LineWidth', 0.8, 'Color', [0.3, 0.3, 0.3], 'Linestyle', '--')

d = 1.5;
tags_analyze = ones([size(omega_h, 1), size(omega_h, 2)]);
NaN_matrix_d = get_NaN_matrix(cat(3, d_rec  (:, :, plot_ind), ...
                                  	 d_dtheta_dp_ma_omega(:, :, plot_ind), ...
                                     d_m2   (:, :, plot_ind), ...
                                     d_Adv  (:, :, plot_ind), ...
                                     d_k2   (:, :, plot_ind), ...
                                     d_l2   (:, :, plot_ind), ...
                                     d_k2_l2(:, :, plot_ind)), d, tags_analyze);

%{
temp_var = sqrt(nanvar([reshape(d_rec(:, :, plot_ind)   .* NaN_matrix_d, [length(NaN_matrix_d(:)), 1]), ...
                        reshape(d_m2(:, :, plot_ind)    .* NaN_matrix_d, [length(NaN_matrix_d(:)), 1]), ...
                        reshape(d_Adv(:, :, plot_ind)   .* NaN_matrix_d, [length(NaN_matrix_d(:)), 1]), ...
                        reshape(d_k2(:, :, plot_ind)    .* NaN_matrix_d, [length(NaN_matrix_d(:)), 1]), ...
                        reshape(d_l2(:, :, plot_ind)    .* NaN_matrix_d, [length(NaN_matrix_d(:)), 1]), ...
                        reshape(d_J(:, :, plot_ind)     .* NaN_matrix_d, [length(NaN_matrix_d(:)), 1]), ...
                        reshape(d_k2_l2(:, :, plot_ind) .* NaN_matrix_d, [length(NaN_matrix_d(:)), 1])]));
[temp_ind_2] = abs(d_rec(:, :, plot_ind)  .*NaN_matrix_d - repmat(nanmean(d_rec(:, :, plot_ind)  .*NaN_matrix_d, 2), 1, size(NaN_matrix_d, 2))) > 3*temp_var(1) | ...
               abs(d_m2(:, :, plot_ind)   .*NaN_matrix_d - repmat(nanmean(d_m2(:, :, plot_ind)   .*NaN_matrix_d, 2), 1, size(NaN_matrix_d, 2))) > 3*temp_var(2) | ...
               abs(d_Adv(:, :, plot_ind)  .*NaN_matrix_d - repmat(nanmean(d_Adv(:, :, plot_ind)  .*NaN_matrix_d, 2), 1, size(NaN_matrix_d, 2))) > 3*temp_var(3) | ...
               abs(d_k2(:, :, plot_ind)   .*NaN_matrix_d - repmat(nanmean(d_k2(:, :, plot_ind)   .*NaN_matrix_d, 2), 1, size(NaN_matrix_d, 2))) > 3*temp_var(4) | ...
               abs(d_l2(:, :, plot_ind)   .*NaN_matrix_d - repmat(nanmean(d_l2(:, :, plot_ind)   .*NaN_matrix_d, 2), 1, size(NaN_matrix_d, 2))) > 3*temp_var(5) | ...
               abs(d_J(:, :, plot_ind)    .*NaN_matrix_d - repmat(nanmean(d_J(:, :, plot_ind)    .*NaN_matrix_d, 2), 1, size(NaN_matrix_d, 2))) > 3*temp_var(6) | ...
               abs(d_k2_l2(:, :, plot_ind).*NaN_matrix_d - repmat(nanmean(d_k2_l2(:, :, plot_ind).*NaN_matrix_d, 2), 1, size(NaN_matrix_d, 2))) > 3*temp_var(7);

NaN_matrix_d(temp_ind_2) = NaN;
%}

omega_QG_h_mean_d = nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_d, 2);

p1 = plot(plot_x, nanmean(-d_omega_QG(:, :, plot_ind).* NaN_matrix_d, 2) ./ omega_QG_h_mean_d);
p2 = plot(plot_x, nanmean(-d_sigma(:, :, plot_ind)   .* NaN_matrix_d, 2) ./ omega_QG_h_mean_d);
p3 = plot(plot_x, nanmean(-d_J(:, :, plot_ind)       .* NaN_matrix_d, 2) ./ omega_QG_h_mean_d);
p4 = plot(plot_x, nanmean(-d_k2_l2(:, :, plot_ind)   .* NaN_matrix_d, 2) ./ omega_QG_h_mean_d);
p5 = plot(plot_x, nanmean(-d_m2(:, :, plot_ind)      .* NaN_matrix_d, 2) ./ omega_QG_h_mean_d);
p6 = plot(plot_x, nanmean(-d_Adv(:, :, plot_ind)     .* NaN_matrix_d, 2) ./ omega_QG_h_mean_d);
p7 = plot(plot_x, nanmean(-d_dtheta_dp_ma_omega(:, :, plot_ind) .* NaN_matrix_d, 2) ./ omega_QG_h_mean_d);
p8 = plot(plot_x, nanmean(-d_rec(:, :, plot_ind)     .* NaN_matrix_d, 2) ./ omega_QG_h_mean_d);

set(p1, 'LineWidth', 1.5, 'Color', 'black');
set(p2, 'LineWidth', 1.5, 'Color', colors(1, :));
set(p3, 'LineWidth', 1.5, 'Color', colors(2, :));
set(p4, 'LineWidth', 1.5, 'Color', colors(3, :));
set(p5, 'LineWidth', 1.5, 'Color', colors(5, :));
set(p6, 'LineWidth', 1.5, 'Color', colors(6, :));
set(p7, 'LineWidth', 1.5, 'Color', colors(7, :), 'Linestyle', '--');
set(p8, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', ':');

lgd = legend( ...
[p1, p2, p3, p4, p5, p6, p7, p8], ...
{'$\Delta\omega_{QG}$' , ...
'$\Delta\omega_{\overline{\sigma}}$' , ...
'$\Delta\omega_{\overline{J}}$' , ...
'$\Delta\omega_{k^2_*, l^2_*}$' , ...
'$\Delta\omega_{m^2_*}$' , ...
'$\Delta\omega_{\overline{Adv}}$' , ...
'$\Delta\omega_{\overline{\frac{T}{\theta}\frac{d\theta}{dp}|_{\theta*}\omega}}$' , ...
'$\Delta\omega_{rec}$' , ...
}, 'location', 'best', ...
'interpreter', 'latex');
%lgd.Position = [0.460, 0.3523, 0.1147, 0.31];
title([plot_level_name, ' Dry Decomposition']);
axis([-70 70 -0.601 0.601]);
xlabel('Latitude');
ylabel('-\Delta\omega_{QG}/\omega_{QG}');
hold off;
saveas(gca, [plot_path, 'zonal_omega_decomposition_', plot_level_name], 'png');

% contour plots

temp = d_sigma(:, :, plot_ind);
plot_title = 'contribution of $\sigma$ (Pa/s)';
plot_filename = 'dry_d_sigma';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = d_k2(:, :, plot_ind);
plot_title = 'contribution of $k^2$ (Pa/s)';
plot_filename = 'dry_d_k2';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = d_k2(:, :, plot_ind) + d_l2(:, :, plot_ind);
plot_title = 'contribution of $k^2$ and $l^2$ (Pa/s)';
plot_filename = 'dry_d_k2_l2';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = d_J(:, :, plot_ind);
plot_title = 'contribution of $J$ (Pa/s)';
plot_filename = 'dry_d_J';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = d_Adv(:, :, plot_ind);
plot_title = 'contribution of $Adv$ (Pa/s)';
plot_filename = 'dry_d_Adv';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = d_m2(:, :, plot_ind);
plot_title = 'contribution of $m^2$ (Pa/s)';
plot_filename = 'dry_d_m2';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = d_rec(:, :, plot_ind) .* NaN_matrix_d;
plot_title = 'sum of all terms';
plot_filename = 'dry_d_rec';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);


%% sixth and last, decomposition of the moist formulation of omega_QG_h_rec_moist
% plot decomposition for the paper
d_omega_QG = (omega_QG_r - omega_QG_h);
[latent_h, latent_r, excessive_h, excessive_r, epsilon_h, epsilon_r, ...
 sigma_m_h, sigma_m_r, k2_m_h, k2_m_r, denom_h, denom_r, ...
 d_sigma_m, d_Adv_m, d_k2_l2_m, d_m2_m, d_rec_m, d_rec_m_full, d_J_res, ...
 d_excessive_1, d_excessive_2, d_epsilon, J_res_h, J_res_r] = moist_decomposition(...
        k2_h, k2_r, l2_h, l2_r, m2_h, m2_r, sigma_h, sigma_r, Plevels, F0, ...
        dtheta_dp_ma_h, dtheta_dp_ma_r, omega_h, omega_r, omega_QG_h_rec, omega_QG_r_rec, ...
        J_h, J_r, Adv_h, Adv_r, sigma_star_h, sigma_star_r, k2_h, k2_r);

% linear expansion
d_denom         = (denom_r - denom_h) .* (Adv_h + kappa ./ Plevels .* l2_h .* (excessive_h + epsilon_h)) ./ denom_h.^2;
d_excessive     = - kappa ./ Plevels .* (l2_r.*excessive_r - l2_h.*excessive_h) ./ denom_h;
d_omega         = omega_r - omega_h;
d_rec_all       = d_denom + d_Adv_m + d_epsilon + d_excessive;

%% Reconstruction of omega_QG in the moist formulation
if SIGMA_2
    sigma_m_h = sigma_h;
    [latent_h(:), excessive_h(:), epsilon_h(:)] = deal(0);
end

if SIGMA_2
    sigma_m_r = sigma_r;
    [latent_r(:), excessive_r(:), epsilon_r(:)] = deal(0);
end

temp = epsilon_h(:, :, plot_ind);
plot_title = 'historical $\epsilon$';
plot_filename = 'epsilon_h';
climit = [-0.2, 0.2];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = epsilon_r(:, :, plot_ind);
plot_title = 'rcp85 $\epsilon$';
plot_filename = 'epsilon_r';
climit = [-0.2, 0.2];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

% reconstruction from the moist formulation:
omega_QG_h_rec_moist = - (Adv_h + kappa ./ Plevels .* l2_h .* (excessive_h + epsilon_h)) ./ denom_h;
omega_QG_r_rec_moist = - (Adv_r + kappa ./ Plevels .* l2_r .* (excessive_r + epsilon_r)) ./ denom_r;

% omega_QG_h_rec_accu, omega_QG_h_rec_rf are exactly the same as omega_QG_h_rec since they are 
% from teh same equation!
omega_QG_h_rec_accu = - Adv_h ./ ((sigma_h + C_h ./ (k2_h .* omega_QG_h_rec)).*k2_h + F0.^2.*m2_h);
omega_QG_r_rec_accu = - Adv_r ./ ((sigma_r + C_r ./ (k2_r .* omega_QG_r_rec)).*k2_r + F0.^2.*m2_r);

% contribution from epsilon
omega_QG_h_rec_rf = - kappa ./ Plevels .* l2_h .* (epsilon_h) ./ denom_h;
omega_QG_r_rec_rf = - kappa ./ Plevels .* l2_r .* (epsilon_r) ./ denom_r;

temp = omega_QG_h_rec_rf(:, :, plot_ind);
plot_title = 'Historical omega QG from $\epsilon$';
plot_filename = 'omega_rec_rf_h';
climit = omega_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = omega_QG_r_rec_rf(:, :, plot_ind);
plot_title = 'RCP85 omega QG from $\epsilon$';
plot_filename = 'omega_rec_rf_r';
climit = omega_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = omega_QG_r_rec_rf(:, :, plot_ind) - omega_QG_h_rec_rf(:, :, plot_ind);
plot_title = '$\Delta\omega_{QG}$ from reduced forcing';
plot_filename = 'omega_rec_rf_change';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = omega_QG_h_rec_moist(:, :, plot_ind);
plot_title = 'historical moist reconstruction';
plot_filename = 'omega_rec_moist_h';
climit = omega_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = omega_QG_r_rec_moist(:, :, plot_ind);
plot_title = 'rcp85 moist reconstruction';
plot_filename = 'omega_rec_moist_r';
climit = omega_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = omega_QG_r_rec_moist(:, :, plot_ind) - omega_QG_h_rec_moist(:, :, plot_ind);
plot_title = '$\Delta\omega_{QG}$ from moist reconstruction (Pa/s)';
plot_filename = 'omega_rec_moist_change';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = sigma_m_r(:, :, plot_ind) - sigma_m_h(:, :, plot_ind);
plot_title = 'change of $\sigma + \frac{RT}{p\theta}\frac{d\theta}{dp}|_{\theta^*}$';
plot_filename = 'sigma_m_change';
climit = [-1.0e-6, 1.0e-6];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

% zonal mean plots
d = 1.5;
NaN_matrix_m = get_NaN_matrix(cat(3, d_omega_QG (:, :, plot_ind), ...
                                  d_Adv_m    (:, :, plot_ind), ...
                                  d_k2_l2_m  (:, :, plot_ind), ...
                                  d_m2_m     (:, :, plot_ind), ...
                                  d_Adv_m    (:, :, plot_ind), ...
                                  d_rec_m    (:, :, plot_ind), ...
                                  d_J_res    (:, :, plot_ind), ...
                                  d_epsilon  (:, :, plot_ind), ...
                                  d_denom    (:, :, plot_ind)), d, tags_analyze);

omega_QG_h_mean_m = nanmean(omega_QG_h_rec(:, :, plot_ind) .* NaN_matrix_m, 2);
omega_h_mean_m = nanmean(omega_h(:, :, plot_ind) .* NaN_matrix_m, 2);


figure('pos',figsize_zonal);
hold on;
grid on;
p0 = plot(plot_x, nanmean(d_omega(:, :, plot_ind)     .* NaN_matrix_m, 2) ./ omega_h_mean_m);
p1 = plot(plot_x, nanmean(d_omega_QG(:, :, plot_ind)  .* NaN_matrix_m, 2) ./ omega_QG_h_mean_m);
p2 = plot(plot_x, nanmean(d_Adv_m(:, :, plot_ind)     .* NaN_matrix_m, 2) ./ omega_QG_h_mean_m);
p3 = plot(plot_x, nanmean(d_denom(:, :, plot_ind)     .* NaN_matrix_m, 2) ./ omega_QG_h_mean_m);
p4 = plot(plot_x, nanmean(d_epsilon(:, :, plot_ind)   .* NaN_matrix_m, 2) ./ omega_QG_h_mean_m);
p5 = plot(plot_x, nanmean(d_excessive(:, :, plot_ind) .* NaN_matrix_m, 2) ./ omega_QG_h_mean_m);
p6 = plot(plot_x, nanmean(d_rec_all(:, :, plot_ind)   .* NaN_matrix_m, 2) ./ omega_QG_h_mean_m);

set(p0, 'LineWidth', 1, 'Color', 'red');
set(p1, 'LineWidth', 1, 'Color', 'black');
set(p2, 'LineWidth', 1, 'Color', 'magenta');
set(p3, 'LineWidth', 1, 'Color', 'cyan');
set(p4, 'LineWidth', 1, 'Color', 'green');
set(p5, 'LineWidth', 1, 'Color', colors(3, :));
set(p6, 'LineWidth', 1, 'Color', 'black', 'Linestyle', '--');

lgd = legend( ...
[p0, p1, p2, p3, p4, p5, p6], ...
{'$\Delta\omega/\omega_{QG}$' , ...
'$\Delta\omega_{QG}/\omega_{QG}$' , ...
'$\Delta\omega_{Adv}/\omega_{QG}$' , ...
'$\Delta\omega_{denom}/\omega_{QG}$' , ...
'$\Delta\omega_{\epsilon}/\omega_{QG}$' , ...
'$\Delta\omega_{exce}/\omega_{QG}$' , ...
'$\Delta\omega_{rec}/\omega_{QG}$' , ...
}, ...
'location', 'best', ...
'interpreter', 'latex');
lgd.Position = [0.460, 0.3523, 0.1147, 0.31];
title([plot_level_name, ' Moist Decomposition']);
axis([-70 70 -0.601 0.601]);
xlabel('Latitude');
ylabel('\Delta\omega_{rec}/\omega_{QG}');
hold off;
saveas(gca, [plot_path, 'zonal_omega_decomposition_moist_', plot_level_name], 'png');


% contour plots

temp = d_denom(:, :, plot_ind);
plot_title = 'contribution of the denominator, $k^2(\sigma + \frac{RT}{p\theta}\frac{d\theta}{dp}|_{\theta^*} + f_0^2m^2)$ (Pa/s)';
plot_filename = 'moist_d_denom';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = d_epsilon(:, :, plot_ind);
plot_title = 'contribution of $\epsilon$ (Pa/s)';
plot_filename = 'moist_d_epsilon';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = d_Adv_m(:, :, plot_ind);
plot_title = 'contribution of RHS (Pa/s)';
plot_filename = 'moist_d_Adv';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = d_excessive(:, :, plot_ind);
plot_title = 'contribution of higher-order heating (Pa/s)';
plot_filename = 'moist_d_higherorder';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = d_rec_all(:, :, plot_ind);
plot_title = 'sum of all terms, moist decomposition (Pa/s)';
plot_filename = 'moist_d_all';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

omega_h_mean_d = nanmean(omega_h(:, :, plot_ind) .* NaN_matrix_d, 2);
plot_path_ensemble = plot_path;
plot_zonal_decomposition_paper

% get variable names except the two large event array, and save into the .mat file
var_names = who;
ind = ones(size(var_names));
for i = 1 : length(var_names)
    if ~isempty(strfind(var_names{i}, 'events_sta_historical')) || ...
       ~isempty(strfind(var_names{i}, 'events_sta_rcp85'))
        ind(i) = 0;
    end
end

save(matfilename, var_names{logical(ind)});

%% diagnostics for VORTICITY
%{
temp = A1_h + A3_h;
plot_title = 'historical term A1 + A3';
plot_filename = 'vorticity_A1A3_h';
climit = [-1e-17, 1e-17];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = A2_h + B_h;
plot_title = 'historical term A2 + B';
plot_filename = 'vorticity_A2B_h';
climit = [-1e-17, 1e-17];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = A1_r + A3_r;
plot_title = 'rcp85 term A1 + A3';
plot_filename = 'vorticity_A1A3_r';
climit = [-1e-17, 1e-17];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = A2_r + B_r;
plot_title = 'rcp85 term A2 + B';
plot_filename = 'vorticity_A2B_r';
climit = [-1e-17, 1e-17];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);
%}




