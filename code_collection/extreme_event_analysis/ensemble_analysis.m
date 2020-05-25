
addpath('/disk7/ziweili/CESM_LENS/source')
addpath('/disk7/ziweili/CESM_LENS/exp/ensemble_plots/')

clear all
v_str = '0.1';
v_str_file = '0.0';
J_OVER_OMEGA = false;
OCEAN_ONLY = false;
OCEAN_DESERT = false;
SIGMA_0 = false;
pwd_str = pwd;
ZONAL_COMPOSITE = false; % use zonal composites for the moist decomposition

if ~isempty(strfind(pwd_str, 'GFDL')) || ...
   ~isempty(strfind(pwd_str, 'daily'))
    EPS = false;
else
    EPS = true;
end

num_threshold = 1;

if ~isempty(strfind(pwd_str, 'CESM')) && ...
    ~isempty(strfind(pwd_str, 'varsig_Q')) && ...
    ~isempty(strfind(pwd_str, 'ocean_desert'))
    
    filenames = {['../001_2D_varsig_Q_ocean_desert/plot_diagnostic.mat'], ...
                 ['../002_2D_varsig_Q_ocean_desert/plot_diagnostic.mat'], ...
                 ['../003_2D_varsig_Q_ocean_desert/plot_diagnostic.mat'], ...
                 ['../004_2D_varsig_Q_ocean_desert/plot_diagnostic.mat'], ...
                 ['../005_2D_varsig_Q_ocean_desert/plot_diagnostic.mat'], ...
                 ['../035_2D_varsig_Q_ocean_desert/plot_diagnostic.mat'], ...
                };
    GFDL = false;
    num_threshold = 1;
    OCEAN_DESERT = true;
    OCEAN_ONLY = true;
    GT_ZERO = true;
    ZONAL_COMPOSITE = true;

elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
    ~isempty(strfind(pwd_str, 'varsig_Q_findcenter_correct_precip_daily'))

    filenames = {['../001_2D_varsig_Q_findcenter_correct_precip_daily/plot_diagnostic.mat'], ...
                 ['../002_2D_varsig_Q_findcenter_correct_precip_daily/plot_diagnostic.mat'], ...
                 ['../003_2D_varsig_Q_findcenter_correct_precip_daily/plot_diagnostic.mat'], ...
                 ['../004_2D_varsig_Q_findcenter_correct_precip_daily/plot_diagnostic.mat'], ...
                 ['../005_2D_varsig_Q_findcenter_correct_precip_daily/plot_diagnostic.mat'], ...
                 ['../035_2D_varsig_Q_findcenter_correct_precip_daily/plot_diagnostic.mat'], ...
                 };
    GFDL = false;
    num_threshold = 5;

elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
    ~isempty(strfind(pwd_str, 'varsig_Q_findcenter')) && ...
    ~isempty(strfind(pwd_str, 'daily'))
    filenames = {['../001_2D_varsig_Q_findcenter_daily/plot_diagnostic.mat'], ...
                 ['../002_2D_varsig_Q_findcenter_daily/plot_diagnostic.mat'], ...
                 ['../003_2D_varsig_Q_findcenter_daily/plot_diagnostic.mat'], ...
                 ['../004_2D_varsig_Q_findcenter_daily/plot_diagnostic.mat'], ...
                 ['../005_2D_varsig_Q_findcenter_daily/plot_diagnostic.mat'], ...
                 ['../035_2D_varsig_Q_findcenter_daily/plot_diagnostic.mat'], ...
                 };
    GFDL = false;
    num_threshold = 5;

elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
    ~isempty(strfind(pwd_str, 'varsig_Q_findcenter_correct_precip'))

    filenames = {['../001_2D_varsig_Q_findcenter_correct_precip/plot_diagnostic.mat'], ...
                 ['../002_2D_varsig_Q_findcenter_correct_precip/plot_diagnostic.mat'], ...
                 ['../003_2D_varsig_Q_findcenter_correct_precip/plot_diagnostic.mat'], ...
                 ['../004_2D_varsig_Q_findcenter_correct_precip/plot_diagnostic.mat'], ...
                 ['../005_2D_varsig_Q_findcenter_correct_precip/plot_diagnostic.mat'], ...
                 ['../035_2D_varsig_Q_findcenter_correct_precip/plot_diagnostic.mat'], ...
                 };
    GFDL = false;
    num_threshold = 30;

elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
    ~isempty(strfind(pwd_str, 'varsig_Q_findcenter'))
    
    % single ensemble member
    %filenames = {['../001_2D_varsig_Q_findcenter/plot_diagnostic.mat']};
    %num_threshold = 0;

    % normal setup
    filenames = {['../001_2D_varsig_Q_findcenter/plot_diagnostic.mat'], ...
                 ['../002_2D_varsig_Q_findcenter/plot_diagnostic.mat'], ...
                 ['../003_2D_varsig_Q_findcenter/plot_diagnostic.mat'], ...
                 ['../004_2D_varsig_Q_findcenter/plot_diagnostic.mat'], ...
                 ['../005_2D_varsig_Q_findcenter/plot_diagnostic.mat'], ...
                 ['../035_2D_varsig_Q_findcenter/plot_diagnostic.mat'], ...
                 };
    GFDL = false;
    num_threshold = 30;

elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
    ~isempty(strfind(pwd_str, 'smoothed_sigma_accu'))
    %{
    filenames = {['../001_2D_map_sigma_smoothed_sigma_accu/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../002_2D_map_sigma_smoothed_sigma_accu/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../003_2D_map_sigma_smoothed_sigma_accu/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../004_2D_map_sigma_smoothed_sigma_accu/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../005_2D_map_sigma_smoothed_sigma_accu/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../035_2D_map_sigma_smoothed_sigma_accu/plot_diagnostic_', v_str_file, '.mat']};
    %}
    filenames = {['../001_2D_map_sigma_smoothed_sigma_accu/plot_diagnostic.mat'], ...
                 ['../002_2D_map_sigma_smoothed_sigma_accu/plot_diagnostic.mat'], ...
                 ['../003_2D_map_sigma_smoothed_sigma_accu/plot_diagnostic.mat'], ...
                 ['../004_2D_map_sigma_smoothed_sigma_accu/plot_diagnostic.mat'], ...
                 ['../005_2D_map_sigma_smoothed_sigma_accu/plot_diagnostic.mat'], ...
                 ['../035_2D_map_sigma_smoothed_sigma_accu/plot_diagnostic.mat']};
    GFDL = false;
    num_threshold = 30;
elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
       ~isempty(strfind(pwd_str, '_traditional'))
    filenames = {['../001_2D_map_sigma_traditional/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../002_2D_map_sigma_traditional/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = false;
elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
       ~isempty(strfind(pwd_str, 'JJA')) && ...
       ~isempty(strfind(pwd_str, 'sigma')) && ...
       ~isempty(strfind(pwd_str, 'traditional'))
    filenames = {['../001_JJA_sigma_traditional/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../002_JJA_sigma_traditional/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../003_JJA_sigma_traditional/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../004_JJA_sigma_traditional/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../005_JJA_sigma_traditional/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../035_JJA_sigma_traditional/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = false;
elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
       ~isempty(strfind(pwd_str, 'JJA')) && ...
       ~isempty(strfind(pwd_str, 'sigma'))
    filenames = {['../001_JJA_sigma/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../002_JJA_sigma/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../003_JJA_sigma/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../004_JJA_sigma/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../005_JJA_sigma/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../035_JJA_sigma/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = false;

elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
       ~isempty(strfind(pwd_str, 'JJA'))
    filenames = {['../001_JJA/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../002_JJA/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../003_JJA/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../004_JJA/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../005_JJA/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../035_JJA/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = false;
elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
       ~isempty(strfind(pwd_str, 'DJF')) && ...
       ~isempty(strfind(pwd_str, 'sigma')) && ...
       ~isempty(strfind(pwd_str, 'traditional'))
    filenames = {['../001_DJF_sigma_traditional/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../002_DJF_sigma_traditional/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../003_DJF_sigma_traditional/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../004_DJF_sigma_traditional/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../005_DJF_sigma_traditional/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../035_DJF_sigma_traditional/plot_diagnostic_', v_str_file, '.mat']};

    GFDL = false;
elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
       ~isempty(strfind(pwd_str, 'DJF')) && ...
       ~isempty(strfind(pwd_str, 'sigma'))
    filenames = {['../001_DJF_sigma/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../002_DJF_sigma/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../003_DJF_sigma/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../004_DJF_sigma/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../005_DJF_sigma/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../035_DJF_sigma/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = false;
elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
       ~isempty(strfind(pwd_str, 'JJA'))
    filenames = {['../001_DJF/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../002_DJF/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../003_DJF/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../004_DJF/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../005_DJF/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../035_DJF/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = false;
elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
       ~isempty(strfind(pwd_str, 'sigma')) && ...
       ~isempty(strfind(pwd_str, 'daily'))
    
    filenames = {['../001_2D_map_sigma_daily/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../002_2D_map_sigma_daily/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../003_2D_map_sigma_daily/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../004_2D_map_sigma_daily/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../005_2D_map_sigma_daily/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../035_2D_map_sigma_daily/plot_diagnostic_', v_str_file, '.mat']};
    
    %{
    filenames = {['../001_2D_map_sigma_daily/plot_diagnostic.mat'], ...
                 ['../002_2D_map_sigma_daily/plot_diagnostic.mat'], ...
                 ['../003_2D_map_sigma_daily/plot_diagnostic.mat'], ...
                 ['../004_2D_map_sigma_daily/plot_diagnostic.mat'], ...
                 ['../005_2D_map_sigma_daily/plot_diagnostic.mat'], ...
                 ['../035_2D_map_sigma_daily/plot_diagnostic.mat']};
    %}
    GFDL = false;
elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
       ~isempty(strfind(pwd_str, 'sigma'))
    filenames = {['../001_2D_map_sigma/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../002_2D_map_sigma/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../003_2D_map_sigma/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../004_2D_map_sigma/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../005_2D_map_sigma/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../035_2D_map_sigma/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = false;
elseif ~isempty(strfind(pwd_str, 'CESM'))

    filenames = {['../001_2D_map/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../002_2D_map/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../003_2D_map/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../004_2D_map/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../005_2D_map/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../035_2D_map/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = false;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'varsig_Q_findcenter_Adv_center'))
    filenames = {['../2D_varsig_Q_findcenter_Adv_center/plot_diagnostic.mat']};
    GFDL = true;
    num_threshold = 15;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'varsig_Q_findcenter_no_smooth2a_filter'))
    filenames = {['../2D_varsig_Q_findcenter_no_smooth2a_filter/plot_diagnostic.mat']};
    GFDL = true;
    num_threshold = 15;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'varsig_Q_findcenter_large_domain'))
    filenames = {['../2D_varsig_Q_findcenter_large_domain/plot_diagnostic.mat']};
    GFDL = true;
    num_threshold = 15;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
        ~isempty(strfind(pwd_str, 'varsig_Q_findcenter_small_domain'))
    filenames = {['../2D_varsig_Q_findcenter_small_domain/plot_diagnostic.mat']};
    GFDL = true;
    num_threshold = 15;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'varsig_Q_findcenter_wb'))
    filenames = {['../2D_varsig_Q_findcenter_wb/plot_diagnostic.mat']};
    GFDL = true;
    num_threshold = 15;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'varsig_tra_wb'))
    filenames = {['../2D_map_varsig_tra_wb_findcenter/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = true;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'varsig_Adv'))
    filenames = {['../2D_varsig_Adv/plot_diagnostic.mat']};
    GFDL = true;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'varsig_Q_findcenter')) && ...
       ~isempty(strfind(pwd_str, 'daily'))
    filenames = {['../2D_varsig_Q_findcenter_daily_99.5/plot_diagnostic.mat']};
    GFDL = true;
    num_threshold = 15;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'varsig_Q_findcenter'))
    filenames = {['../2D_varsig_Q_findcenter/plot_diagnostic.mat']};
    GFDL = true;
    num_threshold = 15;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, '99.5'))
    filenames = {['../2D_grid_99.5/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = true;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'JJA')) && ...
       ~isempty(strfind(pwd_str, 'sigma'))
    filenames = {['../JJA_2D_grid_sigma/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = true;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'DJF')) && ...
       ~isempty(strfind(pwd_str, 'sigma'))
    filenames = {['../DJF_2D_grid_sigma/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = true;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'sigma')) && ...
       ~isempty(strfind(pwd_str, 'daily'))
    filenames = {['../2D_grid_sigma_daily_99.5/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = true;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'sigma'))
    filenames = {['../2D_grid_sigma/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = true;
elseif strfind(pwd_str, 'GFDL')
    filenames = {['../2D_grid/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = true;
end

ensemble_read_data_v1;

if exist('sigma_nc_h') && SIGMA_0
    sigma_h = sigma_nc_h;
    sigma_r = sigma_nc_r;
end

if OCEAN_DESERT
    SMOOTH = false;
    plot_level = 87500;
else
    SMOOTH = true;
    plot_level = 50000;
end

% read in global temperature changes
if GFDL
    load('/disk7/ziweili/test1_GFDL/exp/global_temperature/global_avg_T.mat');
else
    load('/disk7/ziweili/CESM_LENS/exp/global_temperature/global_avg_T.mat');
end
delta_T = mean(T_avg_r) - mean(T_avg_h);

figsize_zonal = [10, 10, 600, 200];

plot_level_name = [num2str(plot_level/100.), 'hPa'];

plevels = [100000, 97500, 95000, 92500, 90000, 87500, 85000, 82500, 80000, 77500, 75000, 70000, 65000, ...
            60000, 55000, 50000, 45000, 40000, 35000, 30000, 25000, 22500, 20000, 17500, 15000]';
%plevels = [85000, 82500, 80000, 77500, 75000, 70000, 65000, ...
%           60000, 55000, 50000, 45000, 40000, 35000, 30000, 25000, 22500, 20000, 17500, 15000]';
plot_ind = plevels == plot_level;
set(0,'DefaultFigureVisible','off');
omega_range = [-1.5, 1.5];
change_range = [-0.5, 0.5];
figsize = [10 10 450 180];

plot_path_ensemble = ['plots_', v_str, '/'];
if ~exist(plot_path_ensemble)
    mkdir(plot_path_ensemble);
end

% if OCEAN_ONLY, get events that are only on the ocean
tags_analyze = ones(size(X));
if OCEAN_ONLY
    tags_land = zeros(size(X));
    temp_X = mod(X + 180, 360) - 180;
    tags_land = double(reshape(landmask(Y(:), temp_X(:)), size(temp_X)));
    [x_, y_] = find(tags_land);
    box_min = 3;
    xx = (box_min - 1) / 2;
    %xx = 0;
    for k = 1 : length(x_)
        tags_land(mod([-xx : xx] + x_(k) - 1, length(tags_land(:, 1))) + 1, ...
                  mod([-xx : xx] + y_(k) - 1, length(tags_land(1, :))) + 1) = 1;
    end
    tags_analyze = tags_analyze & (1 - tags_land);
end

% if OCEAN_DESERT, get events that are in the subtropical ocean deserts

if OCEAN_DESERT
    tags_ocean_desert = zeros(size(X));
    % south Pacific
    tags_ocean_desert(X >= 240 & X <= 290 & Y > -25 & Y < -7) = 1;
    % south Atlantic
    tags_ocean_desert((X >= 330 & Y > -25 & Y < -7) | ...
                      (X <=  10 & Y > -25 & Y < -7)) = 1;

    temp = tags_ocean_desert;
    plot_title = 'ocean desert regions';
    plot_filename = 'ocean_desert';
    climit = [0, 1];
    plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);
    tags_analyze = tags_analyze & tags_ocean_desert;
end

num_event_h_3D = num_event_h;
num_event_r_3D = num_event_r;
tags_analyze_3D = repmat(tags_analyze, [1, 1, length(plevels)]);
num_event_h_3D = num_event_h_3D .* tags_analyze_3D;
num_event_r_3D = num_event_r_3D .* tags_analyze_3D;
num_event_h = num_event_h_3D(:, :, plevels == plot_level);
num_event_r = num_event_r_3D(:, :, plevels == plot_level);
num_event_d = reshape(min([num_event_h(:)'; num_event_r(:)']), size(num_event_h));
if exist('num_event_discarded_h')
    num_event_discarded_h = num_event_discarded_h(:, :, plevels == plot_level);
    num_event_discarded_r = num_event_discarded_r(:, :, plevels == plot_level);
end

omega_QG_h_rec =  - 1 ./ (k2_h .* sigma_h + m2_h .* F0.^2) .* (Adv_h + C_h);
omega_QG_r_rec =  - 1 ./ (k2_r .* sigma_r + m2_r .* F0.^2) .* (Adv_r + C_r);

temp = num_event_h;
plot_title = 'Multi-ensemble historical number of events';
plot_filename = 'num_event_h';
if ~isempty(strfind(pwd_str, '99.5'))
    climit = [0, 100];
elseif GFDL
    climit = [0, 30];
else
    climit = [0, 100];
end
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = num_event_r;
plot_title = 'Multi-ensemble rcp85 number of events';
plot_filename = 'num_event_r';
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

level_ind = plevels == plot_level;
%tags_ocean_h = double(~isnan(omega_h(:, :, level_ind))) .* tags_ocean;
%tags_ocean_h(tags_ocean_h(:) == 0) = NaN;
%tags_ocean_r = double(~isnan(omega_r(:, :, level_ind))) .* tags_ocean;
%tags_ocean_r(tags_ocean_r(:) == 0) = NaN;

%% comparison between omega and omega_QG in the two climates
% the default case
omega_h_mean_2    = nanmean(omega_h(:, :, plot_ind).*tags_analyze, 2); 
omega_r_mean_2    = nanmean(omega_r(:, :, plot_ind).*tags_analyze, 2);
omega_QG_h_mean_2 = nanmean(omega_QG_h(:, :, plot_ind).*tags_analyze, 2);
omega_QG_r_mean_2 = nanmean(omega_QG_r(:, :, plot_ind).*tags_analyze, 2);
% Paul's interest: Latent heating contribution
omega_QG_h_heating =  - 1 ./ (k2_h .* sigma_h + m2_h .* F0.^2) .* (C_h);
omega_QG_r_heating =  - 1 ./ (k2_r .* sigma_r + m2_r .* F0.^2) .* (C_r);
omega_QG_h_heating_mean = nanmean(omega_QG_h_heating(:, :, plot_ind).*tags_analyze, 2);
omega_QG_r_heating_mean = nanmean(omega_QG_r_heating(:, :, plot_ind).*tags_analyze, 2);
% plot zonal mean omega vs. omega_QG
figname = [plot_path_ensemble, 'zonal_omega_comparison_', num2str(plevels(level_ind)/100), 'hPa'];
zonal_omega
% the full-boundary case (even if you compared this case, you don't have enough information for a decomposition)
% You will still have to get a run with full boundaries
omega_h_mean_2    = nanmean(omega_h(:, :, plot_ind).*tags_analyze, 2);
omega_r_mean_2    = nanmean(omega_r(:, :, plot_ind).*tags_analyze, 2);
if exist('omega_QG_b_h')
    omega_QG_h_mean_2 = nanmean(omega_QG_b_h(:, :, plot_ind).*tags_analyze, 2);
    omega_QG_r_mean_2 = nanmean(omega_QG_b_r(:, :, plot_ind).*tags_analyze, 2);
end
figname = [plot_path_ensemble, 'zonal_omega_comparison_full_b_', num2str(plevels(level_ind)/100), 'hPa'];
zonal_omega

% 500hPa omega
temp = omega_h(:, :, plot_ind);
%plot_title = 'Multi-ensemble historical $\omega$';
plot_title = '(a) CESM historical $\omega$';
plot_filename = 'omega_h';
climit = omega_range;
if OCEAN_DESERT
    climit = [-0.5, 0.5];
end
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_h, num_threshold, true);

temp = omega_r(:, :, plot_ind);
plot_title = 'Multi-ensemble rcp85 omega';
plot_filename = 'omega_r';
climit = omega_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_r, num_threshold, false);

% vertically averaged omega
temp = nanmean(omega_r - omega_h, 3);
%temp = nanmean(omega_r - omega_h, 3) ./ nanmean(omega_r, 3);
plot_title = 'Change of vertically averaged omega (Pa/s)';
plot_filename = 'omega_vertical_average_change';
climit = [-0.2, 0.2];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_r, num_threshold, false);


temp = omega_QG_h(:, :, plot_ind);
%plot_title = 'Multi-ensemble historical $\omega_{QG}$';
plot_title = '(b) CESM historical $\omega_{QG}$';
plot_filename = 'omega_QG_h';
climit = omega_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_h, num_threshold, false);

temp = omega_QG_r(:, :, plot_ind);
plot_title = 'Multi-ensemble rcp85 omega QG';
plot_filename = 'omega_QG_r';
climit = omega_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_r, num_threshold, false);

% historical sigma
temp = sigma_h(:, :, plot_ind);
plot_title = 'Historical $\sigma$';
plot_filename = 'sigma_h';
climit = [2e-6, 5e-6];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, false, climit, ones(size(temp)), 0, false);

% historical sigma
temp = sigma_r(:, :, plot_ind);
plot_title = 'RCP8.5 $\sigma$';
plot_filename = 'sigma_r';
climit = [2e-6, 5e-6];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, false, climit, ones(size(temp)), 0, false);

temp = k2_h(:, :, plot_ind);
plot_title = 'Multi-ensemble historical $k^2$';
plot_filename = 'k2_h';
climit = [-1.0e-10, 1.0e-10];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_h, num_threshold, true);

temp = k2_r(:, :, plot_ind);
plot_title = 'Multi-ensemble rcp85 $k^2$';
plot_filename = 'k2_r';
climit = [-1.0e-10, 1.0e-10];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_r, num_threshold, true);

temp = (k2_h(:, :, plot_ind)./k2_r(:, :, plot_ind) - 1) / delta_T;
plot_title = 'Multi-ensemble $k_h^2/k_r^2 - 1$';
plot_filename = 'k2_change';
climit = [-0.13, 0.13];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, false, climit, num_event_d, num_threshold, true);

temp = k2_r(:, :, plot_ind) - k2_h(:, :, plot_ind);
plot_title = 'Multi-ensemble $k_r^2 - k_h^2$';
plot_filename = 'k2_r_minus_k2_h';
climit = [-3.1e-11, 3.1e-11];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, false, climit, num_event_d, num_threshold, true);

temp = l2_r(:, :, plot_ind) - l2_h(:, :, plot_ind);
plot_title = 'Multi-ensemble change in $k_J^2$';
plot_filename = 'l2_r_minus_l2_h';
climit = [-3.1e-11, 3.1e-11];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, false, climit, num_event_d, num_threshold, true);

temp = l2_h(:, :, plot_ind);
plot_title = 'Multi-ensemble historical $l^2$';
plot_filename = 'l2_h';
climit = [-1.0e-10, 1.0e-10];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_h, num_threshold, true);

temp = l2_r(:, :, plot_ind);
plot_title = 'Multi-ensemble rcp85 $l^2$';
plot_filename = 'l2_r';
climit = [-1.0e-10, 1.0e-10];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_r, num_threshold, true);

temp = omega_QG_h_rec(:, :, plot_ind);
plot_title = 'Reconstructed historical omega QG';
plot_filename = 'omega_QG_h_rec';
climit = omega_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_h, num_threshold, false);

temp = omega_QG_r_rec(:, :, plot_ind);
plot_title = 'Reconstructed rcp85 omega QG';
plot_filename = 'omega_QG_r_rec';
climit = omega_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_r, num_threshold, false);

temp = Le_h;
plot_title = 'historical $Le$, (meters)';
plot_filename = 'Le_h';
climit = [0, 6.0e5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, false, climit, num_event_h, num_threshold, false);

temp = Le_r;
plot_title = 'rcp85 $Le$, (meters)';
plot_filename = 'Le_r';
climit = [0, 6.0e5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, false, climit, num_event_r, num_threshold, false);

temp = (Le_r.^2 - Le_h.^2)./Le_h.^2;
plot_title = '$\Delta Le^2/Le^2$';
plot_filename = 'Le_change';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, false, climit, num_event_d, num_threshold, true);

% fractional changes in sigma
temp = (sigma_r(:, :, plot_ind) - sigma_h(:, :, plot_ind)) ./ sigma_h(:, :, plot_ind);
plot_title = '$-\Delta\sigma/\sigma$';
plot_filename = 'sigma_change';
climit = [-0.7, 0.7];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, false, climit, ones(size(temp)), 0, true);

% fractional changes in Q
temp = -(J_r(:, :, plot_ind) - J_h(:, :, plot_ind)) ./ J_h(:, :, plot_ind);
plot_title = '$-\Delta J/J$';
plot_filename = 'J_change';
climit = [-1.0, 1.0];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, false, climit, ones(size(temp)), 0, true);

adv_h_mean = nanmean(Adv_h(:, :, plot_ind), 2);
Adv_h_mean = repmat(adv_h_mean, 1, size(Adv_h, 2));;
temp = (Adv_r(:, :, plot_ind) - Adv_h(:, :, plot_ind)) ./ Adv_h_mean;
plot_title = 'Adv fractional change';
plot_filename = 'Adv_change';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true, true);

% difference between omega and omega_QG
temp = omega_h(:, :, plot_ind) - omega_QG_h(:, :, plot_ind);
plot_title = 'Historical $\omega - \omega_{QG}$';
plot_filename = 'omega_minus_omega_QG_h';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_h, num_threshold, true);

temp = omega_r(:, :, plot_ind) - omega_QG_r(:, :, plot_ind);
plot_title = 'RCP8.5 $\omega - \omega_{QG}$';
plot_filename = 'omega_minus_omega_QG_r';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_r, num_threshold, true);

% test difference between sigma_r and sigma_h, only for ocean deserts
if OCEAN_DESERT
    temp = (sigma_r(:, :, plot_ind) - sigma_h(:, :, plot_ind)) ./ sigma_h(:, :, plot_ind);
    plot_title = '$\Delta\sigma/\sigma$';
    plot_filename = 'sigma_fractional_change';
    climit = change_range;
    plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_r, num_threshold, true);

    temp = sigma_r(:, :, plot_ind) - sigma_h(:, :, plot_ind);
    plot_title = '$\Delta\sigma$';
    plot_filename = 'sigma_absolute_change';
    climit = [-1e-6, 1e-6];
    plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_r, num_threshold, true);

end

%% dry decomposition

d_sigma = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* (sigma_r - sigma_h) .* k2_h .* omega_QG_h_rec;
d_k2    = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* sigma_h .* (k2_r - k2_h) .* omega_QG_h;
d_l2    = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* kappa ./ plot_level .* (l2_r - l2_h) .* J_h;
d_J     = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* kappa ./ plot_level .* k2_h .* (J_r - J_h);
d_m2    = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* F0.^2 .* (m2_r - m2_h) .* omega_QG_h_rec;
d_Adv   = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* (Adv_r - Adv_h);
d_dtheta_dp_ma_omega = - 1./(sigma_h .* k2_h + F0.^2.*m2_h) .* ...
        kappa .* cpd ./ plot_level .* k2_h .* (omega_r .* dtheta_dp_ma_r - omega_h .* dtheta_dp_ma_h);
        %kappa .* cpd ./ plot_level .* k2_h .* (J_r - J_h) / cpd;
d_rec   = d_sigma + d_k2 + d_l2 + d_J + d_m2 + d_Adv;
d_omega_QG = omega_QG_r - omega_QG_h;
d_omega    = omega_r - omega_h;
d_k2_l2 = d_k2 + d_l2;

% zonal mean plots

colors = get(gca,'colororder');
plot_x = lat_series(lat_indices);

NaN_matrix_d = get_NaN_matrix(cat(3, d_omega_QG          (:, :, plot_ind), ...
                                  d_sigma             (:, :, plot_ind), ...
                                  d_J                 (:, :, plot_ind), ...
                                  d_k2_l2             (:, :, plot_ind), ...
                                  d_m2                (:, :, plot_ind), ...
                                  d_Adv               (:, :, plot_ind), ...
                                  d_rec               (:, :, plot_ind), ...
                                  d_dtheta_dp_ma_omega(:, :, plot_ind)), 10.0, tags_analyze);

% calculate zonal average omega and omega_QG
omega_QG_h_mean_d = nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_d, 2);
omega_QG_r_mean_d = nanmean(omega_QG_r(:, :, plot_ind) .* NaN_matrix_d, 2);
Omega_QG_h_mean_d = repmat(omega_QG_h_mean_d, 1, size(omega_QG_h, 2));
omega_h_mean_d    = nanmean(omega_h   (:, :, plot_ind) .* NaN_matrix_d, 2);
omega_r_mean_d    = nanmean(omega_r   (:, :, plot_ind) .* NaN_matrix_d, 2);
Omega_h_mean_d    = repmat(omega_h_mean_d, 1, size(omega_h, 2));

figure('pos', figsize_zonal);
hold on;
plot(-90:90, zeros(size(-90:90)), 'LineWidth', 0.8, 'Color', [0.3, 0.3, 0.3], 'Linestyle', '--')

p1 = plot(plot_x, nanmean(d_omega_QG(:, :, plot_ind) .* NaN_matrix_d, 2) ./ omega_QG_h_mean_d);
p2 = plot(plot_x, nanmean(d_sigma (:, :, plot_ind) .* NaN_matrix_d, 2) ./ omega_QG_h_mean_d);
p3 = plot(plot_x, nanmean(d_J     (:, :, plot_ind) .* NaN_matrix_d, 2) ./ omega_QG_h_mean_d);
p4 = plot(plot_x, nanmean(d_k2_l2 (:, :, plot_ind) .* NaN_matrix_d, 2) ./ omega_QG_h_mean_d);
p5 = plot(plot_x, nanmean(d_m2    (:, :, plot_ind) .* NaN_matrix_d, 2) ./ omega_QG_h_mean_d);
p6 = plot(plot_x, nanmean(d_Adv   (:, :, plot_ind) .* NaN_matrix_d, 2) ./ omega_QG_h_mean_d);
p7 = plot(plot_x, nanmean(d_dtheta_dp_ma_omega(:, :, plot_ind) .* NaN_matrix_d, 2) ./ omega_QG_h_mean_d);
p8 = plot(plot_x, nanmean(d_rec   (:, :, plot_ind) .* NaN_matrix_d, 2) ./ omega_QG_h_mean_d);

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
lgd.Position = [0.460, 0.3523, 0.1147, 0.31];
title([plot_level_name, ' Dry Decomposition']);
axis([-70 70 -0.601 0.601]);
xlabel('Latitude');
ylabel('\Delta\omega_{QG}/\omega_{QG}');
hold off;
saveas(gca, [plot_path_ensemble, 'zonal_omega_decomposition_', plot_level_name], 'png');

% contour plots

temp = (omega_QG_r(:, :, plot_ind) - omega_QG_h(:, :, plot_ind)) ./ Omega_QG_h_mean_d ./ delta_T * 100;
plot_title = '$\Delta\omega_{QG}/\omega_{QG}$';
plot_filename = 'omega_QG_fractional_change';
climit = [-0.15, 0.15] * 100;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true, true);

temp = omega_QG_r(:, :, plot_ind) - omega_QG_h(:, :, plot_ind);
plot_title = '$\Delta\omega_{QG}$';
plot_filename = 'omega_QG_absolute_change';
climit = [-0.2, 0.2];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, false, climit, num_event_d, num_threshold, true);

temp = (omega_r(:, :, plot_ind) - omega_h(:, :, plot_ind)) ./ Omega_h_mean_d ./ delta_T * 100;
plot_title = '$\Delta\omega/\omega$ (\%/K)';
plot_filename = 'omega_change';
climit = [-0.15, 0.15] * 100;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, 0, true, true);

temp = d_sigma(:, :, plot_ind) ./ Omega_QG_h_mean_d;
plot_title = 'contribution of $\sigma$';
plot_filename = 'dry_d_sigma';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true, true);

temp = d_k2(:, :, plot_ind) ./ Omega_QG_h_mean_d;
plot_title = '(a) Contribution of $k^2$';
plot_filename = 'dry_d_k2';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true, true);

temp = d_l2(:, :, plot_ind) ./ Omega_QG_h_mean_d;
plot_title = '(b) Contribution of $k_J^2$';
plot_filename = 'dry_d_l2';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true, true);

temp = d_J(:, :, plot_ind) ./ Omega_QG_h_mean_d;
plot_title = 'contribution of $J$';
plot_filename = 'dry_d_J';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true, true);

temp = d_m2(:, :, plot_ind) ./ Omega_QG_h_mean_d;
plot_title = 'contribution of $m^2$';
plot_filename = 'dry_d_m2';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true, true);

temp = d_Adv(:, :, plot_ind) ./ Omega_QG_h_mean_d;
plot_title = 'contribution of $Adv$';
plot_filename = 'dry_d_Adv';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true, true);

temp = d_rec(:, :, plot_ind) ./ Omega_QG_h_mean_d;
plot_title = 'Dry decompoistion, sum of all terms';
plot_filename = 'dry_d_all';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true, true);

temp = - (k2_r - k2_h) ./ k2_h + (l2_r - l2_h) ./ l2_h - (sigma_r - sigma_h) ./ sigma_h + (J_r - J_h) ./ J_h;
temp = temp(:, :, plot_ind);
plot_title = '$-\frac{dk^2}{k^2} + \frac{dk_J^2}{k_J^2} - \frac{d\sigma}{\sigma} + \frac{dJ}{J}$';
plot_filename = 'dry_corrected_T18';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, false, climit, num_event_d, num_threshold, true, true);

% precipitation
temp = precip_h;
plot_title = 'Historical precipitation (mm/day)';
plot_filename = 'precip_h';
climit = [0, 300];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, false);

temp = precip_r;
plot_title = 'RCP8.5 precipitation (mm/day)';
plot_filename = 'precip_r';
climit = [0, 300];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, false);

% change of precipitation
temp = (precip_r - precip_h) ./ (precip_h);
plot_title = 'Change of extreme precipitation';
plot_filename = 'precip_change';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

%% moist decomposition
if ~exist('sigma_star_h')
    sigma_star_h = - kappa ./ Plevels .* dtheta_dp_ma_h * cpd;
    sigma_star_r = - kappa ./ Plevels .* dtheta_dp_ma_r * cpd;
end
k2_star_h = k2_h;
k2_star_r = k2_r;

[latent_h, latent_r, excessive_h, excessive_r, epsilon_h, epsilon_r, ...
 sigma_m_h, sigma_m_r, k2_m_h, k2_m_r, denom_h, denom_r, ...
 d_sigma_m, d_Adv_m, d_k2_l2_m, d_m2_m, d_rec_m, d_rec_m_full, d_J_res, ...
 d_excessive_1, d_excessive_2, d_epsilon, J_res_h, J_res_r] = moist_decomposition(...
        k2_h, k2_r, l2_h, l2_r, m2_h, m2_r, sigma_h, sigma_r, Plevels, F0, ...
        dtheta_dp_ma_h, dtheta_dp_ma_r, omega_h, omega_r, omega_QG_h_rec, omega_QG_r_rec, ...
        J_h, J_r, Adv_h, Adv_r, sigma_star_h, sigma_star_r, k2_h, k2_r);


NaN_matrix_m = get_NaN_matrix(cat(3, d_omega_QG (:, :, plot_ind), ...
                                  d_Adv_m    (:, :, plot_ind), ...
                                  d_k2_l2_m  (:, :, plot_ind), ...
                                  d_m2_m     (:, :, plot_ind), ...
                                  d_Adv_m    (:, :, plot_ind), ...
                                  d_rec_m    (:, :, plot_ind), ...
                                  d_J_res  (:, :, plot_ind), ...
                                  d_epsilon  (:, :, plot_ind)), 10.0, tags_analyze);


% Paul's question: is J_res generally positive?
temp = J_res_h(:, :, plot_ind) .* NaN_matrix_m;
plot_title = 'Historical $\overline{J}_{res}$ (Wkg$^{-1}$)';
plot_filename = 'J_res_h';
climit = [-0.2, 0.2];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true, false);

temp = J_res_r(:, :, plot_ind) .* NaN_matrix_m;
plot_title = 'RCP8.5 $\overline{J}_{res}$ (Wkg$^{-1}$)';
plot_filename = 'J_res_r';
climit = [-0.2, 0.2];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true, false);

% Reviewer #3's question: how to you compare excessive and epsilon? 
temp = epsilon_h(:, :, plot_ind) .* NaN_matrix_m;
plot_title = 'Historical $\epsilon$ (Wkg$^{-1}$)';
plot_filename = 'epsilon_h';
climit = [-0.2, 0.2];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true, false);

temp = epsilon_r(:, :, plot_ind) .* NaN_matrix_m;
plot_title = 'RCP8.5 $\epsilon$ (Wkg$^{-1}$)';
plot_filename = 'epsilon_r';
climit = [-0.2, 0.2];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true, false);

temp = (epsilon_r(:, :, plot_ind) - epsilon_h(:, :, plot_ind)) .* NaN_matrix_m;
plot_title = '$\epsilon_r - \epsilon_h$ (Wkg$^{-1}$)';
plot_filename = 'epsilon_change';
climit = [-0.05, 0.05];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true, false);

temp = excessive_h(:, :, plot_ind) .* NaN_matrix_m;
plot_title = 'Historical $\overline{J_{res}} - \bar{\epsilon}$ (Wkg$^{-1}$)';
plot_filename = 'excessive_h';
climit = [-0.2, 0.2];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true, false);

temp = excessive_r(:, :, plot_ind) .* NaN_matrix_m;
plot_title = 'Historical $\overline{J_{res}} - \bar{\epsilon}$ (Wkg$^{-1}$)';
plot_filename = 'excessive_r';
climit = [-0.2, 0.2];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true, false);

% Paul's question: How do you compare J_res and Adv in historical? 
temp = kappa ./ plot_level .* l2_h(:, :, plot_ind) .* J_res_h(:, :, plot_ind) .* NaN_matrix_m;
plot_title = 'Historical $\frac{\kappa}{p}k_J^2\overline{J}_{res}$ (mkg$^{-1}$s$^{-1}$)';
plot_filename = 'term_J_res_h';
climit = [-1.0e-16, 1.0e-16];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true, false);

temp = Adv_h(:, :, plot_ind);
plot_title = 'Historical $\overline{Adv}$ (mkg$^{-1}$s$^{-1}$)';
plot_filename = 'term_Adv_h';
climit = [-1.0e-16, 1.0e-16];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true, false);

%%%%%%%%%%%%%%%%%%%% composite event as a function of latitute %%%%%%%%%%%%%%%%%%%%

if ZONAL_COMPOSITE

    omega_QG_h_comp = nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2);
    omega_QG_r_comp = nanmean(omega_QG_r(:, :, plot_ind) .* NaN_matrix_m, 2);
    omega_h_comp = nanmean(omega_h(:, :, plot_ind) .* NaN_matrix_m, 2);
    omega_r_comp = nanmean(omega_r(:, :, plot_ind) .* NaN_matrix_m, 2);
    % k2 and sigma
    k2_sigma_omega_QG_h_comp = nanmean(k2_h(:, :, plot_ind) .* sigma_h(:, :, plot_ind) .* omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2);
    k2_sigma_omega_QG_r_comp = nanmean(k2_r(:, :, plot_ind) .* sigma_r(:, :, plot_ind) .* omega_QG_r(:, :, plot_ind) .* NaN_matrix_m, 2);
    sigma_omega_QG_h_comp = nanmean(sigma_h(:, :, plot_ind) .* omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2);
    sigma_omega_QG_r_comp = nanmean(sigma_r(:, :, plot_ind) .* omega_QG_r(:, :, plot_ind) .* NaN_matrix_m, 2);
    k2_h_comp = k2_sigma_omega_QG_h_comp ./ sigma_omega_QG_h_comp;
    k2_r_comp = k2_sigma_omega_QG_r_comp ./ sigma_omega_QG_r_comp;
    sigma_h_comp = sigma_omega_QG_h_comp ./ omega_QG_h_comp;
    sigma_r_comp = sigma_omega_QG_r_comp ./ omega_QG_r_comp;
    % l2
    l2_J_h_comp = nanmean(l2_h(:, :, plot_ind) .* J_h(:, :, plot_ind) .* NaN_matrix_m, 2);
    l2_J_r_comp = nanmean(l2_r(:, :, plot_ind) .* J_r(:, :, plot_ind) .* NaN_matrix_m, 2);
    J_h_comp = nanmean(J_h(:, :, plot_ind) .* NaN_matrix_m, 2);
    J_r_comp = nanmean(J_r(:, :, plot_ind) .* NaN_matrix_m, 2);
    l2_h_comp = l2_J_h_comp ./ J_h_comp;
    l2_r_comp = l2_J_r_comp ./ J_r_comp;
    % m2
    m2_omega_QG_h_comp = nanmean(m2_h(:, :, plot_ind) .* omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2);
    m2_omega_QG_r_comp = nanmean(m2_r(:, :, plot_ind) .* omega_QG_r(:, :, plot_ind) .* NaN_matrix_m, 2);
    m2_h_comp = m2_omega_QG_h_comp ./ omega_QG_h_comp;
    m2_r_comp = m2_omega_QG_r_comp ./ omega_QG_r_comp;
    % sigma_star and residual heating
    sigma_star_omega_QG_h_comp = nanmean(omega_QG_h(:, :, plot_ind) .* sigma_star_h(:, :, plot_ind) .* NaN_matrix_m, 2);
    sigma_star_omega_QG_r_comp = nanmean(omega_QG_r(:, :, plot_ind) .* sigma_star_r(:, :, plot_ind) .* NaN_matrix_m, 2);
    sigma_star_h_comp = sigma_star_omega_QG_h_comp ./ omega_QG_h_comp;
    sigma_star_r_comp = sigma_star_omega_QG_r_comp ./ omega_QG_r_comp;
    J_res_h_comp = J_h_comp + k2_h_comp ./ l2_h_comp .* plot_level / kappa .* sigma_star_h_comp .* omega_QG_h_comp;
    J_res_r_comp = J_r_comp + k2_r_comp ./ l2_r_comp .* plot_level / kappa .* sigma_star_r_comp .* omega_QG_r_comp;

    sigma_m_h_comp = sigma_h_comp - sigma_star_h_comp;
    sigma_m_r_comp = sigma_r_comp - sigma_star_r_comp;
    denom_m_h_comp = k2_h_comp .* sigma_m_h_comp + f0(lat_indices).^2 .* m2_h_comp;
    denom_m_r_comp = k2_r_comp .* sigma_m_r_comp + f0(lat_indices).^2 .* m2_r_comp;
    Adv_h_comp = nanmean(Adv_h(:, :, plot_ind) .* NaN_matrix_m, 2);
    Adv_r_comp = nanmean(Adv_r(:, :, plot_ind) .* NaN_matrix_m, 2);
       
    %% linear expansion of the composites
    d_sigma_m_comp  = (sigma_m_r_comp - sigma_m_h_comp) .* k2_h_comp .* ...
                        (Adv_h_comp + kappa / plot_level * l2_h_comp .* J_res_h_comp) ./ denom_m_h_comp.^2;
    d_k2_m_comp     = (k2_r_comp - k2_h_comp) .* sigma_m_h_comp .* ...
                        (Adv_h_comp + kappa / plot_level * l2_h_comp .* J_res_h_comp) ./ denom_m_h_comp.^2;
    d_l2_m_comp     = - kappa / plot_level * (l2_r_comp - l2_h_comp) .* J_res_h_comp ./ denom_m_h_comp;
    d_k2_l2_m_comp  = d_k2_m_comp + d_l2_m_comp;
    d_m2_m_comp     = (m2_r_comp - m2_h_comp) .* f0(lat_indices).^2 .* ...
                        (Adv_h_comp + kappa / plot_level * l2_h_comp .* J_res_h_comp) ./ denom_m_h_comp.^2;
    d_denom_comp    = (denom_m_r_comp - denom_m_h_comp) .* ...
                        (Adv_h_comp + kappa / plot_level * l2_h_comp .* J_res_h_comp) ./ denom_m_h_comp.^2;
    d_Adv_m_comp    = - (Adv_r_comp - Adv_h_comp) ./ denom_m_h_comp;
    d_J_res_comp    = - kappa / plot_level * l2_h_comp .* (J_res_r_comp - J_res_h_comp) ./ denom_m_h_comp;
    d_omega_comp    = omega_r_comp - omega_h_comp;
    d_omega_QG_comp = omega_QG_r_comp - omega_QG_h_comp;
    d_rec_m_full_comp = d_sigma_m_comp + d_k2_l2_m_comp + d_m2_m_comp + d_Adv_m_comp + d_J_res_comp;

end

% calculate the temperature advection and tendency
T_adv_QG_h = J_h + sigma_h .* omega_QG_h_rec .* Plevels ./ kappa .* k2_h ./ l2_h;
T_adv_QG_r = J_r + sigma_r .* omega_QG_r_rec .* Plevels ./ kappa .* k2_r ./ l2_r;
if exist('T_adv_h')
    temp = T_adv_h(:, :, plot_ind) - repmat(nanmean(T_adv_h(:, :, plot_ind), 2), 1, size(T_adv_h, 2));
    plot_title = 'Temperature advection in the thermodynamic eq.';
    plot_filename = 'T_adv_h';
    climit = [-0.05, 0.05];
    plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

    temp = T_adv_QG_h(:, :, plot_ind) - repmat(nanmean(T_adv_QG_h(:, :, plot_ind), 2), 1, size(T_adv_h, 2));
    plot_title = '$J - \sigma\omega_{QG}\frac{p}{\kappa}$';
    plot_filename = 'T_adv_QG_h';
    climit = [-0.05, 0.05];
    plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);
end
temp = denom_h(:, :, plot_ind);
plot_title = 'Denominator, $-(k_m^2\sigma_m + m^2f_0^2)$';
plot_filename = 'denom_h';
climit = [-0.5e-16, 0.5e-16];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = denom_r(:, :, plot_ind);
plot_title = 'Denominator, $-(k_m^2\sigma_m + m^2f_0^2)$';
plot_filename = 'denom_r';
climit = [-0.5e-16, 0.5e-16];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

if J_OVER_OMEGA
    excessive_h(:) = 0;
    excessive_r(:) = 0;
    [epsilon_h(:), epsilon_r(:)] = deal(0);
end

% moist reconstruction
omega_QG_h_rec_moist = - (Adv_h + kappa ./ Plevels .* l2_h .* J_res_h) ./ denom_h;
omega_QG_r_rec_moist = - (Adv_r + kappa ./ Plevels .* l2_r .* J_res_r) ./ denom_r;

d_excessive     = d_excessive_1 + d_excessive_2;
d_denom         = d_sigma_m + d_k2_l2_m + d_m2_m;

if J_OVER_OMEGA
    d_sigma_m   = (sigma_m_r - sigma_m_h) .* k2_m_h .* Adv_h ./ denom_h.^2;
    d_k2_l2_m   = sigma_m_h .* (k2_m_r - k2_m_h) .* Adv_h ./ denom_h.^2;
    d_excessive_1 = d_excessive;
    d_excessive_1(:) = 0;
    d_excessive_2(:) = 0;
end

d_rec_m_full = d_sigma_m + d_k2_l2_m + d_m2_m + d_Adv_m + d_J_res;
d_rec_m      = d_sigma_m + d_k2_l2_m + d_m2_m + d_Adv_m;

omega_QG_h_mean_m = nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2);
omega_QG_r_mean_m = nanmean(omega_QG_r(:, :, plot_ind) .* NaN_matrix_m, 2);
Omega_QG_h_mean_m = repmat(omega_QG_h_mean_m, 1, size(omega_QG_h, 2));
omega_h_mean_m    = nanmean(omega_h   (:, :, plot_ind) .* NaN_matrix_m, 2);
omega_r_mean_m    = nanmean(omega_r   (:, :, plot_ind) .* NaN_matrix_m, 2);
Omega_h_mean_m    = repmat(omega_h_mean_m, 1, size(omega_h, 2));

% plot zonally-averaged precipitation
precip_h_mean_d = nanmean(precip_h .* NaN_matrix_d, 2);
precip_r_mean_d = nanmean(precip_r .* NaN_matrix_d, 2);

figure('pos', figsize_zonal);
hold on;
plot(-90:90, zeros(size(-90:90)), 'LineWidth', 0.8, 'Color', [0.3, 0.3, 0.3], 'Linestyle', '--')
p1 = plot(plot_x, nanmean((precip_r - precip_h) .* NaN_matrix_d, 2) ./ precip_h_mean_d / delta_T * 100);
set(p1, 'LineWidth', 1.5, 'Color', colors(1, :));
lgd.Position = [0.460, 0.3523, 0.1147, 0.31];
title('Sensitivity of precipitation extremes');
axis([-70 70 -2.5 25]);
xlabel('Latitude');
ylabel('Sensitivity (%/K)');
hold off;
saveas(gca, [plot_path_ensemble, 'precip_change_zonal'], 'png');


%% plot graphs for the paper

FULL = true;
SMOOTH_2 = true; % tag for smoothing in plot_decomposition_fancy
if OCEAN_DESERT
    SMOOTH_2 = false;
end
plot_decomposition_fancy

% plot zonally averaged decompositions
if ZONAL_COMPOSITE
    plot_zonal_decomposition_paper_zonal_composite
else 
    plot_zonal_decomposition_paper
end

% plot the zonally-averaged magnitude of terms
plot_terms
plot_terms_2D

% plot the 2D map of the moist decomposition
plot_decomposition_2D_moist

% plot magnified response of k and k_J
plot_dry_d_k2_l2

temp = NaN_matrix_m;
plot_title = 'Moist NaN matrix';
plot_filename = 'NaN_matrix_m';
climit = [0, 1];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);


figure('pos',figsize_zonal);
hold on;
grid on;
p0 = plot(plot_x, nanmean(d_omega       (:, :, plot_ind) .* NaN_matrix_m, 2) ./ nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2));
p1 = plot(plot_x, nanmean(d_omega_QG    (:, :, plot_ind) .* NaN_matrix_m, 2) ./ nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2));
p2 = plot(plot_x, nanmean(d_Adv_m       (:, :, plot_ind) .* NaN_matrix_m, 2) ./ nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2));
p3 = plot(plot_x, nanmean(d_denom       (:, :, plot_ind) .* NaN_matrix_m, 2) ./ nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2));
p4 = plot(plot_x, nanmean(d_epsilon     (:, :, plot_ind) .* NaN_matrix_m, 2) ./ nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2));
p5 = plot(plot_x, nanmean(d_J_res     (:, :, plot_ind) .* NaN_matrix_m, 2) ./ nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2));
p6 = plot(plot_x, nanmean(d_rec_m_full  (:, :, plot_ind) .* NaN_matrix_m, 2) ./ nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2));

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
'$\Delta\omega_{J''}/\omega_{QG}$' , ...
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
saveas(gca, [plot_path_ensemble, 'zonal_omega_decomposition_moist_', plot_level_name], 'png');

% contour plots

temp = d_k2_l2_m(:, :, plot_ind) ./ Omega_QG_h_mean_m / delta_T * 100;
plot_title = 'contribution of $k^2$ and $l^2$';
plot_filename = 'moist_d_k2_l2';
climit = [-0.15, 0.15] * 100;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true, true);

temp = d_sigma_m(:, :, plot_ind) ./ Omega_QG_h_mean_m / delta_T * 100;
plot_title = 'contribution of $\sigma_m$';
plot_filename = 'moist_d_sigma_m';
climit = [-0.15, 0.15] * 100;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true, true);

temp = d_m2_m(:, :, plot_ind) ./ Omega_QG_h_mean_m / delta_T * 100;
plot_title = 'contribution of $m^2$';
plot_filename = 'moist_d_m2';
climit = [-0.15, 0.15] * 100;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true, true);

temp = d_denom(:, :, plot_ind) ./ Omega_QG_h_mean_m / delta_T * 100;
plot_title = 'contribution of the denominator, $k^2(\sigma + \frac{RT}{p\theta}\frac{d\theta}{dp}|_{\theta^*} + f_0^2m^2)$';
plot_filename = 'moist_d_denom';
climit = [-0.15, 0.15] * 100;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true, true);

temp = d_epsilon(:, :, plot_ind) ./ Omega_QG_h_mean_m / delta_T * 100;
plot_title = 'contribution of $\epsilon$';
plot_filename = 'moist_d_epsilon';
climit = [-0.15, 0.15] * 100;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true, true);

temp = d_Adv_m(:, :, plot_ind) ./ Omega_QG_h_mean_m / delta_T * 100;
plot_title = 'contribution of $Adv$';
plot_filename = 'moist_d_Adv';
climit = [-0.15, 0.15] * 100;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true, true);

temp = d_J_res(:, :, plot_ind) ./ Omega_QG_h_mean_m / delta_T * 100;
plot_title = 'contribution of $J_{res}$ heating';
plot_filename = 'moist_d_higherorder';
climit = [-0.15, 0.15] * 100;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true, true);

temp = d_rec_m_full(:, :, plot_ind) ./ Omega_QG_h_mean_m / delta_T * 100;
plot_title = 'sum of all terms, moist decomposition';
plot_filename = 'moist_d_all';
climit = [-0.15, 0.15] * 100;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true, true);

temp = (d_k2_l2_m(:, :, plot_ind) + ...
       d_sigma_m(:, :, plot_ind) + ...
       d_m2_m(:, :, plot_ind) + ...
       d_Adv_m(:, :, plot_ind)) ./ Omega_QG_h_mean_m;
plot_title = 'sum of QG terms, moist decomposition';
plot_filename = 'moist_d_QG_theory';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = omega_QG_h_rec_moist(:, :, plot_ind);
plot_title = 'historical moist reconstruction';
plot_filename = 'omega_rec_moist_h';
climit = omega_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_h, num_threshold, false);

temp = omega_QG_r_rec_moist(:, :, plot_ind);
plot_title = 'rcp85 moist reconstruction';
plot_filename = 'omega_rec_moist_r';
climit = omega_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_r, num_threshold, false);

sigma_m_h_mean_d = nanmean(sigma_m_h   (:, :, plot_ind) .* NaN_matrix_d, 2);
Sigma_m_h_mean_d = repmat(sigma_m_h_mean_d, 1, size(sigma_m_h, 2));
temp = (sigma_m_r(:, :, plot_ind) - sigma_m_h(:, :, plot_ind)) ./ Sigma_m_h_mean_d / delta_T * 100;
plot_title = '$\Delta\sigma_m/\sigma_m$ (\%/K)';
plot_filename = 'sigma_m_change';
climit = [-0.15, 0.15] * 100
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, 0, true);


%% term_analysis
term_analysis

%% J term estimation via moist-adiabatic process
fig = figure('pos', [10, 10, 300, 270]);
%ind_lat_scatter = abs(lat) >= 30;
%ind_lat_scatter = (lat >= 30 & lat <= 70) | (lat >= -50 & lat <= -30);
ind_lat_scatter = lat >= -70 & lat <= -60;
ind_lat_scatter = abs(lat) <= 30;
if OCEAN_DESERT
    ind_lat_scatter = (1:length(lat))';
end
temp_1 = J_r(ind_lat_scatter, :, plot_ind);
%temp_2 = omega_r(ind_lat_scatter, :, plot_ind) .* dtheta_dp_ma_r(ind_lat_scatter, :, plot_ind) * cpd;
temp_2 = - omega_r(ind_lat_scatter, :, plot_ind) .* sigma_star_r(ind_lat_scatter, :, plot_ind) * plot_level / kappa;
temp_3 = - omega_r(ind_lat_scatter, :, plot_ind) .* sigma_r(ind_lat_scatter, :, plot_ind) * plot_level / kappa;
temp_4 = omega_QG_r(ind_lat_scatter, :, plot_ind) .* dtheta_dp_ma_r(ind_lat_scatter, :, plot_ind) * cpd;
s1 = scatter(temp_2(:), temp_1(:), 1, 'o', 'filled', 'MarkerFaceColor', colors(2, :));
%s1 = scatter(temp_3(:), temp_1(:), 1, 'o', 'filled', 'MarkerFaceColor', colors(2, :));
hold on
temp_1 = J_h(ind_lat_scatter, :, plot_ind);
%temp_2 = omega_h(ind_lat_scatter, :, plot_ind) .* dtheta_dp_ma_h(ind_lat_scatter, :, plot_ind) * cpd;
temp_2 = - omega_h(ind_lat_scatter, :, plot_ind) .* sigma_star_h(ind_lat_scatter, :, plot_ind) * plot_level / kappa;
temp_3 = - omega_h(ind_lat_scatter, :, plot_ind) .* sigma_h(ind_lat_scatter, :, plot_ind) * plot_level / kappa;
temp_4 = omega_QG_h(ind_lat_scatter, :, plot_ind) .* dtheta_dp_ma_h(ind_lat_scatter, :, plot_ind) * cpd;
s2 = scatter(temp_2(:), temp_1(:), 1, 'o', 'filled', 'MarkerFaceColor', colors(1, :));
%s2 = scatter(temp_3(:), temp_1(:), 1, 'o', 'filled', 'MarkerFaceColor', colors(1, :));
plot([-10, 10], [-10, 10], 'k--');
legend([s2, s1], {'Historical', 'RCP8.5'}, 'location', 'southeast', 'interpreter', 'latex')
legend boxoff
xlim([-0.05, 1.3])
ylim([-0.05, 1.3])
xlabel('$-\frac{p}{\kappa}\overline{\omega}\hat{\sigma}^*$ (J kg$^{-1}$s$^{-1}$)', 'interpreter', 'latex')
ylabel('$\overline{J}$ (J kg$^{-1}$s$^{-1}$)', 'interpreter', 'latex')
%ylabel('$-\frac{p}{\kappa}\omega_{QG}\sigma^*$', 'interpreter', 'latex')
%ylabel('$-\frac{p}{\kappa}\omega\sigma$', 'interpreter', 'latex')
set(gca, 'Position', [0.17, 0.15, 0.75, 0.75])
set(gca, 'TickLabelInterpreter','latex')
set(gca, 'TickDir', 'out');
%figname = [plot_path_ensemble, 'J_estimation_scatter_', plot_level_name, '_tropics'];
figname = [plot_path_ensemble, 'J_estimation_scatter_', plot_level_name];
if EPS
    print(fig, figname, '-depsc', '-r0', '-painters')
else
    saveas(gca, figname, 'png')
end
clf;

% plot for AGU pre
plot_for_pre_12122018

% Paul's idea of decomposition using Adv and m2
%Adv_decomposition

% omega_QG vs omega plot
omega_vs_omega_QG

%% Some diagnostics useful for the paper
ind_lat = abs(lat) > 30;
if OCEAN_DESERT
    ind_lat = (1:length(lat))';
end
Lat = repmat(lat, 1, length(lon));
Lon = repmat(lon', length(lat), 1);
dlambda = (lon(2) - lon(1)) / 180 * pi;
dlat = lat(2) - lat(1);
dS = R * (sin(min(Lat + dlat/2, 90)/180*pi) - sin(max(Lat - dlat/2, -90)/180*pi)) * dlambda;
dS_d = dS .* NaN_matrix_d;
dS_m = dS .* NaN_matrix_m;
% fraction of discarded events
temp_Ind = ((Lat > 45 | Lon < 65 | Lon > 105) & Lat > 30) | ...
             Lat < -30;
N_ratio = nansum(nansum(num_event_discarded_h(temp_Ind(:)))) / ...
          nansum(nansum(num_event_h(temp_Ind(:))));
N_fraction = N_ratio / (1 + N_ratio);
disp(['N_fraction = ', num2str(N_fraction)])
% The other discarded events in num_event_h that are not included in num_event_discarded_h
% is because of the pressure masking

% averaged changes in the extratropics, in %/K
spherical_mean = @(NaN_matrix, dS, var, omega_mean, ind_lat, delta_T, plot_ind) ...
                    nanmean(nanmean(dS(ind_lat, :) .* var(ind_lat, :, plot_ind) .* NaN_matrix(ind_lat, :), 2) ...
                            ./ omega_mean(ind_lat)) / ...
                    nanmean(nanmean(dS(ind_lat, :) .* NaN_matrix(ind_lat, :) .* ...
                    var(ind_lat, :, plot_ind) ./ var(ind_lat, :, plot_ind))) / delta_T;
% Paul's idea of comparing the rates of means, and means of rates
spherical_mean_v2 = @(NaN_matrix, dS, var, omega, ind_lat, delta_T, plot_ind) ...
                    nanmean(nanmean(dS(ind_lat, :) .* var(ind_lat, :, plot_ind) .* NaN_matrix(ind_lat, :) ...
                            ./ omega(ind_lat, :, plot_ind), 2)) / ...
                    nanmean(nanmean(dS(ind_lat, :) .* NaN_matrix(ind_lat, :) .* ...
                    var(ind_lat, :, plot_ind) ./ var(ind_lat, :, plot_ind))) / delta_T;
                    % result, less than 2% difference from the spherical_mean function
d_omega_mean_dry    = spherical_mean(NaN_matrix_d, dS_d, d_omega, omega_h_mean_d, ind_lat, delta_T, plot_ind);
d_omega_mean_dry_v2 = spherical_mean_v2(NaN_matrix_d, dS_d, d_omega, omega_h, ind_lat, delta_T, plot_ind);
d_omega_QG_mean_dry = spherical_mean(NaN_matrix_d, dS_d, d_omega_QG, omega_QG_h_mean_d, ind_lat, delta_T, plot_ind);
d_omega_QG_mean_dry_v2 = spherical_mean_v2(NaN_matrix_d, dS_d, d_omega_QG, omega_QG_h, ind_lat, delta_T, plot_ind);
% comparing omega and omega_QG
omega_h_mean_temp = spherical_mean(NaN_matrix_d, dS_d, omega_h, ones(size(omega_h_mean_d)), ind_lat, 1, plot_ind);
omega_QG_h_mean_temp = spherical_mean(NaN_matrix_d, dS_d, omega_QG_h, ones(size(omega_h_mean_d)), ind_lat, 1, plot_ind);
omega_r_mean_temp = spherical_mean(NaN_matrix_d, dS_d, omega_r, ones(size(omega_r_mean_d)), ind_lat, 1, plot_ind);
omega_QG_r_mean_temp = spherical_mean(NaN_matrix_d, dS_d, omega_QG_r, ones(size(omega_r_mean_d)), ind_lat, 1, plot_ind);
% percentage of omega_QG's underestimation of omega
disp(['In historical, omega_QG underestimates omega by ', num2str(1 - omega_QG_h_mean_temp / omega_h_mean_temp)]) 
disp(['In RCP8.5, omega_QG underestimates omega by ', num2str(1 - omega_QG_r_mean_temp / omega_r_mean_temp)])

% averaged contributions, dry decomposition
d_precip_mean = spherical_mean(NaN_matrix_d, dS_d, precip_r - precip_h, precip_h_mean_d, ind_lat, delta_T, 1);
d_sigma_mean  = spherical_mean(NaN_matrix_d, dS_d, d_sigma, omega_QG_h_mean_d, ind_lat, delta_T, plot_ind);
d_J_mean      = spherical_mean(NaN_matrix_d, dS_d, d_J,     omega_QG_h_mean_d, ind_lat, delta_T, plot_ind);
d_k2_l2_mean  = spherical_mean(NaN_matrix_d, dS_d, d_k2_l2, omega_QG_h_mean_d, ind_lat, delta_T, plot_ind);
d_m2_mean     = spherical_mean(NaN_matrix_d, dS_d, d_m2,    omega_QG_h_mean_d, ind_lat, delta_T, plot_ind);
d_Adv_mean    = spherical_mean(NaN_matrix_d, dS_d, d_Adv,   omega_QG_h_mean_d, ind_lat, delta_T, plot_ind);

d_omega_mean_moist    = spherical_mean(NaN_matrix_m, dS_m, d_omega, omega_h_mean_m, ind_lat, delta_T, plot_ind);
d_omega_QG_mean_moist = spherical_mean(NaN_matrix_m, dS_m, d_omega_QG, omega_QG_h_mean_m, ind_lat, delta_T, plot_ind);

% averaged horizontal contribution, to reply reviewer #3
ind_lat_subtropics = abs(lat) < 40;
d_k2_l2_mean_subtropics_abs = spherical_mean...
    (NaN_matrix_d, dS_d, abs(d_k2_l2), abs(omega_QG_h_mean_d), ind_lat_subtropics, delta_T, plot_ind);
d_k2_mean_subtropics_abs = spherical_mean...
    (NaN_matrix_d, dS_d, abs(d_k2), abs(omega_QG_h_mean_d), ind_lat_subtropics, delta_T, plot_ind);

% averaged contributions, moist decomposition
d_sigma_m_mean = spherical_mean(NaN_matrix_m, dS_m, d_sigma_m, omega_QG_h_mean_m, ind_lat, delta_T, plot_ind);
d_J_res_mean   = spherical_mean(NaN_matrix_m, dS_m, d_J_res,   omega_QG_h_mean_m, ind_lat, delta_T, plot_ind);
d_k2_l2_m_mean = spherical_mean(NaN_matrix_m, dS_m, d_k2_l2_m, omega_QG_h_mean_m, ind_lat, delta_T, plot_ind);
d_m2_m_mean    = spherical_mean(NaN_matrix_m, dS_m, d_m2_m,    omega_QG_h_mean_m, ind_lat, delta_T, plot_ind);
d_Adv_m_mean   = spherical_mean(NaN_matrix_m, dS_m, d_Adv_m,   omega_QG_h_mean_m, ind_lat, delta_T, plot_ind);

% generation of latex code for a table
% save values for the table
matfilename = 'extra_tropical_means';
save([pwd_str, '/', matfilename], 'd_precip_mean', ...
    'd_omega_mean_dry'  , 'd_omega_QG_mean_dry'  , ...
    'd_omega_mean_moist', 'd_omega_QG_mean_moist', ...
    'd_sigma_mean'  , 'd_J_mean'    , 'd_k2_l2_mean'  , 'd_m2_mean'  , 'd_Adv_mean', ...
    'd_sigma_m_mean', 'd_J_res_mean', 'd_k2_l2_m_mean', 'd_m2_m_mean', 'd_Adv_m_mean');


%%%%%%%%%%%%%%%%%%%%%%% calculation and figures to reply to the reviewers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% term in front of dL^2/L^2
alpha = (m2_h .* F0.^2 .* omega_QG_h) .* C_h + Adv_h .* k2_h .* sigma_h .* omega_QG_h ./ ...
        (Adv_h + C_h) ./ (- m2_h .* F0.^2 .* omega_QG_h - k2_h .* sigma_h .* omega_QG_h);
alpha_mean = nanmean(nanmean(alpha(:, :, plevels == plot_level))); % global average
alpha_max = nanmax(nanmax(alpha(:, :, plevels == plot_level))); % global average
disp(['The coefficient in front of delta L^2/L^2 is ', num2str(alpha_mean), ' for global average, and '...
        num2str(alpha_max), ' for maximum.'])

% J vs sigma*omega, lat-lon plot
temp = J_h(:, :, plot_ind);
plot_title = 'historical diabatic heating (J kg$^{-1}$s$^{-1}$)';
plot_filename = 'J_h';
climit = [-0.03, 1.6];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_h, num_threshold, false);

temp = - omega_h(:, :, plot_ind) .* sigma_star_h(:, :, plot_ind) * plot_level / kappa;
plot_title = 'historical moist-adiabatic reconstruction of J (J kg$^{-1}$s$^{-1}$)';
plot_filename = 'J_h_moist_diabatic';
climit = [-0.03, 1.6];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_h, num_threshold, false);

temp = J_r(:, :, plot_ind);
plot_title = 'RCP8.5 diabatic heating (J kg$^{-1}$s$^{-1}$)';
plot_filename = 'J_r';
climit = [-0.03, 1.6];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_r, num_threshold, false);

temp = - omega_r(:, :, plot_ind) .* sigma_star_r(:, :, plot_ind) * plot_level / kappa;
plot_title = 'RCP8.5 moist-adiabatic reconstruction of J (J kg$^{-1}$s$^{-1}$)';
plot_filename = 'J_r_moist_diabatic';
climit = [-0.03, 1.6];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_r, num_threshold, false);

%%%%%%%%%% get averaged vertical profile of omega and omega_QG %%%%%%%%%%%%%
disp('averaged vertical profile of omega')
disp(squeeze(nanmean(nanmean(omega_h .* tags_analyze_3D, 1), 2)))
disp('averaged vertical profile of omega_QG')
disp(squeeze(nanmean(nanmean(omega_r .* tags_analyze_3D, 1), 2)))

% get the RMS error of inversion, reviewer #2
inversion_error = sqrt(spherical_mean(NaN_matrix_d, dS_d, (omega_h - omega_QG_h).^2, ...
                    ones(size(omega_h_mean_d)), ind_lat, 1, plot_ind))

% check lower-level vertical velocities
[js, is] = ind2sub(size(omega_h), find(omega_h(:, :, plot_ind) > 0));
[j_max, i_max] = ind2sub(size(omega_h), find(omega_h(:, :, plot_ind) == max(max(omega_h(:, :, plot_ind)))));
% the vertical structure with largest (weakest) omega at plot_level:
squeeze(omega_h(j_max, i_max, :))


% the fraction of masked points over the subtropics (poleward of 7-40)
temp_lat_ind = (lat < 40 & lat > 3) | (lat > -40 & lat < -3);
mask_fraction = sum(sum(num_event_d(temp_lat_ind, :) < num_threshold)) / sum(temp_lat_ind) / length(lon);
disp(['a fraction of ', num2str(mask_fraction), ' has been masked'])

% plot precipitation quantiles for comparison with precip in CESM-LE
if ~GFDL
    if ~isempty(strfind(pwd_str, 'daily'))
        quantile_file_h = '/net/aimsir/archive1/ziweili/CESM_LENS/output/precip_99.9th_quantile_h_daily_corrected.mat';
    else
        quantile_file_h = '/net/aimsir/archive1/ziweili/CESM_LENS/output/precip_99.9th_quantile_h_corrected.mat';
    end
    temp = load(quantile_file_h, 'Q');
    temp = temp.Q(:, 23:170)';
    plot_title = 'Historical 99.9th percentile precipitation (mm/day)';
    plot_filename = 'precip_h_qualtile';
    climit = [0, 200];
    plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_r, num_threshold, false);
end

