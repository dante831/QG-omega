
addpath('/disk7/ziweili/CESM_LENS/source')
addpath('/disk7/ziweili/CESM_LENS/exp/ensemble_plots/')

v_str = '0.1';
v_str_file = '0.0';
J_OVER_OMEGA = false;
OCEAN_ONLY = false;

pwd_str = pwd;
if ~isempty(strfind(pwd_str, 'CESM')) && ...
   ~isempty(strfind(pwd_str, 'smoothed_sigma_accu'))
    filenames = {['../001_2D_map_sigma_smoothed_sigma_accu/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../002_2D_map_sigma_smoothed_sigma_accu/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../003_2D_map_sigma_smoothed_sigma_accu/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../004_2D_map_sigma_smoothed_sigma_accu/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../005_2D_map_sigma_smoothed_sigma_accu/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../035_2D_map_sigma_smoothed_sigma_accu/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = false;
    num_threshold = 20;
elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
       ~isempty(strfind(pwd_str, '_traditional'))
    filenames = {['../001_2D_map_sigma_traditional/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../002_2D_map_sigma_traditional/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = false;
    num_threshold = 1;
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
    num_threshold = 1;

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
    num_threshold = 1;

elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
       ~isempty(strfind(pwd_str, 'JJA'))
    filenames = {['../001_JJA/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../002_JJA/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../003_JJA/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../004_JJA/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../005_JJA/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../035_JJA/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = false;
    num_threshold = 1;
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
    num_threshold = 1;
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
    num_threshold = 1;
elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
       ~isempty(strfind(pwd_str, 'JJA'))
    filenames = {['../001_DJF/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../002_DJF/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../003_DJF/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../004_DJF/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../005_DJF/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../035_DJF/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = false;
    num_threshold = 1;
elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
       ~isempty(strfind(pwd_str, 'sigma')) && ...
       ~isempty(strfind(pwd_str, 'daily'))
    filenames = {['../001_2D_map_sigma_daily/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../002_2D_map_sigma_daily/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../003_2D_map_sigma_daily/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../004_2D_map_sigma_daily/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../005_2D_map_sigma_daily/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../035_2D_map_sigma_daily/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = false;
    num_threshold = 1;
elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
       ~isempty(strfind(pwd_str, 'sigma'))
    %filenames = {['../001_2D_map_sigma/plot_diagnostic_', v_str_file, '_ugvg.mat'], ...
    %             ['../002_2D_map_sigma/plot_diagnostic_', v_str_file, '_ugvg.mat'], ...
    %             ['../003_2D_map_sigma/plot_diagnostic_', v_str_file, '_ugvg.mat'], ...
    %             ['../004_2D_map_sigma/plot_diagnostic_', v_str_file, '_ugvg.mat'], ...
    %             ['../005_2D_map_sigma/plot_diagnostic_', v_str_file, '_ugvg.mat']};
    filenames = {['../001_2D_map_sigma/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../002_2D_map_sigma/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../003_2D_map_sigma/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../004_2D_map_sigma/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../005_2D_map_sigma/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../035_2D_map_sigma/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = false;
    num_threshold = 1;
elseif ~isempty(strfind(pwd_str, 'CESM'))

    filenames = {['../001_2D_map/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../002_2D_map/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../003_2D_map/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../004_2D_map/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../005_2D_map/plot_diagnostic_', v_str_file, '.mat'], ...
                 ['../035_2D_map/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = false;
    num_threshold = 1;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'varsig_tra_wb'))
    filenames = {['../2D_map_varsig_tra_wb/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = true;
    num_threshold = 1;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, '99.5'))
    filenames = {['../2D_grid_99.5/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = true;
    num_threshold = 1;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'JJA')) && ...
       ~isempty(strfind(pwd_str, 'sigma'))
    filenames = {['../JJA_2D_grid_sigma/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = true;
    num_threshold = 1;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'DJF')) && ...
       ~isempty(strfind(pwd_str, 'sigma'))
    filenames = {['../DJF_2D_grid_sigma/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = true;
    num_threshold = 1;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'sigma')) && ...
       ~isempty(strfind(pwd_str, 'daily'))
    filenames = {['../2D_grid_sigma_daily/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = true;
    num_threshold = 1;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'sigma'))
    filenames = {['../2D_grid_sigma/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = true;
    num_threshold = 1;
elseif strfind(pwd_str, 'GFDL')
    filenames = {['../2D_grid/plot_diagnostic_', v_str_file, '.mat']};
    GFDL = true;
    num_threshold = 1;
end

ensemble_read_data_v1;
figsize_zonal = [10, 10, 600, 200];

if GFDL
    SMOOTH = true;
else
    SMOOTH = true;
end
plot_level = 50000;
plot_level_name = [num2str(plot_level/100.), 'hPa'];

plevels = [100000, 97500, 95000, 92500, 90000, 87500, 85000, 82500, 80000, 77500, 75000, 70000, 65000, ...
            60000, 55000, 50000, 45000, 40000, 35000, 30000, 25000, 22500, 20000, 17500, 15000]';
%    plevels = [85000, 82500, 80000, 77500, 75000, 70000, 65000, ...
%               60000, 55000, 50000, 45000, 40000, 35000, 30000, 25000, 22500, 20000, 17500, 15000]';
plot_ind = plevels == plot_level;
set(0,'DefaultFigureVisible','off');
omega_range = [-1.5, 1.5];
change_range = [-0.5, 0.5];

plot_path_ensemble = ['plots_', v_str, '/'];
if ~exist(plot_path_ensemble)
    mkdir(plot_path_ensemble);
end

plot_terms;
num_event_h = num_event_h(:, :, plevels == plot_level);
num_event_r = num_event_r(:, :, plevels == plot_level);
num_event_d = reshape(min([num_event_h(:)'; num_event_r(:)']), size(num_event_h));

omega_QG_h_rec =  - 1 ./ (k2_h .* sigma_h + m2_h .* F0.^2) .* (Adv_h + C_h);
omega_QG_r_rec =  - 1 ./ (k2_r .* sigma_r + m2_r .* F0.^2) .* (Adv_r + C_r);

% remove events whose omega_QG > 0
num_event_h(omega_QG_h_rec(:, :, plevels == plot_level) >= 0) = 0;
num_event_r(omega_QG_r_rec(:, :, plevels == plot_level) >= 0) = 0;

figsize = [10 10 450 180];
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

% if OCEAN_ONLY, get events that are only on the ocean
tags_land = zeros(size(X));
if OCEAN_ONLY
    temp_X = mod(X + 180, 360) - 180;
    tags_land = double(reshape(landmask(Y(:), temp_X(:)), size(temp_X)));
    [x_, y_] = find(tags_land);
    box_min = 21;
    xx = (box_min - 1) / 2;
    %xx = 0;
    for k = 1 : length(x_)
        tags_land(mod([-xx : xx] + x_(k) - 1, length(tags_land(:, 1))) + 1, ...
                  mod([-xx : xx] + y_(k) - 1, length(tags_land(1, :))) + 1) = 1;
    end
end
tags_analyze = 1 - tags_land;

level_ind = plevels == 50000;
%tags_ocean_h = double(~isnan(omega_h(:, :, level_ind))) .* tags_ocean;
%tags_ocean_h(tags_ocean_h(:) == 0) = NaN;
%tags_ocean_r = double(~isnan(omega_r(:, :, level_ind))) .* tags_ocean;
%tags_ocean_r(tags_ocean_r(:) == 0) = NaN;

% comparison between omega and omega_QG in the two climates
figure('pos',figsize_zonal);
colors = get(gca,'colororder');
hold on;
grid on;

omega_h_mean_2    = nanmean(omega_h(:, :, plot_ind).*tags_analyze, 2); 
omega_r_mean_2    = nanmean(omega_r(:, :, plot_ind).*tags_analyze, 2);
omega_QG_h_mean_2 = nanmean(omega_QG_h(:, :, plot_ind).*tags_analyze, 2);
omega_QG_r_mean_2 = nanmean(omega_QG_r(:, :, plot_ind).*tags_analyze, 2);
% Paul's interest: Latent heating contribution
omega_QG_h_heating =  - 1 ./ (k2_h .* sigma_h + m2_h .* F0.^2) .* (C_h);
omega_QG_r_heating =  - 1 ./ (k2_r .* sigma_r + m2_r .* F0.^2) .* (C_r);
omega_QG_h_heating_mean = nanmean(omega_QG_h_heating(:, :, plot_ind).*tags_analyze, 2);
omega_QG_r_heating_mean = nanmean(omega_QG_r_heating(:, :, plot_ind).*tags_analyze, 2);

p1 = plot(plot_x, - omega_h_mean_2);
p2 = plot(plot_x, - omega_QG_h_mean_2);
p3 = plot(plot_x, - omega_r_mean_2);
p4 = plot(plot_x, - omega_QG_r_mean_2);
%p5 = plot(plot_x, - omega_QG_h_heating_mean);
%p6 = plot(plot_x, - omega_QG_r_heating_mean);

set([p1, p2, p3, p4, p5, p6], 'LineWidth', 1)
set(p1, 'LineWidth', 1, 'Color', colors(1, :), 'Linestyle', '-');
set(p2, 'LineWidth', 1, 'Color', colors(1, :), 'Linestyle', '-.');
set(p3, 'LineWidth', 1, 'Color', colors(2, :), 'Linestyle', '-');
set(p4, 'LineWidth', 1, 'Color', colors(2, :), 'Linestyle', '-.');
%set(p5, 'LineWidth', 1, 'Color', colors(1, :), 'Linestyle', ':');
%set(p6, 'LineWidth', 1, 'Color', colors(2, :), 'Linestyle', ':');

%lgd = legend( ...
%[p1, p2, p3, p4, p5, p6], {...
%'Historical $\omega$', 'Historical $\omega_{QG}$', ...
%'RCP8.5 $\omega$'    , 'RCP8.5 $\omega_{QG}$', ...
%'Historical heating' , 'RCP8.5 heating'}, ...
%'location', 'best', ...
%'interpreter', 'latex');
lgd = legend( ...
[p1, p2, p3, p4], {...
'Historical $\omega$', 'Historical $\omega_{QG}$', ...
'RCP8.5 $\omega$'    , 'RCP8.5 $\omega_{QG}$'}, ...
'location', 'best', ...
'interpreter', 'latex');

title('Comparison of zonal averaged \omega');
axis([-70 70 0.0 1.20]);
xlabel('Latitude');
ylabel('-\omega (Pa/s)')
hold off;

saveas(gca, [plot_path_ensemble, 'zonal_omega_comparison_', num2str(plevels(level_ind)/100), 'hPa'], 'png');
%saveas(gca, [plot_path_ensemble, 'zonal_omega_comparison'], 'png');


temp = omega_h(:, :, plot_ind);
%plot_title = 'Multi-ensemble historical $\omega$';
plot_title = '(a) CESM historical $\omega$';
plot_filename = 'omega_h';
climit = omega_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_h, num_threshold, false);

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



temp = k2_h(:, :, plot_ind);
plot_title = 'Multi-ensemble historical $k^2$';
plot_filename = 'k2_h';
climit = [0.0e-11, 1.0e-10];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_h, num_threshold, false);

temp = k2_r(:, :, plot_ind);
plot_title = 'Multi-ensemble rcp85 $k^2$';
plot_filename = 'k2_r';
climit = [0.0e-11, 1.0e-10];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_r, num_threshold, false);

temp = k2_h(:, :, plot_ind)./k2_r(:, :, plot_ind) - 1;
plot_title = 'Multi-ensemble $k_h^2/k_r^2 - 1$';
plot_filename = 'k2_change';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, false, climit, num_event_d, num_threshold, true);

temp = k2_r(:, :, plot_ind) - k2_h(:, :, plot_ind);
plot_title = 'Multi-ensemble $k_r^2 - k_h^2$';
plot_filename = 'k2_r_minus_k2_h';
climit = [-1.0e-11, 1.0e-11];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_r, num_threshold, false);

temp = l2_r(:, :, plot_ind) - l2_h(:, :, plot_ind);
plot_title = 'Multi-ensemble $l_r^2 - l_h^2$';
plot_filename = 'l2_r_minus_l2_h';
climit = [-1.0e-11, 1.0e-11];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_r, num_threshold, false);

temp = l2_h(:, :, plot_ind);
plot_title = 'Multi-ensemble historical $l^2$';
plot_filename = 'l2_h';
climit = [0.0e-11, 1.0e-10];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_h, num_threshold, false);

temp = l2_r(:, :, plot_ind);
plot_title = 'Multi-ensemble rcp85 $l^2$';
plot_filename = 'l2_r';
climit = [0.0e-11, 1.0e-10];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_r, num_threshold, false);

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
    smooth_window_x, smooth_window_y, false, climit, ones(size(temp)), 0, false);

temp = (Le_r.^2 - Le_h.^2)./Le_h.^2;
plot_title = '$\Delta Le^2/Le^2$';
plot_filename = 'Le_change';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, false, climit, ones(size(temp)), 0, true);

adv_h_mean = nanmean(Adv_h(:, :, plot_ind), 2);
Adv_h_mean = repmat(adv_h_mean, 1, size(Adv_h, 2));;
temp = (Adv_r(:, :, plot_ind) - Adv_h(:, :, plot_ind)) ./ Adv_h_mean;
plot_title = 'Adv fractional change';
plot_filename = 'Adv_change';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true);

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
d_omega_QG = (omega_QG_r - omega_QG_h);
d_omega    = (omega_r - omega_h);
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
                                  d_dtheta_dp_ma_omega(:, :, plot_ind)), 1.5, tags_analyze);
%{
[temp_ind] = abs(d_omega_QG          (:, :, plot_ind)) > 1.5 | ...
             abs(d_sigma             (:, :, plot_ind)) > 1.5 | ...
             abs(d_J                 (:, :, plot_ind)) > 1.5 | ...
             abs(d_k2_l2             (:, :, plot_ind)) > 1.5 | ...
             abs(d_m2                (:, :, plot_ind)) > 1.5 | ...
             abs(d_Adv               (:, :, plot_ind)) > 1.5 | ...
             abs(d_rec               (:, :, plot_ind)) > 1.5 | ...
             abs(d_dtheta_dp_ma_omega(:, :, plot_ind)) > 1.5;
NaN_matrix_d = ones(size(d_omega_QG  (:, :, plot_ind)));
NaN_matrix_d(temp_ind) = NaN;
NaN_matrix_d = NaN_matrix_d .* tags_analyze;
%}

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

%p1 = plot(plot_x, d_omega_QG_comp (:, plot_ind));
%p2 = plot(plot_x, d_sigma_comp  (:, plot_ind));
%p3 = plot(plot_x, d_J_comp      (:, plot_ind));
%p4 = plot(plot_x, d_k2_l2_comp  (:, plot_ind));
%p5 = plot(plot_x, d_m2_comp     (:, plot_ind));
%p6 = plot(plot_x, d_Adv_comp    (:, plot_ind));
%p7 = plot(plot_x, nanmean(d_dtheta_dp_ma_omega(:, :, plot_ind) .* NaN_matrix_d(:, :, plot_ind), 2));
%p8 = plot(plot_x, d_rec_comp         (:, plot_ind));


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

temp = (omega_QG_r(:, :, plot_ind) - omega_QG_h(:, :, plot_ind)) ./ Omega_QG_h_mean_d;
plot_title = '$- \Delta\omega_{QG}/\omega_{QG}$';
plot_filename = 'omega_QG_change';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = (omega_r(:, :, plot_ind) - omega_h(:, :, plot_ind)) ./ Omega_h_mean_d;
plot_title = '$\Delta\omega/\omega$';
plot_filename = 'omega_change';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_sigma(:, :, plot_ind) ./ Omega_QG_h_mean_d;
plot_title = 'contribution of $\sigma$';
plot_filename = 'dry_d_sigma';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_k2(:, :, plot_ind) ./ Omega_QG_h_mean_d;
plot_title = '(a) Contribution of $k^2$';
plot_filename = 'dry_d_k2';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_l2(:, :, plot_ind) ./ Omega_QG_h_mean_d;
plot_title = '(b) Contribution of $k_J^2$';
plot_filename = 'dry_d_l2';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);


temp = d_k2_l2(:, :, plot_ind) ./ Omega_QG_h_mean_d;
plot_title = 'contribution of $k^2$ and $l^2$';
plot_filename = 'dry_d_k2_l2';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_J(:, :, plot_ind) ./ Omega_QG_h_mean_d;
plot_title = 'contribution of $J$';
plot_filename = 'dry_d_J';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_m2(:, :, plot_ind) ./ Omega_QG_h_mean_d;
plot_title = 'contribution of $m^2$';
plot_filename = 'dry_d_m2';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_Adv(:, :, plot_ind) ./ Omega_QG_h_mean_d;
plot_title = 'contribution of $Adv$';
plot_filename = 'dry_d_Adv';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_rec(:, :, plot_ind) ./ Omega_QG_h_mean_d;
plot_title = 'sum of all terms';
plot_filename = 'dry_d_all';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

%% moist decomposition

if strcmp(v_str, '0.0') && J_OVER_OMEGA

    % version A, subtract J from both sides
    sigma_m_h   = sigma_h + kappa ./ Plevels .* J_h ./ omega_QG_h_rec .* l2_h ./ k2_h;
    sigma_m_r   = sigma_r + kappa ./ Plevels .* J_r ./ omega_QG_r_rec .* l2_r ./ k2_r;
    k2_m_h      = (sigma_h .* k2_h + kappa ./ Plevels .* J_h ./ omega_QG_h_rec .* l2_h) ./ sigma_m_h; 
    k2_m_r      = (sigma_r .* k2_r + kappa ./ Plevels .* J_r ./ omega_QG_r_rec .* l2_r) ./ sigma_m_r;
    % version A.1, define sigma_m as sigma + kappa/p*J/omega_QG, problem: this method is not stable in k2
    %sigma_m_h   = sigma_h + kappa ./ Plevels .* J_h ./ omega_QG_h_rec;
    %sigma_m_r   = sigma_r + kappa ./ Plevels .* J_r ./ omega_QG_r_rec;
    % version B
    %sigma_m_h   = max(sigma_h + kappa ./ Plevels .* J_h ./ omega_QG_h_rec .* l2_h ./ k2_h, 0);
    %sigma_m_r   = max(sigma_r + kappa ./ Plevels .* J_r ./ omega_QG_r_rec .* l2_r ./ k2_r, 0);
    %k2_m_h = k2_h;
    %k2_m_r = k2_r;
    latent_h = 0;
    latent_r = 0;
    denom_h     = sigma_m_h .* k2_m_h + F0.^2.*m2_h;
    denom_r     = sigma_m_r .* k2_m_r + F0.^2.*m2_r;
elseif strcmp(v_str, '1.0') || strcmp(v_str, '0.0')
    sigma_m_h   = sigma_h + kappa .* cpd ./ Plevels .* dtheta_dp_ma_h .* l2_h ./ k2_h;
    sigma_m_r   = sigma_r + kappa .* cpd ./ Plevels .* dtheta_dp_ma_r .* l2_r ./ k2_r;
    latent_h    = cpd .* dtheta_dp_ma_h .* omega_QG_h_rec;
    latent_r    = cpd .* dtheta_dp_ma_r .* omega_QG_r_rec;
    %sigma_m_r   = sigma_r + kappa .* cpd ./ Plevels .* dtheta_dp_ma_r;
    %sigma_m_h   = sigma_h + kappa .* cpd ./ Plevels .* dtheta_dp_ma_h;
    %latent_h    = cpd .* dtheta_dp_ma_h .* omega_QG_h_rec .* k2_h./l2_h;
    %latent_r    = cpd .* dtheta_dp_ma_r .* omega_QG_r_rec .* k2_r./l2_r;
    denom_h     = sigma_m_h .* k2_h + F0.^2.*m2_h;
    denom_r     = sigma_m_r .* k2_r + F0.^2.*m2_r;
    k2_m_h      = (sigma_h .* k2_h + kappa ./ Plevels .* J_h ./ omega_QG_h_rec .* l2_h) ./ sigma_m_h;
    k2_m_r      = (sigma_r .* k2_r + kappa ./ Plevels .* J_r ./ omega_QG_r_rec .* l2_r) ./ sigma_m_r;
    k2_m_h = k2_h;
    k2_m_r = k2_r;
elseif strcmp(v_str, '1.1') || strcmp(v_str, '0.1')
    sigma_m_h   = sigma_h + kappa .* cpd ./ Plevels .* dtheta_dp_ma_h;
    sigma_m_r   = sigma_r + kappa .* cpd ./ Plevels .* dtheta_dp_ma_r;
    latent_h    = cpd .* dtheta_dp_ma_h .* omega_QG_h_rec .* k2_h./l2_h;
    latent_r    = cpd .* dtheta_dp_ma_r .* omega_QG_r_rec .* k2_r./l2_r;
    denom_h     = sigma_m_h .* k2_h + F0.^2.*m2_h;
    denom_r     = sigma_m_r .* k2_r + F0.^2.*m2_r;
    k2_m_h = k2_h;
    k2_m_r = k2_r;
end

%%%%%%%%%%%%%%%%%%%% composite event as a function of latitute %%%%%%%%%%%%%%%%%%%%
%{
omega_QG_h_comp = nanmean(omega_QG_h(:, :, plot_ind), 2);
omega_QG_r_comp = nanmean(omega_QG_r(:, :, plot_ind), 2);
omega_h_comp = nanmean(omega_h(:, :, plot_ind), 2);
omega_r_comp = nanmean(omega_r(:, :, plot_ind), 2);
denom_h_comp = nanmean(denom_h(:, :, plot_ind) .* omega_QG_h(:, :, plot_ind), 2) ./ omega_QG_h_comp;
denom_r_comp = nanmean(denom_r(:, :, plot_ind) .* omega_QG_r(:, :, plot_ind), 2) ./ omega_QG_r_comp;
sigma_m_omega_k2_h_comp = nanmean(sigma_m_h(:, :, plot_ind) .* omega_QG_h(:, :, plot_ind) .* k2_m_h(:, :, plot_ind), 2);
sigma_m_omega_k2_r_comp = nanmean(sigma_m_r(:, :, plot_ind) .* omega_QG_r(:, :, plot_ind) .* k2_m_r(:, :, plot_ind), 2);
sigma_m_omega_h_comp = nanmean(sigma_m_h(:, :, plot_ind) .* omega_QG_h(:, :, plot_ind), 2);
sigma_m_omega_r_comp = nanmean(sigma_m_r(:, :, plot_ind) .* omega_QG_r(:, :, plot_ind), 2);
k2_m_h_comp = sigma_m_omega_k2_h_comp ./ sigma_m_omega_h_comp;
k2_m_r_comp = sigma_m_omega_k2_r_comp ./ sigma_m_omega_r_comp;
sigma_m_h_comp = sigma_m_omega_h_comp ./ omega_QG_h_comp;
sigma_m_r_comp = sigma_m_omega_r_comp ./ omega_QG_r_comp;
Adv_h_comp = nanmean(Adv_h(:, :, plot_ind), 2);
Adv_r_comp = nanmean(Adv_r(:, :, plot_ind), 2);
m2_omega_h_comp = nanmean(m2_h(:, :, plot_ind) .* omega_QG_h(:, :, plot_ind), 2);
m2_omega_r_comp = nanmean(m2_r(:, :, plot_ind) .* omega_QG_r(:, :, plot_ind), 2);
m2_h_comp = m2_omega_h_comp ./ omega_QG_h_comp;
m2_r_comp = m2_omega_r_comp ./ omega_QG_r_comp;

%% linear expansion of the composites
d_sigma_m_comp  = (sigma_m_r_comp - sigma_m_h_comp) .* k2_m_h_comp .* Adv_h_comp ./ denom_h_comp.^2;
d_k2_l2_m_comp  = sigma_m_h_comp .* (k2_m_r_comp - k2_m_h_comp) .* Adv_h_comp ./ denom_h_comp.^2;
d_m2_m_comp     = (m2_r_comp - m2_h_comp) .* F0(:, 1).^2 .* Adv_h_comp ./ denom_h_comp.^2;
d_denom_comp    = (denom_r_comp - denom_h_comp) .* Adv_h_comp ./ denom_h_comp.^2;
d_Adv_m_comp    = - (Adv_r_comp - Adv_h_comp) ./ denom_h_comp;
d_omega_comp    = omega_r_comp - omega_h_comp;
d_omega_QG_comp = omega_QG_r_comp - omega_QG_h_comp;
d_rec_m_full_comp = d_sigma_m_comp + d_k2_l2_m_comp + d_m2_m_comp + d_Adv_m_comp;
%}

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

excessive_h = cpd .* dtheta_dp_ma_h .* omega_h - latent_h;
excessive_r = cpd .* dtheta_dp_ma_r .* omega_r - latent_r;
epsilon_h   = J_h - excessive_h - latent_h;
epsilon_r   = J_r - excessive_r - latent_r;
J_prime_h   = J_h - latent_h;
J_prime_r   = J_r - latent_r;

if J_OVER_OMEGA
    excessive_h(:) = 0;
    excessive_r(:) = 0;
    [epsilon_h(:), epsilon_r(:)] = deal(0);
end

% moist reconstruction
omega_QG_h_rec_moist = - (Adv_h + kappa ./ Plevels .* l2_h .* J_prime_h) ./ denom_h;
omega_QG_r_rec_moist = - (Adv_r + kappa ./ Plevels .* l2_r .* J_prime_r) ./ denom_r;

% linear expansion
%{
d_sigma_m       = (sigma_m_r - sigma_m_h) .* k2_m_h .* ...
                    (Adv_h + kappa ./ plot_level .* l2_h .* (epsilon_h + excessive_h)) ./ denom_h.^2;
d_k2_l2_m       = (k2_m_r - k2_m_h) .* sigma_m_h .* ...
                    (Adv_h + kappa ./ plot_level .* l2_h .* (epsilon_h + excessive_h)) ./ denom_h.^2 ...
                    - kappa ./ plot_level .* (l2_r - l2_h) .* excessive_h ./ denom_h;
d_m2_m          = (m2_r - m2_h) .* repmat(F0, [1, 1, length(plevels)]).^2 .* ...
                    (Adv_h + kappa ./ plot_level .* l2_h .* (epsilon_h + excessive_h)) ./ denom_h.^2;
d_denom         = (denom_r - denom_h) .* (Adv_h + kappa ./ plot_level .* l2_h .* (epsilon_h + excessive_h)) ./ denom_h.^2;
d_Adv_m         = - (Adv_r - Adv_h) ./ denom_h;
d_epsilon       = - kappa ./ plot_level .* (l2_r.*epsilon_r - l2_h.*epsilon_h) ./ denom_h;
d_excessive     = - kappa ./ plot_level .* l2_h .* (excessive_r - excessive_h) ./ denom_h;
%}
if J_OVER_OMEGA
    %d_Adv_m = d_Adv_m - kappa ./ plot_level .* l2_h .* (excessive_r - excessive_h) ./ denom_h;
    d_excessive(:) = 0;
end
if strcmp(v_str, '0.1')
    % choice 1
    d_k2_l2_m       = - ((k2_m_r - k2_m_h) .* sigma_m_h .* omega_QG_h_rec + ...
                        kappa ./ Plevels .* (l2_r - l2_h) .* J_prime_h) ./ denom_h;
    d_J_prime       = - kappa ./ Plevels .* l2_h .* (J_prime_r - J_prime_h) ./ denom_h;
    d_excessive_1   = - kappa ./ Plevels .* cpd .* (dtheta_dp_ma_r - dtheta_dp_ma_h) .* (l2_h.*omega_h - k2_h.*omega_QG_h_rec) ./ denom_h;
    d_excessive_2   = - kappa ./ Plevels .* cpd .* dtheta_dp_ma_h .* l2_h.*(omega_r - omega_h + ...
                                                        k2_h./l2_h.*omega_QG_h_rec - k2_r./l2_r.*omega_QG_r_rec) ./ denom_h;

    % choice 2
    %d_k2_l2_m       = - (k2_m_r - k2_m_h) .* sigma_m_h .* omega_QG_h_rec ./ denom_h;
    %d_J_prime       = - kappa ./ Plevels .* (l2_r .* J_prime_r - l2_h .* J_prime_h) ./ denom_h;
    %d_excessive_1   = - kappa ./ Plevels .* cpd .* (dtheta_dp_ma_r - dtheta_dp_ma_h) .* (l2_h.*omega_h - k2_h.*omega_QG_h_rec) ./ denom_h;
    %d_excessive_2   = - kappa ./ Plevels .* cpd .* dtheta_dp_ma_h .* ((l2_r.*omega_r - k2_r.*omega_QG_r_rec) - ...
    %                                                                  (l2_h.*omega_h - k2_h.*omega_QG_h_rec)) ./ denom_h;
    
    
    d_excessive     = d_excessive_1 + d_excessive_2;
    d_sigma_m       = - (sigma_m_r - sigma_m_h) .* k2_m_h .* omega_QG_h_rec ./ denom_h;
    d_m2_m          = - (m2_r - m2_h) .* repmat(F0, [1, 1, length(plevels)]).^2 .* omega_QG_h_rec ./ denom_h;
    %d_denom         = (denom_r - denom_h) .* (Adv_h + kappa ./ Plevels .* l2_h .* J_prime_h) ./ denom_h.^2;
    d_denom         = d_sigma_m + d_k2_l2_m + d_m2_m;
    d_Adv_m         = - (Adv_r - Adv_h) ./ denom_h;
    d_epsilon       = - kappa ./ Plevels .* l2_h .* (epsilon_r - epsilon_h) ./ denom_h;
    
else
    d_excessive_1   = - kappa ./ Plevels .* cpd .* (dtheta_dp_ma_r - dtheta_dp_ma_h) .* l2_h.*(omega_h - omega_QG_h_rec) ./ denom_h;
    d_excessive_2   = - kappa ./ Plevels .* cpd .* l2_h .* ((omega_r - omega_QG_r_rec) - ...
                                                               (omega_h - omega_QG_h_rec)) .* dtheta_dp_ma_h ./ denom_h;
    d_excessive     = d_excessive_1 + d_excessive_2;
    d_J_prime       = d_excessive_1 + d_excessive_2 + d_epsilon;
end
d_omega         = omega_r - omega_h;

if J_OVER_OMEGA
    d_sigma_m   = (sigma_m_r - sigma_m_h) .* k2_m_h .* Adv_h ./ denom_h.^2;
    d_k2_l2_m   = sigma_m_h .* (k2_m_r - k2_m_h) .* Adv_h ./ denom_h.^2;
    d_excessive_1 = d_excessive;
    d_excessive_1(:) = 0;
    d_excessive_2(:) = 0;
end

%d_rec_m_full = (denom_r - denom_h) .* (Adv_h + kappa ./ Plevels .* l2_h .* J_prime_h) ./ denom_h.^2 ...
%    - (Adv_r - Adv_h + kappa ./ Plevels .* (l2_r .* J_prime_r - l2_h .* J_prime_h)) ./ denom_h; % this agrees with the change of omega_QG pretty well
%d_rec_m_full = (denom_r - denom_h) .* (Adv_h + kappa ./ Plevels .* l2_h .* J_prime_h) ./ denom_h.^2 ...
%    - (Adv_r - Adv_h + kappa ./ Plevels .* ((l2_r - l2_h) .* J_prime_h + l2_h .* (J_prime_r - J_prime_h))) ./ denom_h;
d_rec_m_full = d_sigma_m + d_k2_l2_m + d_m2_m + d_Adv_m + d_J_prime;

%omega_QG_r_rec_moist = - (Adv_r + kappa ./ Plevels .* l2_r .* J_prime_r) ./ denom_r;
d_rec_m      = d_sigma_m + d_k2_l2_m + d_m2_m + d_Adv_m;

NaN_matrix_m = get_NaN_matrix(cat(3, d_omega_QG (:, :, plot_ind), ...
                                  d_Adv_m    (:, :, plot_ind), ...
                                  d_k2_l2_m  (:, :, plot_ind), ...
                                  d_m2_m     (:, :, plot_ind), ...
                                  d_Adv_m    (:, :, plot_ind), ...
                                  d_rec_m    (:, :, plot_ind), ...
                                  d_J_prime  (:, :, plot_ind), ...
                                  d_epsilon  (:, :, plot_ind)), 1.5, tags_analyze);
%{
[temp_ind] = abs(d_omega_QG (:, :, plot_ind)) > 1.5 | ...
             abs(d_Adv_m    (:, :, plot_ind)) > 1.5 | ...
             abs(d_k2_l2_m  (:, :, plot_ind)) > 1.5 | ...
             abs(d_m2_m     (:, :, plot_ind)) > 1.5 | ...
             abs(d_Adv_m    (:, :, plot_ind)) > 1.5 | ...
             abs(d_rec_m    (:, :, plot_ind)) > 1.5 | ...
             abs(d_J_prime  (:, :, plot_ind)) > 1.5 | ...
             abs(d_epsilon  (:, :, plot_ind)) > 1.5 | ...
             denom_h(:, :, plot_ind) < 0;
NaN_matrix_m = ones(size(d_omega_QG  (:, :, plot_ind)));
NaN_matrix_m(temp_ind) = NaN;
NaN_matrix_m = NaN_matrix_m .* tags_analyze;
%}

omega_QG_h_mean_m = nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2);
omega_QG_r_mean_m = nanmean(omega_QG_r(:, :, plot_ind) .* NaN_matrix_m, 2);
Omega_QG_h_mean_m = repmat(omega_QG_h_mean_m, 1, size(omega_QG_h, 2));
omega_h_mean_m    = nanmean(omega_h   (:, :, plot_ind) .* NaN_matrix_m, 2);
omega_r_mean_m    = nanmean(omega_r   (:, :, plot_ind) .* NaN_matrix_m, 2);
Omega_h_mean_m    = repmat(omega_h_mean_m, 1, size(omega_h, 2));

temp_var = sqrt(nanvar([reshape(d_denom(:, :, plot_ind)     .* NaN_matrix_m, [length(NaN_matrix_m(:)), 1]), ...
                        reshape(d_Adv_m(:, :, plot_ind)     .* NaN_matrix_m, [length(NaN_matrix_m(:)), 1]), ...
                        reshape(d_rec_m(:, :, plot_ind)     .* NaN_matrix_m, [length(NaN_matrix_m(:)), 1]), ...
                        reshape(d_J_prime(:, :, plot_ind) .* NaN_matrix_m, [length(NaN_matrix_m(:)), 1]), ...
                        reshape(d_epsilon(:, :, plot_ind)   .* NaN_matrix_m, [length(NaN_matrix_m(:)), 1])]));
[temp_ind_2] = abs(d_denom(:, :, plot_ind)     .*NaN_matrix_m - repmat(nanmean(d_denom(:, :, plot_ind)    .*NaN_matrix_m, 2), 1, size(NaN_matrix_m, 2))) > 3*temp_var(1) | ...
               abs(d_Adv_m(:, :, plot_ind)     .*NaN_matrix_m - repmat(nanmean(d_Adv_m(:, :, plot_ind)    .*NaN_matrix_m, 2), 1, size(NaN_matrix_m, 2))) > 3*temp_var(2) | ...
               abs(d_rec_m(:, :, plot_ind)     .*NaN_matrix_m - repmat(nanmean(d_rec_m(:, :, plot_ind)    .*NaN_matrix_m, 2), 1, size(NaN_matrix_m, 2))) > 3*temp_var(3) | ...
               abs(d_J_prime(:, :, plot_ind) .*NaN_matrix_m - repmat(nanmean(d_J_prime(:, :, plot_ind).*NaN_matrix_m, 2), 1, size(NaN_matrix_m, 2))) > 3*temp_var(4) | ...
               abs(d_epsilon(:, :, plot_ind)   .*NaN_matrix_m - repmat(nanmean(d_epsilon(:, :, plot_ind)  .*NaN_matrix_m, 2), 1, size(NaN_matrix_m, 2))) > 3*temp_var(5);
NaN_matrix_m(temp_ind_2) = NaN;

FULL = false;
plot_decomposition_fancy
FULL = true;
plot_decomposition_fancy

plot_zonal_decomposition_paper

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
p5 = plot(plot_x, nanmean(d_J_prime     (:, :, plot_ind) .* NaN_matrix_m, 2) ./ nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2));
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

temp = d_k2_l2_m(:, :, plot_ind) ./ Omega_QG_h_mean_m;
plot_title = 'contribution of $k^2$ and $l^2$';
plot_filename = 'moist_d_k2_l2';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_sigma_m(:, :, plot_ind) ./ Omega_QG_h_mean_m;
plot_title = 'contribution of $\sigma_m$';
plot_filename = 'moist_d_sigma_m';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_m2_m(:, :, plot_ind) ./ Omega_QG_h_mean_m;
plot_title = 'contribution of $m^2$';
plot_filename = 'moist_d_m2';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_denom(:, :, plot_ind) ./ Omega_QG_h_mean_m;
plot_title = 'contribution of the denominator, $k^2(\sigma + \frac{RT}{p\theta}\frac{d\theta}{dp}|_{\theta^*} + f_0^2m^2)$';
plot_filename = 'moist_d_denom';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_epsilon(:, :, plot_ind) ./ Omega_QG_h_mean_m;
plot_title = 'contribution of $\epsilon$';
plot_filename = 'moist_d_epsilon';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_Adv_m(:, :, plot_ind) ./ Omega_QG_h_mean_m;
plot_title = 'contribution of RHS';
plot_filename = 'moist_d_Adv';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_J_prime(:, :, plot_ind) ./ Omega_QG_h_mean_m;
plot_title = 'contribution of J'' heating';
plot_filename = 'moist_d_higherorder';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_rec_m_full(:, :, plot_ind) ./ Omega_QG_h_mean_m;
plot_title = 'sum of all terms, moist decomposition';
plot_filename = 'moist_d_all';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

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

temp = (sigma_m_r(:, :, plot_ind) - sigma_m_h(:, :, plot_ind)) ./ sigma_m_h(:, :, plot_ind);
plot_title = '$\Delta\sigma_m/\sigma_m$';
plot_filename = 'sigma_m_change';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);


%% term_analysis
term_analysis

%% J term estimation via moist-adiabatic process
figure('pos', [10, 10, 300, 270])
temp_1 = J_r(:, :, plot_ind);
temp_2 = omega_r(:, :, plot_ind) .* dtheta_dp_ma_r(:, :, plot_ind) * cpd;
temp_3 = - omega_r(:, :, plot_ind) .* sigma_r(:, :, plot_ind) * plot_level / kappa;
temp_4 = omega_QG_r(:, :, plot_ind) .* dtheta_dp_ma_r(:, :, plot_ind) * cpd;
s1 = scatter(temp_1(:), temp_2(:), 1, 'o', 'filled', 'MarkerFaceColor', colors(2, :));
%s1 = scatter(temp_1(:), temp_4(:), 1, 'o', 'filled', 'MarkerFaceColor', colors(2, :));
hold on
temp_1 = J_h(:, :, plot_ind);
temp_2 = omega_h(:, :, plot_ind) .* dtheta_dp_ma_h(:, :, plot_ind) * cpd;
temp_3 = - omega_h(:, :, plot_ind) .* sigma_h(:, :, plot_ind) * plot_level / kappa;
temp_4 = omega_QG_h(:, :, plot_ind) .* dtheta_dp_ma_h(:, :, plot_ind) * cpd;
s2 = scatter(temp_1(:), temp_2(:), 1.5, 'o', 'filled', 'MarkerFaceColor', colors(1, :));
%s2 = scatter(temp_1(:), temp_4(:), 1, 'o', 'filled', 'MarkerFaceColor', colors(1, :));
plot([-10, 10], [-10, 10], 'k--');
legend([s2, s1], {'Historical', 'RCP85'}, 'location', 'best')
legend boxoff
xlim([-0.05, 1.2])
ylim([-0.05, 1.2])
xlabel('$J$ (W/kg)', 'interpreter', 'latex')
ylabel('$c_p\frac{T}{\theta}\frac{d\theta}{dp}\big|_{\theta^*}\omega$', 'interpreter', 'latex')
%ylabel('$c_p\frac{T}{\theta}\frac{d\theta}{dp}\big|_{\theta^*}\omega_{QG}$', 'interpreter', 'latex')
%ylabel('$\sigma\omega$', 'interpreter', 'latex')
ylb = get(gca,'ylabel');
ylb.Position = [0.12, 0.64 -1] ;
ylb.Rotation = 0;
%set(ylb, 'Rotation', 0)
set(gca, 'Position', [0.1000, 0.1500, 0.750, 0.750])
set(gca, 'TickLabelInterpreter','latex')
%title('Diabatic heating vs. $c_p\frac{T}{\theta}\frac{d\theta}{dp}|_\theta^*\omega$', 'interpreter', 'latex')
saveas(gca, [plot_path_ensemble, 'J_estimation_scatter_', plot_level_name], 'png')
clf;

% plot for AGU pre
plot_for_pre_12122018

% Paul's idea of decomposition using Adv and m2
%Adv_decomposition
