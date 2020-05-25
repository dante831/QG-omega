
addpath('../../source/');
addpath('/disk7/ziweili/CESM_LENS/exp/seasonal_on_land/')
v_str = '1.0';
pwd_str = pwd;
if ~isempty(strfind(pwd_str, 'CESM')) && ...
   ~isempty(strfind(pwd_str, 'JJA')) && ...
   ~isempty(strfind(pwd_str, 'sigma'))
    filenames = {['../001_JJA_sigma/plot_diagnostic_', v_str, '.mat'], ...
                 ['../002_JJA_sigma/plot_diagnostic_', v_str, '.mat'], ...
                 ['../003_JJA_sigma/plot_diagnostic_', v_str, '.mat'], ...
                 ['../004_JJA_sigma/plot_diagnostic_', v_str, '.mat'], ...
                 ['../005_JJA_sigma/plot_diagnostic_', v_str, '.mat'], ...
                 ['../035_JJA_sigma/plot_diagnostic_', v_str, '.mat']};
    GFDL = false;
    num_threshold = 1;
elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
       ~isempty(strfind(pwd_str, 'JJA'))
    filenames = {['../001_JJA/plot_diagnostic_', v_str, '.mat'], ...
                 ['../002_JJA/plot_diagnostic_', v_str, '.mat'], ...
                 ['../003_JJA/plot_diagnostic_', v_str, '.mat'], ...
                 ['../004_JJA/plot_diagnostic_', v_str, '.mat'], ...
                 ['../005_JJA/plot_diagnostic_', v_str, '.mat'], ...
                 ['../035_JJA/plot_diagnostic_', v_str, '.mat']};
    GFDL = false;
    num_threshold = 1;
elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
       ~isempty(strfind(pwd_str, 'DJF')) && ...
       ~isempty(strfind(pwd_str, 'sigma'))
    filenames = {['../001_DJF_sigma/plot_diagnostic_', v_str, '.mat'], ...
                 ['../002_DJF_sigma/plot_diagnostic_', v_str, '.mat'], ...
                 ['../003_DJF_sigma/plot_diagnostic_', v_str, '.mat'], ...
                 ['../004_DJF_sigma/plot_diagnostic_', v_str, '.mat'], ...
                 ['../005_DJF_sigma/plot_diagnostic_', v_str, '.mat'], ...
                 ['../035_DJF_sigma/plot_diagnostic_', v_str, '.mat']};
    GFDL = false;
    num_threshold = 1;
elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
       ~isempty(strfind(pwd_str, 'JJA'))
    filenames = {['../001_DJF/plot_diagnostic_', v_str, '.mat'], ...
                 ['../002_DJF/plot_diagnostic_', v_str, '.mat'], ...
                 ['../003_DJF/plot_diagnostic_', v_str, '.mat'], ...
                 ['../004_DJF/plot_diagnostic_', v_str, '.mat'], ...
                 ['../005_DJF/plot_diagnostic_', v_str, '.mat'], ...
                 ['../035_DJF/plot_diagnostic_', v_str, '.mat']};
    GFDL = false;
    num_threshold = 1;
elseif ~isempty(strfind(pwd_str, 'CESM'))

    filenames = {['../001_2D_map/plot_diagnostic_', v_str, '.mat'], ...
                 ['../002_2D_map/plot_diagnostic_', v_str, '.mat'], ...
                 ['../003_2D_map/plot_diagnostic_', v_str, '.mat'], ...
                 ['../004_2D_map/plot_diagnostic_', v_str, '.mat'], ...
                 ['../005_2D_map/plot_diagnostic_', v_str, '.mat'], ...
                 ['../035_2D_map/plot_diagnostic_', v_str, '.mat']};
    GFDL = false;
    num_threshold = 30;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, '99.5'))
    filenames = {['../2D_grid_99.5/plot_diagnostic_', v_str, '.mat']};
    GFDL = true;
    num_threshold = 50;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'JJA')) && ...
       ~isempty(strfind(pwd_str, 'sigma'))
    filenames = {['../JJA_2D_grid_sigma/plot_diagnostic_', v_str, '.mat']};
    GFDL = true;
    num_threshold = 1;
elseif strfind(pwd_str, 'GFDL')
    filenames = {['../2D_grid/plot_diagnostic_', v_str, '.mat']};
    GFDL = true;
    num_threshold = 10;
end

ensemble_read_data_v1;
figsize_zonal = [10, 10, 600, 200];

if GFDL
    SMOOTH = false;
else
    SMOOTH = true;
end
plot_level = 50000;
plevels = [85000, 82500, 80000, 77500, 75000, 70000, 65000, ...
           60000, 55000, 50000, 45000, 40000, 35000, 30000, 25000, 22500, 20000, 17500, 15000]';
plot_ind = plevels == plot_level;
set(0,'DefaultFigureVisible','off');
ABSOLUTE = true;
omega_range = [-1.5, 1.5];
change_range = [-0.5, 0.5];

plot_path_ensemble = ['plots_', v_str, '/'];
if ~exist(plot_path_ensemble)
    mkdir(plot_path_ensemble);
end

plot_error;
num_event_h = num_event_h(:, :, plevels == plot_level);
num_event_r = num_event_r(:, :, plevels == plot_level);
num_event_d = reshape(min([num_event_h(:)'; num_event_r(:)']), size(num_event_h));

omega_QG_h_rec =  - 1 ./ (k2_h .* sigma_h + m2_h .* F0.^2) .* (Adv_h + C_h);
omega_QG_r_rec =  - 1 ./ (k2_r .* sigma_r + m2_r .* F0.^2) .* (Adv_r + C_r);

% remove events whose omega_QG > 0
num_event_h(omega_QG_h_rec(:, :, plevels == plot_level) >= 0) = 0;
num_event_r(omega_QG_r_rec(:, :, plevels == plot_level) >= 0) = 0;

% calculate zonal average omega and omega_QG
omega_QG_h_mean = nanmean(omega_QG_h(:, :, plot_ind), 2);
Omega_QG_h_mean = repmat(omega_QG_h_mean, 1, size(omega_QG_h, 2));
omega_h_mean = nanmean(omega_h(:, :, plot_ind), 2);
Omega_h_mean = repmat(omega_h_mean, 1, size(omega_h, 2));


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



temp = omega_h(:, :, plot_ind);
plot_title = 'Multi-ensemble historical omega';
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
plot_title = 'Multi-ensemble historical omega QG';
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
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);


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


temp = (omega_QG_r(:, :, plot_ind) - omega_QG_h(:, :, plot_ind)) ./ Omega_QG_h_mean;
plot_title = '$- \Delta\omega_{QG}/\omega_{QG}$';
plot_filename = 'omega_QG_change';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = (omega_r(:, :, plot_ind) - omega_h(:, :, plot_ind)) ./ Omega_h_mean;
plot_title = '$\Delta\omega/\omega$';
plot_filename = 'omega_change';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

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
    smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_h, num_threshold, false);

temp = Le_r;
plot_title = 'rcp85 $Le$, (meters)';
plot_filename = 'Le_r';
climit = [0, 6.0e5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = (Le_r.^2 - Le_h.^2)./Le_h.^2;
plot_title = '$\Delta Le^2/Le^2$';
plot_filename = 'Le_change';
climit = [-0.5, 0.5];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, ...
    smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true);



%% moist version


d_sigma = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* (sigma_r - sigma_h) .* k2_h .* omega_QG_h_rec;
d_k2    = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* sigma_h .* (k2_r - k2_h) .* omega_QG_h;
d_l2    = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* kappa ./ plot_level .* (l2_r - l2_h) .* J_h;
d_J     = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* kappa ./ plot_level .* k2_h .* (J_r - J_h);
d_m2    = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* F0.^2 .* (m2_r - m2_h) .* omega_QG_h_rec;
d_Adv   = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* (Adv_r - Adv_h);
d_dtheta_dp_ma_omega = - 1./(sigma_h .* k2_h + F0.^2.*m2_h) .* ...
        kappa .* cpd ./ plot_level .* l2_h .* (omega_r .* dtheta_dp_ma_r - omega_h .* dtheta_dp_ma_h);
d_rec   = d_sigma + d_k2 + d_l2 + d_J + d_m2 + d_Adv;
d_omega_QG = (omega_QG_r - omega_QG_h);
d_omega    = (omega_r - omega_h);
d_k2_l2 = d_k2 + d_l2;

% zonal mean plots

colors = get(gca,'colororder');
plot_x = lat_series(lat_indices);

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

%{
[d_sigma_comp, d_k2_comp, d_l2_comp, d_J_comp, d_m2_comp, ...
 d_Adv_comp, d_rec_comp, d_omega_QG_comp] = deal(zeros([size(sigma_h, 1), size(sigma_h, 3)]));
[X_h_dry, X_r_dry] = deal(zeros([size(sigma_h, 1), 7, size(sigma_h, 3)]));

for j = 1 : size(sigma_h, 1)
[d_sigma_comp(j, :), d_k2_comp(j, :), d_l2_comp(j, :), d_J_comp(j, :), d_m2_comp(j, :), ...
 d_Adv_comp(j, :), d_rec_comp(j, :), d_omega_QG_comp(j, :), X_h_dry(j, :, :), X_r_dry(j, :, :)] = ...
        composite_changes_dry(squeeze(sigma_h(j, :, :) .* NaN_matrix(j, :, :)), ...
                              squeeze(k2_h(j, :, :) .* NaN_matrix(j, :, :)), ...
                              squeeze(l2_h(j, :, :) .* NaN_matrix(j, :, :)), ...
                              squeeze(m2_h(j, :, :).*F0(j,1).^2 .* NaN_matrix(j, :, :)), ...
                              squeeze(Adv_h(j, :, :) .* NaN_matrix(j, :, :)), ...
                              squeeze(J_h(j, :, :) .* NaN_matrix(j, :, :)), ...
                              squeeze(omega_QG_h_rec(j, :, :) .* NaN_matrix(j, :, :)), ...
                              squeeze(sigma_r(j, :, :) .* NaN_matrix(j, :, :)), ...
                              squeeze(k2_r(j, :, :) .* NaN_matrix(j, :, :)), ...
                              squeeze(l2_r(j, :, :) .* NaN_matrix(j, :, :)), ...
                              squeeze(m2_r(j, :, :).*F0(j,1).^2 .* NaN_matrix(j, :, :)), ...
                              squeeze(Adv_r(j, :, :) .* NaN_matrix(j, :, :)), ...
                              squeeze(J_r(j, :, :) .* NaN_matrix(j, :, :)), ...
                              squeeze(omega_QG_r_rec(j, :, :) .* NaN_matrix(j, :, :)), ...
                              plevels);
end
d_k2_l2_comp = d_k2_comp + d_l2_comp;
%}
%{
figure('pos', figsize_zonal);
hold on;
plot(-90:90, zeros(size(-90:90)), 'LineWidth', 0.8, 'Color', [0.3, 0.3, 0.3], 'Linestyle', '--')

p1 = plot(plot_x, nanmean(d_omega_QG(:, :, plot_ind) .* NaN_matrix_d, 2));
p2 = plot(plot_x, nanmean(d_sigma (:, :, plot_ind) .* NaN_matrix_d, 2));
p3 = plot(plot_x, nanmean(d_J     (:, :, plot_ind) .* NaN_matrix_d, 2));
p4 = plot(plot_x, nanmean(d_k2_l2 (:, :, plot_ind) .* NaN_matrix_d, 2));
p5 = plot(plot_x, nanmean(d_m2    (:, :, plot_ind) .* NaN_matrix_d, 2));
p6 = plot(plot_x, nanmean(d_Adv   (:, :, plot_ind) .* NaN_matrix_d, 2));
p7 = plot(plot_x, nanmean(d_dtheta_dp_ma_omega(:, :, plot_ind) .* NaN_matrix_d(:, :, plot_ind), 2));
p8 = plot(plot_x, nanmean(d_rec   (:, :, plot_ind) .* NaN_matrix_d(:, :, plot_ind), 2));

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
axis([-70 70 -0.301 0.301]);
xlabel('Latitude');
ylabel('\Delta\omega_{QG}/\omega_{QG}');
hold off;
saveas(gca, [plot_path_ensemble, 'zonal_omega_decomposition_', plot_level_name], 'png');
%}
% contour plots

temp = d_sigma(:, :, plot_ind) ./ Omega_QG_h_mean;
plot_title = 'contribution of $\sigma$';
plot_filename = 'dry_d_sigma';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_k2(:, :, plot_ind) ./ Omega_QG_h_mean;
plot_title = 'contribution of $k^2$';
plot_filename = 'dry_d_k2';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_k2_l2(:, :, plot_ind) ./ Omega_QG_h_mean;
plot_title = 'contribution of $k^2$ and $l^2$';
plot_filename = 'dry_d_k2_l2';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_J(:, :, plot_ind) ./ Omega_QG_h_mean;
plot_title = 'contribution of $J$';
plot_filename = 'dry_d_J';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_m2(:, :, plot_ind) ./ Omega_QG_h_mean;
plot_title = 'contribution of $m^2$';
plot_filename = 'dry_d_m2';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_Adv(:, :, plot_ind) ./ Omega_QG_h_mean;
plot_title = 'contribution of $Adv$';
plot_filename = 'dry_d_Adv';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_rec(:, :, plot_ind) ./ Omega_QG_h_mean;
plot_title = 'sum of all terms';
plot_filename = 'dry_d_all';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

%% plot for the paper, 1st figure
%{
figsize_multipanel = [10, 10, 600, 375];
fig = figure('pos', figsize_multipanel);
hold on
row = 3;
col = 2;
climit = change_range;
ax = gca;
temp = d_omega(:, :, plot_ind) ./ Omega_h_mean;
plot_title = '$\Delta\omega/\omega$';
subplot_2D_map(ax, row, col, 1, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_omega_QG(:, :, plot_ind) ./ Omega_QG_h_mean;
plot_title = '$\Delta\omega_{QG}/\omega_{QG}$';;
subplot_2D_map(ax, row, col, 2, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_rec(:, :, plot_ind) ./ Omega_QG_h_mean;
plot_title = 'Sum of all terms';
subplot_2D_map(ax, row, col, 3, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_J(:, :, plot_ind) ./ Omega_QG_h_mean;
plot_title = 'Contribution of J';
subplot_2D_map(ax, row, col, 4, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_k2_l2(:, :, plot_ind) ./ Omega_QG_h_mean;
plot_title = 'Contribution of $k^2$ and $l^2$';
subplot_2D_map(ax, row, col, 5, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_sigma(:, :, plot_ind) ./ Omega_QG_h_mean;
plot_title = 'Contribution of $\sigma$';
subplot_2D_map(ax, row, col, 6, figsize_multipanel, ...
    temp, X, Y, plot_title, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

saveas(gca, [plot_path_ensemble, 'dry_decomposition'], 'png')
%saveas(gca, [plot_path_ensemble, 'dry_decomposition'], 'pdf')
%}

%% moist decomposition

% exessive may be abandoned since the contribution from epsilon and excessive largely cancel each other
if strcmp(v_str, '1.0') || strcmp(v_str, '0.0')
    sigma_m_h   = sigma_h + kappa .* cpd ./ plot_level .* dtheta_dp_ma_h .* l2_h ./ k2_h;
    sigma_m_r   = sigma_r + kappa .* cpd ./ plot_level .* dtheta_dp_ma_r .* l2_r ./ k2_r;
    latent_h    = cpd .* dtheta_dp_ma_h .* omega_QG_h_rec;
    latent_r    = cpd .* dtheta_dp_ma_r .* omega_QG_r_rec;
elseif strcmp(v_str, '1.1')
    sigma_m_h   = sigma_h + kappa .* cpd ./ plot_level .* dtheta_dp_ma_h;
    sigma_m_r   = sigma_r + kappa .* cpd ./ plot_level .* dtheta_dp_ma_r;
    latent_h    = cpd .* dtheta_dp_ma_h .* omega_QG_h_rec .* k2_h./l2_h;
    latent_r    = cpd .* dtheta_dp_ma_r .* omega_QG_r_rec .* k2_r./l2_r;
end
denom_h     = sigma_m_h .* k2_h + F0.^2.*m2_h;

%excessive_h = cpd .* dtheta_dp_ma_h .* (omega_h - k2_h./l2_h .* omega_QG_h_rec);
%excessive_r = cpd .* dtheta_dp_ma_r .* (omega_r - k2_r./l2_r .* omega_QG_r_rec);
excessive_h = cpd .* dtheta_dp_ma_h .* omega_h - latent_h;
epsilon_h   = J_h - excessive_h - latent_h;
excessive_r = cpd .* dtheta_dp_ma_r .* omega_r - latent_r;
epsilon_r   = J_r - excessive_r - latent_r;
denom_r     = sigma_m_r .* k2_r + F0.^2.*m2_r;

% reconstruction
omega_QG_h_rec_moist = - (Adv_h + kappa / plot_level * l2_h .* (excessive_h + epsilon_h)) ./ denom_h;
omega_QG_r_rec_moist = - (Adv_r + kappa / plot_level * l2_r .* (excessive_r + epsilon_r)) ./ denom_r;

% linear expansion
d_sigma_m       = ((sigma_r - sigma_h) .* k2_h + Ra ./ Plevels .* l2_h .* (dtheta_dp_ma_r - dtheta_dp_ma_h)) .* ...
                    (Adv_h + kappa ./ plot_level .* l2_h .* (epsilon_h + excessive_h)) ./ denom_h.^2;
d_k2_l2_m       = ((k2_r - k2_h) .* sigma_h + Ra ./ Plevels .* (l2_r - l2_h) .* dtheta_dp_ma_h) .* ...
                    (Adv_h + kappa ./ plot_level .* l2_h .* (epsilon_h + excessive_h)) ./ denom_h.^2;
d_m2_m          = (m2_r - m2_h) .* repmat(F0, [1, 1, length(plevels)]).^2 .* ...
                    (Adv_h + kappa ./ plot_level .* l2_h .* (epsilon_h + excessive_h)) ./ denom_h.^2;
d_denom         = (denom_r - denom_h) .* (Adv_h + kappa ./ plot_level .* l2_h .* (epsilon_h + excessive_h)) ./ denom_h.^2;
d_Adv_m         = - (Adv_r - Adv_h) ./ denom_h;
d_epsilon       = - kappa ./ plot_level .* (l2_r.*epsilon_r - l2_h.*epsilon_h) ./ denom_h;
d_excessive     = - kappa ./ plot_level .* (l2_r.*excessive_r - l2_h.*excessive_h) ./ denom_h;
d_excessive_1   = - kappa ./ plot_level .* cpd .* (dtheta_dp_ma_r - dtheta_dp_ma_h) .* k2_h.*(omega_h - omega_QG_h_rec) ./ denom_h;
d_excessive_2   = - kappa ./ plot_level .* cpd .* (l2_r .* (omega_r - omega_QG_r_rec) - ...
                                                   l2_h .* (omega_h - omega_QG_h_rec)) .* dtheta_dp_ma_h ./ denom_h;
d_omega         = omega_r - omega_h;
d_rec_m         = d_sigma_m + d_k2_l2_m + d_m2_m + d_Adv_m;
d_rec_m_full      = d_sigma_m + d_k2_l2_m + d_m2_m + d_Adv_m + d_epsilon + d_excessive;

% zonal mean plots
%{
colors = get(gca,'colororder');
plot_x = lat_series(lat_indices);
[temp_ind] = abs(d_omega_QG)> 1.5 | ...
             abs(d_Adv_m)   > 1.5 | ...
             abs(d_k2_l2_m) > 1.5 | ...
             abs(d_m2_m)    > 1.5 | ...
             abs(d_Adv_m)   > 1.5 | ...
             abs(d_rec_m)     > 1.5 | ...
             abs(d_excessive(:, :, plot_ind)) > 1.5 | ...
             abs(d_epsilon  (:, :, plot_ind)) > 1.5 | ...
             denom_h(:, :, plot_ind) < 0;
NaN_matrix = ones(size(d_omega_QG));
NaN_matrix(temp_ind) = NaN;
%}


[temp_ind] = abs(d_omega_QG (:, :, plot_ind)) > 1.5 | ...
             abs(d_Adv_m    (:, :, plot_ind)) > 1.5 | ...
             abs(d_k2_l2_m  (:, :, plot_ind)) > 1.5 | ...
             abs(d_m2_m     (:, :, plot_ind)) > 1.5 | ...
             abs(d_Adv_m    (:, :, plot_ind)) > 1.5 | ...
             abs(d_rec_m    (:, :, plot_ind)) > 1.5 | ...
             abs(d_excessive(:, :, plot_ind)) > 1.5 | ...
             abs(d_epsilon  (:, :, plot_ind)) > 1.5 | ...
             denom_h(:, :, plot_ind) < 0;
NaN_matrix_m = ones(size(d_omega_QG  (:, :, plot_ind)));
NaN_matrix_m(temp_ind) = NaN;

temp_var = sqrt(nanvar([reshape(d_denom(:, :, plot_ind)     .* NaN_matrix_m, [length(NaN_matrix_m(:)), 1]), ...
                        reshape(d_Adv_m(:, :, plot_ind)     .* NaN_matrix_m, [length(NaN_matrix_m(:)), 1]), ...
                        reshape(d_rec_m(:, :, plot_ind)     .* NaN_matrix_m, [length(NaN_matrix_m(:)), 1]), ...
                        reshape(d_excessive(:, :, plot_ind) .* NaN_matrix_m, [length(NaN_matrix_m(:)), 1]), ...
                        reshape(d_epsilon(:, :, plot_ind)   .* NaN_matrix_m, [length(NaN_matrix_m(:)), 1])]));
[temp_ind_2] = abs(d_denom(:, :, plot_ind)     .*NaN_matrix_m - repmat(nanmean(d_denom(:, :, plot_ind)    .*NaN_matrix_m, 2), 1, size(NaN_matrix_m, 2))) > 3*temp_var(1) | ...
               abs(d_Adv_m(:, :, plot_ind)     .*NaN_matrix_m - repmat(nanmean(d_Adv_m(:, :, plot_ind)    .*NaN_matrix_m, 2), 1, size(NaN_matrix_m, 2))) > 3*temp_var(2) | ...
               abs(d_rec_m(:, :, plot_ind)     .*NaN_matrix_m - repmat(nanmean(d_rec_m(:, :, plot_ind)    .*NaN_matrix_m, 2), 1, size(NaN_matrix_m, 2))) > 3*temp_var(3) | ...
               abs(d_excessive(:, :, plot_ind) .*NaN_matrix_m - repmat(nanmean(d_excessive(:, :, plot_ind).*NaN_matrix_m, 2), 1, size(NaN_matrix_m, 2))) > 3*temp_var(4) | ...
               abs(d_epsilon(:, :, plot_ind)   .*NaN_matrix_m - repmat(nanmean(d_epsilon(:, :, plot_ind)  .*NaN_matrix_m, 2), 1, size(NaN_matrix_m, 2))) > 3*temp_var(5);
NaN_matrix_m(temp_ind_2 | temp_ind) = NaN;

FULL = false;
plot_decomposition_fancy
FULL = true;
plot_decomposition_fancy

temp = NaN_matrix_m;
plot_title = 'Moist NaN matrix';
plot_filename = 'NaN_matrix_m';
climit = [0, 1];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);


d_rec_m_full       = d_denom + d_Adv_m + d_epsilon + d_excessive;
figure('pos',figsize_zonal);
hold on;
grid on;
p0 = plot(plot_x, nanmean(d_omega         (:, :, plot_ind) .* NaN_matrix_m, 2) ./ nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2));
p1 = plot(plot_x, nanmean(d_omega_QG      (:, :, plot_ind) .* NaN_matrix_m, 2) ./ nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2));
p2 = plot(plot_x, nanmean(d_Adv_m       (:, :, plot_ind) .* NaN_matrix_m, 2) ./ nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2));
p3 = plot(plot_x, nanmean(d_denom       (:, :, plot_ind) .* NaN_matrix_m, 2) ./ nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2));
p4 = plot(plot_x, nanmean(d_epsilon     (:, :, plot_ind) .* NaN_matrix_m, 2) ./ nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2));
p5 = plot(plot_x, nanmean(d_excessive   (:, :, plot_ind) .* NaN_matrix_m, 2) ./ nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2));
p6 = plot(plot_x, nanmean(d_rec_m_full     (:, :, plot_ind) .* NaN_matrix_m, 2) ./ nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2));

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
if ~ABSOLUTE
    ylabel('-\Delta\omega_{rec}/\omega_{QG}');
else
    ylabel('\Delta\omega_{rec}');
end
hold off;
saveas(gca, [plot_path_ensemble, 'zonal_omega_decomposition_moist_', plot_level_name], 'png');

% contour plots

temp = d_k2_l2_m(:, :, plot_ind) ./ Omega_QG_h_mean;
plot_title = 'contribution of $k^2$ and $l^2$';
plot_filename = 'moist_d_k2_l2';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_sigma_m(:, :, plot_ind) ./ Omega_QG_h_mean;
plot_title = 'contribution of $\sigma_m$';
plot_filename = 'moist_d_sigma_m';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_m2_m(:, :, plot_ind) ./ Omega_QG_h_mean;
plot_title = 'contribution of $m^2$';
plot_filename = 'moist_d_m2';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_denom(:, :, plot_ind) ./ Omega_QG_h_mean;
plot_title = 'contribution of the denominator, $k^2(\sigma + \frac{RT}{p\theta}\frac{d\theta}{dp}|_{\theta^*} + f_0^2m^2)$';
plot_filename = 'moist_d_denom';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_epsilon(:, :, plot_ind) ./ Omega_QG_h_mean;
plot_title = 'contribution of $\epsilon$';
plot_filename = 'moist_d_epsilon';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_Adv_m(:, :, plot_ind) ./ Omega_QG_h_mean;
plot_title = 'contribution of RHS';
plot_filename = 'moist_d_Adv';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_excessive(:, :, plot_ind) ./ Omega_QG_h_mean;
plot_title = 'contribution of higher-order heating';
plot_filename = 'moist_d_higherorder';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = d_rec_m_full(:, :, plot_ind) ./ Omega_QG_h_mean;
plot_title = 'sum of all terms, moist decomposition';
plot_filename = 'moist_d_all';
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_ensemble, smooth_window_x, smooth_window_y, SMOOTH, climit, num_event_d, num_threshold, true);

temp = (d_k2_l2_m(:, :, plot_ind) + ...
       d_sigma_m(:, :, plot_ind) + ...
       d_m2_m(:, :, plot_ind) + ...
       d_Adv_m(:, :, plot_ind)) ./ Omega_QG_h_mean;
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
figure('pos', [10, 10, 400, 350])
temp_1 = J_r(:, :, plot_ind);
temp_2 = omega_r(:, :, plot_ind) .* dtheta_dp_ma_r(:, :, plot_ind) * cpd;
s1 = scatter(temp_1(:), temp_2(:), 1, 'o', 'filled', 'MarkerFaceColor', colors(2, :));
hold on
temp_1 = J_h(:, :, plot_ind);
temp_2 = omega_h(:, :, plot_ind) .* dtheta_dp_ma_h(:, :, plot_ind) * cpd;
s2 = scatter(temp_1(:), temp_2(:), 1.5, 'o', 'filled', 'MarkerFaceColor', colors(1, :));
plot([-10, 10], [-10, 10], 'k--');
legend([s2, s1], {'Historical', 'RCP85'}, 'location', 'best')
xlim([-0.05, 1.2])
ylim([-0.05, 1.2])
xlabel('Diabatic heating (W/kg)')
ylabel('$c_p\frac{T}{\theta}\frac{d\theta}{dp}\Big|_{\theta^*}\omega$', 'interpreter', 'latex')
ylb = get(gca,'ylabel');
ylb.Position = [-0.25, 0.5750 -1];
ylb.Rotation = 0;
%set(ylb, 'Rotation', 0)
set(gca, 'Position', [0.2000, 0.1500, 0.7500, 0.7500])
%title('Diabatic heating vs. $c_p\frac{T}{\theta}\frac{d\theta}{dp}|_\theta^*\omega$', 'interpreter', 'latex')
saveas(gca, [plot_path_ensemble, 'J_estimation_scatter'], 'png')
clf;


