
addpath('/disk7/ziweili/CESM_LENS/source')
addpath('/disk7/ziweili/CESM_LENS/exp/ensemble_plots/')

clear all
pwd_str = pwd;

if ~isempty(strfind(pwd_str, 'CESM')) && ...
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
    num_threshold = 1;
elseif ~isempty(strfind(pwd_str, 'CESM')) && ...
    ~isempty(strfind(pwd_str, 'varsig_Q_findcenter'))
    filenames = {['../001_2D_varsig_Q_findcenter/plot_diagnostic.mat'], ...
                 ['../002_2D_varsig_Q_findcenter/plot_diagnostic.mat'], ...
                 ['../003_2D_varsig_Q_findcenter/plot_diagnostic.mat'], ...
                 ['../004_2D_varsig_Q_findcenter/plot_diagnostic.mat'], ...
                 ['../005_2D_varsig_Q_findcenter/plot_diagnostic.mat'], ...
                 ['../035_2D_varsig_Q_findcenter/plot_diagnostic.mat'], ...
                 };
    GFDL = false;
    num_threshold = 40;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'varsig_Q_findcenter_wb'))
    filenames = {['../2D_varsig_Q_findcenter_wb/plot_diagnostic.mat']};
    GFDL = true;
    num_threshold = 1;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'varsig_Q_findcenter')) && ...
       ~isempty(strfind(pwd_str, 'daily'))
    filenames = {['../2D_varsig_Q_findcenter_daily_99.5/plot_diagnostic.mat']};
    GFDL = true;
    num_threshold = 1;
elseif ~isempty(strfind(pwd_str, 'GFDL')) && ...
       ~isempty(strfind(pwd_str, 'varsig_Q_findcenter'))
    filenames = {['../2D_varsig_Q_findcenter/plot_diagnostic.mat']};
    GFDL = true;
    num_threshold = 1;
end

ensemble_read_data_v1;

figsize_zonal = [10, 10, 600, 200];

if GFDL
    SMOOTH = false;
else
    SMOOTH = true;
end
plot_level = 50000;
plot_level_name = [num2str(plot_level/100.), 'hPa'];

plevels = [100000, 97500, 95000, 92500, 90000, 87500, 85000, 82500, 80000, 77500, 75000, 70000, 65000, ...
            60000, 55000, 50000, 45000, 40000, 35000, 30000, 25000, 22500, 20000, 17500, 15000]';
plot_ind = plevels == plot_level;
set(0,'DefaultFigureVisible','off');
omega_range = [-1.5, 1.5];
change_range = [-0.5, 0.5];

plot_path_ensemble = ['plots/'];
if ~exist(plot_path_ensemble)
    mkdir(plot_path_ensemble);
end

num_event_h = num_event_h(:, :, plevels == plot_level);
num_event_r = num_event_r(:, :, plevels == plot_level);
num_event_d = reshape(min([num_event_h(:)'; num_event_r(:)']), size(num_event_h));
if exist('num_event_discarded_h')
    num_event_discarded_h = num_event_discarded_h(:, :, plevels == plot_level);
    num_event_discarded_r = num_event_discarded_r(:, :, plevels == plot_level);
end

omega_QG_h_rec =  - 1 ./ (k2_h .* sigma_h + m2_h .* F0.^2) .* (Adv_h + C_h);
omega_QG_r_rec =  - 1 ./ (k2_r .* sigma_r + m2_r .* F0.^2) .* (Adv_r + C_r);

level_ind = plevels == 50000;

%% comparison between omega and omega_QG in the two climates
% the default case
omega_h_mean_2    = nanmean(omega_h(:, :, plot_ind), 2); 
omega_r_mean_2    = nanmean(omega_r(:, :, plot_ind), 2);
omega_QG_h_mean_2 = nanmean(omega_QG_h(:, :, plot_ind), 2);
omega_QG_r_mean_2 = nanmean(omega_QG_r(:, :, plot_ind), 2);
% Paul's interest: Latent heating contribution
omega_QG_h_heating =  - 1 ./ (k2_h .* sigma_h + m2_h .* F0.^2) .* (C_h);
omega_QG_r_heating =  - 1 ./ (k2_r .* sigma_r + m2_r .* F0.^2) .* (C_r);
omega_QG_h_heating_mean = nanmean(omega_QG_h_heating(:, :, plot_ind), 2);
omega_QG_r_heating_mean = nanmean(omega_QG_r_heating(:, :, plot_ind), 2);
% plot zonal mean omega vs. omega_QG
figname = [plot_path_ensemble, 'zonal_omega_comparison_', num2str(plevels(level_ind)/100), 'hPa'];
zonal_omega
% the full-boundary case (even if you compared this case, you don't have enough information for a decomposition)
% You will still have to get a run with full boundaries
omega_h_mean_2    = nanmean(omega_h(:, :, plot_ind), 2);
omega_r_mean_2    = nanmean(omega_r(:, :, plot_ind), 2);
if exist('omega_QG_b_h')
    omega_QG_h_mean_2 = nanmean(omega_QG_b_h(:, :, plot_ind), 2);
    omega_QG_r_mean_2 = nanmean(omega_QG_b_r(:, :, plot_ind), 2);
end
figname = [plot_path_ensemble, 'zonal_omega_comparison_full_b_', num2str(plevels(level_ind)/100), 'hPa'];
zonal_omega

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
                                  d_dtheta_dp_ma_omega(:, :, plot_ind)), 1.5);

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

% moist reconstruction
omega_QG_h_rec_moist = - (Adv_h + kappa ./ Plevels .* l2_h .* J_res_h) ./ denom_h;
omega_QG_r_rec_moist = - (Adv_r + kappa ./ Plevels .* l2_r .* J_res_r) ./ denom_r;

d_excessive     = d_excessive_1 + d_excessive_2;
d_denom         = d_sigma_m + d_k2_l2_m + d_m2_m;
d_rec_m_full = d_sigma_m + d_k2_l2_m + d_m2_m + d_Adv_m + d_J_res;
d_rec_m      = d_sigma_m + d_k2_l2_m + d_m2_m + d_Adv_m;

NaN_matrix_m = get_NaN_matrix(cat(3, d_omega_QG (:, :, plot_ind), ...
                                  d_Adv_m    (:, :, plot_ind), ...
                                  d_k2_l2_m  (:, :, plot_ind), ...
                                  d_m2_m     (:, :, plot_ind), ...
                                  d_Adv_m    (:, :, plot_ind), ...
                                  d_rec_m    (:, :, plot_ind), ...
                                  d_J_res  (:, :, plot_ind), ...
                                  d_epsilon  (:, :, plot_ind)), 1.5);

omega_QG_h_mean_m = nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_m, 2);
omega_QG_r_mean_m = nanmean(omega_QG_r(:, :, plot_ind) .* NaN_matrix_m, 2);
Omega_QG_h_mean_m = repmat(omega_QG_h_mean_m, 1, size(omega_QG_h, 2));
omega_h_mean_m    = nanmean(omega_h   (:, :, plot_ind) .* NaN_matrix_m, 2);
omega_r_mean_m    = nanmean(omega_r   (:, :, plot_ind) .* NaN_matrix_m, 2);
Omega_h_mean_m    = repmat(omega_h_mean_m, 1, size(omega_h, 2));


% read in global temperature changes
if GFDL
    load('/disk7/ziweili/test1_GFDL/exp/global_temperature/global_avg_T.mat');
else
    load('/disk7/ziweili/CESM_LENS/exp/global_temperature/global_avg_T.mat');
end
delta_T = mean(T_avg_r) - mean(T_avg_h);

% Figure 1
FULL = true;
plot_decomposition_fancy

% Figure 2
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

%% J term estimation via moist-adiabatic process
figure('pos', [10, 10, 300, 270])
ind_lat = abs(lat) >= 30;
temp_1 = J_r(ind_lat, :, plot_ind);
temp_2 = - omega_r(ind_lat, :, plot_ind) .* sigma_star_r(ind_lat, :, plot_ind) * plot_level / kappa;
temp_3 = - omega_r(ind_lat, :, plot_ind) .* sigma_r(ind_lat, :, plot_ind) * plot_level / kappa;
temp_4 = omega_QG_r(ind_lat, :, plot_ind) .* dtheta_dp_ma_r(ind_lat, :, plot_ind) * cpd;
s1 = scatter(temp_2(:), temp_1(:), 1, 'o', 'filled', 'MarkerFaceColor', colors(2, :));
hold on
temp_1 = J_h(ind_lat, :, plot_ind);
temp_2 = - omega_h(ind_lat, :, plot_ind) .* sigma_star_h(ind_lat, :, plot_ind) * plot_level / kappa;
temp_3 = - omega_h(ind_lat, :, plot_ind) .* sigma_h(ind_lat, :, plot_ind) * plot_level / kappa;
temp_4 = omega_QG_h(ind_lat, :, plot_ind) .* dtheta_dp_ma_h(ind_lat, :, plot_ind) * cpd;
s2 = scatter(temp_2(:), temp_1(:), 1.5, 'o', 'filled', 'MarkerFaceColor', colors(1, :));
plot([-10, 10], [-10, 10], 'k--');
legend([s2, s1], {'Historical', 'RCP8.5'}, 'location', 'southeast', 'interpreter', 'latex')
legend boxoff
xlim([-0.05, 1.2])
ylim([-0.05, 1.2])
xlabel('$-\frac{p}{\kappa}\omega\sigma^*$ (J kg$^{-1}$s$^{-1}$)', 'interpreter', 'latex')
ylabel('$J$ (J kg$^{-1}$s$^{-1}$)', 'interpreter', 'latex')
set(gca, 'Position', [0.17, 0.15, 0.75, 0.75])
set(gca, 'TickLabelInterpreter','latex')
set(gca, 'TickDir', 'out');
saveas(gca, [plot_path_ensemble, 'J_estimation_scatter_', plot_level_name], 'png')
clf;

% plot for AGU pre
plot_for_pre_12122018

% omega_QG vs omega plot
omega_vs_omega_QG

%% Some diagnostics useful for the paper
ind_lat = abs(lat) > 30;
Lat = repmat(lat, 1, length(lon));
Lon = repmat(lon', length(lat), 1);
dlambda = (lon(2) - lon(1)) / 180 * pi;
dlat = lat(2) - lat(1);
dS = R * (sin(min(Lat + dlat/2, 90)/180*pi) - sin(max(Lat - dlat/2, -90)/180*pi)) * dlambda;
dS_d = dS(ind_lat, :) .* NaN_matrix_d(ind_lat, :);
dS_m = dS(ind_lat, :) .* NaN_matrix_m(ind_lat, :);
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
                    nanmean(nanmean(dS .* var(ind_lat, :, plot_ind) .* NaN_matrix(ind_lat, :), 2) ...
                            ./ omega_mean(ind_lat)) / ...
                    nanmean(nanmean(dS)) / delta_T;
% Paul's idea of comparing the rates of means, and means of rates
spherical_mean_v2 = @(NaN_matrix, dS, var, omega, ind_lat, delta_T, plot_ind) ...
                    nanmean(nanmean(dS .* var(ind_lat, :, plot_ind) .* NaN_matrix(ind_lat, :) ...
                            ./ omega(ind_lat, :, plot_ind), 2)) / ...
                    nanmean(nanmean(dS)) / delta_T;
                    % result, less than 2% difference from the spherical_mean function
d_omega_mean_dry    = spherical_mean(NaN_matrix_d, dS_d, d_omega, omega_h_mean_d, ind_lat, delta_T, plot_ind);
d_omega_mean_dry_v2 = spherical_mean_v2(NaN_matrix_d, dS_d, d_omega, omega_h, ind_lat, delta_T, plot_ind);
d_omega_QG_mean_dry = spherical_mean(NaN_matrix_d, dS_d, d_omega_QG, omega_QG_h_mean_d, ind_lat, delta_T, plot_ind);
d_omega_QG_mean_dry_v2 = spherical_mean_v2(NaN_matrix_d, dS_d, d_omega_QG, omega_QG_h, ind_lat, delta_T, plot_ind);
% comparing omega and omega_QG
omega_mean_temp = spherical_mean(NaN_matrix_d, dS_d, omega_h, ones(size(omega_h_mean_d)), ind_lat, 1, plot_ind);
omega_QG_mean_temp = spherical_mean(NaN_matrix_d, dS_d, omega_QG_h, ones(size(omega_h_mean_d)), ind_lat, 1, plot_ind);
disp(['omega_QG underestimates omega by ', num2str(1 - omega_QG_mean_temp / omega_mean_temp)]) 
% percentage of omega_QG's underestimation of omega

% averaged contributions, dry decomposition
d_precip     = spherical_mean(NaN_matrix_d, dS_d, precip_r - precip_h, precip_h_mean_d, ...
                    ind_lat, delta_T, 1);
d_sigma_mean = spherical_mean(NaN_matrix_d, dS_d, d_sigma, omega_QG_h_mean_d, ind_lat, delta_T, plot_ind);
d_J_mean     = spherical_mean(NaN_matrix_d, dS_d, d_J,     omega_QG_h_mean_d, ind_lat, delta_T, plot_ind);
d_k2_l2_mean = spherical_mean(NaN_matrix_d, dS_d, d_k2_l2, omega_QG_h_mean_d, ind_lat, delta_T, plot_ind);
d_m2_mean    = spherical_mean(NaN_matrix_d, dS_d, d_m2,    omega_QG_h_mean_d, ind_lat, delta_T, plot_ind);
d_Adv_mean   = spherical_mean(NaN_matrix_d, dS_d, d_Adv,   omega_QG_h_mean_d, ind_lat, delta_T, plot_ind);

d_omega_mean_moist    = spherical_mean(NaN_matrix_m, dS_m, d_omega, omega_h_mean_m, ind_lat, delta_T, plot_ind);
d_omega_QG_mean_moist = spherical_mean(NaN_matrix_m, dS_m, d_omega_QG, omega_QG_h_mean_m, ind_lat, delta_T, plot_ind);

% averaged contributions, moist decomposition
d_sigma_m_mean = spherical_mean(NaN_matrix_m, dS_m, d_sigma_m, omega_QG_h_mean_m, ind_lat, delta_T, plot_ind);
d_J_res_mean   = spherical_mean(NaN_matrix_m, dS_m, d_J_res,   omega_QG_h_mean_m, ind_lat, delta_T, plot_ind);
d_k2_l2_m_mean = spherical_mean(NaN_matrix_m, dS_m, d_k2_l2_m, omega_QG_h_mean_m, ind_lat, delta_T, plot_ind);
d_m2_m_mean    = spherical_mean(NaN_matrix_m, dS_m, d_m2_m,    omega_QG_h_mean_m, ind_lat, delta_T, plot_ind);
d_Adv_m_mean   = spherical_mean(NaN_matrix_m, dS_m, d_Adv_m,   omega_QG_h_mean_m, ind_lat, delta_T, plot_ind);

% generation of latex code for a table
% save values for the table
matfilename = 'extra_tropical_means';
save([pwd_str, '/', matfilename], 'd_precip', ...
    'd_omega_mean_dry'  , 'd_omega_QG_mean_dry'  , ...
    'd_omega_mean_moist', 'd_omega_QG_mean_moist', ...
    'd_sigma_mean'  , 'd_J_mean'    , 'd_k2_l2_mean'  , 'd_m2_mean'  , 'd_Adv_mean', ...
    'd_sigma_m_mean', 'd_J_res_mean', 'd_k2_l2_m_mean', 'd_m2_m_mean', 'd_Adv_m_mean');

