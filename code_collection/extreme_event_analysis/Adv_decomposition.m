% Paul's idea on Mar. 26th, 2019. In the email he wrote:

%   Lastly, could you check if the approximation omega = adv/(f^2m^2)
%   would work well for climate change (e.g. lower panels of fig 2 and fig
%   4)? This means assuming sigma_m=J'=0 for both the control climate and
%   climate change.

%% Adv decomposition
% new form according to Paul:
d_m2_2    = - 1./(F0.^2.*m2_h) .* F0.^2 .* (m2_r - m2_h) .* omega_QG_h_rec;
d_Adv_2   = - 1./(F0.^2.*m2_h) .* (Adv_r - Adv_h);
% correct dry form:
% d_m2_2    = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* F0.^2 .* (m2_r - m2_h) .* omega_QG_h_rec;
% d_Adv_2   = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* (Adv_r - Adv_h);
d_rec_2   = d_m2_2 + d_Adv_2;
d_omega_QG = (omega_QG_r - omega_QG_h);
d_omega    = (omega_r - omega_h);

colors = get(gca,'colororder');
plot_x = lat_series(lat_indices);

[temp_ind] = abs(d_omega             (:, :, plot_ind)) > 1.5 | ...
             abs(d_omega_QG          (:, :, plot_ind)) > 1.5 | ...
             abs(d_m2_2                (:, :, plot_ind)) > 1.5 | ...
             abs(d_Adv_2               (:, :, plot_ind)) > 1.5;
NaN_matrix_2 = ones(size(d_omega_QG(:, :, plot_ind)));
NaN_matrix_2(temp_ind) = NaN;

figure('pos',figsize_zonal);
hold on;
plot(-90:90, zeros(size(-90:90)), 'LineWidth', 0.8, 'Color', [0.3, 0.3, 0.3], 'Linestyle', '--')
%omega_h_mean_2    = nanmean(omega_h   (:, :, plot_ind) .* NaN_matrix_2, 2);
%omega_QG_h_mean_2 = nanmean(omega_QG_h(:, :, plot_ind) .* NaN_matrix_2, 2);

p0 = plot(plot_x, nanmean(d_omega   (:, :, plot_ind).*NaN_matrix_2, 2) ./ omega_h_mean);
p1 = plot(plot_x, nanmean(d_omega_QG(:, :, plot_ind).*NaN_matrix_2, 2) ./ omega_QG_h_mean);
p5 = plot(plot_x, nanmean(d_m2_2    (:, :, plot_ind).*NaN_matrix_2, 2) ./ omega_QG_h_mean);
p6 = plot(plot_x, nanmean(d_Adv_2   (:, :, plot_ind).*NaN_matrix_2, 2) ./ omega_QG_h_mean);
p8 = plot(plot_x, nanmean(d_rec_2   (:, :, plot_ind).*NaN_matrix_2, 2) ./ omega_QG_h_mean);

set(p0, 'LineWidth', 1.5, 'Color', 'red');
set(p1, 'LineWidth', 1.5, 'Color', 'black');
set(p5, 'LineWidth', 1.5, 'Color', colors(5, :));
set(p6, 'LineWidth', 1.5, 'Color', colors(6, :));
set(p8, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', ':');


lgd = legend([p0, p1, p2, p4, p6], ...
{'$\Delta\omega/\omega_{QG}$' , ...
'$\Delta\omega_{QG}/\omega_{QG}$' , ...
'$\Delta\omega_{Adv}/\omega_{QG}$' , ...
'$\Delta\omega_{m^2}/\omega_{QG}$' , ...
'$\Delta\omega_{rec}/\omega_{QG}$' , ...
}, ...
'location', 'best', ...
'interpreter', 'latex');
lgd.Position = [0.460, 0.3523, 0.1147, 0.31];
title([plot_level_name, ' Adv Decomposition']);
axis([-70 70 -0.601 0.601]);
xlabel('Latitude');
ylabel('\Delta\omega_{QG}/\omega_{QG}');

hold off;
saveas(gca, [plot_path_ensemble, 'zonal_omega_decomposition_Adv_', plot_level_name], 'png');




