%% moist decomposition
%{
% linear expansion
d_denom         = (denom_r - denom_h) .* (AB_h_mean + kappa ./ plot_level .* l2_h_comp .* (epsilon_h + excessive_h)) ./ denom_h.^2;
d_Adv           = - (AB_r_mean - AB_h_mean) ./ denom_h;
%d_epsilon       = - kappa ./ plot_level .* (epsilon_r - epsilon_h) .* l2_h_comp ./ denom_h;
%d_excessive     = - kappa ./ plot_level .* (excessive_r - excessive_h) .* l2_h_comp ./ denom_h;
d_epsilon       = - kappa ./ plot_level .* (l2_r_comp.*epsilon_r - l2_h_comp.*epsilon_h) ./ denom_h;
d_excessive     = - kappa ./ plot_level .* (l2_r_comp.*excessive_r - l2_h_comp.*excessive_h) ./ denom_h;

% accurate version:
%d_denom         = (denom_r - denom_h) .* (AB_h_mean + kappa ./ plot_level .* l2_h_comp .* (epsilon_h + excessive_h)) ./ denom_h ./ denom_r;
%d_Adv           = - (AB_r_mean - AB_h_mean) ./ denom_r;
%d_epsilon       = - kappa ./ plot_level .* (l2_r_comp.*epsilon_r - l2_h_comp.*epsilon_h) ./ denom_r;
%d_excessive     = - kappa ./ plot_level .* (l2_r_comp.*excessive_r - l2_h_comp.*excessive_h) ./ denom_r;

d_rec_all       = d_denom + d_Adv + d_epsilon + d_excessive;

d_denom         = d_denom           ./ omega_QG_h_rec_0;
d_epsilon       = d_epsilon         ./ omega_QG_h_rec_0;
d_Adv           = d_Adv             ./ omega_QG_h_rec_0;
d_excessive     = d_excessive       ./ omega_QG_h_rec_0;
d_rec_all       = d_rec_all         ./ omega_QG_h_rec_0;
%}

indN = plot_x > 40;
indS = plot_x < -40;
plot_x_N = plot_x(indN);
plot_x_S = plot_x(indS);
colors = get(gca,'colororder');
bias = 30;

figure('pos', [10 10 500 240]);
hold on;
%grid on;
plot(-70:70, zeros(size(-70:70)), 'LineWidth', 0.8, 'Color', [0.3, 0.3, 0.3], 'Linestyle', '--')
%clear('alpha')
%alpha(0.5)

%p0 = plot(plot_x_S + bias, d_all      (indS));
%{
p1 = plot(plot_x_S + bias, d_all_QG   (indS));
p2 = plot(plot_x_S + bias, d_Adv      (indS));
p3 = plot(plot_x_S + bias, d_denom    (indS));
p4 = plot(plot_x_S + bias, d_excessive(indS));
p5 = plot(plot_x_S + bias, d_epsilon  (indS));
p6 = plot(plot_x_S + bias, d_rec_all  (indS));
%}
%if ABSOLUTE
%    omega_QG_h_mean = ones(size(omega_QG_h_rec, 1));
%else
    omega_QG_h_mean = - nanmean(omega_QG_h_rec(:, :, plot_ind) .* NaN_matrix, 2);
%end
p1 = plot(plot_x_S + bias, nanmean(d_all_QG   (indS, :, plot_ind) .* NaN_matrix(indS, :), 2) ./ omega_QG_h_mean(indS));
p2 = plot(plot_x_S + bias, nanmean(d_Adv      (indS, :, plot_ind) .* NaN_matrix(indS, :), 2) ./ omega_QG_h_mean(indS));
%p3 = plot(plot_x_S + bias, nanmean(d_denom    (indS, :, plot_ind) .* NaN_matrix(indS, :), 2) ./ omega_QG_h_mean(indS));
p3_1 = plot(plot_x_S + bias, nanmean(d_sigma_m(indS, :, plot_ind) .* NaN_matrix(indS, :), 2) ./ omega_QG_h_mean(indS));
p3_2 = plot(plot_x_S + bias, nanmean(d_k2_l2  (indS, :, plot_ind) .* NaN_matrix(indS, :), 2) ./ omega_QG_h_mean(indS));
p3_3 = plot(plot_x_S + bias, nanmean(d_m2     (indS, :, plot_ind) .* NaN_matrix(indS, :), 2) ./ omega_QG_h_mean(indS));
%p4 = plot(plot_x_S + bias, nanmean(d_excessive(indS, :, plot_ind) .* NaN_matrix(indS, :), 2) ./ omega_QG_h_mean(indS));
p4_1 = plot(plot_x_S + bias, nanmean(d_excessive_1(indS, :, plot_ind) .* NaN_matrix(indS, :), 2) ./ omega_QG_h_mean(indS));
p4_2 = plot(plot_x_S + bias, nanmean(d_excessive_2(indS, :, plot_ind) .* NaN_matrix(indS, :), 2) ./ omega_QG_h_mean(indS));
p5 = plot(plot_x_S + bias, nanmean(d_epsilon  (indS, :, plot_ind) .* NaN_matrix(indS, :), 2) ./ omega_QG_h_mean(indS));
p6 = plot(plot_x_S + bias, nanmean(d_rec_all  (indS, :, plot_ind) .* NaN_matrix(indS, :), 2) ./ omega_QG_h_mean(indS));


%set(p0, 'LineWidth', 1, 'Color', 'red');
set(p1, 'LineWidth', 1.5, 'Color', 'black');
set(p2, 'LineWidth', 1.5, 'Color', colors(6, :));
%set(p3, 'LineWidth', 1.5, 'Color', 'red');
set(p3_1, 'LineWidth', 1.5, 'Color', colors(1, :));
set(p3_2, 'LineWidth', 1.5, 'Color', colors(3, :));
set(p3_3, 'LineWidth', 1.5, 'Color', colors(5, :));
%set(p4, 'LineWidth', 1.5, 'Color', colors(4, :));
set(p4_1, 'LineWidth', 1.5, 'Color', colors(4, :));
set(p4_2, 'LineWidth', 1.5, 'Color', colors(2, :));
set(p5, 'LineWidth', 1.5, 'Color', colors(7, :));
set(p6, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', ':');

%p0 = plot(plot_x_N - bias, d_all      (indN));
p1 = plot(plot_x_N - bias, nanmean(d_all_QG   (indN, :, plot_ind) .* NaN_matrix(indN, :), 2) ./ omega_QG_h_mean(indN));
p2 = plot(plot_x_N - bias, nanmean(d_Adv      (indN, :, plot_ind) .* NaN_matrix(indN, :), 2) ./ omega_QG_h_mean(indN));
%p3 = plot(plot_x_N - bias, nanmean(d_denom    (indN, :, plot_ind) .* NaN_matrix(indN, :), 2) ./ omega_QG_h_mean(indN));
p3_1 = plot(plot_x_N - bias, nanmean(d_sigma_m(indN, :, plot_ind) .* NaN_matrix(indN, :), 2) ./ omega_QG_h_mean(indN));
p3_2 = plot(plot_x_N - bias, nanmean(d_k2_l2  (indN, :, plot_ind) .* NaN_matrix(indN, :), 2) ./ omega_QG_h_mean(indN));
p3_3 = plot(plot_x_N - bias, nanmean(d_m2     (indN, :, plot_ind) .* NaN_matrix(indN, :), 2) ./ omega_QG_h_mean(indN));
%p4 = plot(plot_x_N - bias, nanmean(d_excessive(indN, :, plot_ind) .* NaN_matrix(indN, :), 2) ./ omega_QG_h_mean(indN));
p4_1 = plot(plot_x_N - bias, nanmean(d_excessive_1(indN, :, plot_ind) .* NaN_matrix(indN, :), 2) ./ omega_QG_h_mean(indN));
p4_2 = plot(plot_x_N - bias, nanmean(d_excessive_2(indN, :, plot_ind) .* NaN_matrix(indN, :), 2) ./ omega_QG_h_mean(indN));
p5 = plot(plot_x_N - bias, nanmean(d_epsilon  (indN, :, plot_ind) .* NaN_matrix(indN, :), 2) ./ omega_QG_h_mean(indN));
p6 = plot(plot_x_N - bias, nanmean(d_rec_all  (indN, :, plot_ind) .* NaN_matrix(indN, :), 2) ./ omega_QG_h_mean(indN));

set(p1, 'LineWidth', 1.5, 'Color', 'black');
set(p2, 'LineWidth', 1.5, 'Color', colors(6, :));
%set(p3, 'LineWidth', 1.5, 'Color', 'red');
set(p3_1, 'LineWidth', 1.5, 'Color', colors(1, :));
set(p3_2, 'LineWidth', 1.5, 'Color', colors(3, :));
set(p3_3, 'LineWidth', 1.5, 'Color', colors(5, :));
%set(p4, 'LineWidth', 1.5, 'Color', colors(4, :));
set(p4_1, 'LineWidth', 1.5, 'Color', colors(4, :));
set(p4_2, 'LineWidth', 1.5, 'Color', colors(2, :));
set(p5, 'LineWidth', 1.5, 'Color', colors(7, :));
set(p6, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', ':');




axis([-70 + bias 70 - bias -0.300 0.301]);
xticks([-60.0 + bias, -40.0 + bias, 0.0, 40.0 - bias, 60.0 - bias])
xticklabels({'-60', '-40', '0', '40', '60'})

lgd = legend( ...
[p1, p2, p3_1, p3_2, p3_3, p4_1, p4_2, p5, p6], ...
{...'$\Delta\omega$' , ...
'$\Delta\omega_{QG}$' , ...
'$\Delta\omega_{Adv}$' , ...
...'$\Delta\omega_{denom}$' , ...
'$\Delta\omega_{\sigma_m}$', ...
'$\Delta\omega_{k^2,l^2}$', ...
'$\Delta\omega_{m^2}$', ...
'$\Delta\omega_{hi-order1}$' , ...
'$\Delta\omega_{hi-order2}$' , ...
'$\Delta\omega_{\epsilon}$' , ...
'$\Delta\omega_{rec}$' , ...
}, ...
'location', 'best', ...
'interpreter', 'latex');
lgd.Position = [0.460, 0.3523, 0.1147, 0.31];
title([plot_level_name, ' Moist Decomposition']);
xlabel('Latitude');
%if ~ABSOLUTE    
    ylabel('-\Delta\omega/\omega');
%else
%    ylabel('\Delta\omega (Pa/s)');
%end
hold off;
break_axis('position', -5.0, 'axis', 'x', 'length', 0.2)
break_axis('position',  5.0, 'axis', 'x', 'length', 0.2)

saveas(gca, [plot_path_ensemble, 'omega_decomposition_moist_fancy_', plot_level_name, '_full'], 'png');


