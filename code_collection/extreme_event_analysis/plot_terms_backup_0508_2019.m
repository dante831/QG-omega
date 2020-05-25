
%climate_error{1}.rhs_mean1

set(0,'DefaultFigureVisible','off');


error_h = - C_h + sigma_h .* (- k2_h) .* omega_QG_h + F0.^2 .* (- m2_h) .* omega_QG_h - Adv_h;
error_r = - C_r + sigma_r .* (- k2_r) .* omega_QG_r + F0.^2 .* (- m2_r) .* omega_QG_r - Adv_r;
term1_h = - C_h + sigma_h .* (- k2_h) .* omega_QG_h;
term1_r = - C_r + sigma_r .* (- k2_r) .* omega_QG_r;
term2_h = Adv_h;
term2_r = Adv_r;
term3_h = F0.^2 .* (- m2_h) .* omega_QG_h;
term3_r = F0.^2 .* (- m2_r) .* omega_QG_r;
term4_h = sigma_h .* (- k2_h) .* omega_QG_h;
term4_r = sigma_r .* (- k2_r) .* omega_QG_r;
term5_h = - C_h;
term5_r = - C_r;

error_h = nanmean(error_h(:, :, plot_ind), 2);
error_r = nanmean(error_r(:, :, plot_ind), 2);
term1_h = nanmean(term1_h(:, :, plot_ind), 2);
term1_r = nanmean(term1_r(:, :, plot_ind), 2);
term2_h = nanmean(term2_h(:, :, plot_ind), 2);
term2_r = nanmean(term2_r(:, :, plot_ind), 2);
term3_h = nanmean(term3_h(:, :, plot_ind), 2);
term3_r = nanmean(term3_r(:, :, plot_ind), 2);
term4_h = nanmean(term4_h(:, :, plot_ind), 2);
term4_r = nanmean(term4_r(:, :, plot_ind), 2);
term5_h = nanmean(term5_h(:, :, plot_ind), 2);
term5_r = nanmean(term5_r(:, :, plot_ind), 2);

% segregate the data sequence to plot filled graph
ind = plot_x >= 0 | plot_x <= -0;
if ind(1) == 1
    starts = [1, find(diff(ind) == 1) + 1];
else
    starts = find(diff(ind) == 1) + 1;
end
if ind(end) == 1
    ends = [find(diff(ind) == -1), length(plot_x)];
else
    ends = find(diff(ind) == -1);
end

error_h(~ind) = NaN;
error_r(~ind) = NaN;
term1_h(~ind) = NaN;
term1_r(~ind) = NaN;
term2_h(~ind) = NaN;
term2_r(~ind) = NaN;
term3_h(~ind) = NaN;
term3_r(~ind) = NaN;
term4_h(~ind) = NaN;
term4_r(~ind) = NaN;
term5_h(~ind) = NaN;
term5_r(~ind) = NaN;

figure('pos',[10 10 500 200]);
hold on;
%grid on;
plot(-90:90, zeros(size(-90:90)), 'LineWidth', 0.8, 'Color', [0.3, 0.3, 0.3], 'Linestyle', ':')
plot_x = lat_series(lat_indices);
p11 = plot(plot_x, error_h, 'LineWidth', 1, 'Color', 'blue');
p12 = plot(plot_x, error_r, 'LineWidth', 1, 'Color', 'blue', 'LineStyle', '--');
%p21 = plot(plot_x, term1_h, 'LineWidth', 1, 'Color', 'magenta');
%p22 = plot(plot_x, term1_r, 'LineWidth', 1, 'Color', 'magenta', 'LineStyle', '--');
p31 = plot(plot_x, term2_h, 'LineWidth', 1, 'Color', 'cyan');
p32 = plot(plot_x, term2_r, 'LineWidth', 1, 'Color', 'cyan', 'LineStyle', '--');
p41 = plot(plot_x, term3_h, 'LineWidth', 1, 'Color', 'red');
p42 = plot(plot_x, term3_r, 'LineWidth', 1, 'Color', 'red', 'LineStyle', '--');
p51 = plot(plot_x, term4_h, 'LineWidth', 1, 'Color', 'black');
p52 = plot(plot_x, term4_r, 'LineWidth', 1, 'Color', 'black', 'LineStyle', '--');
p61 = plot(plot_x, term5_h, 'LineWidth', 1, 'Color', 'green');
p62 = plot(plot_x, term5_r, 'LineWidth', 1, 'Color', 'green', 'LineStyle', '--');
clear('alpha')
for i = 1 : length(starts)
    ind2 = starts(i) : ends(i);
    fill([plot_x(ind2)' fliplr(plot_x(ind2)')], [error_h(ind2)' fliplr(error_r(ind2)')], 'blue', 'edgealpha', 0);
    %fill([plot_x(ind2)' fliplr(plot_x(ind2)')], [term1_h(ind2)' fliplr(term1_r(ind2)')], 'magenta', 'edgealpha', 0);
    %fill([plot_x(ind2)' fliplr(plot_x(ind2)')], [term2_h(ind2)' fliplr(term2_r(ind2)')], 'yellow', 'edgealpha', 0);
    fill([plot_x(ind2)' fliplr(plot_x(ind2)')], [term2_h(ind2)' fliplr(term2_r(ind2)')], 'cyan', 'edgealpha', 0);
    fill([plot_x(ind2)' fliplr(plot_x(ind2)')], [term3_h(ind2)' fliplr(term3_r(ind2)')], 'red', 'edgealpha', 0);
    fill([plot_x(ind2)' fliplr(plot_x(ind2)')], [term4_h(ind2)' fliplr(term4_r(ind2)')], 'black', 'edgealpha', 0);
    fill([plot_x(ind2)' fliplr(plot_x(ind2)')], [term5_h(ind2)' fliplr(term5_r(ind2)')], 'green', 'edgealpha', 0);
    
    alpha(.15);
end

lgd = legend( ...
[p11, p31, p41, p51, p61], ...
{'numerical error' , ...
'$Adv$', ...
'$f_0^2\frac{\partial^2}{\partial p^2}\omega_{QG}$', ...
'$\nabla^2(\sigma\omega_{QG})$', ...
'$\frac{\kappa}{p}\nabla^2J$'}, ...
'location', 'best', ...
'Interpreter', 'latex');
lgd.Position = [0.460, 0.3823, 0.11, 0.28];

if ~isempty(strfind(pwd_str, 'CESM'))
    axis([-80 80 -2.5e-16 2.5e-16]);
else
    axis([-80 80 -0.5e-16 0.5e-16]);
end
ylabel('$m/(kg\cdot s)$', 'Interpreter', 'latex');
xticks([-60:20:60])
xticklabels({'$60^\circ$S', '$40^\circ$S', '$20^\circ$S', '$0^\circ$', '$20^\circ$N', '$40^\circ$N', '$60^\circ$N'})
set(gca, 'TickLabelInterpreter', 'latex')
hold off;
saveas(gca, [plot_path_ensemble, 'terms_', plot_level_name], 'png');



