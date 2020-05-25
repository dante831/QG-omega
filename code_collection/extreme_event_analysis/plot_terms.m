
set(0,'DefaultFigureVisible','off');

error_h = - C_h + sigma_h .* (- k2_h) .* omega_QG_h + F0.^2 .* (- m2_h) .* omega_QG_h - Adv_h;
error_r = - C_r + sigma_r .* (- k2_r) .* omega_QG_r + F0.^2 .* (- m2_r) .* omega_QG_r - Adv_r;
term2_h = - Adv_h;
term2_r = - Adv_r;
term3_h = F0.^2 .* (- m2_h) .* omega_QG_h;
term3_r = F0.^2 .* (- m2_r) .* omega_QG_r;
term4_h = sigma_h .* (- k2_h) .* omega_QG_h;
term4_r = sigma_r .* (- k2_r) .* omega_QG_r;
term5_h = - C_h;
term5_r = - C_r;

error_h = nanmean(error_h(:, :, plot_ind), 2);
error_r = nanmean(error_r(:, :, plot_ind), 2);
term2_h = nanmean(term2_h(:, :, plot_ind), 2);
term2_r = nanmean(term2_r(:, :, plot_ind), 2);
term3_h = nanmean(term3_h(:, :, plot_ind), 2);
term3_r = nanmean(term3_r(:, :, plot_ind), 2);
term4_h = nanmean(term4_h(:, :, plot_ind), 2);
term4_r = nanmean(term4_r(:, :, plot_ind), 2);
term5_h = nanmean(term5_h(:, :, plot_ind), 2);
term5_r = nanmean(term5_r(:, :, plot_ind), 2);

if OCEAN_DESERT
    latmin_N = 5;
    latmax_S = -5;
    bias = 0;
else
    bias = 25;
    latmin_N = 30;
    latmax_S = -30;
end
% segregate the data sequence to plot filled graph
ind = plot_x >= latmin_N | plot_x <= latmax_S;
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
term2_h(~ind) = NaN;
term2_r(~ind) = NaN;
term3_h(~ind) = NaN;
term3_r(~ind) = NaN;
term4_h(~ind) = NaN;
term4_r(~ind) = NaN;
term5_h(~ind) = NaN;
term5_r(~ind) = NaN;

indN = plot_x > latmin_N;
indS = plot_x < latmax_S;
plot_x_N = plot_x(indN);
plot_x_S = plot_x(indS);

figure('pos',[10 10 500 200]);
hold on;
plot(-90:90, zeros(size(-90:90)), 'LineWidth', 0.8, 'Color', [0.3, 0.3, 0.3], 'Linestyle', ':')
plot_x = lat_series(lat_indices);

%{
colors_2 = [[170, 68,153]
          [ 68,170,153]; ...
          [ 17,119, 51]; ...
          [153,153, 51]; ...
          [204,102,119]; ...
          [221,204,119]; ...
          [136, 34, 85]; ...
          [ 51, 34,136]; ...
          [136,204,238]; ...
          [170, 68,153]];
colors_2 = colors_2 / 256;
%}
colors_2 = get(gca,'colororder');
colors_2(7, 1) = colors_2(7, 1) + 0.2;
colors_2(7, 3) = colors_2(7, 3) - 0.05;
colors_2 = colors_2([1:3, 5:7], :);
%{
p11_S = plot(plot_x_S + bias, error_h(indS), 'LineWidth', 1.5, 'Color', colors_2(1, :));
p12_S = plot(plot_x_S + bias, error_r(indS), 'LineWidth', 1.5, 'Color', colors_2(1, :), 'LineStyle', '--');
p31_S = plot(plot_x_S + bias, term2_h(indS), 'LineWidth', 1.5, 'Color', 'cyan');
p32_S = plot(plot_x_S + bias, term2_r(indS), 'LineWidth', 1.5, 'Color', 'cyan', 'LineStyle', '--');
p41_S = plot(plot_x_S + bias, term3_h(indS), 'LineWidth', 1.5, 'Color', 'red');
p42_S = plot(plot_x_S + bias, term3_r(indS), 'LineWidth', 1.5, 'Color', 'red', 'LineStyle', '--');
p51_S = plot(plot_x_S + bias, term4_h(indS), 'LineWidth', 1.5, 'Color', 'black');
p52_S = plot(plot_x_S + bias, term4_r(indS), 'LineWidth', 1.5, 'Color', 'black', 'LineStyle', '--');
p61_S = plot(plot_x_S + bias, term5_h(indS), 'LineWidth', 1.5, 'Color', 'green');
p62_S = plot(plot_x_S + bias, term5_r(indS), 'LineWidth', 1.5, 'Color', 'green', 'LineStyle', '--');

p11_N = plot(plot_x_N - bias, error_h(indN), 'LineWidth', 1.5, 'Color', colors_2(1, :));
p12_N = plot(plot_x_N - bias, error_r(indN), 'LineWidth', 1.5, 'Color', colors_2(1, :), 'LineStyle', '--');
p31_N = plot(plot_x_N - bias, term2_h(indN), 'LineWidth', 1.5, 'Color', 'cyan');
p32_N = plot(plot_x_N - bias, term2_r(indN), 'LineWidth', 1.5, 'Color', 'cyan', 'LineStyle', '--');
p41_N = plot(plot_x_N - bias, term3_h(indN), 'LineWidth', 1.5, 'Color', 'red');
p42_N = plot(plot_x_N - bias, term3_r(indN), 'LineWidth', 1.5, 'Color', 'red', 'LineStyle', '--');
p51_N = plot(plot_x_N - bias, term4_h(indN), 'LineWidth', 1.5, 'Color', 'black');
p52_N = plot(plot_x_N - bias, term4_r(indN), 'LineWidth', 1.5, 'Color', 'black', 'LineStyle', '--');
p61_N = plot(plot_x_N - bias, term5_h(indN), 'LineWidth', 1.5, 'Color', 'green');
p62_N = plot(plot_x_N - bias, term5_r(indN), 'LineWidth', 1.5, 'Color', 'green', 'LineStyle', '--');
%}
p11_S = plot(plot_x_S + bias, error_h(indS), 'LineWidth', 1.5, 'Color', colors_2(1, :));
p12_S = plot(plot_x_S + bias, error_r(indS), 'LineWidth', 1.5, 'Color', colors_2(1, :), 'LineStyle', '--');
p31_S = plot(plot_x_S + bias, term2_h(indS), 'LineWidth', 1.5, 'Color', colors_2(2, :));
p32_S = plot(plot_x_S + bias, term2_r(indS), 'LineWidth', 1.5, 'Color', colors_2(2, :), 'LineStyle', '--');
p41_S = plot(plot_x_S + bias, term3_h(indS), 'LineWidth', 1.5, 'Color', colors_2(3, :));
p42_S = plot(plot_x_S + bias, term3_r(indS), 'LineWidth', 1.5, 'Color', colors_2(3, :), 'LineStyle', '--');
p51_S = plot(plot_x_S + bias, term4_h(indS), 'LineWidth', 1.5, 'Color', colors_2(4, :));
p52_S = plot(plot_x_S + bias, term4_r(indS), 'LineWidth', 1.5, 'Color', colors_2(4, :), 'LineStyle', '--');
p61_S = plot(plot_x_S + bias, term5_h(indS), 'LineWidth', 1.5, 'Color', colors_2(5, :));
p62_S = plot(plot_x_S + bias, term5_r(indS), 'LineWidth', 1.5, 'Color', colors_2(5, :), 'LineStyle', '--');

p11_N = plot(plot_x_N - bias, error_h(indN), 'LineWidth', 1.5, 'Color', colors_2(1, :));
p12_N = plot(plot_x_N - bias, error_r(indN), 'LineWidth', 1.5, 'Color', colors_2(1, :), 'LineStyle', '--');
p31_N = plot(plot_x_N - bias, term2_h(indN), 'LineWidth', 1.5, 'Color', colors_2(2, :));
p32_N = plot(plot_x_N - bias, term2_r(indN), 'LineWidth', 1.5, 'Color', colors_2(2, :), 'LineStyle', '--');
p41_N = plot(plot_x_N - bias, term3_h(indN), 'LineWidth', 1.5, 'Color', colors_2(3, :));
p42_N = plot(plot_x_N - bias, term3_r(indN), 'LineWidth', 1.5, 'Color', colors_2(3, :), 'LineStyle', '--');
p51_N = plot(plot_x_N - bias, term4_h(indN), 'LineWidth', 1.5, 'Color', colors_2(4, :));
p52_N = plot(plot_x_N - bias, term4_r(indN), 'LineWidth', 1.5, 'Color', colors_2(4, :), 'LineStyle', '--');
p61_N = plot(plot_x_N - bias, term5_h(indN), 'LineWidth', 1.5, 'Color', colors_2(5, :));
p62_N = plot(plot_x_N - bias, term5_r(indN), 'LineWidth', 1.5, 'Color', colors_2(5, :), 'LineStyle', '--');

clear('alpha')

ind2 = starts(1) : ends(1);
fill([plot_x(ind2)' + bias, fliplr(plot_x(ind2)' + bias)], [error_h(ind2)' fliplr(error_r(ind2)')], colors_2(1, :), 'edgealpha', 0);
fill([plot_x(ind2)' + bias, fliplr(plot_x(ind2)' + bias)], [term2_h(ind2)' fliplr(term2_r(ind2)')], colors_2(2, :), 'edgealpha', 0);
fill([plot_x(ind2)' + bias, fliplr(plot_x(ind2)' + bias)], [term3_h(ind2)' fliplr(term3_r(ind2)')], colors_2(3, :), 'edgealpha', 0);
fill([plot_x(ind2)' + bias, fliplr(plot_x(ind2)' + bias)], [term4_h(ind2)' fliplr(term4_r(ind2)')], colors_2(4, :), 'edgealpha', 0);
fill([plot_x(ind2)' + bias, fliplr(plot_x(ind2)' + bias)], [term5_h(ind2)' fliplr(term5_r(ind2)')], colors_2(5, :), 'edgealpha', 0);

ind2 = starts(2) : ends(2);
fill([plot_x(ind2)' - bias, fliplr(plot_x(ind2)' - bias)], [error_h(ind2)' fliplr(error_r(ind2)')], colors_2(1, :), 'edgealpha', 0);
fill([plot_x(ind2)' - bias, fliplr(plot_x(ind2)' - bias)], [term2_h(ind2)' fliplr(term2_r(ind2)')], colors_2(2, :), 'edgealpha', 0);
fill([plot_x(ind2)' - bias, fliplr(plot_x(ind2)' - bias)], [term3_h(ind2)' fliplr(term3_r(ind2)')], colors_2(3, :), 'edgealpha', 0);
fill([plot_x(ind2)' - bias, fliplr(plot_x(ind2)' - bias)], [term4_h(ind2)' fliplr(term4_r(ind2)')], colors_2(4, :), 'edgealpha', 0); 
fill([plot_x(ind2)' - bias, fliplr(plot_x(ind2)' - bias)], [term5_h(ind2)' fliplr(term5_r(ind2)')], colors_2(5, :), 'edgealpha', 0); 
    
alpha(.15);

lgd = legend( ...
[p51_S, p41_S, p11_S, p31_S, p61_S], ...
{'$\overline{\nabla^2(\sigma\omega_{QG})}$', ...
 '$f_0^2\overline{\partial_p^2\omega_{QG}}$', ...
 'numerical error' , ...
 '$-\overline{Adv}$', ...
 '$\frac{\kappa}{p}\overline{\nabla^2J}$'}, ...
'location', 'best', 'Interpreter', 'latex');
lgd.Position = [0.84, 0.185, 0.11, 0.7];
set(lgd,'Box','off')

if ~isempty(strfind(pwd_str, 'CESM'))
    if OCEAN_DESERT
        axis([-33 + bias, 33 - bias, -1.5e-16, 1.5e-16]);
    else
        axis([-73 + bias, 73 - bias, -2.5e-16, 2.5e-16]);
    end
else
    axis([-73 + bias, 73 - bias, -0.6e-16, 0.6e-16]);
    yticks([-0.5e-16, -0.25e-16, 0, 0.25e-16, 0.5e-16])
end
set(gca, 'TickDir', 'out');

ax = gca;
set(ax, 'Position', [0.09, 0.17, 0.70, 0.73]);
if ~OCEAN_DESERT
    break_axis('handle', ax, 'position', 0.0, 'axis', 'x', 'length', 0.40)
end

xlabel('Latitude', 'interpreter', 'latex')
ylabel('m kg$^{-1}$s$^{-1}$', 'Interpreter', 'latex');
if OCEAN_DESERT
    xticks([-70 + bias, -50 + bias, -30 + bias, -10 + bias, 10 - bias, 30 - bias, 50 - bias, 70 - bias])
    xticklabels({'$70^\circ$S', '$50^\circ$S', '$30^\circ$S', '$10^\circ$S', '$10^\circ$N', '$30^\circ$N', '$50^\circ$N', '$70^\circ$N'})
else
    xticks([-70 + bias, -50 + bias, -30 + bias, 30 - bias, 50 - bias, 70 - bias])
    xticklabels({'$70^\circ$S', '$50^\circ$S', '$30^\circ$S', '$30^\circ$N', '$50^\circ$N', '$70^\circ$N'})
end
set(gca, 'TickLabelInterpreter', 'latex')

hold off;
%if EPS
%    print(gcf, [plot_path_ensemble, 'terms_', plot_level_name], '-depsc', '-r0', '-painters')
%else
    print(gcf, [plot_path_ensemble, 'terms_', plot_level_name], '-dpng', '-r300')
%end

