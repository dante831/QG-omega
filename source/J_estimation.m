
% one-to-one plot of the approximation of diabatic heating using latent heating

figure('pos', [10, 10, 300, 270])
ind_lat = abs(lat) >= 30;
temp_1 = J_r(ind_lat, :);
temp_2 = - omega_r(ind_lat, :) .* sigma_star_r(ind_lat, :) * plot_level / kappa;
s1 = scatter(temp_2(:), temp_1(:), 1, 'o', 'filled', 'MarkerFaceColor', colors(2, :));
hold on
temp_1 = J_h(ind_lat, :);
temp_2 = - omega_h(ind_lat, :) .* sigma_star_h(ind_lat, :) * plot_level / kappa;
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



