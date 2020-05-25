
pwd_str = pwd;
if ~isempty(strfind(pwd_str, 'GFDL')) || ...
   ~isempty(strfind(pwd_str, 'daily'))
    EPS = false;
else
    EPS = true;    
end


% read in global temperature changes
if strfind(pwd_str, 'CESM')
    load('/disk7/ziweili/CESM_LENS/exp/global_temperature/global_avg_T.mat');
elseif strfind(pwd_str, 'GFDL')
    load('/disk7/ziweili/test1_GFDL/exp/global_temperature/global_avg_T.mat');
end
delta_T = mean(T_avg_r) - mean(T_avg_h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dry decomposition, zonal average
figsize_decomp = [10, 10, 500, 250];
left_margin     = 0.06;
right_margin    = 0.22;
down_margin     = 0.13;
up_margin       = 0.08;
distance        = 0.10;
width           = 1 - left_margin - right_margin;
height          = (1 - down_margin - up_margin - distance) / 2;

position_1 = [left_margin, down_margin + height + distance, width, height];
position_2 = [left_margin, down_margin                    , width, height];

omega_QG_h_mean = nanmean(omega_QG_h_rec(:, :, plot_ind) .* NaN_matrix_m, 2);

fig = figure('pos', figsize_decomp);

lat_1 = -30;
lat_2 = 30;
bias = 25;
if exist('OCEAN_DESERT') && OCEAN_DESERT
    lat_1 = -5;
    lat_2 = 5;
    bias = 0;
end
indN = plot_x > lat_2;
indS = plot_x < lat_1;
plot_x_N = plot_x(indN);
plot_x_S = plot_x(indS);
colors = get(gca,'colororder');
colors(7, 1) = colors(7, 1) + 0.2;
colors(7, 3) = colors(7, 3) - 0.05;

ax1 = subplot(2, 1, 1);
plot(-90:90, zeros(size(-90:90)), 'LineWidth', 0.5, 'Color', [0.3, 0.3, 0.3], 'Linestyle', '-.')
hold on;

pwd_str = pwd;
if ~isempty(strfind(pwd_str, 'GFDL')) || OCEAN_DESERT
    f = @(x) x;
else
    f = @one_two_one; % define the smoothing function
end

p0_d = plot(ax1, plot_x_S + bias, f(nanmean(d_omega             (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2)) ./ omega_h_mean_d(indS));
p1_d = plot(ax1, plot_x_S + bias, f(nanmean(d_omega_QG          (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p2_d = plot(ax1, plot_x_S + bias, f(nanmean(d_sigma             (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p3_d = plot(ax1, plot_x_S + bias, f(nanmean(d_J                 (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p4_d = plot(ax1, plot_x_S + bias, f(nanmean(d_k2_l2             (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p5_d = plot(ax1, plot_x_S + bias, f(nanmean(d_m2                (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p6_d = plot(ax1, plot_x_S + bias, f(nanmean(d_Adv               (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p7_d = plot(ax1, plot_x_S + bias, f(nanmean(d_dtheta_dp_ma_omega(indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));
p8_d = plot(ax1, plot_x_S + bias, f(nanmean(d_rec               (indS, :, plot_ind) .* NaN_matrix_d(indS, :), 2)) ./ omega_QG_h_mean_d(indS));

set(p0_d, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', '--');
set(p1_d, 'LineWidth', 1.5, 'Color', 'black');
set(p2_d, 'LineWidth', 0.8, 'Color', colors(1, :), 'Marker', 's', 'MarkerSize', 2);
set(p3_d, 'LineWidth', 1.5, 'Color', colors(7, :));
set(p4_d, 'LineWidth', 1.5, 'Color', colors(3, :));
set(p5_d, 'LineWidth', 1.5, 'Color', colors(5, :));
set(p6_d, 'LineWidth', 1.5, 'Color', colors(6, :));
set(p7_d, 'LineWidth', 1.5, 'Color', colors(7, :), 'Linestyle', '-.');
set(p8_d, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', ':');

p0_d = plot(ax1, plot_x_N - bias, f(nanmean(d_omega             (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2)) ./ omega_h_mean_d(indN));
p1_d = plot(ax1, plot_x_N - bias, f(nanmean(d_omega_QG          (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));
p2_d = plot(ax1, plot_x_N - bias, f(nanmean(d_sigma             (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));
p3_d = plot(ax1, plot_x_N - bias, f(nanmean(d_J                 (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));
p4_d = plot(ax1, plot_x_N - bias, f(nanmean(d_k2_l2             (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));
p5_d = plot(ax1, plot_x_N - bias, f(nanmean(d_m2                (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));
p6_d = plot(ax1, plot_x_N - bias, f(nanmean(d_Adv               (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));
p7_d = plot(ax1, plot_x_N - bias, f(nanmean(d_dtheta_dp_ma_omega(indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));
p8_d = plot(ax1, plot_x_N - bias, f(nanmean(d_rec               (indN, :, plot_ind) .* NaN_matrix_d(indN, :), 2)) ./ omega_QG_h_mean_d(indN));

set(p0_d, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', '--');
set(p1_d, 'LineWidth', 1.5, 'Color', 'black');
set(p2_d, 'LineWidth', 0.8, 'Color', colors(1, :), 'Marker', 's', 'MarkerSize', 2);
set(p3_d, 'LineWidth', 1.5, 'Color', colors(7, :));
set(p4_d, 'LineWidth', 1.5, 'Color', colors(3, :));
set(p5_d, 'LineWidth', 1.5, 'Color', colors(5, :));
set(p6_d, 'LineWidth', 1.5, 'Color', colors(6, :));
set(p7_d, 'LineWidth', 1.5, 'Color', colors(7, :), 'Linestyle', '-.');
set(p8_d, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', ':');

title(['(a) Dry Decomposition'], 'FontWeight', 'normal', 'interpreter', 'latex')
ylabel('\%K$^{-1}$', 'Interpreter', 'latex', 'FontSize', 9);
ylab1 = get(ax1, 'ylabel');
set(ylab1, 'rotation', 0)
set(ylab1,'Units','normalized');
ylab1.Position(1) = ylab1.Position(1) + 0.06;
ylab1.Position(2) = ylab1.Position(2) / height * (2.1 * height);


if strfind(pwd_str, 'GFDL')
    axis([-73 + bias 73 - bias -0.40 0.65]);
    yticks([-0.15:0.05:0.15]*delta_T);
    yticklabels({'-15', '-10', '-5', '0', '5', '10', '15'})
else
    axis([-73 + bias 73 - bias -0.20 0.30]);
    yticks([-0.06:0.03:0.06]*delta_T);
    yticklabels({'-6', '-3', '0', '3', '6'})
end

if exist('OCEAN_ONLY') && OCEAN_ONLY
    axis([-33 + bias 33 - bias -0.6 0.6]);
    yticks([-0.16:0.08:0.16]*delta_T);
    yticklabels({'-16', '-8', '0', '8', '16'})
end

xticks([-70.0 + bias, -50 + bias, lat_1 + bias, ...
        lat_2 - bias,  50 - bias, 70.0 - bias])
xticklabels([])

set([ax1], 'Position', position_1);
set(ax1, 'TickLabelInterpreter', 'latex')
break_axis('handle', ax1, 'position', 0.0, 'axis', 'x', 'length', 0.20)
set(ax1,'TickDir','out');
%ax1.XTicklabels = [];
box off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% moist decomposition, zonal average
ax2 = subplot(2, 1, 2);

d_rec_m_full_paper    = d_sigma_m + d_k2_l2_m + d_m2_m + d_Adv_m + d_J_res;

plot(-90:90, zeros(size(-90:90)), 'LineWidth', 0.5, 'Color', [0.3, 0.3, 0.3], 'Linestyle', '-.')
hold on;

p1_m   = plot(plot_x_S + bias, f(nanmean(d_omega_QG (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS));
p2_m   = plot(plot_x_S + bias, f(nanmean(d_Adv_m    (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS));
p3_1_m = plot(plot_x_S + bias, f(nanmean(d_sigma_m  (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS));
p3_2_m = plot(plot_x_S + bias, f(nanmean(d_k2_l2_m  (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS));
p3_3_m = plot(plot_x_S + bias, f(nanmean(d_m2_m     (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS));
p4_1_m = plot(plot_x_S + bias, f(nanmean(d_J_res    (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS));
%p4_2_m = plot(plot_x_S + bias, f(nanmean(d_excessive_2(indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS), 'o');
%p5_m   = plot(plot_x_S + bias, f(nanmean(d_epsilon    (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS));
p6_m   = plot(plot_x_S + bias, f(nanmean(d_rec_m_full_paper(indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS));

set(p1_m, 'LineWidth', 1.5, 'Color', 'black');
set(p2_m, 'LineWidth', 1.5, 'Color', colors(6, :));
set(p3_1_m, 'LineWidth', 0.8, 'Color', colors(1, :), 'Marker', 's', 'MarkerSize', 2);
set(p3_2_m, 'LineWidth', 1.5, 'Color', colors(3, :));
set(p3_3_m, 'LineWidth', 1.5, 'Color', colors(5, :));
set(p4_1_m, 'LineWidth', 1.5, 'Color', colors(7, :), 'Linestyle', '--', 'MarkerSize', 2);
%set(p4_2_m, 'LineWidth', 0.8, 'Color', colors(2, :), 'Linestyle', '--', 'MarkerSize', 2);
%set(p5_m, 'LineWidth', 1.5, 'Color', colors(7, :));
set(p6_m, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', ':');

p1_m   = plot(plot_x_N - bias, f(nanmean(d_omega_QG(indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN));
p2_m   = plot(plot_x_N - bias, f(nanmean(d_Adv_m   (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN));
p3_1_m = plot(plot_x_N - bias, f(nanmean(d_sigma_m (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN));
p3_2_m = plot(plot_x_N - bias, f(nanmean(d_k2_l2_m (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN));
p3_3_m = plot(plot_x_N - bias, f(nanmean(d_m2_m    (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN));
p4_1_m = plot(plot_x_N - bias, f(nanmean(d_J_res   (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN));
%p4_2_m = plot(plot_x_N - bias, f(nanmean(d_excessive_2(indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN), 'o');
%p5_m   = plot(plot_x_N - bias, f(nanmean(d_epsilon    (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN));
p6_m   = plot(plot_x_N - bias, f(nanmean(d_rec_m_full_paper(indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN));

set(p1_m, 'LineWidth', 1.5, 'Color', 'black');
set(p2_m, 'LineWidth', 1.5, 'Color', colors(6, :));
set(p3_1_m, 'LineWidth', 0.8, 'Color', colors(1, :), 'Marker', 's', 'MarkerSize', 2);
set(p3_2_m, 'LineWidth', 1.5, 'Color', colors(3, :));
set(p3_3_m, 'LineWidth', 1.5, 'Color', colors(5, :));
set(p4_1_m, 'LineWidth', 1.5, 'Color', colors(7, :), 'Linestyle', '--');
%set(p4_2_m, 'LineWidth', 0.8, 'Color', colors(2, :), 'Linestyle', '--', 'MarkerSize', 2);
%set(p5_m, 'LineWidth', 1.5, 'Color', colors(7, :));
set(p6_m, 'LineWidth', 1.5, 'Color', 'black', 'Linestyle', ':');



if strfind(pwd_str, 'GFDL')
    axis([-73 + bias 73 - bias -0.25 0.38]);
    yticks([-0.075:0.025:0.075]*delta_T);
    yticklabels({'-7.5', '-5', '-2.5', '0', '2.5', '5', '7.5'})
else
    axis([-73 + bias 73 - bias -0.09 0.15]);
    yticks([-0.06:0.015:0.06]*delta_T);
    yticklabels({'-6', '-4.5', '-3', '-1.5', '0', '1.5', '3', '4.5', '6'})
end

if exist('OCEAN_ONLY') && OCEAN_ONLY
    axis([-33 + bias 33 - bias -1.2 1.2]);
    yticks([-0.24:0.12:0.24]*delta_T);
    yticklabels({'-24', '-12', '0', '12', '24'})
end 

xticks([-70.0 + bias, -50 + bias, lat_1 + bias, ...
        lat_2 - bias, 50 - bias, 70.0 - bias]);
xticklabels({'$70^\circ$S', '$50^\circ$S', ['$', num2str(abs(lat_1)), '^\circ$S'], ...
             ['$', num2str(lat_2), '^\circ$N'], '$50^\circ$N', '$70^\circ$N'})
if exist('OCEAN_ONLY') && OCEAN_ONLY
    xticks([-30 + bias, lat_1 + bias, lat_2 - bias, 30 - bias]);
    xticklabels({'$30^\circ$S', ['$', num2str(abs(lat_1)), '^\circ$S'], ...
    ['$', num2str(lat_2), '^\circ$N'], '$30^\circ$N'})
end
set(ax2, 'TickLabelInterpreter', 'latex')

    
lgd = legend( ...
[p7_d, p3_d, p0_d, p1_m, p6_m, p3_3_m, p2_m, p3_2_m, p3_1_m, p4_1_m], ...
{'$-\frac{p}{\kappa}\overline{\sigma^*\omega}$', ...
 '$\overline{J}$'...
 ...%'${c_p{\frac{T}{\theta}\frac{d\theta}{dp}\big|_{\theta*}\omega}}$' , ...
 ...'(a) $\frac{p}{\kappa}\sigma^*$; (b) $\overline{Res}$', ...
 '$\overline{\omega}$' , ...
 '$\overline{\omega}_{QG}$' , ...
 '${\omega_{QG-sum}}$', ...
 '$m$', ...
 '$\overline{Adv}$' , ...
 '$k,k_J$', ...
 '(a) $\hat{\sigma}$, (b) $\hat{\sigma}_m$', ...
 '$\overline{J}_{res}$'}, ...
 ...'$J''_{\frac{d\theta}{dp}\big|_{\theta*}}$' , ...
 ...'$J''_{\omega-\omega_{QG}}$'},  ...
 'location', 'best', ...
 'interpreter', 'latex');
lgd.Position = [left_margin + width + 0.07, 0.10, 0.08, 0.80];
set(lgd,'Box','off')

%{
lgd = legend( ...
[p7_d, p3_d, p0_d, p1_m, p6_m, p3_3_m, p2_m, p3_2_m, p3_1_m, p4_1_m], ...
{'$\frac{p}{\kappa}\sigma^*$', ...
 '${J}$'...
 ...%'${c_p{\frac{T}{\theta}\frac{d\theta}{dp}\big|_{\theta*}\omega}}$' , ...
 ...'(a) $\frac{p}{\kappa}\sigma^*$; (b) $\overline{Res}$', ...
 '${\omega}$' , ...
 '${\omega}_{QG}$' , ...
 '${\omega_{QG-rec}}$', ...
 '$m$', ...
 '${Adv}$' , ...
 '$k,k_J$', ...
 '(a) ${\sigma}$; (b) ${\sigma}_m$', ...
 '${J}_{Res}$'}, ...
 ...'$J''_{\frac{d\theta}{dp}\big|_{\theta*}}$' , ...
 ...'$J''_{\omega-\omega_{QG}}$'},  ...
 'location', 'best', ...
 'interpreter', 'latex');
lgd.Position = [left_margin + width + 0.07, 0.10, 0.08, 0.80];
set(lgd,'Box','off')
%}

title(['(b) Moist Decomposition'], 'FontWeight', 'normal', 'Interpreter', 'latex');
xlabel('Latitude', 'Interpreter', 'latex');
ylabel('\%K$^{-1}$', 'Interpreter', 'latex', 'FontSize', 9);
ylab2 = get(ax2, 'ylabel');
set(ylab2, 'rotation', 0)
set(ylab2,'Units','normalized');
ylab2.Position(1) = ylab2.Position(1) + 0.05;
ylab2.Position(2) = ylab2.Position(2) / height * (2.1 * height);

set([ax2], 'Position', position_2);
set(ax2, 'TickLabelInterpreter', 'latex')
break_axis('handle', ax2, 'position', 0.0, 'axis', 'x', 'length', 0.20)
set(ax2,'TickDir','out');
hold off;
box off;


figname = [plot_path_ensemble, 'omega_decomposition_moist_paper_', plot_level_name];

if EPS
    print(fig, figname, '-depsc', '-r0', '-painters')
else
    saveas(fig, figname, 'png')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% decomposition of the higher-order heating
fig2 = figure('pos', [10, 10, 500, 150]);

left_margin     = 0.08;
right_margin    = 0.20;
down_margin     = 0.13;
up_margin       = 0.08;
distance        = 0.10;
width           = 1 - left_margin - right_margin;
height          = 1 - down_margin - up_margin - distance;

position = [left_margin, down_margin, width, height];

ax = subplot(1, 1, 1);

plot(-90:90, zeros(size(-90:90)), 'LineWidth', 0.5, 'Color', [0.3, 0.3, 0.3], 'Linestyle', '-.')
hold on;

p1 = plot(plot_x_S + bias, f(nanmean(d_J_res      (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS), 's');
p2 = plot(plot_x_S + bias, f(nanmean(d_excessive_1(indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS));
p3 = plot(plot_x_S + bias, f(nanmean(d_excessive_2(indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS));
p4 = plot(plot_x_S + bias, f(nanmean(d_epsilon    (indS, :, plot_ind) .* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS));
%p5 = plot(plot_x_S + bias, f(nanmean((d_excessive_1(indS, :, plot_ind) + ...
%                                      d_excessive_2(indS, :, plot_ind) + ...
%                                      d_epsilon    (indS, :, plot_ind)).* NaN_matrix_m(indS, :), 2)) ./ omega_QG_h_mean_m(indS));


set(p1, 'LineWidth', 0.8, 'Color', colors(4, :), 'Linestyle', '--', 'MarkerSize', 2);
set(p2, 'LineWidth', 1.5, 'Color', colors(1, :));
set(p3, 'LineWidth', 1.5, 'Color', colors(3, :));
set(p4, 'LineWidth', 1.5, 'Color', colors(2, :));

p1 = plot(plot_x_N - bias, f(nanmean(d_J_res      (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN), 's');
p2 = plot(plot_x_N - bias, f(nanmean(d_excessive_1(indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN));
p3 = plot(plot_x_N - bias, f(nanmean(d_excessive_2(indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN));
p4 = plot(plot_x_N - bias, f(nanmean(d_epsilon    (indN, :, plot_ind) .* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN));
%p5 = plot(plot_x_N - bias, f(nanmean((d_excessive_1(indN, :, plot_ind) + ...
%                                      d_excessive_2(indN, :, plot_ind) + ...
%                                      d_epsilon    (indN, :, plot_ind)).* NaN_matrix_m(indN, :), 2)) ./ omega_QG_h_mean_m(indN));


set(p1, 'LineWidth', 0.8, 'Color', colors(4, :), 'Linestyle', '--', 'MarkerSize', 2);
set(p2, 'LineWidth', 1.5, 'Color', colors(1, :));
set(p3, 'LineWidth', 1.5, 'Color', colors(3, :));
set(p4, 'LineWidth', 1.5, 'Color', colors(2, :));

if strfind(pwd_str, 'GFDL')
    axis([-73 + bias 73 - bias -0.20 0.35]);
    yticks([-0.075:0.025:0.075]*delta_T);
    yticklabels({'-7.5', '-5', '-2.5', '0', '2.5', '5', '7.5'})
else
    axis([-73 + bias 73 - bias -0.12 0.15]);
    yticks([-0.06:0.015:0.06]*delta_T);
    yticklabels({'-6', '-4.5', '-3', '-1.5', '0', '1.5', '3', '4.5', '6'})
end
xticks([-70.0 + bias, -50 + bias, lat_1 + bias, ...
        lat_2 - bias, 50 - bias, 70.0 - bias]);
xticklabels({'$70^\circ$S', '$50^\circ$S', ['$', num2str(abs(lat_1)), '^\circ$S'], ...
             ['$', num2str(lat_2), '^\circ$N'], '$50^\circ$N', '$70^\circ$N'})
set(ax, 'TickLabelInterpreter', 'latex')


lgd = legend( ...
[p1, p2, p3, p4], ...
{'$J''$', ...
 '$J''_{\frac{T}{\theta}\frac{d\theta}{dp}\big|_{\theta*}}$' , ...
 '$J''_{k_J^2\omega-k_*^2\omega_{QG}}$', ...
 '$J''_{\epsilon}$'}, ...
 'location', 'best', ...
 'interpreter', 'latex');
lgd.Position = [left_margin + width + 0.07, 0.15, 0.08, 0.60];
set(lgd,'Box','off')

xlabel('Latitude', 'Interpreter', 'latex');
ylabel('\%/K$^{-1}$', 'Interpreter', 'latex', 'FontSize', 9);
ylab = get(ax,'ylabel');
set(ylab, 'rotation', 0)
set(ylab,'Units','normalized');
ylab.Position(1) = ylab.Position(1) + 0.02;
ylab.Position(2) = ylab.Position(2) / height * (4 * height + 2.2 * distance);

set([ax], 'Position', position);
set(ax, 'TickLabelInterpreter', 'latex')
break_axis('handle', ax, 'position', 0.0, 'axis', 'x', 'length', 0.20)
set(ax,'TickDir','out');
hold off;
box off;

figname = [plot_path_ensemble, 'unbalanced_decomposition_moist_paper_', plot_level_name];
if EPS
    print(fig2, figname, '-depsc', '-r0', '-painters')
else
    saveas(fig2, figname, 'png')
end




