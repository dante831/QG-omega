
files = {};

files{1} = './data/extra_tropical_means_CESM_6hourly.mat';
files{2} = './data/extra_tropical_means_GFDL_6hourly.mat';
files{3} = './data/extra_tropical_means_CESM_daily.mat';
files{4} = './data/extra_tropical_means_GFDL_daily.mat';

[d_precip, ...
 d_omega_dry, d_omega_QG_dry, ...
 d_omega_moist, d_omega_QG_moist, ...
 d_sigma  , d_J    , d_k2_l2  , d_m2  , d_Adv, ...
 d_sigma_m, d_J_res, d_k2_l2_m, d_m2_m, d_Adv_m] = deal(zeros(length(files), 1));
[data_d, data_m] = deal(zeros(length(files), 7));


for i = 1 : length(files)
    load(files{i})
    data_d(i, :) = [d_omega_mean_dry, d_omega_QG_mean_dry, d_J_mean, ...
                    d_sigma_mean, d_Adv_mean, d_k2_l2_mean, d_m2_mean];
    data_m(i, :) = [d_omega_mean_moist, d_omega_QG_mean_moist, d_J_res_mean, ...
                    d_sigma_m_mean, d_Adv_m_mean, d_k2_l2_m_mean, d_m2_m_mean];
end

% round to the nearest 0.1%, with 0.05% rounded away from zero, 
% and transform to percentage
round_ = @(x) abs(round(x*100, 1)).*sign(x)/100;
data_d = round_(data_d)*100; 
data_m = round_(data_m)*100;

output_file = 'code.out';
f_out = fopen(output_file, 'w');

var_name_d = {'$\overline{\omega}$', '$\overline{\omega}_{QG}$', ...
              '$\overline{J}$', '$\hat{\sigma}$', '$\overline{Adv}$', '$k, k_J$', '$m$'};
var_name_m = {'$\overline{\omega}$', '$\overline{\omega}_{QG}$', ...
              '$\overline{J}_{res}$', '$\hat{\sigma}_m$', '$\overline{Adv}$', '$k, k_J$', '$m$'};

fprintf(f_out, '%dry decomposition\n');
formatSpec = ['%1$-26s & %2$5.1f & %3$5.1f && %4$10.1f & %5$10.1f \\\\[0.15cm]\n'];
for j = 1 : length(var_name_d)
    fprintf(f_out, formatSpec, var_name_d{j}, data_d(1, j), data_d(2, j), data_d(3, j), data_d(4, j));
end
fprintf(f_out, '\n');
for j = 1 : length(var_name_m)
    fprintf(f_out, formatSpec, var_name_m{j}, data_m(1, j), data_m(2, j), data_m(3, j), data_m(4, j));
end
fclose(f_out);

