plot_path_term = [plot_path_ensemble, 'terms/'];
if ~exist(plot_path_term)
    mkdir(plot_path_term);
end



% term comparisons
% 1. sigma*k2*omega_QG
temp = - sigma_h(:, :, plot_ind) .* k2_h(:, :, plot_ind) .* omega_QG_h_rec(:, :, plot_ind);
plot_title = 'historical term $\nabla^2(\sigma\omega_{QG})$';
plot_filename = 'term_sigmak2omega_h';
climit = [-3e-16, 3e-16];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_term, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true);

temp = - sigma_r(:, :, plot_ind) .* k2_r(:, :, plot_ind) .* omega_QG_r_rec(:, :, plot_ind);
plot_title = 'rcp85 term $\nabla^2(\sigma\omega_{QG})$';
plot_filename = 'term_sigmak2omega_r';
climit = [-3e-16, 3e-16];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_term, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true);

% 2. m2*f0^2*omega_QG
temp = - m2_h(:, :, plot_ind) .* F0(:, :, plot_ind).^2 .* omega_QG_h_rec(:, :, plot_ind);
plot_title = 'historical term $f_0^2\frac{\partial^2\omega_QG}{\partial p^2}$';
plot_filename = 'term_f2m2omega_h';
climit = [-0.5e-16, 0.5e-16];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_term, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true);

temp = - m2_r(:, :, plot_ind) .* F0(:, :, plot_ind).^2 .* omega_QG_r_rec(:, :, plot_ind);
plot_title = 'rcp85 term $f_0^2\frac{\partial^2\omega_QG}{\partial p^2}$';
plot_filename = 'term_f2m2omega_r';
climit = [-0.5e-16, 0.5e-16];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_term, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true);

% 3. Adv
temp = Adv_h(:, :, plot_ind);
plot_title = 'historical term $Adv$';
plot_filename = 'term_Adv_h';
climit = [-0.5e-16, 0.5e-16];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_term, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true);

temp = Adv_r(:, :, plot_ind);
plot_title = 'rcp85 term $Adv$';
plot_filename = 'term_Adv_r';
climit = [-0.5e-16, 0.5e-16];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_term, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true);

temp = (Adv_r(:, :, plot_ind) - Adv_h(:, :, plot_ind)) ./ Adv_h(:, :, plot_ind);
plot_title = 'Fractional change of term $Adv$';
plot_filename = 'term_Adv_change';
%climit = [-0.2e-16, 0.2e-16];
climit = change_range;
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_term, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true);


% 4. sigma*k2*omega_QG - kappa/p*l2*J
temp = - sigma_h(:, :, plot_ind) .* k2_h(:, :, plot_ind) .* omega_QG_h_rec(:, :, plot_ind) - C_h(:, :, plot_ind);
plot_title = 'historical term $\nabla^2(\sigma\omega_{QG} + \frac{\kappa}{p}J)$';
plot_filename = 'term_sigmak2omega-C_h';
climit = [-0.5e-16, 0.5e-16];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_term, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true);

temp = - sigma_r(:, :, plot_ind) .* k2_r(:, :, plot_ind) .* omega_QG_r_rec(:, :, plot_ind) - C_r(:, :, plot_ind);
plot_title = 'rcp85 term $\nabla^2(\sigma\omega_{QG} + \frac{\kappa}{p}J)$';
plot_filename = 'term_sigmak2omega-C_r';
climit = [-0.5e-16, 0.5e-16];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_term, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true);

% 5. C
temp = C_h(:, :, plot_ind);
plot_title = 'historical term $\frac{\kappa}{p}\nabla^2J$';
plot_filename = 'term_C_h';
climit = [-3e-16, 3e-16];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_term, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true);

% 6. change of m2*f0^2*omega_QG
temp = - m2_h(:, :, plot_ind) .* F0(:, :, plot_ind).^2 .* omega_QG_h_rec(:, :, plot_ind) - ...
      (- m2_r(:, :, plot_ind) .* F0(:, :, plot_ind).^2 .* omega_QG_r_rec(:, :, plot_ind));
plot_title = 'change of term $f_0^2\frac{\partial^2\omega_QG}{\partial p^2}$';
plot_filename = 'term_f2m2omega_change';
climit = [-0.1e-16, 0.1e-16];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_term, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, true);




temp = sigma_h(:, :, plot_ind);
plot_title = 'historical term $\sigma$';
plot_filename = 'term_sigma_h';
climit = [0, 5e-6];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_term, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = sigma_r(:, :, plot_ind);
plot_title = 'rcp85 term $\sigma$';
plot_filename = 'term_sigma_r';
climit = [0, 5e-6];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_term, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = - Ra ./ plot_level * dtheta_dp_ma_h(:, :, plot_ind);
plot_title = 'historical term $-\frac{R}{p}\frac{d\theta}{dp}|_{\theta}^*$';
plot_filename = 'term_dtheta_dp_ma_h';
climit = [0, 5e-6];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_term, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = - Ra ./ plot_level * dtheta_dp_ma_r(:, :, plot_ind);
plot_title = 'rcp85 term $-\frac{R}{p}\frac{d\theta}{dp}|_{\theta}^*$';
plot_filename = 'term_dtheta_dp_ma_r';
climit = [0, 5e-6];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_term, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

temp = denom_h(:, :, plot_ind);
plot_title = 'historical denominator $k^2(\sigma + \frac{RT}{p\theta}\frac{d\theta}{dp}|_{\theta^*} + f_0^2m^2)$';
plot_filename = 'term_denom_h';
climit = [0, 8e-17];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_term, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);


temp = denom_r(:, :, plot_ind);
plot_title = 'rcp85 denominator $k^2(\sigma + \frac{RT}{p\theta}\frac{d\theta}{dp}|_{\theta^*} + f_0^2m^2)$';
plot_filename = 'term_denom_r';
climit = [0, 8e-17];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_term, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);



temp = sigma_h(:, :, plot_ind) .* k2_r(:, :, plot_ind);
plot_title = 'historical term $\sigma k^2$';
plot_filename = 'term_sigmak2_h';
climit = [1e-16, 2e-16];
plot_2D_map_shifted(figsize, temp, X, Y, plot_filename, plot_title, plot_path_term, smooth_window_x, smooth_window_y, SMOOTH, climit, ones(size(temp)), 0, false);

