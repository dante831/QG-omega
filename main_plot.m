
addpath('./source');

pwd_str = pwd;

% invert the omega equation of an example extreme-precipitation event, and plot figure 1
inversion_example

mat_filenames = {'CESM_6hourly.mat', ...
                 'CESM_daily.mat', ...
                 'GFDL_6hourly.mat', ...
                 'GFDL_6hourly_full_b.mat', ...
                 'GFDL_daily.mat'};

for f_ind = 1 : length(mat_filenames)
    filename = ['./data/', mat_filenames{f_ind}];
    
    % load prepared 3D fields from data files
    load(filename)
        % Each indivdual file contains the following fields. 
        % Suffix '_h' denotes values averaged across events in the 
        % historical climate, and '_r' denotes values for RCP8.5.
        % 'omega_h' and 'omega_r': pressure velocity
        % 'omega_QG_h' and 'omega_QG_r': pressure velocity inverted from teh QG-omega equation
        % 'precip_h' and 'precip_r': precipitation
        % 'Adv_h' and 'Adv_r': advection forcing in the Q-vector form

    % define the threshold of number of events for masking
    if strfind(filename, 'CESM_6hourly.mat')
        num_threshold = 30;
    elseif strfind(filename, 'GFDL_6hourly.mat') | ...
           strfind(filename, 'GFDL_daily.mat') | ...
           strfind(filename, 'GFDL_6hourly_full_b.mat')
        num_threshold = 15;
    elseif strfind(filename, 'CESM_daily.mat')
        num_threshold = 5;
    end

    if ~isempty(strfind(filename, 'GFDL'))
        GFDL = true;
    else
        GFDL = false;
    end

    % define zonal-mean figure size
    figsize_zonal = [10, 10, 600, 200];

    % whether smooth 2D fields or not, depending on the model
    if GFDL
        SMOOTH = false;
    else
        SMOOTH = true;
    end

    % some constants
    Ra = 287.04;        % specific gas constant for dry air, in J/(kg*K)
    R = 6371000.0;      % Earth's radius, in meters
    cpd= 1005.7;        % specific heat capacity under constant pressure, in J/(kg*K)
    kappa = Ra / cpd;   
    plot_level = 50000; % level at question, in Pascals
    plot_level_name = [num2str(plot_level/100.), 'hPa'];

    [X, Y] = meshgrid(lon, lat);

    % plot settings
    set(0,'DefaultFigureVisible','off');
    omega_range = [-1.5, 1.5];
    change_range = [-0.5, 0.5];
    smooth_window_x = 1;
    smooth_window_y = 1;

    plot_path_ensemble = ['plots/', mat_filenames{f_ind}(1:end-4), '/'];
    if ~exist(plot_path_ensemble)
        mkdir(plot_path_ensemble);
    end

    num_event_d = reshape(min([num_event_h(:)'; num_event_r(:)']), size(num_event_h));
    omega_QG_h_rec =  - 1 ./ (k2_h .* sigma_h + m2_h .* F0.^2) .* (Adv_h + C_h);
    omega_QG_r_rec =  - 1 ./ (k2_r .* sigma_r + m2_r .* F0.^2) .* (Adv_r + C_r);

    %% dry decomposition

    d_sigma = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* (sigma_r - sigma_h) .* k2_h .* omega_QG_h_rec;
    d_k2    = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* sigma_h .* (k2_r - k2_h) .* omega_QG_h;
    d_l2    = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* kappa ./ plot_level .* (l2_r - l2_h) .* J_h;
    d_J     = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* kappa ./ plot_level .* k2_h .* (J_r - J_h);
    d_m2    = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* F0.^2 .* (m2_r - m2_h) .* omega_QG_h_rec;
    d_Adv   = - 1./(sigma_h.*k2_h + F0.^2.*m2_h) .* (Adv_r - Adv_h);
    d_dtheta_dp_ma_omega = - 1./(sigma_h .* k2_h + F0.^2.*m2_h) .* ...
            kappa .* cpd ./ plot_level .* k2_h .* (omega_r .* dtheta_dp_ma_r - omega_h .* dtheta_dp_ma_h);
            % note that dtheta_dp_ma_h already has a T/theta factor
    d_rec   = d_sigma + d_k2 + d_l2 + d_J + d_m2 + d_Adv;
    d_omega_QG = omega_QG_r - omega_QG_h;
    d_omega    = omega_r - omega_h;
    d_k2_l2 = d_k2 + d_l2;

    % zonal mean plots

    colors = get(gca,'colororder');
    plot_x = lat_series(lat_indices);

    % obtain NaN_matrix for the dry case. Points that have d_omega_QG, d_sigma, etc. larger than 1000% are not displayed 
    % or enter the following calculation. 
    NaN_matrix_d = get_NaN_matrix(cat(3, d_omega_QG, d_sigma, d_J, d_k2_l2, d_m2, ...
            d_Adv, d_rec, d_dtheta_dp_ma_omega), 10.0);

    % calculate zonal average omega and omega_QG
    omega_QG_h_mean_d = nanmean(omega_QG_h .* NaN_matrix_d, 2);
    omega_QG_r_mean_d = nanmean(omega_QG_r .* NaN_matrix_d, 2);
    omega_h_mean_d    = nanmean(omega_h    .* NaN_matrix_d, 2);
    omega_r_mean_d    = nanmean(omega_r    .* NaN_matrix_d, 2);

    % plot zonally-averaged precipitation
    precip_h_mean_d = nanmean(precip_h .* NaN_matrix_d, 2);
    precip_r_mean_d = nanmean(precip_r .* NaN_matrix_d, 2);

    %% moist decomposition
    k2_star_h = k2_h;
    k2_star_r = k2_r;

    [latent_h, latent_r, excessive_h, excessive_r, epsilon_h, epsilon_r, ...
     sigma_m_h, sigma_m_r, k2_m_h, k2_m_r, denom_h, denom_r, ...
     d_sigma_m, d_Adv_m, d_k2_l2_m, d_m2_m, d_rec_m, d_rec_m_full, d_J_res, ...
     d_excessive_1, d_excessive_2, d_epsilon, J_res_h, J_res_r] = moist_decomposition(...
            k2_h, k2_r, l2_h, l2_r, m2_h, m2_r, sigma_h, sigma_r, Plevels, F0, ...
            dtheta_dp_ma_h, dtheta_dp_ma_r, omega_h, omega_r, omega_QG_h_rec, omega_QG_r_rec, ...
            J_h, J_r, Adv_h, Adv_r, sigma_star_h, sigma_star_r, k2_h, k2_r);

    % moist reconstruction
    omega_QG_h_rec_moist = - (Adv_h + kappa ./ Plevels .* l2_h .* J_res_h) ./ denom_h;
    omega_QG_r_rec_moist = - (Adv_r + kappa ./ Plevels .* l2_r .* J_res_r) ./ denom_r;

    d_excessive     = d_excessive_1 + d_excessive_2;
    d_denom         = d_sigma_m + d_k2_l2_m + d_m2_m;
    d_rec_m_full = d_sigma_m + d_k2_l2_m + d_m2_m + d_Adv_m + d_J_res;
    d_rec_m      = d_sigma_m + d_k2_l2_m + d_m2_m + d_Adv_m;

    % obtain NaN_matrix for the moist case.
    NaN_matrix_m = get_NaN_matrix(cat(3, d_omega_QG, d_Adv_m, d_k2_l2_m, d_m2_m, ...
            d_Adv_m, d_rec_m, d_J_res, d_epsilon), 10.0);

    omega_QG_h_mean_m = nanmean(omega_QG_h .* NaN_matrix_m, 2);
    omega_QG_r_mean_m = nanmean(omega_QG_r .* NaN_matrix_m, 2);
    Omega_QG_h_mean_m = repmat(omega_QG_h_mean_m, 1, size(omega_QG_h, 2));
    omega_h_mean_m    = nanmean(omega_h    .* NaN_matrix_m, 2);
    omega_r_mean_m    = nanmean(omega_r    .* NaN_matrix_m, 2);
    Omega_h_mean_m    = repmat(omega_h_mean_m, 1, size(omega_h, 2));

    % read in global temperature changes
    if GFDL
        load('/disk7/ziweili/test1_GFDL/exp/global_temperature/global_avg_T.mat');
    else
        load('/disk7/ziweili/CESM_LENS/exp/global_temperature/global_avg_T.mat');
    end
    delta_T = mean(T_avg_r) - mean(T_avg_h);

    Omega_QG_h_mean_d = repmat(omega_QG_h_mean_d, 1, size(omega_QG_h, 2));
    Omega_h_mean_d = repmat(omega_h_mean_d, 1, size(omega_h, 2));

    %% plot figure 2, lat-lon map of terms in the QG-omega equation
    figure_2

    %% figures 3, 6, and S2
    omega_vs_omega_QG

    %% plot figure 4, lat-lon map of dry decomposition
    figure_4

    %% plot figure 5, zonal averaged dry and moist decompositions
    figure_5

    %% figure 7, estimation of J term via moist-adiabatic process
    figure_7

    %% figure 8, lat-lon map of moist decomposition
    figure_8

    %% plot the magnitude of zonally-averaged terms
    plot_terms

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

    % averaged changes in the extratropics, in %/K
    spherical_mean = @(NaN_matrix, dS, var, omega_mean, ind_lat, delta_T) ...
                        nanmean(nanmean(dS .* var(ind_lat, :) .* NaN_matrix(ind_lat, :), 2) ...
                                ./ omega_mean(ind_lat)) / ...
                        nanmean(nanmean(dS)) / delta_T;
    d_omega_mean_dry    = spherical_mean(NaN_matrix_d, dS_d, d_omega, omega_h_mean_d, ind_lat, delta_T);
    d_omega_QG_mean_dry = spherical_mean(NaN_matrix_d, dS_d, d_omega_QG, omega_QG_h_mean_d, ind_lat, delta_T);
    % comparing omega and omega_QG
    omega_mean_temp = spherical_mean(NaN_matrix_d, dS_d, omega_h, ones(size(omega_h_mean_d)), ind_lat, 1);
    omega_QG_mean_temp = spherical_mean(NaN_matrix_d, dS_d, omega_QG_h, ones(size(omega_h_mean_d)), ind_lat, 1);
    disp(['omega_QG underestimates omega by ', num2str(1 - omega_QG_mean_temp / omega_mean_temp)]) 
    % percentage of omega_QG's underestimation of omega

    % averaged contributions, dry decomposition
    d_precip     = spherical_mean(NaN_matrix_d, dS_d, precip_r - precip_h, precip_h_mean_d, ...
                        ind_lat, delta_T);
    d_sigma_mean = spherical_mean(NaN_matrix_d, dS_d, d_sigma, omega_QG_h_mean_d, ind_lat, delta_T);
    d_J_mean     = spherical_mean(NaN_matrix_d, dS_d, d_J,     omega_QG_h_mean_d, ind_lat, delta_T);
    d_k2_l2_mean = spherical_mean(NaN_matrix_d, dS_d, d_k2_l2, omega_QG_h_mean_d, ind_lat, delta_T);
    d_m2_mean    = spherical_mean(NaN_matrix_d, dS_d, d_m2,    omega_QG_h_mean_d, ind_lat, delta_T);
    d_Adv_mean   = spherical_mean(NaN_matrix_d, dS_d, d_Adv,   omega_QG_h_mean_d, ind_lat, delta_T);

    d_omega_mean_moist    = spherical_mean(NaN_matrix_m, dS_m, d_omega, omega_h_mean_m, ind_lat, delta_T);
    d_omega_QG_mean_moist = spherical_mean(NaN_matrix_m, dS_m, d_omega_QG, omega_QG_h_mean_m, ind_lat, delta_T);

    % averaged contributions, moist decomposition
    d_sigma_m_mean = spherical_mean(NaN_matrix_m, dS_m, d_sigma_m, omega_QG_h_mean_m, ind_lat, delta_T);
    d_J_res_mean   = spherical_mean(NaN_matrix_m, dS_m, d_J_res,   omega_QG_h_mean_m, ind_lat, delta_T);
    d_k2_l2_m_mean = spherical_mean(NaN_matrix_m, dS_m, d_k2_l2_m, omega_QG_h_mean_m, ind_lat, delta_T);
    d_m2_m_mean    = spherical_mean(NaN_matrix_m, dS_m, d_m2_m,    omega_QG_h_mean_m, ind_lat, delta_T);
    d_Adv_m_mean   = spherical_mean(NaN_matrix_m, dS_m, d_Adv_m,   omega_QG_h_mean_m, ind_lat, delta_T);

    % generation of latex code for a table
    % save values for the table
    matfilename = ['./data/extra_tropical_means_', mat_filenames{f_ind}(1:end-4)];
    save([pwd_str, '/', matfilename], 'd_precip', ...
        'd_omega_mean_dry'  , 'd_omega_QG_mean_dry'  , ...
        'd_omega_mean_moist', 'd_omega_QG_mean_moist', ...
        'd_sigma_mean'  , 'd_J_mean'    , 'd_k2_l2_mean'  , 'd_m2_mean'  , 'd_Adv_mean', ...
        'd_sigma_m_mean', 'd_J_res_mean', 'd_k2_l2_m_mean', 'd_m2_m_mean', 'd_Adv_m_mean');

end

% generate latex code for table in the supplementary information
latex_code_for_table

