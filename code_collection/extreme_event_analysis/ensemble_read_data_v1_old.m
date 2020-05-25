
clear('data')
data = {};
n_data = length(filenames);

for f = 1 : n_data

    load(filenames{f});

    data{f}.J_l2_h_full = J_l2_h_full;
    data{f}.J_l2_r_full = J_l2_r_full;
    data{f}.J_h_full = J_h_full;
    data{f}.J_r_full = J_r_full;
    data{f}.omega_sigma_h_full = omega_sigma_h_full;
    data{f}.omega_sigma_r_full = omega_sigma_r_full;
    data{f}.omega_k2_sigma_h_full = omega_k2_sigma_h_full;
    data{f}.omega_k2_sigma_r_full = omega_k2_sigma_r_full;
    data{f}.m2_h_full = m2_h_full;
    data{f}.m2_r_full = m2_r_full;
    data{f}.d2_dp2_omega_QG_grid_h_full = d2_dp2_omega_QG_grid_h_full;
    data{f}.d2_dp2_omega_QG_grid_r_full = d2_dp2_omega_QG_grid_r_full;
    data{f}.sigma_h_full = sigma_h_full;
    data{f}.sigma_r_full = sigma_r_full;
    if exist('sigma_star_h')
        data{f}.sigma_star_h = sigma_star_h;
        data{f}.sigma_star_r = sigma_star_r;
    end
    data{f}.Adv_h_full = Adv_h_full;
    data{f}.Adv_r_full = Adv_r_full;
    data{f}.C_h_full = C_h_full;
    data{f}.C_r_full = C_r_full;
    data{f}.dtheta_dp_ma_h_full = dtheta_dp_ma_h_full;
    data{f}.dtheta_dp_ma_r_full = dtheta_dp_ma_r_full;
    data{f}.Le_h = Le_h;
    data{f}.Le_r = Le_r;
    %if exist('dtheta_dp_ma_avg_h_full')
    %    data{f}.dtheta_dp_ma_avg_h_full = dtheta_dp_ma_avg_h_full;
    %    data{f}.dtheta_dp_ma_avg_r_full = dtheta_dp_ma_avg_r_full;
    %end
    if exist('T_adv_h_full')
        data{f}.T_adv_h_full = T_adv_h_full;
        data{f}.T_adv_r_full = T_adv_r_full;
    end

    data{f}.omega_QG_h_full      = omega_QG_h_full;
    data{f}.omega_QG_r_full      = omega_QG_r_full;
    data{f}.omega_h_full         = omega_h_full;
    data{f}.omega_r_full         = omega_r_full;
    data{f}.num_event_h          = repmat(num_event_historical.num_event, [1, 1, size(omega_h_full, 3)]);
    data{f}.num_event_r          = repmat(num_event_rcp85.num_event, [1, 1, size(omega_h_full, 3)]);

    % remove omega_QG > 0
    %temp_ind_h = omega_QG_h_full(:, :, plevels == 50000) >= 0;
    data{f}.num_event_h(omega_QG_h_full >= 0) = 0;
    data{f}.num_event_r(omega_QG_r_full >= 0) = 0;
    % remove J < 0
    %data{f}.num_event_h(J_h_full <= 0) = 0;
    %data{f}.num_event_r(J_r_full <= 0) = 0;
    % remove k2 < 0
    data{f}.num_event_h(k2_h_full <= 0) = 0;
    data{f}.num_event_r(k2_r_full <= 0) = 0;
    % remove l2 < 0
    %data{f}.num_event_h(J_l2_h_full <= 0) = 0;
    %data{f}.num_event_r(J_l2_r_full <= 0) = 0;
    % remove m2 < 0
    %data{f}.num_event_h(d2_dp2_omega_QG_grid_h_full <= 0) = 0;
    %data{f}.num_event_r(d2_dp2_omega_QG_grid_r_full <= 0) = 0;
    
    % set all nan to zero
    names = fieldnames(data{f});
    
    names(ismember(names, 'Le_h')) = [];
    names(ismember(names, 'Le_r')) = [];
    for n = 1 : length(names)
        temp = getfield(data{f}, names{n});
        nan_ind = isnan(temp);
        temp(nan_ind) = 0;
        if strfind(names{n}, '_h')
            data{f}.num_event_h(nan_ind) = 0;
        elseif strfind(names{n}, '_r')
            data{f}.num_event_r(nan_ind) = 0;
        end
        data{f}.(names{n}) = temp;
    end
    for n = 1 : length(names)
        temp = getfield(data{f}, names{n});
        zero_ind_h = data{f}.num_event_h == 0;
        zero_ind_r = data{f}.num_event_r == 0;
        if strfind(names{n}, '_h')
            temp(zero_ind_h) = 0;
        elseif strfind(names{n}, '_r')
            temp(zero_ind_r) = 0;
        end
        data{f}.(names{n}) = temp;
    end
end

[sigma_omega_h_temp, sigma_omega_r_temp, ...
 Adv_h_temp, Adv_r_temp, C_h_temp, C_r_temp, ...
 dtheta_dp_ma_h_temp, dtheta_dp_ma_r_temp, ...
 sigma_star_omega_QG_h_temp, sigma_star_omega_QG_r_temp, ...
 dtheta_dp_ma_omega_h_temp, dtheta_dp_ma_omega_r_temp, ...
 dtheta_dp_ma_avg_h_temp, dtheta_dp_ma_avg_r_temp, ...
 omega_QG_h_temp, omega_QG_r_temp, omega_h_temp, omega_r_temp, ...
 omega_k2_sigma_h_temp, omega_sigma_h_temp, ...
 omega_k2_sigma_r_temp, omega_sigma_r_temp, ...
 J_l2_h_temp, J_h_temp, ...
 J_l2_r_temp, J_r_temp, ...
 ...%d2_dp2_omega_QG_grid_h_temp, d2_dp2_omega_QG_grid_r_temp, ...
 sigma_omega_QG_h_temp, sigma_omega_QG_r_temp, ...
 omega_QG_m2_h_temp, omega_QG_m2_r_temp, ...
 num_event_h, num_event_r, ...
 omega_QG_h_temp, ...
 T_adv_h_temp, T_adv_r_temp] = deal(zeros(size(data{1}.J_h_full)));
[Le_h_temp, Le_r_temp] = deal(zeros(size(data{1}.J_h_full(:, :, 10))));
for i = 1 : n_data
    % deal with AB
    Adv_h_temp               = Adv_h_temp                 + ...
        data{i}.Adv_h_full                .* data{i}.num_event_h;
    Adv_r_temp               = Adv_r_temp                 + ...
        data{i}.Adv_r_full                .* data{i}.num_event_r;
    % deal with C
    C_h_temp                = C_h_temp                  + ...
        data{i}.C_h_full                 .* data{i}.num_event_h;
    C_r_temp                = C_r_temp                  + ...
        data{i}.C_r_full                 .* data{i}.num_event_r;
    % deal with sigma
    sigma_omega_QG_h_temp            = sigma_omega_QG_h_temp              + ...
        data{i}.sigma_h_full .* data{i}.omega_QG_h_full .* data{i}.num_event_h;
    sigma_omega_QG_r_temp            = sigma_omega_QG_r_temp              + ...
        data{i}.sigma_r_full .* data{i}.omega_QG_r_full .* data{i}.num_event_r;
    dtheta_dp_ma_h_temp     = dtheta_dp_ma_h_temp       + ...
        data{i}.dtheta_dp_ma_h_full      .* data{i}.num_event_h;
    dtheta_dp_ma_omega_h_temp = dtheta_dp_ma_omega_h_temp + ...
        data{i}.dtheta_dp_ma_h_full .* data{i}.omega_h_full .* data{i}.num_event_h;
    dtheta_dp_ma_omega_r_temp = dtheta_dp_ma_omega_r_temp + ...
        data{i}.dtheta_dp_ma_r_full .* data{i}.omega_r_full .* data{i}.num_event_r;
    dtheta_dp_ma_r_temp     = dtheta_dp_ma_r_temp       + ...
        data{i}.dtheta_dp_ma_r_full      .* data{i}.num_event_r;
    if exist('sigma_star_h')
        sigma_star_omega_QG_h_temp = sigma_star_omega_QG_h_temp + ...
            data{f}.sigma_star_h .* data{i}.omega_QG_h_full .* data{i}.num_event_h;
        sigma_star_omega_QG_r_temp = sigma_star_omega_QG_r_temp + ...
            data{f}.sigma_star_r .* data{i}.omega_QG_r_full .* data{i}.num_event_r;
    end
    %if exist('dtheta_dp_ma_avg_h_full')
    %    size(dtheta_dp_ma_avg_h_temp)
    %    size(data{i}.dtheta_dp_ma_avg_h_full)
    %    size(data{i}.num_event_h)
    %    dtheta_dp_ma_avg_h_temp = dtheta_dp_ma_avg_h_temp + ...
    %        data{i}.dtheta_dp_ma_avg_h_full .* data{i}.num_event_h;
    %    dtheta_dp_ma_avg_r_temp = dtheta_dp_ma_avg_r_temp + ...
    %        data{i}.dtheta_dp_ma_avg_r_full .* data{i}.num_event_r;
    %end
    % deal with k2, l2, m2
    omega_k2_sigma_h_temp   = omega_k2_sigma_h_temp     + ...
        data{i}.omega_k2_sigma_h_full    .* data{i}.num_event_h;
    omega_k2_sigma_r_temp   = omega_k2_sigma_r_temp     + ...
        data{i}.omega_k2_sigma_r_full    .* data{i}.num_event_r;
    omega_sigma_h_temp      = omega_sigma_h_temp        + ...
        data{i}.omega_sigma_h_full       .* data{i}.num_event_h;
    omega_sigma_r_temp      = omega_sigma_r_temp        + ...
        data{i}.omega_sigma_r_full       .* data{i}.num_event_r;
    J_l2_h_temp             = J_l2_h_temp               + ...
        data{i}.J_l2_h_full              .* data{i}.num_event_h;
    J_l2_r_temp             = J_l2_r_temp               + ...
        data{i}.J_l2_r_full              .* data{i}.num_event_r;
    J_h_temp                = J_h_temp + ...
        data{i}.J_h_full .* data{i}.num_event_h;
    J_r_temp                = J_r_temp + ...
        data{i}.J_r_full .* data{i}.num_event_r;
    omega_QG_m2_h_temp      = omega_QG_m2_h_temp + ...
        data{i}.m2_h_full .* data{i}.omega_QG_h_full .* data{i}.num_event_h;
    omega_QG_m2_r_temp      = omega_QG_m2_r_temp + ...
        data{i}.m2_r_full .* data{i}.omega_QG_r_full .* data{i}.num_event_r;
    %d2_dp2_omega_QG_grid_h_temp = d2_dp2_omega_QG_grid_h_temp + ...
    %    data{i}.d2_dp2_omega_QG_grid_h_full .* data{i}.num_event_h;
    %d2_dp2_omega_QG_grid_r_temp = d2_dp2_omega_QG_grid_r_temp + ...
    %    data{i}.d2_dp2_omega_QG_grid_r_full .* data{i}.num_event_r;
    if isfield(data{i}, 'T_adv_h_full')
        T_adv_h_temp = T_adv_h_temp + data{i}.T_adv_h_full .* data{i}.num_event_h;
        T_adv_r_temp = T_adv_r_temp + data{i}.T_adv_r_full .* data{i}.num_event_r;
    end
    % total events at each grid point
    num_event_h = num_event_h + data{i}.num_event_h;
    num_event_r = num_event_r + data{i}.num_event_r;
    % omega
    omega_h_temp = omega_h_temp + data{i}.omega_h_full .* data{i}.num_event_h;
    omega_r_temp = omega_r_temp + data{i}.omega_r_full .* data{i}.num_event_r;
    % omega_QG
    omega_QG_h_temp = omega_QG_h_temp + data{i}.omega_QG_h_full .* data{i}.num_event_h;
    omega_QG_r_temp = omega_QG_r_temp + data{i}.omega_QG_r_full .* data{i}.num_event_r;
    % e-folding distance
    nan_ind_Le_h = isnan(data{i}.Le_h);
    nan_ind_Le_r = isnan(data{i}.Le_r);
    temp_num_h = data{i}.num_event_h(:, :, 10) .* double(~nan_ind_Le_h);
    temp_num_r = data{i}.num_event_r(:, :, 10) .* double(~nan_ind_Le_r);
    data{i}.Le_h(nan_ind_Le_h) = 0;
    data{i}.Le_r(nan_ind_Le_r) = 0;
    Le_h_temp = Le_h_temp + data{i}.Le_h .* temp_num_h;
    Le_r_temp = Le_r_temp + data{i}.Le_r .* temp_num_r;
end

omega_h     = omega_h_temp      ./ num_event_h;
omega_r     = omega_r_temp      ./ num_event_r;
omega_QG_h  = omega_QG_h_temp   ./ num_event_h;
omega_QG_r  = omega_QG_r_temp   ./ num_event_r;
Adv_h       = Adv_h_temp         ./ num_event_h;
Adv_r       = Adv_r_temp         ./ num_event_r;
C_h         = C_h_temp          ./ num_event_h;
C_r         = C_r_temp          ./ num_event_r;
J_h         = J_h_temp          ./ num_event_h;
J_r         = J_r_temp          ./ num_event_r;
omega_QG_m2_h_temp = omega_QG_m2_h_temp ./ num_event_h;
omega_QG_m2_r_temp = omega_QG_m2_r_temp ./ num_event_r;
omega_sigma_h_temp = omega_sigma_h_temp ./ num_event_h;
omega_sigma_r_temp = omega_sigma_r_temp ./ num_event_r;
omega_k2_sigma_h_temp = omega_k2_sigma_h_temp ./ num_event_h;
omega_k2_sigma_r_temp = omega_k2_sigma_r_temp ./ num_event_r;
%sigma_omega_QG_h_temp = sigma_omega_QG_h_temp ./ num_event_h;
%sigma_omega_QG_r_temp = sigma_omega_QG_r_temp ./ num_event_r;
dtheta_dp_ma_h = dtheta_dp_ma_h_temp ./ num_event_h;
dtheta_dp_ma_r = dtheta_dp_ma_r_temp ./ num_event_r;
if exist('dtheta_dp_ma_omega_full')
    dtheta_dp_ma_h = dtheta_dp_ma_omega_h_temp ./ num_event_h ./ omega_h;
    dtheta_dp_ma_r = dtheta_dp_ma_omega_r_temp ./ num_event_r ./ omega_r;
end
if exist('sigma_h_star')
    sigma_star_omega_h = sigma_star_omega_QG_h_temp ./ num_event_h ./ omega_QG_h;
    sigma_star_omega_r = sigma_star_omega_QG_r_temp ./ num_event_r ./ omega_QG_r;
end
dtheta_dp_ma_avg_h = dtheta_dp_ma_avg_h_temp ./ num_event_h;
dtheta_dp_ma_avg_r = dtheta_dp_ma_avg_r_temp ./ num_event_r;
if exist('T_adv_h_full')
    T_adv_h     = T_adv_h_temp ./ num_event_h;
    T_adv_r     = T_adv_r_temp ./ num_event_r;
end
sigma_h     = omega_sigma_h_temp ./ omega_QG_h;
sigma_r     = omega_sigma_r_temp ./ omega_QG_r;
k2_h        = omega_k2_sigma_h_temp ./ omega_sigma_h_temp;
k2_r        = omega_k2_sigma_r_temp ./ omega_sigma_r_temp;
Plevels = repmat(reshape(plevels', 1, 1, length(plevels)), size(Adv_h, 1), size(Adv_h, 2), 1);
l2_h        = C_h_temp ./ (kappa ./ Plevels .* J_h_temp);
l2_r        = C_r_temp ./ (kappa ./ Plevels .* J_r_temp);
%l2_h        = J_l2_h_temp ./ J_h_temp;
%l2_r        = J_l2_r_temp ./ J_r_temp;
m2_h        = omega_QG_m2_h_temp ./ omega_QG_h;
m2_r        = omega_QG_m2_r_temp ./ omega_QG_r;
%m2_h        = - d2_dp2_omega_QG_grid_h_temp    ./ omega_QG_h_temp;
%m2_r        = - d2_dp2_omega_QG_grid_r_temp    ./ omega_QG_r_temp;
Le_h        = Le_h_temp ./ num_event_h(:, :, 10);
Le_r        = Le_r_temp ./ num_event_r(:, :, 10);

% set all zeros to nan
names = {'omega_h', 'omega_r', 'omega_QG_h', 'omega_QG_r', ...
         'Adv_h', 'Adv_r', 'C_h', 'C_r', ...
         'sigma_h', 'sigma_r', 'J_h', 'J_r', ...
         'dtheta_dp_ma_h', 'dtheta_dp_ma_r', ...
         'dtheta_dp_ma_avg_h', 'dtheta_dp_ma_avg_r', ...
         'k2_h', 'k2_r', 'l2_h', 'l2_r', ...
         'm2_h', 'm2_r', 'T_adv_h_temp', 'T_adv_r_temp'};
if exist('sigma_star_omega_h')
    names = {names{:}, 'sigma_star_omega_h', 'sigma_star_omega_r'};
end
for n = 1 : length(names)
    temp = eval(names{n});
    zero_ind = temp == 0;
    temp(zero_ind) = NaN;
    eval([names{n}, ' = temp;'])
end
F0 = repmat(F0, [1, 1, length(plevels)]);

