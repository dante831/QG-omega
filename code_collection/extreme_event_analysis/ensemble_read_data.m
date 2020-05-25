
data = {};
n_data = length(filenames);
for f = 1 : n_data

    load(filenames{f});

    data{f}.J_l2_h = J_l2_h;
    data{f}.J_l2_r = J_l2_r;
    data{f}.J_h = J_h;
    data{f}.J_r = J_r;
    data{f}.omega_sigma_h = omega_sigma_h;
    data{f}.omega_sigma_r = omega_sigma_r;
    data{f}.omega_k2_sigma_h = omega_k2_sigma_h;
    data{f}.omega_k2_sigma_r = omega_k2_sigma_r;
    data{f}.d2_dp2_omega_QG_grid_h = d2_dp2_omega_QG_grid_h;
    data{f}.d2_dp2_omega_QG_grid_r = d2_dp2_omega_QG_grid_r;
    data{f}.sigma_h = sigma_h;
    data{f}.sigma_r = sigma_r;
    data{f}.AB_h = AB_h;
    data{f}.AB_r = AB_r;
    data{f}.C_h = C_h;
    data{f}.C_r = C_r;
    data{f}.dtheta_dp_ma_h = dtheta_dp_ma_h;
    data{f}.dtheta_dp_ma_r = dtheta_dp_ma_r;

    data{f}.omega_QG_h      = omega_QG_h;
    data{f}.omega_QG_r      = omega_QG_r;
    data{f}.omega_h         = omega_h;
    data{f}.omega_r         = omega_r;
    data{f}.num_event_h     = num_event_historical.num_event;
    data{f}.num_event_r     = num_event_rcp85.num_event;

    % set all nan to zero
    names = fieldnames(data{f});
    for n = 1 : length(names)
        temp = getfield(data{f}, names{n});
        nan_ind = isnan(temp);
        temp(nan_ind) = 0;
        data{f}.(names{n}) = temp;
    end
end

[sigma_h_temp, sigma_r_temp, ...
 AB_h_temp, AB_r_temp, C_h_temp, C_r_temp, ...
 dtheta_dp_ma_h_temp, dtheta_dp_ma_r_temp, ...
 omega_QG_h_temp, omega_QG_r_temp, omega_h_temp, omega_r_temp, ...
 omega_k2_sigma_h_temp, omega_sigma_h_temp, ...
 omega_k2_sigma_r_temp, omega_sigma_r_temp, ...
 J_l2_h_temp, J_h_temp, ...
 J_l2_r_temp, J_r_temp, ...
 d2_dp2_omega_QG_grid_h_temp, d2_dp2_omega_QG_grid_r_temp, ...
 num_event_h, num_event_r, ...
 omega_QG_h_temp] = deal(zeros(size(data{1}.J_h)));
for i = 1 : n_data
    % deal with AB
    AB_h_temp               = AB_h_temp                 + ...
        data{i}.AB_h                .* data{i}.num_event_h;
    AB_r_temp               = AB_r_temp                 + ...
        data{i}.AB_r                .* data{i}.num_event_r;
    % deal with C
    C_h_temp                = C_h_temp                  + ...
        data{i}.C_h                 .* data{i}.num_event_h;
    C_r_temp                = C_r_temp                  + ...
        data{i}.C_r                 .* data{i}.num_event_r;
    % deal with sigma
    sigma_h_temp            = sigma_h_temp              + ...
        data{i}.sigma_h             .* data{i}.num_event_h;
    sigma_r_temp            = sigma_r_temp              + ...
        data{i}.sigma_r             .* data{i}.num_event_r;
    dtheta_dp_ma_h_temp     = dtheta_dp_ma_h_temp       + ...
        data{i}.dtheta_dp_ma_h      .* data{i}.num_event_h;
    dtheta_dp_ma_r_temp     = dtheta_dp_ma_r_temp       + ...
        data{i}.dtheta_dp_ma_r      .* data{i}.num_event_r;
    % deal with k2, l2, m2
    omega_k2_sigma_h_temp   = omega_k2_sigma_h_temp     + ...
        data{i}.omega_k2_sigma_h    .* data{i}.num_event_h;
    omega_k2_sigma_r_temp   = omega_k2_sigma_r_temp     + ...
        data{i}.omega_k2_sigma_r    .* data{i}.num_event_r;
    omega_sigma_h_temp      = omega_sigma_h_temp        + ...
        data{i}.omega_sigma_h       .* data{i}.num_event_h;
    omega_sigma_r_temp      = omega_sigma_r_temp        + ...
        data{i}.omega_sigma_r       .* data{i}.num_event_r;
    J_l2_h_temp             = J_l2_h_temp               + ...
        data{i}.J_l2_h              .* data{i}.num_event_h;
    J_l2_r_temp             = J_l2_r_temp               + ...
        data{i}.J_l2_r              .* data{i}.num_event_r;
    J_h_temp                = J_h_temp + ...
        data{i}.J_h .* data{i}.num_event_h;
    J_r_temp                = J_r_temp + ...
        data{i}.J_r .* data{i}.num_event_r;
    d2_dp2_omega_QG_grid_h_temp = d2_dp2_omega_QG_grid_h_temp + ...
        data{i}.d2_dp2_omega_QG_grid_h .* data{i}.num_event_h;
    d2_dp2_omega_QG_grid_r_temp = d2_dp2_omega_QG_grid_r_temp + ...
        data{i}.d2_dp2_omega_QG_grid_r .* data{i}.num_event_r;
    % total events at each grid point
    num_event_h = num_event_h + data{i}.num_event_h;
    num_event_r = num_event_r + data{i}.num_event_r;
    % omega
    omega_h_temp = omega_h_temp + data{i}.omega_h .* data{i}.num_event_h;
    omega_r_temp = omega_r_temp + data{i}.omega_r .* data{i}.num_event_r;
    % omega_QG
    omega_QG_h_temp = omega_QG_h_temp + data{i}.omega_QG_h .* data{i}.num_event_h;
    omega_QG_r_temp = omega_QG_r_temp + data{i}.omega_QG_r .* data{i}.num_event_r;
end


omega_h     = omega_h_temp      ./ num_event_h;
omega_r     = omega_r_temp      ./ num_event_r;
omega_QG_h  = omega_QG_h_temp   ./ num_event_h;
omega_QG_r  = omega_QG_r_temp   ./ num_event_r;
AB_h        = AB_h_temp         ./ num_event_h;
AB_r        = AB_r_temp         ./ num_event_r;
C_h         = C_h_temp          ./ num_event_h;
C_r         = C_r_temp          ./ num_event_r;
sigma_h     = sigma_h_temp      ./ num_event_h;
sigma_r     = sigma_r_temp      ./ num_event_r;
J_h         = J_h_temp          ./ num_event_h;
J_r         = J_r_temp          ./ num_event_r;
dtheta_dp_ma_h = dtheta_dp_ma_h_temp ./ num_event_h;
dtheta_dp_ma_r = dtheta_dp_ma_r_temp ./ num_event_r;
k2_h = omega_k2_sigma_h_temp ./ omega_sigma_h_temp;
k2_r = omega_k2_sigma_r_temp ./ omega_sigma_r_temp;
l2_h = J_l2_h_temp ./ J_h_temp;
l2_r = J_l2_r_temp ./ J_r_temp;
m2_h = - d2_dp2_omega_QG_grid_h_temp    ./ omega_QG_h_temp;
m2_r = - d2_dp2_omega_QG_grid_r_temp    ./ omega_QG_r_temp;


