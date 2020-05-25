
clear('data')
data = {};
n_data = length(filenames);

for ind_f = 1 : n_data

    load(filenames{ind_f});

    data{ind_f}.J_l2_h = J_l2_h;
    data{ind_f}.J_l2_r = J_l2_r;
    data{ind_f}.J_h = J_h;
    data{ind_f}.J_r = J_r;
    data{ind_f}.omega_sigma_h = omega_sigma_h;
    data{ind_f}.omega_sigma_r = omega_sigma_r;
    data{ind_f}.omega_k2_sigma_h = omega_k2_sigma_h;
    data{ind_f}.omega_k2_sigma_r = omega_k2_sigma_r;
    data{ind_f}.m2_h = m2_h;
    data{ind_f}.m2_r = m2_r;
    data{ind_f}.d2_dp2_omega_QG_grid_h = d2_dp2_omega_QG_grid_h;
    data{ind_f}.d2_dp2_omega_QG_grid_r = d2_dp2_omega_QG_grid_r;
    data{ind_f}.sigma_h = sigma_h;
    data{ind_f}.sigma_r = sigma_r;
    if exist('sigma_star_h')
        data{ind_f}.sigma_star_h = sigma_star_h;
        data{ind_f}.sigma_star_r = sigma_star_r;
    end
    if exist('sigma_nc_h')
        data{ind_f}.sigma_nc_h = sigma_nc_h;
        data{ind_f}.sigma_nc_r = sigma_nc_r;
    end
    data{ind_f}.Adv_h = Adv_h;
    data{ind_f}.Adv_r = Adv_r;
    data{ind_f}.C_h = C_h;
    data{ind_f}.C_r = C_r;
    data{ind_f}.dtheta_dp_ma_h = dtheta_dp_ma_h;
    data{ind_f}.dtheta_dp_ma_r = dtheta_dp_ma_r;
    data{ind_f}.Le_h = Le_h;
    data{ind_f}.Le_r = Le_r;
    if exist('omega_500hPa_composite_grid_h')
        data{ind_f}.omega_500hPa_composite_grid_h = omega_500hPa_composite_grid_h;
        data{ind_f}.omega_500hPa_composite_grid_r = omega_500hPa_composite_grid_r;
    end
    %if exist('dtheta_dp_ma_avg_h')
    %    data{ind_f}.dtheta_dp_ma_avg_h = dtheta_dp_ma_avg_h;
    %    data{ind_f}.dtheta_dp_ma_avg_r = dtheta_dp_ma_avg_r;
    %end
    if exist('omega_QG_b_h')
        data{ind_f}.omega_QG_b_h = omega_QG_b_h;
        data{ind_f}.omega_QG_b_r = omega_QG_b_r;
    end
    if exist('T_adv_h')
        data{ind_f}.T_adv_h = T_adv_h;
        data{ind_f}.T_adv_r = T_adv_r;
    end
    if exist('num_event_discarded_h')
        data{ind_f}.num_event_discarded_h = repmat(num_event_discarded_h.num_event_discarded, [1, 1, size(omega_h, 3)]);
        data{ind_f}.num_event_discarded_r = repmat(num_event_discarded_r.num_event_discarded, [1, 1, size(omega_h, 3)]);
    end
    if exist('omega_500hPa_std_h')
        data{ind_f}.omega_500hPa_std_h    = omega_500hPa_std_h;
        data{ind_f}.omega_500hPa_std_r    = omega_500hPa_std_r;
        data{ind_f}.omega_QG_500hPa_std_h = omega_QG_500hPa_std_h;
        data{ind_f}.omega_QG_500hPa_std_r = omega_QG_500hPa_std_r;
    end

    data{ind_f}.omega_QG_h      = omega_QG_h;
    data{ind_f}.omega_QG_r      = omega_QG_r;
    data{ind_f}.omega_h         = omega_h;
    data{ind_f}.omega_r         = omega_r;
    data{ind_f}.precip_h        = precip_h;
    data{ind_f}.precip_r        = precip_r;
    data{ind_f}.num_event_h     = repmat(num_event_historical.num_event, [1, 1, size(omega_h, 3)]);
    data{ind_f}.num_event_r     = repmat(num_event_rcp85.num_event, [1, 1, size(omega_h, 3)]);

    % remove omega_QG > 0 at 500hPa
    if ~exist('GT_ZERO') || GT_ZERO
        temp_ind_h = data{ind_f}.omega_QG_h(:, :, plevels == plot_level) >= 0;
        temp_ind_h = repmat(temp_ind_h, 1, 1, length(plevels));
        temp_ind_r = data{ind_f}.omega_QG_r(:, :, plevels == plot_level) >= 0;
        temp_ind_r = repmat(temp_ind_r, 1, 1, length(plevels));
        if exist('num_event_discarded_h');
            data{ind_f}.num_event_discarded_h(temp_ind_h) = data{ind_f}.num_event_discarded_h(temp_ind_h) + ...
                                                            data{ind_f}.num_event_h(temp_ind_h);
            data{ind_f}.num_event_discarded_r(temp_ind_r) = data{ind_f}.num_event_discarded_r(temp_ind_r) + ...
                                                            data{ind_f}.num_event_r(temp_ind_r);
        end
        data{ind_f}.num_event_h(temp_ind_h) = 0;
        data{ind_f}.num_event_r(temp_ind_r) = 0;
    end
    
    % set all nan entries to zero, and set num_event properly
    names = fieldnames(data{ind_f});
    names(ismember(names, 'omega_500hPa_composite_grid_h')) = [];
    names(ismember(names, 'omega_500hPa_composite_grid_r')) = [];
    names(ismember(names, 'Le_h')) = [];
    names(ismember(names, 'Le_r')) = []; % remove this because Le_r does not have z direction
    [names_h, names_r, names_h_2D, names_r_2D] = deal({});
    for n = 1 : length(names)
        if strfind(names{n}, '_h')
            if length(size(getfield(data{ind_f}, names{n}))) == 2
                names_h_2D{end + 1} = names{n};
            else
                names_h{end + 1} = names{n};
            end
        elseif strfind(names{n}, '_r')
            if length(size(getfield(data{ind_f}, names{n}))) == 2
                names_r_2D{end + 1} = names{n};
            else
                names_r{end + 1} = names{n};
            end
        end
    end
    
    % get nan locations in both historical and rcp85
    [nan_ind_h, nan_ind_r, nan_ind_h_2D, nan_ind_r_2D] = deal(zeros(size(data{ind_f}.omega_h)));
    for n = 1 : length(names_h)
        temp = getfield(data{ind_f}, names_h{n});
        nan_ind_h = nan_ind_h | isnan(temp);
    end
    for n = 1 : length(names_h)
        temp = getfield(data{ind_f}, names_h{n});
        temp(nan_ind_h) = 0;
        data{ind_f}.(names_h{n}) = temp;
    end
    for n = 1 : length(names_r)
        temp = getfield(data{ind_f}, names_r{n});
        nan_ind_r = nan_ind_r | isnan(temp);
    end
    for n = 1 : length(names_r)
        temp = getfield(data{ind_f}, names_r{n});
        temp(nan_ind_r) = 0;
        data{ind_f}.(names_r{n}) = temp;
    end
    for n = 1 : length(names_h)
        temp = getfield(data{ind_f}, names_h{n});
        nan_ind_h_2D = nan_ind_h_2D | isnan(temp);
    end
    for n = 1 : length(names_h)
        temp = getfield(data{ind_f}, names_h{n});
        temp(nan_ind_h_2D) = 0;
        data{ind_f}.(names_h{n}) = temp;
    end
    for n = 1 : length(names_r)
        temp = getfield(data{ind_f}, names_r{n});
        nan_ind_r_2D = nan_ind_r_2D | isnan(temp);
    end
    for n = 1 : length(names_r)
        temp = getfield(data{ind_f}, names_r{n});
        temp(nan_ind_r_2D) = 0;
        data{ind_f}.(names_r{n}) = temp;
    end


    for n = 1 : length(names)
        if strfind(names{n}, '_h')
            if exist('num_event_discarded_h');
                data{ind_f}.num_event_discarded_h(nan_ind_h) = data{ind_f}.num_event_discarded_h(nan_ind_h) + ...
                                                               data{ind_f}.num_event_h(nan_ind_h);
            end
            % this does not have overcounting issues because data{ind_f}.num_event_h(nan_ind_h) is set to zero
            data{ind_f}.num_event_h(nan_ind_h) = 0;
        elseif strfind(names{n}, '_r')
            if exist('num_event_discarded_r');
                data{ind_f}.num_event_discarded_r(nan_ind_r) = data{ind_f}.num_event_discarded_r(nan_ind_r) + ...
                                                               data{ind_f}.num_event_r(nan_ind_r);
            end
            data{ind_f}.num_event_r(nan_ind_r) = 0;
        end
    end

    % set all entries that correspond to zero num_event to zero as well
    names(ismember(names, 'num_event_discarded_h')) = []; % don't modify these entries
    names(ismember(names, 'num_event_discarded_r')) = [];
    for n = 1 : length(names)
        temp = getfield(data{ind_f}, names{n});
        if length(size(temp)) == 2
            zero_ind_h = data{ind_f}.num_event_h(:, :, plevels == plot_level) == 0;
            zero_ind_r = data{ind_f}.num_event_r(:, :, plevels == plot_level) == 0;
        else
            zero_ind_h = data{ind_f}.num_event_h == 0;
            zero_ind_r = data{ind_f}.num_event_r == 0;
        end
        if strfind(names{n}, '_h')
            temp(zero_ind_h) = 0;
        elseif strfind(names{n}, '_r')
            temp(zero_ind_r) = 0;
        end
        data{ind_f}.(names{n}) = temp;
    end
end

[precip_h_temp, precip_r_temp] = deal(zeros([size(data{1}.J_h, 1), size(data{1}.J_h, 2)]));
[sigma_omega_h_temp, sigma_omega_r_temp, ...
 Adv_h_temp, Adv_r_temp, C_h_temp, C_r_temp, ...
 dtheta_dp_ma_h_temp, dtheta_dp_ma_r_temp, ...
 sigma_star_omega_QG_h_temp, sigma_star_omega_QG_r_temp, ...
 dtheta_dp_ma_omega_h_temp, dtheta_dp_ma_omega_r_temp, ...
 dtheta_dp_ma_avg_h_temp, dtheta_dp_ma_avg_r_temp, ...
 omega_QG_h_temp, omega_QG_r_temp, omega_h_temp, omega_r_temp, ...
 omega_QG_b_h_temp, omega_QG_b_r_temp, ...
 omega_k2_sigma_h_temp, omega_sigma_h_temp, ...
 omega_k2_sigma_r_temp, omega_sigma_r_temp, ...
 sigma_nc_h_temp, sigma_nc_r_temp, ...
 J_l2_h_temp, J_h_temp, ...
 J_l2_r_temp, J_r_temp, ...
 ...%d2_dp2_omega_QG_grid_h_temp, d2_dp2_omega_QG_grid_r_temp, ...
 sigma_omega_QG_h_temp, sigma_omega_QG_r_temp, ...
 omega_QG_m2_h_temp, omega_QG_m2_r_temp, ...
 num_event_h, num_event_r, ...
 num_event_discarded_h, num_event_discarded_r, ...
 T_adv_h_temp, T_adv_r_temp] = deal(zeros(size(data{1}.J_h)));
[Le_h_temp, Le_r_temp] = deal(zeros(size(data{1}.J_h(:, :, 1))));
for i = 1 : n_data
    % deal with AB
    Adv_h_temp               = Adv_h_temp                 + ...
        data{i}.Adv_h                .* data{i}.num_event_h;
    Adv_r_temp               = Adv_r_temp                 + ...
        data{i}.Adv_r                .* data{i}.num_event_r;
    % deal with C
    C_h_temp                = C_h_temp                  + ...
        data{i}.C_h                 .* data{i}.num_event_h;
    C_r_temp                = C_r_temp                  + ...
        data{i}.C_r                 .* data{i}.num_event_r;
    % deal with sigma
    if exist('sigma_nc_h')
        sigma_nc_h_temp     = sigma_nc_h_temp + data{i}.sigma_nc_h .* data{i}.num_event_h;
        sigma_nc_r_temp     = sigma_nc_r_temp + data{i}.sigma_nc_r .* data{i}.num_event_r;
    end
    sigma_omega_QG_h_temp            = sigma_omega_QG_h_temp              + ...
        data{i}.sigma_h .* data{i}.omega_QG_h .* data{i}.num_event_h;
    sigma_omega_QG_r_temp            = sigma_omega_QG_r_temp              + ...
        data{i}.sigma_r .* data{i}.omega_QG_r .* data{i}.num_event_r;
    dtheta_dp_ma_h_temp     = dtheta_dp_ma_h_temp       + ...
        data{i}.dtheta_dp_ma_h      .* data{i}.num_event_h;
    dtheta_dp_ma_omega_h_temp = dtheta_dp_ma_omega_h_temp + ...
        data{i}.dtheta_dp_ma_h .* data{i}.omega_h .* data{i}.num_event_h;
    dtheta_dp_ma_omega_r_temp = dtheta_dp_ma_omega_r_temp + ...
        data{i}.dtheta_dp_ma_r .* data{i}.omega_r .* data{i}.num_event_r;
    dtheta_dp_ma_r_temp     = dtheta_dp_ma_r_temp       + ...
        data{i}.dtheta_dp_ma_r      .* data{i}.num_event_r;
    if exist('sigma_star_h')
        sigma_star_omega_QG_h_temp = sigma_star_omega_QG_h_temp + ...
            data{i}.sigma_star_h .* data{i}.omega_QG_h .* data{i}.num_event_h;
        sigma_star_omega_QG_r_temp = sigma_star_omega_QG_r_temp + ...
            data{i}.sigma_star_r .* data{i}.omega_QG_r .* data{i}.num_event_r;
    end
    %if exist('dtheta_dp_ma_avg_h')
    %    size(dtheta_dp_ma_avg_h_temp)
    %    size(data{i}.dtheta_dp_ma_avg_h)
    %    size(data{i}.num_event_h)
    %    dtheta_dp_ma_avg_h_temp = dtheta_dp_ma_avg_h_temp + ...
    %        data{i}.dtheta_dp_ma_avg_h .* data{i}.num_event_h;
    %    dtheta_dp_ma_avg_r_temp = dtheta_dp_ma_avg_r_temp + ...
    %        data{i}.dtheta_dp_ma_avg_r .* data{i}.num_event_r;
    %end
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
    omega_QG_m2_h_temp      = omega_QG_m2_h_temp + ...
        data{i}.m2_h .* data{i}.omega_QG_h .* data{i}.num_event_h;
    omega_QG_m2_r_temp      = omega_QG_m2_r_temp + ...
        data{i}.m2_r .* data{i}.omega_QG_r .* data{i}.num_event_r;
    %d2_dp2_omega_QG_grid_h_temp = d2_dp2_omega_QG_grid_h_temp + ...
    %    data{i}.d2_dp2_omega_QG_grid_h .* data{i}.num_event_h;
    %d2_dp2_omega_QG_grid_r_temp = d2_dp2_omega_QG_grid_r_temp + ...
    %    data{i}.d2_dp2_omega_QG_grid_r .* data{i}.num_event_r;
    if isfield(data{i}, 'T_adv_h')
        T_adv_h_temp = T_adv_h_temp + data{i}.T_adv_h .* data{i}.num_event_h;
        T_adv_r_temp = T_adv_r_temp + data{i}.T_adv_r .* data{i}.num_event_r;
    end
    % total events at each grid point
    num_event_h = num_event_h + data{i}.num_event_h;
    num_event_r = num_event_r + data{i}.num_event_r;
    % total discarded events
    if isfield(data{i}, 'num_event_discarded_h')
        num_event_discarded_h = num_event_discarded_h + data{i}.num_event_discarded_h;
        num_event_discarded_r = num_event_discarded_r + data{i}.num_event_discarded_r;
    end
    % omega
    omega_h_temp = omega_h_temp + data{i}.omega_h .* data{i}.num_event_h;
    omega_r_temp = omega_r_temp + data{i}.omega_r .* data{i}.num_event_r;
    % omega_QG
    omega_QG_h_temp = omega_QG_h_temp + data{i}.omega_QG_h .* data{i}.num_event_h;
    omega_QG_r_temp = omega_QG_r_temp + data{i}.omega_QG_r .* data{i}.num_event_r;
    % omega_QG with full boundaries
    if isfield(data{i}, 'omega_QG_b_h')
        omega_QG_b_h_temp = omega_QG_b_h_temp + data{i}.omega_QG_b_h .* data{i}.num_event_h;
        omega_QG_b_r_temp = omega_QG_b_r_temp + data{i}.omega_QG_b_r .* data{i}.num_event_r;
    end
    % precipitation values
    precip_h_temp = precip_h_temp + data{i}.precip_h .* data{i}.num_event_h(:, :, plevels == plot_level);
    precip_r_temp = precip_r_temp + data{i}.precip_r .* data{i}.num_event_r(:, :, plevels == plot_level);
    % get composite 2D fields at 500hPa to calculate e-folding distance later
    if exist('omega_500hPa_composite_grid_h')
        if i == 1
            % this is to partially define the structure of omega_500hPa_composite_grid_h for later use
            omega_500hPa_composite_grid_h_full = data{i}.omega_500hPa_composite_grid_h;
            omega_500hPa_composite_grid_r_full = data{i}.omega_500hPa_composite_grid_r;
        else
            temp_ind = find(~cellfun(@isempty, data{i}.omega_500hPa_composite_grid_h) | ...
                            ~cellfun(@isempty, data{i}.omega_500hPa_composite_grid_r));
            for ind = 1 : length(temp_ind)
                j = temp_ind(ind);
                [lat_j, lon_i] = ind2sub(size(omega_500hPa_composite_grid_h_full), j);
                if ~isempty(data{i}.omega_500hPa_composite_grid_h{j})
                    if isempty(omega_500hPa_composite_grid_h_full{j})
                        omega_500hPa_composite_grid_h_full{j} = ...
                                data{i}.omega_500hPa_composite_grid_h{j} * data{i}.num_event_h(lat_j, lon_i, plevels == plot_level);
                    else
                        omega_500hPa_composite_grid_h_full{j} = omega_500hPa_composite_grid_h_full{j} + ...
                                data{i}.omega_500hPa_composite_grid_h{j} * data{i}.num_event_h(lat_j, lon_i, plevels == plot_level);
                    end
                end
                if ~isempty(data{i}.omega_500hPa_composite_grid_r{j})
                    if isempty(omega_500hPa_composite_grid_r_full{j})
                        omega_500hPa_composite_grid_r_full{j} = ...
                                data{i}.omega_500hPa_composite_grid_r{j} * data{i}.num_event_r(lat_j, lon_i, plevels == plot_level);
                    else
                        omega_500hPa_composite_grid_r_full{j} = omega_500hPa_composite_grid_r_full{j} + ...
                                data{i}.omega_500hPa_composite_grid_r{j} * data{i}.num_event_r(lat_j, lon_i, plevels == plot_level);
                    end
                end
            end
        end
    end
    nan_ind_Le_h = isnan(data{i}.Le_h);
    nan_ind_Le_r = isnan(data{i}.Le_r);
    temp_num_h = data{i}.num_event_h(:, :, plevels == plot_level) .* double(~nan_ind_Le_h);
    temp_num_r = data{i}.num_event_r(:, :, plevels == plot_level) .* double(~nan_ind_Le_r);
    data{i}.Le_h(nan_ind_Le_h) = 0;
    data{i}.Le_r(nan_ind_Le_r) = 0;
    Le_h_temp = Le_h_temp + data{i}.Le_h .* temp_num_h;
    Le_r_temp = Le_r_temp + data{i}.Le_r .* temp_num_r;
end

omega_h     = omega_h_temp      ./ num_event_h;
omega_r     = omega_r_temp      ./ num_event_r;
omega_QG_h  = omega_QG_h_temp   ./ num_event_h;
omega_QG_r  = omega_QG_r_temp   ./ num_event_r;
if isfield(data{i}, 'omega_QG_b_h')
    omega_QG_b_h = omega_QG_b_h_temp ./ num_event_h;
    omega_QG_b_r = omega_QG_b_r_temp ./ num_event_r;
end
precip_h    = precip_h_temp     ./ num_event_h(:, :, plevels == plot_level);
precip_r    = precip_r_temp     ./ num_event_r(:, :, plevels == plot_level);
Adv_h       = Adv_h_temp        ./ num_event_h;
Adv_r       = Adv_r_temp        ./ num_event_r;
C_h         = C_h_temp          ./ num_event_h;
C_r         = C_r_temp          ./ num_event_r;
J_h         = J_h_temp          ./ num_event_h;
J_r         = J_r_temp          ./ num_event_r;
if exist('sigma_nc_h')
    sigma_nc_h  = sigma_nc_h_temp   ./ num_event_h;
    sigma_nc_r  = sigma_nc_r_temp   ./ num_event_r;
end
omega_QG_m2_h_temp = omega_QG_m2_h_temp ./ num_event_h;
omega_QG_m2_r_temp = omega_QG_m2_r_temp ./ num_event_r;
omega_sigma_h_temp = omega_sigma_h_temp ./ num_event_h;
omega_sigma_r_temp = omega_sigma_r_temp ./ num_event_r;
omega_k2_sigma_h_temp = omega_k2_sigma_h_temp ./ num_event_h;
omega_k2_sigma_r_temp = omega_k2_sigma_r_temp ./ num_event_r;
sigma_omega_QG_h_temp = sigma_omega_QG_h_temp ./ num_event_h;
sigma_omega_QG_r_temp = sigma_omega_QG_r_temp ./ num_event_r;
dtheta_dp_ma_h = dtheta_dp_ma_h_temp ./ num_event_h;
dtheta_dp_ma_r = dtheta_dp_ma_r_temp ./ num_event_r;
if exist('dtheta_dp_ma_omega')
    dtheta_dp_ma_h = dtheta_dp_ma_omega_h_temp ./ num_event_h ./ omega_h;
    dtheta_dp_ma_r = dtheta_dp_ma_omega_r_temp ./ num_event_r ./ omega_r;
end
if exist('sigma_star_h')
    sigma_star_h = sigma_star_omega_QG_h_temp ./ num_event_h ./ omega_QG_h;
    sigma_star_r = sigma_star_omega_QG_r_temp ./ num_event_r ./ omega_QG_r;
end
dtheta_dp_ma_avg_h = dtheta_dp_ma_avg_h_temp ./ num_event_h;
dtheta_dp_ma_avg_r = dtheta_dp_ma_avg_r_temp ./ num_event_r;
if exist('T_adv_h')
    T_adv_h     = T_adv_h_temp ./ num_event_h;
    T_adv_r     = T_adv_r_temp ./ num_event_r;
end
sigma_h     = sigma_omega_QG_h_temp ./ omega_QG_h;
sigma_r     = sigma_omega_QG_r_temp ./ omega_QG_r;
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
if exist('omega_500hPa_composite_grid_h')
    [Le_h, Le_r] = deal(zeros(size(omega_500hPa_composite_grid_h_full)));
    [latmax, latmin, lonmax, lonmin] = get_lat_lon_limits('/disk7/ziweili/CESM_LENS/exp/001_2D_varsig_Q_ocean_desert/');
    [lat_indices, lon_indices] = latlonindices(lat_series, lon_series, latmin, latmax, lonmin, lonmax);
    lat = lat_series(lat_indices);
    
    temp_num_event_h = num_event_h(:, :, plevels == plot_level);
    temp_ind = find(~cellfun(@isempty, omega_500hPa_composite_grid_h_full));
    for ind = 1 : length(temp_ind)
        j = temp_ind(ind);
        [lat_j, lon_i] = ind2sub(size(omega_500hPa_composite_grid_h_full), j);
        omega_composite = omega_500hPa_composite_grid_h_full{j} / temp_num_event_h(j);
        n_lon = size(omega_composite, 2);
        n_lat = size(omega_composite, 1);
        temp_lon_indices = mod((lon_indices(lon_i) - (n_lon - 1) / 2 : lon_indices(lon_i) + (n_lon - 1) / 2) - 1, length(lon)) + 1;
        temp_lat_indices = (lat_indices(lat_j) - (n_lat - 1) / 2 : lat_indices(lat_j) + (n_lat - 1) / 2);
        event_omega_avg = omega_avg_h(temp_lon_indices, temp_lat_indices, plevels_2 == 50000)';
        Le_h(lat_j, lon_i) = e_fold_v1(omega_composite, zeros(size(event_omega_avg)), dphi, dlambda, lat(lat_j)/180*pi, R);
        omega_500hPa_composite_grid_h_full{j} = omega_composite;
    end
    
    temp_num_event_r = num_event_r(:, :, plevels == plot_level);
    temp_ind = find(~cellfun(@isempty, omega_500hPa_composite_grid_r_full));
    for ind = 1 : length(temp_ind)    
        j = temp_ind(ind);
        [lat_j, lon_i] = ind2sub(size(omega_500hPa_composite_grid_r_full), j);
        omega_composite = omega_500hPa_composite_grid_r_full{j} / temp_num_event_r(j);
        n_lon = size(omega_composite, 2);        
        n_lat = size(omega_composite, 1);
        temp_lon_indices = mod((lon_indices(lon_i) - (n_lon - 1) / 2 : lon_indices(lon_i) + (n_lon - 1) / 2) - 1, length(lon)) + 1;
        temp_lat_indices = (lat_indices(lat_j) - (n_lat - 1) / 2 : lat_indices(lat_j) + (n_lat - 1) / 2);
        event_omega_avg = omega_avg_r(temp_lon_indices, temp_lat_indices, plevels_2 == 50000)';
        Le_r(lat_j, lon_i) = e_fold_v1(omega_composite, zeros(size(event_omega_avg)), dphi, dlambda, lat(lat_j)/180*pi, R);
        omega_500hPa_composite_grid_r_full{j} = omega_composite;
    end
else
    Le_h        = Le_h_temp ./ num_event_h(:, :, plevels == plot_level);
    Le_r        = Le_r_temp ./ num_event_r(:, :, plevels == plot_level);
end

% get standard deviations
if exist('omega_500hPa_std_h')
    [omega_var_h, omega_var_r, omega_QG_var_h, omega_QG_var_r] = deal(zeros([size(data{1}.J_h(:, :, 1)), n_data]));
    for i = 1 : n_data
        omega_var_h   (:, :, i) = data{i}.omega_500hPa_std_h.^2    .* data{i}.num_event_h(:, :, plevels == plot_level);
        omega_var_r   (:, :, i) = data{i}.omega_500hPa_std_r.^2    .* data{i}.num_event_r(:, :, plevels == plot_level);
        omega_QG_var_h(:, :, i) = data{i}.omega_QG_500hPa_std_h.^2 .* data{i}.num_event_h(:, :, plevels == plot_level);
        omega_QG_var_r(:, :, i) = data{i}.omega_QG_500hPa_std_r.^2 .* data{i}.num_event_r(:, :, plevels == plot_level);
    end
    omega_500hPa_std_h = sqrt(squeeze(sum(omega_var_h, 3))./num_event_h(:, :, plevels == plot_level));
    omega_500hPa_std_r = sqrt(squeeze(sum(omega_var_r, 3))./num_event_r(:, :, plevels == plot_level));
    omega_QG_500hPa_std_h = sqrt(squeeze(sum(omega_QG_var_h, 3))./num_event_h(:, :, plevels == plot_level));
    omega_QG_500hPa_std_r = sqrt(squeeze(sum(omega_QG_var_r, 3))./num_event_r(:, :, plevels == plot_level));
end

% set all zeros to nan
names = {'omega_h', 'omega_r', 'omega_QG_h', 'omega_QG_r', ...
         'omega_QG_b_h', 'omega_QG_b_r', ...
         'precip_h', 'precip_r', ...
         'Adv_h', 'Adv_r', 'C_h', 'C_r', ...
         'sigma_h', 'sigma_r', 'J_h', 'J_r', ...
         'dtheta_dp_ma_h', 'dtheta_dp_ma_r', ...
         'dtheta_dp_ma_avg_h', 'dtheta_dp_ma_avg_r', ...
         'k2_h', 'k2_r', 'l2_h', 'l2_r', ...
         'm2_h', 'm2_r', 'T_adv_h_temp', 'T_adv_r_temp'};
if exist('sigma_star_h')
    names = {names{:}, 'sigma_star_h', 'sigma_star_r'};
end
for n = 1 : length(names)
    if exist(names{n})
        temp = eval(names{n});
        zero_ind = temp == 0;
        temp(zero_ind) = NaN;
        eval([names{n}, ' = temp;'])
    end
end


