
    % add functions from ../../source/
    addpath('/net/aimsir/archive1/ziweili/CESM_LENS/source/')

    % some constants

    latmax = 47; % degrees
    latmin = 44;
    lonmax = 183; % degrees
    lonmin = 177;
    n_precip = 1; % n_precip = 1 for 6hourly, n_precip = 4 for daily in CESM-LENS data
    ocean_only = false;
    MASK = true;
    box_min = 15;
    box_max = 29;
    h_max = 500;
    p_min = 10000; % upper limit
    p_max = 55000; % lower limit, events should at least have full non-NaN values at this level
    p_mask = 55000; % lower limit for masking
    NETCDF = true; % whether output event omega to NetCDF files
    NETCDF_interval = 1;
    n_ensemble = 5;
    quantile_file_h = '/net/aimsir/archive1/ziweili/CESM_LENS/output/precip_99.9th_quantile_h_corrected.mat';
    percentile = 0.999;

    FIND_CENTER = true;
    JJA_DJF = 0;
    ACCU_SIGMA = true;
    SIGMA_LOCAL = false;

    set(0,'DefaultFigureVisible','off');

    %set and create paths
    data_path = '/net/aimsir/archive1/ziweili/CESM_LENS/data/';
    input_path = '/net/aimsir/archive1/ziweili/CESM_LENS/output/';
    output_path = '/net/aimsir/archive1/ziweili/CESM_LENS/paper_code_github/QG-omega/data/';
    if ~exist(output_path)
        mkdir(output_path)
    end

    historical_files = {[data_path, ...
        'b.e11.B20TRC5CNBDRD.f09_g16.', num2str(n_ensemble, '%.3d'), '.cam.h2.PRECT.1990010100Z-2005123118Z.nc']};

    lat_series_original = ncread(historical_files{1}, 'lat');
    lon_series_original = ncread(historical_files{1}, 'lon');
    [lat_indices, lon_indices] = latlonindices(lat_series_original, lon_series_original, latmin, latmax, lonmin, lonmax);


    %% extreme_precip part
    
    disp('calculate 99.9th percentile of historical precipitation');
    if ~exist('historical_precip')
        historical_precip = precip_with_mask(n_precip, lon_indices, lat_indices, historical_files, MASK, h_max, ocean_only, p_mask);
    end
    temp = load(quantile_file_h, 'Q');
    historical_quantile = temp.Q(lon_indices, lat_indices); clear('temp');
    [days_historical, event_precip_historical] = deal(cell(length(lat_indices), length(lon_indices)));

    % calculate QG_omega

    for j = 1 : length(lat_indices)
        for i = 1 : length(lon_indices)

            % get lnogitude and day of one event. days are converted to 'days after 1990'
            days_historical{j, i} = find(historical_precip(i, j, :) > historical_quantile(i, j));
            ind = sub2ind(size(historical_precip), ...
                    repmat(i, length(days_historical{j, i}), 1), ...
                    repmat(j, length(days_historical{j, i}), 1), ...
                    days_historical{j, i});
            event_precip_historical{j, i} = historical_precip(ind);

            % convert to real time values
            days_historical{j, i} = days_historical{j, i} * 0.25;

            % get rid of time boundary values
            days_select_1 = days_historical{j, i} >= 365 * (2006 - 1990) - 1; % remove the last day
            days_select_2 = days_historical{j, i} <= 0.25; % remove the first day, 6:00am
            days_select_3 = days_historical{j, i} >= 365 * (2001 - 1990); % remove years 2001 - 2005
            days_select_4 = days_historical{j, i} <= 365 * (1991 - 1990) + 0.25; % remove year 1990
            days_historical{j, i}(days_select_1 | days_select_2 | ...
                                  days_select_3 | days_select_4) = [];

        end
    end
    matfilename = 'precip_events_historical.mat';

    %% event_analysis part

    days_ = days_historical;
    precip_ = event_precip_historical;

    % definitions of useful constants and arrays

    Omega = 7.2921e-5;  % Earth's angular speed in rad/s 
    R = 6371000.0;      % average radius from https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html, in meters
    Ra = 287.04;        % specific gas constant, in J/(kg*K)
    cp = 1005.0;        % heat capacity at 300K, in J/(kg*K)
    kappa = Ra/cp;
    dt = 3600.0 * 6.0;  % time difference, in seconds
    window_x = 1;       % smoothing window. This translates to a (2*window_x+1)-by-(2*window_y+1) smoothing filter
    window_y = 1;
    NetCDF_label = 0;   % NetCDF file counter

    lon_series = lon_series_original(lon_indices); % longitude series that covers the regions in the analysis
    lat_series = lat_series_original(lat_indices); % latitude series that covers the regions in the analysis
    dlambda    = (lon_series(2) - lon_series(1)) / 180.0 * 3.1415926; % x angular increment
    dphi       = (lat_series(2) - lat_series(1)) / 180.0 * 3.1415926; % y angular increment

    [num_event, num_event_discarded, num_event_discarded_misalign] = deal(zeros(length(lat_indices), length(lon_indices))); % define array for number of events
    events_sta = deal(cell(length(lat_indices), length(lon_indices))); % define cell array for "event statistics"
    pwd_str = pwd; % get current path

    jj = 1;
    ii = 1;
    m = 1;

    event_day = days_{jj, ii}(m);
    event_lon = lon_indices(ii);
    event_lat = lat_indices(jj);

    box_max_x = box_max;
    box_max_y = box_max;
    computation_y = box_max_y;
    computation_x = box_max_x;

    event_year = floor((event_day - 0.25) / 365) + 1990;
    string_1 = num2str(n_ensemble, '%.3d'); % string for the ensemble number 
        
    % get whether it is historical or rcp85
    if event_year >= 1990 && event_year <= 2005
        string_2 = 'historical';
        tag = true;
    elseif event_year >= 2071 && event_year <= 2080
        string_2 = 'rcp85';
        tag = true;
    else
        disp(['year ', num2str(event_year), ' does not exist in the data, abandoned']);
        tag = false;
    end
    
    % calculate lat and lon spans for the event
    [event_latspan, event_lonspan, center_y, computation_y] = ...
            event_span_Wrapper(lat_series_original, lon_series_original, event_lat, event_lon, ...
            box_max_y, computation_x, computation_y);

    plevels = double(ncread([input_path, 'ta_', string_1, '_', num2str(event_year), '.nc'], 'level'));
    p_start_0 = 1; % lowest level
    p_start = find(plevels == p_max); % lower boundary start point
    p_end = find(plevels == p_min); % highest level
    level_indices = p_start_0 : p_end;
    level = plevels(level_indices);

    % read in omega field
    [omega, tag] = omega_read_Wrapper(computation_y, computation_x, length(level_indices), ...
            event_year, input_path, event_latspan, event_lonspan, event_day, level_indices, string_1);

    % read in event data
    [T, ua, va, ug, vg, omega, omega_b, q, tag] = event_read_Wrapper_uava_q(computation_y, computation_x, length(level_indices), ...
            event_year, input_path, event_latspan, event_lonspan, event_day, level_indices, string_1, JJA_DJF);
    
    % define NetCDF file name, and check if an old file exists
    nc_filename = [output_path, 'example_event.nc'];
    if exist(nc_filename, 'file') == 2
        delete(nc_filename);
    end

    % define some other useful arrays
    lat = lat_series_original(event_latspan);
    lon = lon_series_original(event_lonspan);
    phi = lat / 180.0 * 3.1415926;
    lambda = lon / 180.0 * 3.1415926;
    event_timespan = event_day + [- 0.25, 0.00, 0.25]; % convert to 'days since 1990'

    var_names = {'T', 'ua', 'va', 'ug', 'vg', 'omega', 'omega_b'};
    for n = 1 : length(var_names)
        eval([var_names{n}, '_in = zeros(length(lon), length(lat), length(level), length(event_timespan));'])
    end
    for n = 1 : length(var_names)
        for k = 1 : length(level)
            for t = 1 : length(event_timespan)
                if strcmp(var_names{n}, 'omega_b')
                    eval([var_names{n}, '_in(:, :, k, t) = ', var_names{n}, '(:, :, k)'';'])
                else
                    eval([var_names{n}, '_in(:, :, k, t) = ', var_names{n}, '(:, :, k, t)'';'])
                end
            end
        end
        writeNetCDF_v2(nc_filename, var_names{n}, eval([var_names{n}, '_in']), ...
                lat, lon, event_timespan, level);
    end

    %% read in fields for the synoptic map
    box_max_x_s = 21;
    box_max_y_s = 21;
    level_indices_0 = p_start_0 : p_end;
    [event_latspan_s, event_lonspan_s, center_y, ~] = ...
            event_span_Wrapper(lat_series_original, lon_series_original, lat_indices(1), lon_indices(1), ...        
            box_max_y_s, box_max_x_s, box_max_y_s);
    lat_s = lat_series_original(event_latspan_s);
    lon_s = lon_series_original(event_lonspan_s);
    [T_s, ua_s, va_s, ug_s, vg_s, omega_s, omega_b_s, ~, ~]= event_read_Wrapper_uava_q(box_max_y_s, box_max_x_s, length(level_indices_0), ...
            event_year, input_path, event_latspan_s, event_lonspan_s, event_day, level_indices_0, string_1, JJA_DJF);
    for n = 1 : length(var_names)
        eval([var_names{n}, '_s_in = zeros(length(lon_s), length(lat_s), length(level), length(event_timespan));'])
    end
    % define NetCDF file name for synoptic map, and check if an old file exists
    nc_filename_s = [output_path, 'example_event_synoptic_map.nc'];
    if exist(nc_filename_s, 'file') == 2
        delete(nc_filename_s);
    end
    % write into the file for the synoptic map
    for n = 1 : length(var_names)
        for k = 1 : length(level)
            for t = 1 : length(event_timespan)
                if strcmp(var_names{n}, 'omega_b')
                    eval([var_names{n}, '_s_in(:, :, k, t) = ', var_names{n}, '_s(:, :, k)'';'])
                else
                	eval([var_names{n}, '_s_in(:, :, k, t) = ', var_names{n}, '_s(:, :, k, t)'';'])
                end
            end
        end
        writeNetCDF_v2(nc_filename_s, var_names{n}, eval([var_names{n}, '_s_in']), ...
                lat_s, lon_s, event_timespan, level);
    end

    %% prepare precipitation time series and fields
	event_precip_full = precip_with_mask(n_precip, event_lonspan, event_latspan, historical_files, MASK, h_max, ocean_only, p_mask);
    event_precip = event_precip_full(:, :, 1:max(event_timespan)/0.25+100);
    save([output_path, 'event_precip.mat'], 'event_precip')

    % get surface pressure for the synoptic plot
    historical_files_ps = {[data_path, ...
        'b.e11.B20TRC5CNBDRD.f09_g16.', num2str(n_ensemble, '%.3d'), '.cam.h2.PS.1990010100Z-2005123118Z.nc']};
    PS_full = ncread(historical_files_ps{1}, 'PS');
    PS = PS_full(event_lonspan_s, event_latspan_s, event_timespan/0.25);
    PS_mean = mean(PS_full(event_lonspan_s, event_latspan_s, :), 3);
    PS = PS - PS_mean; % get anomaly
    save([output_path, 'example_event.mat'], 'PS', 'event_precip', ...
            'lat_series_original', 'lon_series_original', 'event_lonspan', 'event_latspan')


