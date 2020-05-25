
    % compute prescribed percentile of precipitation extremes
    % final edition, before submission. edited on Jan. 16th, 2019 -- Ziwei

    % add functions from ../../source/
    addpath('../../source/');

    % some constants
    
    latmax = 70; % degrees
    latmin = -70;
    lonmax = 360; % degrees
    lonmin = 0;
    HISTORICAL = false;
    RCP85 = true;
    n = 1; % n = 1 for 6hourly, n = 4 for daily in CESM-LENS data
    ocean_only = false;
    MASK = true;
    box_min = 15;
    box_max = 29;
    h_max = 500;
    p_min = 10000; % upper limit
    p_max = 55000; % lower limit, events should at least have full non-NaN values at this level
    p_mask = 55000; % lower limit for masking
    NETCDF = true; % whether output event omega to NetCDF files
    rng(831); % set random seed for NetCDF file generation
    NETCDF_interval = 200;
    n_ensemble = 1;
    quantile_file_h = '/net/aimsir/archive1/ziweili/CESM_LENS/output/precip_99.9th_quantile_h.mat';
    quantile_file_r = '/net/aimsir/archive1/ziweili/CESM_LENS/output/precip_99.9th_quantile_r.mat';
    percentile = 0.999;
    FIND_CENTER = true;
    JJA_DJF = 0;
    ACCU_SIGMA = true;
    TRADITIONAL_ADV = false;
    
    %set and create paths
    data_path = '/net/aimsir/archive1/ziweili/CESM_LENS/data/';
    input_path = '/net/aimsir/archive1/ziweili/CESM_LENS/output/';
    output_path = './event_output/';
    plot_path = './event_plots/';
    if ~exist(output_path)
        mkdir(output_path)
    end
    if ~exist(plot_path)
        mkdir(plot_path)
    end

    historical_files = {[data_path, ...
        'b.e11.B20TRC5CNBDRD.f09_g16.', num2str(n_ensemble, '%.3d'), '.cam.h2.PRECT.1990010100Z-2005123118Z.nc']};
    rcp85_files = {[data_path, ...
        'b.e11.BRCP85C5CNBDRD.f09_g16.', num2str(n_ensemble, '%.3d'), '.cam.h2.PRECT.2071010100Z-2080123118Z.nc']};

    lat_series_original = ncread(historical_files{1}, 'lat');
    lon_series_original = ncread(historical_files{1}, 'lon');
    [lat_indices, lon_indices] = latlonindices(lat_series_original, lon_series_original, latmin, latmax, lonmin, lonmax);
    if HISTORICAL
    
        disp('calculate 99.9th percentile of historical precipitation');
        historical_precip = precip_with_mask(n, lon_indices, lat_indices, historical_files, MASK, h_max, ocean_only, p_mask);
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
        clear('historical_precip');
        
        % compute historical extreme precipitation events
        matfilename = 'precip_events_historical.mat';
        event_analysis_v1(box_min, box_max, p_min, p_max, lat_series_original, lon_series_original, lat_indices, lon_indices, ...
                days_historical, event_precip_historical, matfilename, input_path, output_path, NETCDF, NETCDF_interval, n_ensemble, ...
                FIND_CENTER, JJA_DJF, ACCU_SIGMA, TRADITIONAL_ADV);
     
    end
    

    
    if RCP85
        
        disp('calculate 99.9th percentile of rcp85 precipitation');
        rcp85_precip = precip_with_mask(n, lon_indices, lat_indices, rcp85_files, MASK, h_max, ocean_only, p_mask);
        temp = load(quantile_file_r, 'Q');
        rcp85_quantile = temp.Q(lon_indices, lat_indices); clear('temp');
        [days_rcp85, event_precip_rcp85] = deal(cell(length(lat_indices), length(lon_indices)));
        
        % calculate QG_omega
    
        for j = 1 : length(lat_indices)
            for i = 1 : length(lon_indices)
                
                % get lnogitude and day of one event. days are converted to 'days after 1990'
                days_rcp85{j, i} = find(rcp85_precip(i, j, :) > rcp85_quantile(i, j));
                ind = sub2ind(size(rcp85_precip), ...
                        repmat(i, length(days_rcp85{j, i}), 1), ...
                        repmat(j, length(days_rcp85{j, i}), 1), ...
                        days_rcp85{j, i});
                event_precip_rcp85{j, i} = rcp85_precip(ind);

                days_rcp85{j, i} = days_rcp85{j, i} * 0.25 + 365 * (2071 - 1990);

                days_select_1 = days_rcp85{j, i} >= 365 * (2081 - 1990) - 1; % remove 20801231
                days_select_2 = days_rcp85{j, i} <= 0.25 + 365 * (2071 - 1990); % remove 20710101 6:00
                days_rcp85{j, i}(days_select_1 | days_select_2) = [];

            end
        
        end
        
        clear('rcp85_precip');
        
        % compute rcp85 extreme precipitation events
        matfilename = 'precip_events_rcp85.mat';
        
        event_analysis_v1(box_min, box_max, p_min, p_max, lat_series_original, lon_series_original, lat_indices, lon_indices, ...
                days_rcp85, event_precip_rcp85, matfilename, input_path, output_path, NETCDF, NETCDF_interval, n_ensemble, ...
                FIND_CENTER, JJA_DJF, ACCU_SIGMA, TRADITIONAL_ADV);
        
    end
    
    
    
    
    
            
            
    

