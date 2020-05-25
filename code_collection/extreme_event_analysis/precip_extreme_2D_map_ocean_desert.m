
    % compute prescribed percentile of precipitation extremes
    % final edition, before submission. edited on May. 22nd, 2019 -- Ziwei
    % combined the 6-hourly and daily precipitation extremes for CESM and GFDL

    % add functions from ../../source/
    addpath('../../source/');

    % get information from info.m
    info

    % some constants
    pwd_str = pwd;
    
    latmax = 70; % degrees
    latmin = -70;
    lonmax = 360; % degrees
    lonmin = 0;
    ocean_only = false;
    MASK = true;
    
    % determine daily precip extremes or 6-hourly
    if ~exist('DAILY')
        if strfind(pwd_str, 'daily')
            DAILY = true;
        else
            DAILY = false;
        end
    end
 
    % get model-specific values
    if strfind(pwd_str, 'CESM')
        if DAILY
            quantile_file_h = '/net/aimsir/archive1/ziweili/CESM_LENS/output/precip_99.9th_quantile_h_daily.mat';
            quantile_file_r = '/net/aimsir/archive1/ziweili/CESM_LENS/output/precip_99.9th_quantile_r_daily.mat';
            % set number of precip points
            n = 4; % n = 1 for 6hourly, n = 4 for daily in CESM-LENS data
        else
            quantile_file_h = '/net/aimsir/archive1/ziweili/CESM_LENS/output/precip_99.9th_quantile_h.mat';
            quantile_file_r = '/net/aimsir/archive1/ziweili/CESM_LENS/output/precip_99.9th_quantile_r.mat';
            n = 1;
        end
        % determine the maximum and minimum of event box
        box_min = 15;
        box_max = 29;
        data_path = '/net/aimsir/archive1/ziweili/CESM_LENS/data/';
        input_path = '/net/aimsir/archive1/ziweili/CESM_LENS/output/';
        historical_files = {[data_path, ...
            'b.e11.B20TRC5CNBDRD.f09_g16.', num2str(n_ensemble, '%.3d'), '.cam.h2.PRECT.1990010100Z-2005123118Z.nc']};
        rcp85_files = {[data_path, ...
            'b.e11.BRCP85C5CNBDRD.f09_g16.', num2str(n_ensemble, '%.3d'), '.cam.h2.PRECT.2071010100Z-2080123118Z.nc']};
    else
        if DAILY
            n = 8; 
        else
            n = 2;
        end
        box_min = 9;
        box_max = 15;
        datapath = '/net/aimsir/archive1/ziweili/test1_GFDL/data/';
        input_path = datapath;
        historical_files = {...
                [datapath, 'pr_3hr_GFDL-CM3_historical_r1i1p1_1980010100-1984123123.nc'], ...
                [datapath, 'pr_3hr_GFDL-CM3_historical_r1i1p1_1985010100-1989123123.nc'], ...
                [datapath, 'pr_3hr_GFDL-CM3_historical_r1i1p1_1990010100-1994123123.nc'], ...
                [datapath, 'pr_3hr_GFDL-CM3_historical_r1i1p1_1995010100-1999123123.nc']};
        rcp85_files = {...
                [datapath, 'pr_3hr_GFDL-CM3_rcp85_r1i1p1_2081010100-2085123123.nc'], ...
                [datapath, 'pr_3hr_GFDL-CM3_rcp85_r1i1p1_2086010100-2090123123.nc'], ...
                [datapath, 'pr_3hr_GFDL-CM3_rcp85_r1i1p1_2091010100-2095123123.nc'], ...
                [datapath, 'pr_3hr_GFDL-CM3_rcp85_r1i1p1_2096010100-2100123123.nc']};
    end

    h_max = 500;
    p_min = 10000; % upper limit
    p_max = 55000; % lower limit, events should at least have full non-NaN values at this level
    p_mask = 55000; % lower limit for masking
    NETCDF = true; % whether output event omega to NetCDF files
    rng(831); % set random seed for NetCDF file generation
    NETCDF_interval = 200;

    percentile = 0.999;
    % if the percentile is set to 99.5 as indicated by the folder name
    if strfind(pwd_str, '99.5') 
        percentile = 0.995;
    end

    JJA_DJF = 0;
    ACCU_SIGMA = true;
    TRADITIONAL_ADV = false;
    FIND_CENTER = false;
    OCEAN_DESERT = true;
    
    %set and create paths
    output_path = './event_output/';
    plot_path = './event_plots/';
    if ~exist(output_path)
        mkdir(output_path)
    end
    if ~exist(plot_path)
        mkdir(plot_path)
    end

    lat_series_original = ncread(historical_files{1}, 'lat');
    lon_series_original = ncread(historical_files{1}, 'lon');
    [lat_indices, lon_indices] = latlonindices(lat_series_original, lon_series_original, latmin, latmax, lonmin, lonmax);
    
    % get index matrix for ocean deserts
	if OCEAN_DESERT
        [X, Y] = meshgrid(lon_series_original(lon_indices), lat_series_original(lat_indices));
		tags_ocean_desert = zeros(size(X));
		% south Pacific
		tags_ocean_desert(X >= 240 & X <= 290 & Y > -30 & Y < 0) = 1;
		% south Atlantic
		tags_ocean_desert((X >= 330 & Y > -30 & Y < 0) | ...
						  (X <=  10 & Y > -30 & Y < 0)) = 1;
	end
    
    if HISTORICAL
    
        disp('calculate 99.9th percentile of historical precipitation');
        if strfind(pwd_str, 'CESM')
            historical_precip = precip_with_mask(n, lon_indices, lat_indices, historical_files, MASK, h_max, ocean_only, p_mask);
            temp = load(quantile_file_h, 'Q');
            historical_quantile = temp.Q(lon_indices, lat_indices); clear('temp');
        else
            [historical_precip, historical_quantile] = percentile_grid_with_mask(n, lon_indices, lat_indices, ...
                    historical_files, percentile, MASK, h_max, ocean_only, p_mask);
        end

        [days_historical, event_precip_historical] = deal(cell(length(lat_indices), length(lon_indices)));
        
        % get time of extreme events
        for j = 1 : length(lat_indices)
            for i = 1 : length(lon_indices)

                % get lnogitude and day of one event. days are converted to 'days after 1990'
                days_historical{j, i} = find(historical_precip(i, j, :) > historical_quantile(i, j));
                ind = sub2ind(size(historical_precip), ...
                        repmat(i, length(days_historical{j, i}), 1), ...
                        repmat(j, length(days_historical{j, i}), 1), ...
                        days_historical{j, i});
                event_precip_historical{j, i} = historical_precip(ind);

                if OCEAN_DESERT
                    if ~tags_ocean_desert(j, i)
                        days_historical{j, i} = [];
                        event_precip_historical{j, i} = [];
                    end
                end

                % convert to real time values
                if strfind(pwd_str, 'CESM')
                    days_historical{j, i} = days_historical{j, i} * n / 4; % days after 1990
                    % get rid of time boundary values
                    days_select_1 = days_historical{j, i} >= 365 * (2006 - 1990) - 1; % remove the last day
                    days_select_2 = days_historical{j, i} <= 0.25; % remove the first day, 6:00am
                    days_select_3 = days_historical{j, i} >= 365 * (2001 - 1990); % remove years 2001 - 2005
                    days_select_4 = days_historical{j, i} <= 365 * (1991 - 1990) + 0.25; % remove year 1990
                    days_historical{j, i}(days_select_1 | days_select_2 | ...
                                          days_select_3 | days_select_4) = [];
                else
                    days_historical{j, i} = days_historical{j, i} * n / 8; % days after 1980
                    % get rid of time boundary values
                    days_select = days_historical{j, i} <= 0.25; % remove 19800101
                    days_historical{j, i}(days_select) = [];
                end
            end
        end
        clear('historical_precip');
        
        % compute historical extreme precipitation events
        matfilename = 'precip_events_historical.mat';
        event_analysis_v2(box_min, box_max, p_min, p_max, lat_series_original, lon_series_original, lat_indices, lon_indices, ...
                days_historical, event_precip_historical, matfilename, input_path, output_path, NETCDF, NETCDF_interval, n_ensemble, ...
                FIND_CENTER, JJA_DJF, ACCU_SIGMA, TRADITIONAL_ADV, DAILY, ACCU_B);
     
    end
    
    if RCP85
        
        disp('calculate 99.9th percentile of rcp85 precipitation');
        if strfind(pwd_str, 'CESM')
            rcp85_precip = precip_with_mask(n, lon_indices, lat_indices, rcp85_files, MASK, h_max, ocean_only, p_mask);
            temp = load(quantile_file_r, 'Q');
            rcp85_quantile = temp.Q(lon_indices, lat_indices); clear('temp');
        else
            [rcp85_precip, rcp85_quantile] = percentile_grid_with_mask(n, lon_indices, lat_indices, ...
                    rcp85_files, percentile, MASK, h_max, ocean_only, p_mask);
        end
        
        [days_rcp85, event_precip_rcp85] = deal(cell(length(lat_indices), length(lon_indices)));
        
        % get time of extreme events
        for j = 1 : length(lat_indices)
            for i = 1 : length(lon_indices)
                % get lnogitude and day of one event.
                days_rcp85{j, i} = find(rcp85_precip(i, j, :) > rcp85_quantile(i, j));
                ind = sub2ind(size(rcp85_precip), ...
                        repmat(i, length(days_rcp85{j, i}), 1), ...
                        repmat(j, length(days_rcp85{j, i}), 1), ...
                        days_rcp85{j, i});
                event_precip_rcp85{j, i} = rcp85_precip(ind);

                if OCEAN_DESERT
                    if ~tags_ocean_desert(j, i)
                        days_rcp85{j, i} = [];
                        event_precip_rcp85{j, i} = [];
                    end
                end

                if strfind(pwd_str, 'CESM')
                    days_rcp85{j, i} = days_rcp85{j, i} * n / 4 + 365 * (2071 - 1990); % convert to days after 1990
                    days_select_1 = days_rcp85{j, i} >= 365 * (2081 - 1990) - 1; % remove boundary, 20801231
                    days_select_2 = days_rcp85{j, i} <= 0.25 + 365 * (2071 - 1990); % remove 20710101 6:00
                    days_rcp85{j, i}(days_select_1 | days_select_2) = [];
                else
                    days_rcp85{j, i} = days_rcp85{j, i} * n / 8 + 36500 + 365; 
                        % convert to days after 1980; this 365 is important because precipitation files starts from 2081
                    days_select = days_rcp85{j, i} >= 365 * 120 + 365 * 5 * (length(rcp85_files) - 4);
                    days_rcp85{j, i}(days_select) = []; % remove days from 21000101 on and 20991231, 
                                                        % the case with 20991231 is that
                                                        % it is time boundary, unable to
                                                        % calculate QG omega
                end
            end
        end
        clear('rcp85_precip');
        
        % compute rcp85 extreme precipitation events
        matfilename = 'precip_events_rcp85.mat';
        event_analysis_v2(box_min, box_max, p_min, p_max, lat_series_original, lon_series_original, lat_indices, lon_indices, ...
                days_rcp85, event_precip_rcp85, matfilename, input_path, output_path, NETCDF, NETCDF_interval, n_ensemble, ...
                FIND_CENTER, JJA_DJF, ACCU_SIGMA, TRADITIONAL_ADV, DAILY, ACCU_B);
        
    end
    
    
    
    
    

