    
    
    %years = 1991:1995;
    %years = 1996:2000;
    %years = 2071:2075; 
    years = 2076:2080; 
    
    latmax = 90; % degrees
    latmin = -90;
    lonmax = 360; % degrees
    lonmin = 0;
    plevels = [100000, 97500, 95000, 92500, 90000, 87500, 85000, 82500, 80000, 77500, 75000, 70000, 65000, ...
                60000, 55000, 50000, 45000, 40000, 35000, 30000, 25000, 22500, 20000, 17500, 15000, 12500, ...
                10000,  7000,  5000,  3000,  2000,  1000,   700,   500,   300,   200,   100];
    input_path = '/archive1/ziweili/CESM_LENS/data/';
    output_path = '/archive1/ziweili/CESM_LENS/output/';
    n_ensemble = 5;

    for f = 1 : length(years)

        analysis(n_ensemble, years(f), output_path, input_path, latmax, latmin, lonmax, lonmin, plevels);

    end
