function [latmax, latmin, lonmax, lonmin] = get_lat_lon_limits(path_str)

    % search run file to find latmax and latmin
    run_file = dir([path_str, '/precip_extreme_2D_map*.m']);
    c = textread(run_file.name,'%s','delimiter','\n');
    latmax_string = c{find(~cellfun(@isempty,strfind(c,'latmax = ')))};
    latmax = str2double(regexp(latmax_string, '\d+\.?\d*|-\d+\.?\d*|\.?\d+|-\.?\d+', 'match'));
    latmin_string = c{find(~cellfun(@isempty,strfind(c,'latmin = ')))};
    latmin = str2double(regexp(latmin_string, '\d+\.?\d*|-\d+\.?\d*|\.?\d+|-\.?\d+', 'match'));
    lonmax_string = c{find(~cellfun(@isempty,strfind(c,'lonmax = ')))};
    lonmax = str2double(regexp(lonmax_string, '\d+\.?\d*|-\d+\.?\d*|\.?\d+|-\.?\d+', 'match'));
    lonmin_string = c{find(~cellfun(@isempty,strfind(c,'lonmin = ')))};
    lonmin = str2double(regexp(lonmin_string, '\d+\.?\d*|-\d+\.?\d*|\.?\d+|-\.?\d+', 'match'));

return
