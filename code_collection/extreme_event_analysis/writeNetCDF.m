function writeNetCDF(nc_filename, varname, var, lat, lon, time, plevels)

    % the entire solution can be seen at: 
    % http://code.guillaumemaze.org/tips/howtocreateacleannetcdffilefrommatlabusingthebuiltintoolbox

    fillValue = 1e+20;

    if exist(nc_filename, 'file') == 2
        delete(nc_filename);    
    end
    
    % define some global attributes
    ncid = netcdf.create(nc_filename, 'NETCDF4');
    %varid = netcdf.getConstant('GLOBAL');
    NC_GLOBAL = netcdf.getConstant('NC_GLOBAL');
    netcdf.putAtt(ncid, NC_GLOBAL, 'creation_date', datestr(now));

    % define dimensions

    dims = size(var);
    dimlon_id = netcdf.defDim(ncid, 'longitude', dims(1));
    dimlat_id = netcdf.defDim(ncid, 'latitude', dims(2));
    dimp_id = netcdf.defDim(ncid, 'level', dims(3));
    %if length(time) ~= 1
        dimtime_id = netcdf.defDim(ncid, 'time', length(time));
    %end

    % define longitude axis
    
    lon_id = netcdf.defVar(ncid, 'longitude', 'double', dimlon_id);
    netcdf.putAtt(ncid, lon_id, 'long_name', 'longitude');
    netcdf.putAtt(ncid, lon_id, 'units', 'degrees_east');
    netcdf.putVar(ncid, lon_id, lon);
    
    % define latitude axis
    
    lat_id = netcdf.defVar(ncid, 'latitude', 'double', dimlat_id);
    netcdf.putAtt(ncid, lat_id, 'long_name', 'latitude');
    netcdf.putAtt(ncid, lat_id, 'units', 'degrees_north');
    netcdf.putVar(ncid, lat_id, lat);
    
    % define pressure axis
    
    p_id = netcdf.defVar(ncid, 'level', 'double', dimp_id); 
    netcdf.putAtt(ncid, p_id, 'long_name', 'pressure');
    netcdf.putAtt(ncid, p_id, 'units', 'Pa');
    netcdf.putVar(ncid, p_id, plevels);
    
    % define time axis
    
    %if length(time) ~= 1
        time_id = netcdf.defVar(ncid, 'time', 'int', dimtime_id);
        netcdf.putAtt(ncid, time_id, 'long_name', 'time');
        netcdf.putAtt(ncid, time_id, 'calendar', 'noleap');
        %netcdf.putAtt(ncid,var_ti,'calendar','proleptic_gregorian');
        %netcdf.putAtt(ncid, time_id, 'units', 'days since 1860-01-01 00:00:00');
        netcdf.putVar(ncid, time_id, time);
    %end

    % write variable field
    %if length(time) ~= 1
        var_id = netcdf.defVar(ncid, varname, 'double', [dimlon_id, dimlat_id, dimp_id, dimtime_id]);
    %else
    %    var_id = netcdf.defVar(ncid, varname, 'double', [dimlon_id, dimlat_id, dimp_id]);
    %end
    netcdf.defVarFill(ncid, var_id, false, fillValue);
    if ~isempty(var(isnan(var)))
        disp('Warning: NaN in the data, revert to fillValue')
        var(isnan(var)) = fillValue;
    end
    netcdf.putVar(ncid, var_id, var);
    
    % close netcdf file
    
    netcdf.close(ncid);

end

