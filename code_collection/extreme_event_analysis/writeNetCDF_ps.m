function writeNetCDF_ps(nc_filename, varname, var, lat, lon)

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
    
    % write variable field
    
    var_id = netcdf.defVar(ncid, varname, 'double', [dimlon_id, dimlat_id]);
    netcdf.defVarFill(ncid, var_id, false, fillValue);
    if ~isempty(var(isnan(var)))
        disp('Warning: NaN in the data, revert to fillValue')
        var(isnan(var)) = fillValue;
    end
    netcdf.putVar(ncid, var_id, var);
    
    % close netcdf file
    
    netcdf.close(ncid);

end

