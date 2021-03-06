
% combine several files to obtain a common quantile of a desired extremeness

data_path = '/net/aimsir/archive1/ziweili/CESM_LENS/data/';
output_path = '/archive1/ziweili/CESM_LENS/output/';
percentage = 0.995;

historical_files = {...
                    [data_path, 'b.e11.B20TRC5CNBDRD.f09_g16.001.cam.h2.PRECT.1990010100Z-2005123118Z.nc'], ...
                    [data_path, 'b.e11.B20TRC5CNBDRD.f09_g16.002.cam.h2.PRECT.1990010100Z-2005123118Z.nc'], ...
                    [data_path, 'b.e11.B20TRC5CNBDRD.f09_g16.003.cam.h2.PRECT.1990010100Z-2005123118Z.nc'], ...
                    [data_path, 'b.e11.B20TRC5CNBDRD.f09_g16.004.cam.h2.PRECT.1990010100Z-2005123118Z.nc'], ...
                    [data_path, 'b.e11.B20TRC5CNBDRD.f09_g16.005.cam.h2.PRECT.1990010100Z-2005123118Z.nc'], ...
                    [data_path, 'b.e11.B20TRC5CNBDRD.f09_g16.035.cam.h2.PRECT.1990010100Z-2005123118Z.nc']};
rcp85_files = {...
                    [data_path, 'b.e11.BRCP85C5CNBDRD.f09_g16.001.cam.h2.PRECT.2071010100Z-2080123118Z.nc'], ...
                    [data_path, 'b.e11.BRCP85C5CNBDRD.f09_g16.002.cam.h2.PRECT.2071010100Z-2080123118Z.nc'], ...
                    [data_path, 'b.e11.BRCP85C5CNBDRD.f09_g16.003.cam.h2.PRECT.2071010100Z-2080123118Z.nc'], ...
                    [data_path, 'b.e11.BRCP85C5CNBDRD.f09_g16.004.cam.h2.PRECT.2071010100Z-2080123118Z.nc'], ...
                    [data_path, 'b.e11.BRCP85C5CNBDRD.f09_g16.005.cam.h2.PRECT.2071010100Z-2080123118Z.nc'], ...
                    [data_path, 'b.e11.BRCP85C5CNBDRD.f09_g16.035.cam.h2.PRECT.2071010100Z-2080123118Z.nc']};

months = [12, 1, 2];

files = historical_files;
years = 1991:2000;
start_year = 1990;
output_name = [output_path, 'precip_99.5th_quantile_h_DJF.mat'];
multifile_quantile_month(files, years, percentage, output_name, start_year, months);

files = rcp85_files;
years = 2071:2080;
start_year = 2071;
output_name = [output_path, 'precip_99.5th_quantile_r_DJF.mat'];
multifile_quantile_month(files, years, percentage, output_name, start_year, months);


