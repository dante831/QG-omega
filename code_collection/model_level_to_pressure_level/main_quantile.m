
% combine several files to obtain a common quantile of a desired extremeness

data_path = '/net/aimsir/archive1/ziweili/CESM_LENS/data/';
output_path = '/archive1/ziweili/CESM_LENS/output/';
percentage = 0.999;

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

multifile_quantile(historical_files, 1991:2000, percentage, [output_path, 'precip_99.9th_quantile_h_corrected.mat'], 1990);
multifile_quantile(rcp85_files     , 2071:2080, percentage, [output_path, 'precip_99.9th_quantile_r_corrected.mat'], 2071);


