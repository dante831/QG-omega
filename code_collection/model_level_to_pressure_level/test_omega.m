
days_normal = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
month = 1;
%start = [1, 1, 1, days_normal(1 : (month - 1)) + 1];
start = [1, 1, 1, 1];
count = [Inf, Inf, Inf, days_normal(month) * 4];

file_name_1 = '/archive1/ziweili/CESM_LENS/output/omega_1990.nc';
omega = ncread(file_name_1, 'omega_plevel', start, count);
plevels = ncread(file_name_1, 'level');
omega_month_1 = mean(omega, 4);

start2 = [1, 1, 1, (1990 - 1850) * 12 + month];
count2 = [Inf, Inf, Inf, 1];
%omega_month_2 = ncread('/archive1/ziweili/CESM_LENS/data/b.e11.B20TRC5CNBDRD.f09_g16.001.cam.h0.OMEGA.185001-200512.nc', ...
%            'OMEGA', start2, count2);
omega_month_2 = ncread('/archive1/ziweili/CESM_LENS/output/omega_monthly_mean_1990_02_.nc', 'omega_plevel');

figure
contour(omega_month_1(:, :, 16))
title('omega from continuity')
caxis([-0.5, 0.5])
colorbar

figure
contour(omega_month_2(:, :, 16))
title('omega from model output')
caxis([-0.5, 0.5])
colorbar

figure
contourf(omega_month_2(:, :, 16) ./ omega_month_1(:, :, 16))
caxis([-1., 1.])
colorbar


