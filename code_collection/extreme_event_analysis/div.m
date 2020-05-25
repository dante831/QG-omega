function div_A = div(Ax, Ay, Phi, dphi, dlambda)

    r = 6371000.0;

    div_A = 1 ./ (r * cos(Phi)) .* d_dlambda(Ax, dlambda)...
            + 1 ./ (r * cos(Phi)) .* d_dphi(Ay .* cos(Phi), dphi);

    return

%{
% test case one, Cartesian coordinate
load wind
k = 4;
Lambda = x(:,:,k); 
Phi = y(:,:,k); 
ug = u(:,:,k); 
vg = v(:,:,k); 
div_1 = divergence(Lambda, Phi, ug, vg);
contourf(Lambda, Phi, div_1, 100, 'Linestyle', 'None'); 
shading interp
hold on
title('divergence function in Matlab')
%quiver(Lambda, Phi, ug, vg, 'y');
hold off
colormap('copper');
colorbar

div_2 = div(ug, vg, zeros(size(y(:, :, k))), y(2, 1, k) - y(1, 1, k), x(1, 2, k) - x(1, 1, k)) * r;
figure
contourf(Lambda, Phi, div_2, 100, 'Linestyle', 'None'); 
shading interp
hold on
title('Custom-defined div()')
%quiver(Lambda, Phi, ug, vg, 'y');
hold off
colormap('copper');
colorbar
%}

%{
% test case two, spherical coordinate

lambda = [1.9417, 1.9853, 2.0289, 2.0726, 2.1162, 2.1598, 2.2035, 2.2471, 2.2907, 2.3344, 2.3780];
phi =    [0.6807, 0.7156, 0.7505, 0.7854, 0.8203, 0.8552, 0.8901, 0.9250, 0.9599, 0.9948, 1.0297];
[Lambda, Phi] = meshgrid(lambda, phi);
% the following field should be divergence-free in spherical coordinates
temp_ug = Lambda.*Phi;
temp_vg = -1/2 ./ cos(Phi) .* Phi.^2;
div_3 = div(temp_ug, temp_vg, Phi, dphi, dlambda);
figure
contourf(Lambda, Phi, div_3, 100, 'Linestyle', 'None');
hold on
quiver(Lambda, Phi, temp_ug, temp_vg, 'y')
title('Spherically non-divergent velocity field')
hold off
colormap('copper');
colorbar

figure
div_4 = div(temp_ug, temp_vg, zeros(size(Phi)), dphi, dlambda);
contourf(Lambda, Phi, div_4, 100, 'Linestyle', 'None');
hold on
quiver(Lambda, Phi, temp_ug, temp_vg, 'y')
title('Cartisian case')
hold off
colormap('copper');
colorbar

%}

