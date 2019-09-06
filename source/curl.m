function rot_A = curl_spherical(phi, Ax, Ay, dphi, dlambda)

    r = 6371000.0;
    sizeA = size(Ax);
    Ni = sizeA(2);
    Nj = sizeA(1);
    Omega = 7.2921e-5;
    [dAy_dx, dAx_dy] = deal(zeros(sizeA));

    % x (longitude) boundary condition
    dAy_dx(:, 1)  = 1 ./ (r * cos(phi(:))) .* ( - 3/2 * Ay(:, 1)  + 2 * Ay(:, 2)      - 1/2 * Ay(:, 3))      / dlambda;
    dAy_dx(:, Ni) = 1 ./ (r * cos(phi(:))) .* (   3/2 * Ay(:, Ni) - 2 * Ay(:, Ni - 1) + 1/2 * Ay(:, Ni - 2)) / dlambda;

    % y (latitude) boundary condition
    dAx_dy(1, :)  = 1 ./ (r * cos(phi(1))) .* ( - 3/2 * Ax(1, :) *cos(phi(1)) + 2 * Ax(2, :)     *cos(phi(2))      - 1/2 * Ax(3, :)     *cos(phi(3)))      / dphi;
    dAx_dy(Nj, :) = 1 ./ (r * cos(phi(Nj))).* (   3/2 * Ax(Nj, :)*cos(phi(Nj))- 2 * Ax(Nj - 1, :)*cos(phi(Nj - 1)) + 1/2 * Ax(Nj - 2, :)*cos(phi(Nj - 2))) / dphi; 
    
    % differentiate in x direction, vg and delf_x
    
    for i = 2 : Ni - 1
        for j = 1 : Nj
            dAy_dx(j, i) = dAy_dx(j, i) + 1 / (r * cos(phi(j))) * (Ay(j, i + 1) - Ay(j, i - 1)) / (2 * dlambda);
        end
    end
    
    
    for i = 1 : Ni
        for j = 2 : Nj - 1
            dAx_dy(j, i) = dAx_dy(j, i) + 1 / (r * cos(phi(j))) * (Ax(j + 1, i)*cos(phi(j + 1)) - Ax(j - 1, i)*cos(phi(j - 1))) / (2 * dphi);
        end
    end
    
    rot_A = dAy_dx - dAx_dy;

    return

    
%{
% test case one, Cartesian coordinate

r = 6371000.0;
load wind
k = 4;
Lambda = x(:,:,k); 
Phi = y(:,:,k); 
ug = u(:,:,k); 
vg = v(:,:,k); 
cav = curl(Lambda, Phi, ug, vg);
pcolor(Lambda, Phi, cav); 
shading interp
hold on
title('Curl function in Matlab')
quiver(Lambda, Phi, ug, vg, 'y');
hold off
colormap('copper');

rot_A = curl_spherical(ones(size(y(:, 1))), ug, vg, y(2, 1) - y(1, 1), x(1, 2) - x(1, 1)) * r;

figure
pcolor(Lambda, Phi, rot_A); 
shading interp
hold on
title('Custom-defined curl\_spherical')
quiver(Lambda, Phi, ug, vg, 'y');
hold off
colormap('copper');

%}

% test case two, spherical coordinate

%{
dlambda = 0.0436;
dphi = 0.0349;
lambda = [1.9417, 1.9853, 2.0289, 2.0726, 2.1162, 2.1598, 2.2035, 2.2471, 2.2907, 2.3344, 2.3780];
phi =    [0.6807, 0.7156, 0.7505, 0.7854, 0.8203, 0.8552, 0.8901, 0.9250, 0.9599, 0.9948, 1.0297];
[Lambda, Phi] = meshgrid(lambda, phi);
% the following field should be vorticity-free in spherical coordinates
temp_ug = Lambda.*Phi ./ cos(Phi);
temp_vg = 1/2 .* Lambda.^2;
vort_1 = curl_spherical(phi, temp_ug, temp_vg, dphi, dlambda);
figure
contourf(Lambda, Phi, vort_1, 100, 'Linestyle', 'None');
hold on
quiver(Lambda, Phi, temp_ug, temp_vg, 'y')
title('Spherically rotation-free velocity field')
hold off
colormap('copper');
colorbar

figure
vort_2 = curl_spherical(zeros(size(phi)), temp_ug, temp_vg, dphi, dlambda);
contourf(Lambda, Phi, vort_2, 100, 'Linestyle', 'None');
hold on
quiver(Lambda, Phi, temp_ug, temp_vg, 'y')
title('Cartisian case')
hold off
colormap('copper');
colorbar

%}
