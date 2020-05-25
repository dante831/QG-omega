% test if e_fold makes sense...
addpath('source/');
dx = 0.2;
dy = 0.24;
kx = 1.1;
ky = 0.95;
x = 0 : dx : 20 * dx;
y = 0 : dy : 22 * dy;

[X, Y] = meshgrid(x, y);

omega = cos(kx * (X - x((end+1)/2))) .* cos(ky * (Y - y((end+1)/2)));

[Le, Lx, Ly] = e_fold(omega, dy, dx, 0, 1);

temp = - 4*del2(omega, dx, dy) ./ omega; % result from del2 should be multiplied by 4

k2 = temp((end + 1) / 2, (end + 1) / 2); 
k2_0 = kx^2 + ky^2;
lambda_x = 2*pi / kx;
lambda_y = 2*pi / ky;
lambda = 2*pi / sqrt(k2);

% one thing is certain: combined wavenumber is larger than either kx or ky, but the combined 
% e-folding distance is larger as well. The wavenumber k^2 corresponds to the real magnitude
% of the horizontal Laplacian, but the largeness of Le^2 (smallness of 1/Le^2) effectively 
% makes the horizontal term smaller. 

% the effect of this mechanism: 
disp(['The horizontal term in Tandon is ', num2str(1/Le^2 / k2), ' times the real value']);



