function lapz = spherical_laplacian(z, phi, dphi, dlambda)

    R = 6371000.0;
    zsize = size(z);
    lapz = zeros(zsize);
    Ni = zsize(2); % longitude
    Nj = zsize(1); % latitude

    % x boundary condition: implement second order accuracy

    %lapz(1, :)  = lapz(1, :)  + 1 ./ (R^2 * cos(phi(:)').^2) .* ...
    %              (2 * z(1, :)  - 5 * z(2, :)      + 4 * z(3, :)      - z(4, :))      / dlambda^2;
    %lapz(Ni, :) = lapz(Ni, :) + 1 ./ (R^2 * cos(phi(:)').^2) .* ...
    %              (2 * z(Ni, :) - 5 * z(Ni - 1, :) + 4 * z(Ni - 2, :) - z(Ni - 3, :)) / dlambda^2;

    lapz(:, 1)  = lapz(:, 1)  + 1 ./ (R^2 * cos(phi(:)).^2) .* ...
                    (2 * z(:, 1)  - 5 * z(:, 2)      + 4 * z(:, 3)      - z(:, 4))      / dlambda^2;
    lapz(:, Ni) = lapz(:, Ni) + 1 ./ (R^2 * cos(phi(:)).^2) .* ...
                    (2 * z(:, Ni) - 5 * z(:, Ni - 1) + 4 * z(:, Ni - 2) - z(:, Ni - 3)) / dlambda^2;

    % y boundary condition (note the first derivative has a opposite sign in Nj against 1)

    %lapz(:, 1)  = lapz(:, 1)  + (1 / R^2) .* (2 * z(:, 1)  - 5 * z(:, 2)      + 4 * z(:, 3)      - z(:, 4))      / dphi^2 - ...
    %              (tan(phi(1))  / R^2) .* ( - 3/2 * z(:, 1)  + 2 * z(:, 2)      - 1/2 * z(:, 3))      / dphi;
    %lapz(:, Nj) = lapz(:, Nj) + (1 / R^2) .* (2 * z(:, Nj) - 5 * z(:, Nj - 1) + 4 * z(:, Nj - 2) - z(:, Nj - 3)) / dphi^2 - ...
    %              (tan(phi(Nj)) / R^2) .* (   3/2 * z(:, Nj) - 2 * z(:, Nj - 1) + 1/2 * z(:, Nj - 2)) / dphi;
    lapz(1, :)  = lapz(1, :)  + (1 / R^2) .* (2 * z(1, :)  - 5 * z(2, :)      + 4 * z(3, :)      - z(4, :))      / dphi^2 - ...
                    (tan(phi(1))  / R^2) .* ( - 3/2 * z(1, :)  + 2 * z(2, :)      - 1/2 * z(3, :))      / dphi;
    lapz(Nj, :) = lapz(Nj, :) + (1 / R^2) .* (2 * z(Nj, :) - 5 * z(Nj - 1, :) + 4 * z(Nj - 2, :) - z(Nj - 3, :)) / dphi^2 - ...
                    (tan(phi(Nj)) / R^2) .* (   3/2 * z(Nj, :) - 2 * z(Nj - 1, :) + 1/2 * z(Nj - 2, :)) / dphi;

    % in x direction

    for i = 2 : Ni - 1
        for j = 1 : Nj
            lapz(j, i) = lapz(j, i) + 1 / (R^2 * cos(phi(j))^2) * (z(j, i - 1) + z(j, i + 1) - 2 * z(j, i)) / dlambda^2;
            %lapz(i, j) = lapz(i, j) + 1 / (R^2 * cos(phi(j))^2) * (z(i - 1, j) + z(i + 1, j) - 2 * z(i, j)) / dlambda^2;

        end
    end
    
    % in y direction

    for i = 1 : Ni
        for j = 2 : Nj - 1
            lapz(j, i) = lapz(j, i) + (1 / R^2) * (z(j - 1, i) + z(j + 1, i) - 2 * z(j, i)) / dphi^2 - ...
                        (tan(phi(j)) / R^2) * (z(j + 1, i) - z(j - 1, i)) / (2 * dphi);
            %lapz(i, j) = lapz(i, j) + (1 / R^2) * (z(i, j - 1) + z(i, j + 1) - 2 * z(i, j)) / dphi^2 - ...
            %            (tan(phi(j)) / R^2) * (z(i, j + 1) - z(i, j - 1)) / (2 * dphi);
        end
    end

    return

    % testing code: z = sin(phi) * cos(lambda)
    %               lapz = - 1 / R^2 * ( 2 * sin(phi) .* cos(lambda) + cos(lambda) .* sin(phi) / cos(phi).^2) 

    %{
    
    R = 6371000.0;
    dphi = 0.75 / 180 * 3.1415926;
    dlambda = 0.75 / 180.0 * 3.1415926;
    lon = ncread('ERAi_pakistan.nc', 'longitude');
    lat_temp = ncread('ERAi_pakistan.nc', 'latitude');
    latmax = 71;
    latmin = 1;
    lat = lat_temp(latmax : -1 : latmin);
    phi = lat / 180.0 * 3.1415926;
    lambda = lon / 180.0 * 3.1415926;
    [Phi, Lambda] = meshgrid(phi, lambda); % [Y, X] = meshgrid(y, x)
    
    lambda = (1:144)*dlambda;
    phi = ((1:90) - 45.5)*2*dphi;
    [Lambda, Phi] = meshgrid(lambda, phi);
    z_test = sin(Phi) .* cos(Lambda);
    lapz_n = spherical_laplacian(z_test, phi, dphi, dlambda);
    lapz_a = - 1 / R^2 * (2 * sin(Phi) .* cos(Lambda) + cos(Lambda) .* sin(Phi) ./ cos(Phi).^2);

    max(max(abs(lapz_n - lapz_a)))

    %}

