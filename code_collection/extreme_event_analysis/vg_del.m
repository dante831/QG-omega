function [vg_delf, ug, vg] = vg_del(phi, z, f, dphi, dlambda)

    % note, the z that goes into this function should have unit of m^2/s^2, 
    % which is 'height' times g (9.8 m/s^2)

    r = 6371000.0;
    sizez = size(z);
    Ni = sizez(2);
    Nj = sizez(1);
    Omega = 7.2921e-5;
    [ug, vg, delf_x, delf_y] = deal(zeros(sizez));

    % x (longitude) boundary condition
    
    %vg(1, :)  = 1 ./ (2 * Omega * sin(phi(:)')) .* ...
    %            ( - 3/2 * z(1, :)  + 2 * z(2, :)      - 1/2 * z(3, :))      ./ (r * cos(phi(:)') * dlambda);
    %vg(Ni, :) = 1 ./ (2 * Omega * sin(phi(:)')) .* ...
    %            (   3/2 * z(Ni, :) - 2 * z(Ni - 1, :) + 1/2 * z(Ni - 2, :)) ./ (r * cos(phi(:)') * dlambda);
    %delf_x(1, :)  = 1 ./ (r * cos(phi(:)')) .* ( - 3/2 * f(1, :)  + 2 * f(2, :)      - 1/2 * f(3, :))      / dlambda;
    %delf_x(Ni, :) = 1 ./ (r * cos(phi(:)')) .* (   3/2 * f(Ni, :) - 2 * f(Ni - 1, :) + 1/2 * f(Ni - 2, :)) / dlambda;
    
    vg(:, 1)  = 1 ./ (2 * Omega * sin(phi(:))) .* ...
                ( - 3/2 * z(:, 1)  + 2 * z(:, 2)      - 1/2 * z(:, 3))      ./ (r * cos(phi(:)) * dlambda);
    vg(:, Ni) = 1 ./ (2 * Omega * sin(phi(:))) .* ...
                (   3/2 * z(:, Ni) - 2 * z(:, Ni - 1) + 1/2 * z(:, Ni - 2)) ./ (r * cos(phi(:)) * dlambda);
    delf_x(:, 1)  = 1 ./ (r * cos(phi(:))) .* ( - 3/2 * f(:, 1)  + 2 * f(:, 2)      - 1/2 * f(:, 3))      / dlambda;
    delf_x(:, Ni) = 1 ./ (r * cos(phi(:))) .* (   3/2 * f(:, Ni) - 2 * f(:, Ni - 1) + 1/2 * f(:, Ni - 2)) / dlambda;

    % y (latitude) boundary condition
    
    %ug(:, 1)  = - 1 / (2 * Omega * sin(phi(1)))  * ...
    %            ( - 3/2 * z(:, 1)  + 2 * z(:, 2)      - 1/2 * z(:, 3))      / (r * dphi);
    %ug(:, Nj) = - 1 / (2 * Omega * sin(phi(Nj))) * ...
    %            (   3/2 * z(:, Nj) - 2 * z(:, Nj - 1) + 1/2 * z(:, Nj - 2)) / (r * dphi);
    %delf_y(:, 1)  = 1 / r * ( - 3/2 * f(:, 1)  + 2 * f(:, 2)      - 1/2 * f(:, 3))      / dphi;
    %delf_y(:, Nj) = 1 / r * (   3/2 * f(:, Nj) - 2 * f(:, Nj - 1) + 1/2 * f(:, Nj - 2)) / dphi;

    ug(1, :)  = - 1 / (2 * Omega * sin(phi(1)))  * ...
                ( - 3/2 * z(1, :)  + 2 * z(2, :)      - 1/2 * z(3, :))      / (r * dphi);
    ug(Nj, :) = - 1 / (2 * Omega * sin(phi(Nj))) * ...
                (   3/2 * z(Nj, :) - 2 * z(Nj - 1, :) + 1/2 * z(Nj - 2, :)) / (r * dphi);
    delf_y(1, :)  = 1 / r * ( - 3/2 * f(1, :)  + 2 * f(2, :)      - 1/2 * f(3, :))      / dphi;
    delf_y(Nj, :) = 1 / r * (   3/2 * f(Nj, :) - 2 * f(Nj - 1, :) + 1/2 * f(Nj - 2, :)) / dphi;
    
    % differentiate in x direction, vg and delf_x
    
    for i = 2 : Ni - 1
        for j = 1 : Nj
            vg(j, i) = vg(j, i) + 1 / (2 * Omega * sin(phi(j))) * (z(j, i + 1) - z(j, i - 1)) / (2 * r * cos(phi(j)) * dlambda);
            delf_x(j, i) = delf_x(j, i) + 1 / (r * cos(phi(j))) * (f(j, i + 1) - f(j, i - 1)) / (2 * dlambda);
        end
    end
    
    
    for i = 1 : Ni
        for j = 2 : Nj - 1
            ug(j, i)     = ug(j, i) - 1 / (2 * Omega * sin(phi(j))) * (z(j + 1, i) - z(j - 1, i)) / (2 * r * dphi);
            delf_y(j, i) = delf_y(j, i) + 1 / r * (f(j + 1, i) - f(j - 1, i)) / (2 * dphi);
        end
    end
    
    vg_delf = ug .* delf_x + vg .* delf_y;

    return
    
    % testing code: z = r * cos(phi) * sin(lambda), 
    % then ug = sin(phi) * sin(lambda), vg = cos(lambda), 
    %    delz = (cos(lambda) ,-sin(phi)*sin(lambda))
    
    %{
    
    Omega = 7.2921e-5;
    r = 6371000.0;
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
    
    
    
    % 1st case:  z = r * cos(phi) * sin(lambda), f = r * exp(phi * lambda)
    % then ug = sin(phi) * sin(lambda) / (2*Omega*sin(Phi)), vg = cos(lambda) / (2*Omega*sin(Phi)),
    %    delf = (1 / (r .* cos(Phi)) .* Phi .* exp(Lambda .* Phi) , 1 / r .* Lambda .* exp(Lambda .* Phi))
    
    f = r * exp(Phi .* Lambda);
    z = r * cos(Phi) .* sin(Lambda);
    vg_delf_a = sin(Phi) .* sin(Lambda) ./ (2*Omega*sin(Phi)) ./ cos(Phi) .* Phi .* exp(Lambda .* Phi) + ...
                cos(Lambda) ./ (2*Omega*sin(Phi)) .* Lambda .* exp(Lambda .* Phi);
    vg_delf_n = vg_del(phi, z, f, dphi, dlambda);
    max(max(abs(vg_delf_a - vg_delf_n)))


    % 2nd case: z = r * exp(Phi) .* cos(Lambda)
    % vg_delf_a = 0

    z = r * exp(Phi) .* cos(Lambda);
    vg_delf_a = - 1 ./ cos(Phi) .* exp(Phi) .* sin(Lambda) .* ( - exp(Phi) .* cos(Lambda)) + ...
                exp(Phi) .* cos(Lambda) .* ( -1 ./ cos(Phi) .* exp(Phi) .* sin(Lambda));
    vg_delf_n = vg_del(phi, z, z, dphi, dlambda);
    max(max(abs(vg_delf_a - vg_delf_n)))



    %}

