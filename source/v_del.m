function v_delf = v_del(phi, u, v, f, dphi, dlambda)

    r = 6371000.0;
    sizef = size(f);
    Ni = sizef(2);
    Nj = sizef(1);
    Omega = 7.2921e-5;
    [delf_x, delf_y] = deal(zeros(sizef));

    % x (longitude) boundary condition
    
    %delf_x(1, :)  = 1 ./ (r * cos(phi(:)')) .* ( - 3/2 * f(1, :)  + 2 * f(2, :)      - 1/2 * f(3, :))      / dlambda;
    %delf_x(Ni, :) = 1 ./ (r * cos(phi(:)')) .* (   3/2 * f(Ni, :) - 2 * f(Ni - 1, :) + 1/2 * f(Ni - 2, :)) / dlambda;
    
    delf_x(:, 1)  = 1 ./ (r * cos(phi(:))) .* ( - 3/2 * f(:, 1)  + 2 * f(:, 2)      - 1/2 * f(:, 3))      / dlambda;
    delf_x(:, Ni) = 1 ./ (r * cos(phi(:))) .* (   3/2 * f(:, Ni) - 2 * f(:, Ni - 1) + 1/2 * f(:, Ni - 2)) / dlambda;

    % y (latitude) boundary condition
    
    %delf_y(:, 1)  = 1 / r * ( - 3/2 * f(:, 1)  + 2 * f(:, 2)      - 1/2 * f(:, 3))      / dphi;
    %delf_y(:, Nj) = 1 / r * (   3/2 * f(:, Nj) - 2 * f(:, Nj - 1) + 1/2 * f(:, Nj - 2)) / dphi;
    delf_y(1, :)  = 1 / r * ( - 3/2 * f(1, :)  + 2 * f(2, :)      - 1/2 * f(3, :))      / dphi;
    delf_y(Nj, :) = 1 / r * (   3/2 * f(Nj, :) - 2 * f(Nj - 1, :) + 1/2 * f(Nj - 2, :)) / dphi;

    % differentiate in x direction, vg and delf_x
    
    for i = 2 : Ni - 1
        for j = 1 : Nj
            delf_x(j, i) = delf_x(j, i) + 1 / (r * cos(phi(j))) * (f(j, i + 1) - f(j, i - 1)) / (2 * dlambda);
        end
    end
    
    
    for i = 1 : Ni
        for j = 2 : Nj - 1
            delf_y(j, i) = delf_y(j, i) + 1 / r * (f(j + 1, i) - f(j - 1, i)) / (2 * dphi);
        end
    end
    
    v_delf = u .* delf_x + v .* delf_y;

    return
    

