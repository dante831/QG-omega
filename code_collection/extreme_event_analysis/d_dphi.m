function output = d_dphi(f, dphi)

    dims = size(f);
    output = nan(size(f));

    % second-order differentiation

    if sum(dims) == max(dims) + 1

        output(2 : end - 1) = (f(3 : end) - f(1 : end - 2)) / (2 * dphi);
        output(1)   = ( - 3/2 * f(1)   + 2 * f(2)       - 1/2 * f(3))      / dphi;
        output(end) = (   3/2 * f(end) - 2 * f(end - 1) + 1/2 * f(end - 2)) / dphi;

    elseif length(dims) == 2

        %output(:, 2 : end - 1) = (f(:, 3 : end) - f(:, 1 : end - 2)) / (2 * dphi);
        %output(:, 1)   = ( - 3/2 * f(:, 1)   + 2 * f(:, 2)       - 1/2 * f(:, 3))      / dphi;
        %output(:, end) = (   3/2 * f(:, end) - 2 * f(:, end - 1) + 1/2 * f(:, end - 2)) / dphi;
        output(2 : end - 1, :) = (f(3 : end, :) - f(1 : end - 2, :)) / (2 * dphi);
        output(1, :)   = ( - 3/2 * f(1, :)   + 2 * f(2, :)       - 1/2 * f(3, :))      / dphi;
        output(end, :) = (   3/2 * f(end, :) - 2 * f(end - 1, :) + 1/2 * f(end - 2, :)) / dphi;

    else

        disp('dimension error!');

    end

    return


