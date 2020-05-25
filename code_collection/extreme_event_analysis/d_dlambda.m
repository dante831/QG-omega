function output = d_dlambda(f, dlambda)

    dims = size(f);
    output = nan(size(f));
    
    % second-order differentiation

    if sum(dims) == max(dims) + 1
        
        output(2 : end - 1) = (f(3 : end) - f(1 : end - 2)) / (2 * dlambda);
        output(1)   = ( - 3/2 * f(1)   + 2 * f(2)       - 1/2 * f(3))      / dlambda;
        output(end) = (   3/2 * f(end) - 2 * f(end - 1) + 1/2 * f(end - 2)) / dlambda;
    
    elseif length(dims) == 2
        
        %output(2 : end - 1, :) = (f(3 : end, :) - f(1 : end - 2, :)) / (2 * dlambda);
        %output(1, :)   = ( - 3/2 * f(1, :)   + 2 * f(2, :)       - 1/2 * f(3, :))      / dlambda;
        %output(end, :) = (   3/2 * f(end, :) - 2 * f(end - 1, :) + 1/2 * f(end - 2, :)) / dlambda;
        output(:, 2 : end - 1) = (f(:, 3 : end) - f(:, 1 : end - 2)) / (2 * dlambda);
        output(:, 1)   = ( - 3/2 * f(:, 1)   + 2 * f(:, 2)       - 1/2 * f(:, 3)) / dlambda;
        output(:, end) = (   3/2 * f(:, end) - 2 * f(:, end - 1) + 1/2 * f(:, end - 2)) / dlambda;
 
    else

        disp('dimension error!');
    
    end
    
    return
