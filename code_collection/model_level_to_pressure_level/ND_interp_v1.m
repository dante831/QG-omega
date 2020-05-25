function output = ND_interp_v1(data, p, plevels, ps, string)

    % multiple dimensionalized interpolation, 
    % with nth dimension as interpolation dimension
    if nargin == 5
        switch(string)
        case{'linear'}
            interp = @interp1;
        case{'pchip'}
            interp = @pchip;
        otherwise
            disp('invalid interpolation method');
            return
        end
    elseif nargin == 4
        interp = @interp1;
    else
        disp('incorrect number of input variables!')
        return
    end
    

    dims = size(data);
    dim_output = zeros(0, 0);
    for h = 1 : length(dims)
        if h == 3
            dim_output = [dim_output, length(plevels)];
        else
            dim_output = [dim_output, dims(h)];
        end
    end

    output = zeros(dim_output);
    if length(dims) == 3
        plevels_2 = zeros(size(plevels));
        new_p    = permute(p, [3, 1, 2]);
        new_data = permute(data, [3, 1, 2]);
        output = zeros(dim_output([3, 1, 2]));
        for j = 1 : dims(2)
            for i = 1 : dims(1)
                %plevels_2(:) = plevels;
                ind = plevels > ps(i, j);
                %plevels_2(plevels_2 > ps(i, j)) = NaN;
                %output(:, i, j) = interp(squeeze(new_p(:, i, j)), ...
                %        squeeze(new_data(:, i, j)), plevels_2);
                output(:, i, j) = interp(new_p(:, i, j), new_data(:, i, j), plevels);
                output(ind, i, j) = NaN;
                %output(i, j, :) = interp(squeeze(p(i, j, :)), ...
                %        squeeze(data(i, j, :)), plevels_2);
            end
        end
        output = permute(output, [2, 3, 1]);
    elseif length(dims) == 4
        for t = 1 : dims(4)
            for j = 1 : dims(2)
                for i = 1 : dims(1)
                    plevels_2 = plevels;
                    plevels_2(plevels_2 > ps(i, j, t)) = NaN;
                    output(i, j, :, t) = interp(squeeze(p(i, j, :, t)), ...
                            squeeze(data(i, j, :, t)), plevels_2);
                end
            end
        end
    else
        disp('data is not 3 dimentional!');
        return
    end

    return
   
