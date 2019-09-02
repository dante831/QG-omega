function [x_ind, y_ind, computation_x, computation_y, p_start, NaN_tag] = ...
                shrink_box_v1(omega, box_min, p_max, computation_x, computation_y, p_start, p_end, plevels)

    % improved shink box routine. The function starts from the lowest lower boundary, shrink horizontally, and then
    % finds the highest lower boundary that contains no NaN. 
    % update on 11/08/2018: this method works well with the vorticity equation. It still needs to be tested in the
    % full QG omega equation

    %computation_x_0 = computation_x;
    %computation_y_0 = computation_y;
    [computation_y_0, computation_x_0, ~] = size(omega);
    shrink_box = true;
    NaN_tag = false;
    p_start = 1;
    while shrink_box
        x_ind = (computation_x_0 - computation_x) / 2 + 1 : (computation_x_0 + computation_x) / 2;
        y_ind = (computation_y_0 - computation_y) / 2 + 1 : (computation_y_0 + computation_y) / 2;
        temp = omega(y_ind, x_ind, p_start : p_end, :);
        if any(isnan(temp(:))) && (computation_x > box_min && computation_y > box_min)
            computation_x = computation_x - 2;
            computation_y = computation_y - 2;
        elseif any(isnan(temp(:)))
            while any(isnan(temp(:)))
                p_start = p_start + 1;
                if p_start == (find(plevels == p_max) + 1)
                    NaN_tag = true;
                    shrink_box = false;
                    break;
                end
                temp = omega(y_ind, x_ind, p_start : p_end, :);
            end
            shrink_box = false;
        else
            shrink_box = false;
        end    
    end



