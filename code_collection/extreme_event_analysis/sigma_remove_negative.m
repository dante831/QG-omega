function [sigma_out, tag_NaN] = sigma_remove_negative(sigma_accu, sigma)

    sigma_out = sigma_accu;
    alpha = 0.2;
    tag = squeeze(min(min(sigma_accu))) < sigma * alpha;

    if sum(sigma_accu(:) < 0) >= 0.1*length(sigma_accu(:)) % if negative fraction of sigma_accu is higher than 10%, return
        tag_NaN = false;
        return
    end
    s = size(sigma_accu);
    if length(s) == 3 % if only has one time step, set the time dimension to 1
        s(4) = 1;
    end
    while any(tag(:))
        temp_sigma = sigma_out;
        for t = 1 : s(4)
            for k = 1 : s(3)
                [indy, indx] = find(sigma_out(:, :, k, t) < sigma(k, t)*alpha);
                for i = 1 : length(indy)
                    [indy_2, indx_2] = meshgrid(indy(i)-1:indy(i)+1, indx(i)-1:indx(i)+1); % get a 3 by 3 matrix
                    indy_2 = indy_2([1:4, 6:9]); % exclude the center point
                    indx_2 = indx_2([1:4, 6:9]); 
                    remove_ind = indy_2 < 1 | indy_2 > s(1) | indx_2 < 1 | indx_2 > s(2); % if points in the matrix exceeds the event field, neglect them
                    indy_2 = indy_2(~remove_ind);
                    indx_2 = indx_2(~remove_ind);
                    temp_ind = sub2ind(s, indy_2, indx_2, repmat(k, size(indy_2)), repmat(t, size(indy_2)));
                    % pick the mean value of the adjacent points, or sigma(k, t)*alpha, whichever is larger
                    temp_sigma(indy(i), indx(i), k, t) = max([mean(squeeze(sigma_out(temp_ind))), sigma(k, t)*alpha]); 
                end
            end
        end
        sigma_out = temp_sigma;
        tag = squeeze(min(min(sigma_out))) < sigma * alpha;
    end
    
    tag_NaN = true;
    return

