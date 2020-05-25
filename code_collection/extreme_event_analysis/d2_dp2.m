function d2_dp2z = d2_dp2(z, level)

    sizez = size(z);
    Nk = length(level);
    d2_dp2z = zeros(sizez);
    level = double(level);

    if (length(sizez) == 2) && ((sizez(1) == 1) ||(sizez(2) == 1))

        % upper and lower boundary conditions

        %d_dpz(1)  = (z(2)  - z(1))      / (level(2)  - level(1));
        %d_dpz(Nk) = (z(Nk) - z(Nk - 1)) / (level(Nk) - level(Nk - 1));
        
        % inner points

        for k = 2 : Nk - 1
            dp2 = level(k + 1) - level(k);
            dp1 = level(k) - level(k - 1);
            %d_dpz(k) = z(k) * (dp2 - dp1) / (dp1 * dp2) + z(k + 1) * dp1 / (dp2 * (dp1 + dp2)) - ...
            %           z(k - 1) * dp2 / (dp1 * (dp1 + dp2));
            d2_dp2z(k) = z(k) * (-2) / (dp1 * dp2) + ...
                         z(k - 1) * 2 / (dp1 * (dp1 + dp2)) + ...
                         z(k + 1) * 2 / (dp2 * (dp1 + dp2));
        end

        % upper and lower boundaries
        d2_dp2z(1)  = d2_dp2z(2);
        d2_dp2z(Nk) = d2_dp2z(Nk - 1);

    else if length(sizez) == 3

        % upper and lower boundary conditions

        %d_dpz(: ,:, 1)  = (z(: ,:, 2)  - z(: ,:, 1))      / (level(2)  - level(1));
        %d_dpz(: ,:, Nk) = (z(:, :, Nk) - z(:, :, Nk - 1)) / (level(Nk) - level(Nk - 1));

        % inner points

        for k = 2 : Nk - 1
            dp2 = level(k + 1) - level(k);
            dp1 = level(k) - level(k - 1);
            %d_dpz(: ,:, k) = z(: ,:, k) * (dp2 - dp1) / (dp1 * dp2) + z(: ,:, k + 1) * dp1 / (dp2 * (dp1 + dp2)) - ...
            %                 z(: ,:, k - 1) * dp2 / (dp1 * (dp1 + dp2));
            d2_dp2z(: ,:, k) = z(: ,:, k) * (-2) / (dp1 * dp2) + ...
                               z(: ,:, k - 1) * 2 / (dp1 * (dp1 + dp2)) + ...
                               z(:, :, k + 1) * 2 / (dp2 * (dp1 + dp2));
        end

        % upper and lower boundaries
        d2_dp2z(:, :, 1)  = d2_dp2z(:, :, 2);
        d2_dp2z(:, :, Nk) = d2_dp2z(:, :, Nk - 1);

    else

        disp('error with dimension in function d_dp()');
        
    end

    end

    return

    % testing code: z = - 9.8 * 10000 * ln(p / 1000)

    %{

    level = ncread('ERAi_pakistan.nc', 'level');
    d_dpz_a = - 9.8 * 10000 ./ level;
    d_dpz_n = d_dp(- 9.8 * 10000 * ln(level / 1000), level);

    max(max(abs(d_dpz_a - d_dpz_n)))

    %}


