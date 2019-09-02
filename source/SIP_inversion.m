function omega = SIP_diagnostic_v1(r, f0, Ni, Nj, Nk, phi, lambda, level, sigma, dphi, dlambda, rhs, omega_b, LOWER_BOUNDARY)
    
    % SIP method adapted from: 
    % ftp://ftp.springer.de/pub/technik/peric/solvers/lapl3d.f
    %%% some parameters %%%

    Nijk = Ni * Nj * Nk;
    % SIP parameters
    maxit = 300; % changed from 150 to 300 to increase accuracy
    resmax = 1e-14;
    alpha = 0.93;
    %%%%%%%%%%%%%%%%%%%%%%%

    % initialize the arrays and matrixes

    Li = zeros(1, Ni);
    Lk = zeros(1, Nk);

    [AE, AW, AN, AS, AT, AB, AP, Q] = deal(zeros(1, Nijk));
    T = zeros(1, Nijk);
    [UE, UN, UT, RES, LB, LS, LW, LPR] = deal(zeros(1, Nijk)); 

    Nim = Ni - 1;
    Njm = Nj - 1;
    Nkm = Nk - 1;
    Nij = Ni * Nj;

    for i = 1 : Ni
        Li(i) = (i - 1) * Nj;
    end
    for k = 1 : Nk
        Lk(k) = (k - 1) * Nij;
    end

    IJK = repmat(reshape(Lk  ,  1,  1, Nk), Ni, Nj,  1) + ...
          repmat(reshape(Li  , Ni,  1,  1),  1, Nj, Nk) + ...
          repmat(reshape(1:Nj,  1, Nj,  1), Ni,  1, Nk);
    %{
    ijk_2 = zeros(Ni, Nj, Nk);
    for k = 1 : Nk
        for j = 1 : Nj
            for i = 1 : Ni
                ijk_2(i, j, k) = Lk(k) + Li(i) + j;
            end
        end
    end
    %}

    % set up boundary values where non-zero
    % --- note on Mar 12: 
    % The lower and upper boundary is set to zero, whereas lateral boundaries 
    % are set to be the horizontal average of omega
    % --- note on 2017 Jan 03:
    % According to Paul, the horizontal boundary is set to the climatological mean
    % --- note on 2019 Apr 11:
    % introduce LOWER_BOUNDARY, and setting upper boundary to zero
    omega_b_column = zeros(1, Nijk);
    for k = 1 : Nk
        for j = 1 : Nj
            for i = 1 : Ni
                omega_b_column(IJK(i, j, k)) = omega_b(j, i, k);
            end
        end
    end

    %if ACCU_BOUNDARY
    %    T(IJK( 1 : Ni,  1 : Nj, [1, Nk])) = omega_b_column(IJK( 1 : Ni, 1 : Nj,  [1, Nk])); % top and bottom boundary
    %else
    %    T(IJK( 1 : Ni,  1 : Nj, [1, Nk])) = 0; % top and bottom boundary
    %end

    if LOWER_BOUNDARY % if use accurate lower boundary
        T(IJK( 1 : Ni,  1 : Nj, 1)) = omega_b_column(IJK( 1 : Ni, 1 : Nj,  1));
    else % else set lower boundary to zero
        T(IJK( 1 : Ni,  1 : Nj, 1)) = 0;
    end
    T(IJK( 1 : Ni,  1 : Nj, Nk)) = 0; % upper boundary
    T(IJK( 1 : Ni, [1, Nj],  1 : Nk)) = omega_b_column(IJK( 1 : Ni, [1, Nj],  1 : Nk)); % north and south boundary
    T(IJK([1, Ni],  1 : Nj,  1 : Nk)) = omega_b_column(IJK([1, Ni],  1 : Nj,  1 : Nk)); % east and west boundary
    
    for k = 2 : Nkm
        %del2_sigma = spherical_laplacian(sigma(:, :, k), phi, dphi, dlambda);
        dp2 = level(k + 1) - level(k);
        dp1 = level(k) - level(k - 1);
        for j = 2 : Njm
            %for i = 2 : Nim
				ijk = Lk(k) + Li(2 : Nim) + j;
                %AP(ijk) = - 2 / (r^2 * cos(phi(j))^2 * dlambda^2) * sigma(k) ...
                %          - 2 / (r^2 * dphi^2) * sigma(k) ...  
                %          - 2 * f0^2 / (dp1 * dp2);
                AP(ijk) = - 2 / (r^2 * cos(phi(j))^2 * dlambda^2) * sigma(j, 2 : Nim, k) ...
                          - 2 / (r^2 * dphi^2) * sigma(j, 2 : Nim, k) ... 
                          - 2 * f0^2 / (dp1 * dp2); 
                AW(ijk) = (1 / (r^2 * cos(phi(j))^2 * dlambda^2)) * sigma(j, 1 : Nim - 1, k); 
                AS(ijk) = (1 / (r^2 * dphi^2) + tan(phi(j)) / (2 * r^2 * dphi)) * sigma(j - 1, 2 : Nim, k);
                AN(ijk) = (1 / (r^2 * dphi^2) - tan(phi(j)) / (2 * r^2 * dphi)) * sigma(j + 1, 2 : Nim, k);
                AE(ijk) = (1 / (r^2 * cos(phi(j))^2 * dlambda^2)) * sigma(j, 3 : Nim + 1, k);
                AB(ijk) = 2 * f0^2 / (dp1 * (dp1 + dp2));    
                AT(ijk) = 2 * f0^2 / (dp2 * (dp1 + dp2));
                
                Q(ijk) = rhs(j, 2 : Nim, k);
            %end
        end
    end

    % implement boundary conditions (dirichlet)

    % south and north boundary
    temp_ijk = IJK(2 : Nim, 2, 2 : Nkm);
    Q(temp_ijk) = Q(temp_ijk) - AS(temp_ijk) .* T(temp_ijk - 1);
    AS(temp_ijk) = 0.0;
    temp_ijk = IJK(2 : Nim, Njm, 2 : Nkm);
    Q(temp_ijk) = Q(temp_ijk) - AN(temp_ijk) .* T(temp_ijk + 1);
    AN(temp_ijk) = 0.0;

    % west and east boundary
    temp_ijk = IJK(2, 2 : Njm, 2 : Nkm);
    Q(temp_ijk) = Q(temp_ijk) - AW(temp_ijk) .* T(temp_ijk - Nj);
    AW(temp_ijk) = 0.0;
    temp_ijk = IJK(Nim, 2 : Njm, 2 : Nkm);
    Q(temp_ijk) = Q(temp_ijk) - AE(temp_ijk) .* T(temp_ijk + Nj);
    AE(temp_ijk) = 0.0;

    % upper and lower boundary
    temp_ijk = IJK(2 : Nim, 2 : Njm, 2);
    Q(temp_ijk) = Q(temp_ijk) - AB(temp_ijk) .* T(temp_ijk - Nij);
    AB(temp_ijk) = 0.0;
    temp_ijk = IJK(2 : Nim, 2 : Njm, Nkm);
    Q(temp_ijk) = Q(temp_ijk) - AT(temp_ijk) .* T(temp_ijk + Nij);
    AT(temp_ijk) = 0.0;
       
    for k = 2 : Nkm
        for i = 2 : Nim
            for j = 2 : Njm
                ijk = Lk(k) + Li(i) + j;
                LB(ijk) = AB(ijk) / (1.0 + alpha * (UN(ijk - Nij) + UE(ijk - Nij)));
                LW(ijk) = AW(ijk) / (1.0 + alpha * (UN(ijk - Nj) + UT(ijk - Nj)));
                LS(ijk) = AS(ijk) / (1.0 + alpha * (UE(ijk - 1) + UT(ijk - 1)));
                p1 = alpha * (LB(ijk) * UN(ijk - Nij) + LW(ijk) * UN(ijk - Nj));
                p2 = alpha * (LB(ijk) * UE(ijk - Nij) + LS(ijk) * UE(ijk - 1));
                p3 = alpha * (LW(ijk) * UT(ijk - Nj) + LS(ijk) * UT(ijk - 1));
                LPR(ijk) = 1.0 / (AP(ijk) + p1 + p2 + p3 - LB(ijk) * UT(ijk - Nij) - ...
                        LW(ijk) * UE(ijk - Nj) - LS(ijk) * UN(ijk - 1) + 1.E-20);
                UN(ijk) = (AN(ijk) - p1) * LPR(ijk);
                UE(ijk) = (AE(ijk) - p2) * LPR(ijk);
                UT(ijk) = (AT(ijk) - p3) * LPR(ijk);
            end
        end
    end

    for l = 1 : maxit
        resl = 0.0;
        for k = 2 : Nkm
            for i = 2 : Nim
                for j = 2 : Njm
					ijk = Lk(k) + Li(i) + j;
                    RES(ijk) = Q(ijk) - AE(ijk) * T(ijk + Nj) - AW(ijk) * T(ijk - Nj) - ...
                            AN(ijk) * T(ijk + 1) - AS(ijk) * T(ijk - 1) - ...
                            AT(ijk) * T(ijk + Nij) - AB(ijk) * T(ijk - Nij) - ...
                            AP(ijk) * T(ijk);
                    resl = resl + abs(RES(ijk));
                    RES(ijk) = (RES(ijk) - LB(ijk) * RES(ijk - Nij) - ...
                            LW(ijk) * RES(ijk - Nj) - LS(ijk) * RES(ijk - 1)) * LPR(ijk);
                end
            end
        end
    
        if (l == 1)
            res0 = resl;
        end
        rsm = resl / res0;
    
        for k = Nkm : -1 : 2
            for i = Nim : -1 : 2
                for j = Njm : -1 : 2
                    ijk = Lk(k) + Li(i) + j;
                    RES(ijk) = RES(ijk) - UN(ijk) * RES(ijk + 1) - UE(ijk) * RES(ijk + Nj) - ...
                            UT(ijk) * RES(ijk + Nij);
                    T(ijk) = T(ijk) + RES(ijk);
                end
            end
        end
   
        %disp(['iteration = ', num2str(l), ' resl = ', num2str(resl), ' rsm = ', num2str(rsm)]);
        if (rsm < resmax || l == maxit)
            %disp(['iteration = ', num2str(l)]);
            omega = zeros(size(rhs));
            for k = 1 : Nk
                for j = 1 : Nj
                    for i = 1 : Ni
                        omega(j, i, k) = T(IJK(i, j, k));
                    end
                end 
            end
            return
        end
    end
        
