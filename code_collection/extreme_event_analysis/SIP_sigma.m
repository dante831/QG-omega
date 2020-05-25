function omega = SIP_sigma(r, f0, Ni, Nj, Nk, phi, lambda, level, sigma, dphi, dlambda, rhs, omega_b)

    %{
    Ni = 41;
    Nj = 41;
    Nk = 41;
    lx = 1.0; % domain size
    ly = 1.0;
    lz = 1.0;
    dx = lx / (Ni - 1);
    dy = ly / (Nj - 1);
    dz = lz / (Nk - 1);
    %}

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

    %X = linspace(0.0, 1.0, Ni);
    %Y = linspace(0.0, 1.0, Nj);
    %Z = linspace(0.0, 1.0, Nk);
    %dxr2 = -1.0 / dx ^ 2;
    %dyr2 = -1.0 / dy ^ 2;
    %dzr2 = -1.0 / dz ^ 2;
    %apc = -2.0 * (dxr2 + dyr2 + dzr2);
    
    % set up boundary values where non-zero
    % --- note on Mar 12: 
    % The bottom and upper boundary is set to zero, whereas lateral boundaries 
    % are set to be the horizontal average of omega
    
    mean_omega_b = reshape(mean(mean(omega_b, 1), 2), [Nk, 1]);
    
    % top and bottom boundary
    for j = 1 : Nj
        for i = 1 : Ni
            ijk = Lk(1) + Li(i) + j;
            T(ijk) = 0;
            %T(ijk) = omega_b(i, j, 1);
            ijk = Lk(Nk) + Li(i) + j;
            T(ijk) = 0;
            %T(ijk) = omega_b(i, j, Nk);
        end
    end

    % north and south boundary
    for k = 1 : Nk
        for i = 1 : Ni
            ijk = Lk(k) + Li(i) + 1;
            %T(ijk) = 0;

            %T(ijk) = omega_b(i, 1, k);
            T(ijk) = omega_b(1, i, k);

            ijk = Lk(k) + Li(i) + Nj;
            %T(ijk) = 0;

            %T(ijk) = omega_b(i, Nj, k);
            T(ijk) = omega_b(Nj, i, k);

        end
    end

    % east and west boundary
    for k = 1 : Nk
        for j = 1 : Nj
            ijk = Lk(k) + Li(1) + j;
            %T(ijk) = 0;

            %T(ijk) = omega_b(1, j, k);
            T(ijk) = omega_b(j, 1, k);

            ijk = Lk(k) + Li(Ni) + j;
            %T(ijk) = 0;

            %T(ijk) = omega_b(Ni, j, k);
            T(ijk) = omega_b(j, Ni, k);

        end
    end
    
    for k = 2 : Nkm
        del2_sigma = spherical_laplacian(sigma(:, :, k), phi, dphi, dlambda);
        for j = 2 : Njm
            for i = 2 : Nim
                ijk = Lk(k) + Li(i) + j;
                dp2 = level(k + 1) - level(k);
                dp1 = level(k) - level(k - 1);
                %AP(ijk) = apc;
                %AP(ijk) = - 2 / (r^2 * cos(phi(j))^2 * dlambda^2) * sigma(k) ...
                %          - 2 / (r^2 * dphi^2) * sigma(k) ...  
                %          - 2 * f0^2 / (dp1 * dp2);
                AP(ijk) = - 2 / (r^2 * cos(phi(j))^2 * dlambda^2) * sigma(j, i, k) ...
                          - 2 / (r^2 * dphi^2) * sigma(j, i, k) ... 
                          - 2 * f0^2 / (dp1 * dp2) + ...
                          + del2_sigma(j, i);
            
                %AW(ijk) = dxr2;
                AW(ijk) = (1 / (r^2 * cos(phi(j))^2 * dlambda^2)) * sigma(j, i, k);
                %AS(ijk) = dyr2;
                AS(ijk) = (1 / (r^2 * dphi^2) + tan(phi(j)) / (2 * r^2 * dphi)) * sigma(j, i, k);
                %AN(ijk) = dyr2; 
                AN(ijk) = (1 / (r^2 * dphi^2) - tan(phi(j)) / (2 * r^2 * dphi)) * sigma(j, i, k);
                %AE(ijk) = dxr2; 
                AE(ijk) = AW(ijk);
                %AB(ijk) = dzr2;
                AB(ijk) = 2 * f0^2 / (dp1 * (dp1 + dp2));    
                %AT(ijk) = dzr2; 
                AT(ijk) = 2 * f0^2 / (dp2 * (dp1 + dp2));
                %Q(ijk) = 0.0;
                Q(ijk) = rhs(j, i, k);
            end
        end
    end

    % implement boundary conditions (dirichlet)

    % south and north boundary
    for k = 2 : Nkm
        for i = 2 : Nim
            ijk = Lk(k) + Li(i) + 2;
            Q(ijk) = Q(ijk) - AS(ijk) * T(ijk - 1);
            AS(ijk) = 0.0;
        
            ijk = Lk(k) + Li(i) + Njm;
            Q(ijk) = Q(ijk) - AN(ijk) * T(ijk + 1);
            AN(ijk) = 0.0;
        end
    end

    % west and east boundary

    for k = 2 : Nkm
        for j = 2 : Njm
            ijk = Lk(k) + Li(2) + j;
            Q(ijk) = Q(ijk) - AW(ijk) * T(ijk - Nj);
            AW(ijk) = 0.0;
        
            ijk = Lk(k) + Li(Nim) + j;
            Q(ijk) = Q(ijk) - AE(ijk) * T(ijk + Nj);
            AE(ijk) = 0.0;
        end
    end

    % upper and lower boundary

    for i = 2 : Nim
        for j = 2 : Njm
            ijk = Lk(2) + Li(i) + j;
            Q(ijk) = Q(ijk) - AB(ijk) * T(ijk - Nij);
            AB(ijk) = 0.0;

            ijk = Lk(Nkm) + Li(i) + j;
            Q(ijk) = Q(ijk) - AT(ijk) * T(ijk + Nij);
            AT(ijk) = 0.0;
        end
    end

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
            omega = zeros(size(rhs));
            for k = 1 : Nk
                for j = 1 : Nj
                    for i = 1 : Ni
                        ijk = Lk(k) + Li(i) + j;
                        omega(j, i, k) = T(ijk);
                    end
                end 
            end
            return
        end
    end
        
