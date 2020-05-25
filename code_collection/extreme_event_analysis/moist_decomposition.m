function [latent_h, latent_r, excessive_h, excessive_r, epsilon_h, epsilon_r, ...
        sigma_m_h, sigma_m_r, k2_m_h, k2_m_r, denom_h, denom_r, ...
        d_sigma_m, d_Adv_m, d_k2_l2_m, d_m2_m, d_rec_m, d_rec_m_full, d_J_res, ...
        d_excessive_1, d_excessive_2, d_epsilon, J_res_h, J_res_r] = moist_decomposition(...
            k2_h, k2_r, l2_h, l2_r, m2_h, m2_r, sigma_h, sigma_r, Plevels, F0, ...
            dtheta_dp_ma_h, dtheta_dp_ma_r, omega_h, omega_r, omega_QG_h_rec, omega_QG_r_rec, ...
            J_h, J_r, Adv_h, Adv_r, sigma_star_h, sigma_star_r, k2_star_h, k2_star_r)

Ra = 287.04;
R = 6371000.0;
cpd= 1005.7;
kappa = Ra / cpd;

omega_QG_h_rec =  - 1 ./ (k2_h .* sigma_h + m2_h .* F0.^2) .* (Adv_h + l2_h.*kappa./Plevels.*J_h);
omega_QG_r_rec =  - 1 ./ (k2_r .* sigma_r + m2_r .* F0.^2) .* (Adv_r + l2_r.*kappa./Plevels.*J_r);

% the original choice
% dtheta_dp_ma is defined using omega, and sigma_star is defined using omega_QG
sigma_star_h_omega = - kappa .* cpd ./ Plevels .* dtheta_dp_ma_h;
sigma_star_r_omega = - kappa .* cpd ./ Plevels .* dtheta_dp_ma_r;

sigma_m_h   = sigma_h - sigma_star_h;
sigma_m_r   = sigma_r - sigma_star_r;
latent_h    = - sigma_star_h .* Plevels ./ kappa .* omega_QG_h_rec .* k2_star_h./l2_h;
latent_r    = - sigma_star_r .* Plevels ./ kappa .* omega_QG_r_rec .* k2_star_r./l2_r;

k2_m_h = (k2_h .* sigma_h - k2_star_h .* sigma_star_h) ./ sigma_m_h;
k2_m_r = (k2_r .* sigma_r - k2_star_r .* sigma_star_r) ./ sigma_m_r;
%k2_m_h = k2_h;
%k2_m_r = k2_r;
denom_h     = sigma_m_h .* k2_m_h + F0.^2.*m2_h;
denom_r     = sigma_m_r .* k2_m_r + F0.^2.*m2_r;

%% get a smoother version of moist decomposition
% no need for this because this does not change the noisiness of the moist decomposition
%{
denom_h_mean = nanmean(denom_h, 2);
Denom_h_mean = repmat(denom_h_mean, 1, size(denom_h, 2), 1);
denom_r_mean = nanmean(denom_r, 2);
Denom_r_mean = repmat(denom_r_mean, 1, size(denom_h, 2), 1);
denom_h = Denom_h_mean; denom_r = Denom_r_mean;
%}

excessive_h = - sigma_star_h_omega .* Plevels ./ kappa .* omega_h - latent_h;
excessive_r = - sigma_star_r_omega .* Plevels ./ kappa .* omega_r - latent_r;
epsilon_h   = J_h - excessive_h - latent_h;
epsilon_r   = J_r - excessive_r - latent_r;
J_res_h   = J_h - latent_h;
J_res_r   = J_r - latent_r;
 
omega_QG_h_rec_moist = - (Adv_h + kappa ./ Plevels .* l2_h .* J_res_h) ./ denom_h;
omega_QG_r_rec_moist = - (Adv_r + kappa ./ Plevels .* l2_r .* J_res_r) ./ denom_r;

d_excessive_1   = (sigma_star_r - sigma_star_h) .* (l2_h.*omega_h - k2_h.*omega_QG_h_rec) ./ denom_h;
d_excessive_2   = sigma_star_h .* l2_h.*(omega_r - omega_h + ...
                    k2_h./l2_h.*omega_QG_h_rec - k2_r./l2_r.*omega_QG_r_rec) ./ denom_h;                    
d_epsilon       = - kappa ./ Plevels .* l2_h.*(epsilon_r - epsilon_h) ./ denom_h;

%d_J_res       = - kappa ./ Plevels .* l2_h .* (...
%                           (- Plevels./kappa.*sigma_star_r .* (omega_r - k2_r./l2_r.*omega_QG_r_rec) + epsilon_r) - ...
%                           (- Plevels./kappa.*sigma_star_h .* (omega_h - k2_h./l2_h.*omega_QG_h_rec) + epsilon_h)) ...
%                         ./ denom_h;

%J_res_r = J_r + sigma_star_r .* Plevels ./ kappa .* omega_QG_r_rec .* k2_star_r./l2_r;
%J_res_r = J_r + Plevels./kappa.*sigma_star_r .* ( - omega_r + k2_r./l2_r.*omega_QG_r_rec) + ...
%                  Plevels./kappa.*sigma_star_r .* omega_r;


%d_J_res       = - kappa ./ Plevels .* l2_h .* (...
%                           J_res_r - ...
%                           J_res_h)...
%                         ./ denom_h;

% choice 1: Paul agreed that this is the preferable choice
...%{
d_k2_l2_m       = - ((k2_m_r - k2_m_h) .* sigma_m_h .* omega_QG_h_rec + ...
                    kappa ./ Plevels .* (l2_r - l2_h) .* J_res_h) ./ denom_h;
d_J_res       = - kappa ./ Plevels .* l2_h .* (J_res_r - J_res_h) ./ denom_h;
d_excessive_1   = (sigma_star_r - sigma_star_h) .* (l2_h.*omega_h - k2_h.*omega_QG_h_rec) ./ denom_h;
d_excessive_2   = sigma_star_h .* l2_h.*(omega_r - omega_h + ...
                                                    k2_h./l2_h.*omega_QG_h_rec - k2_r./l2_r.*omega_QG_r_rec) ./ denom_h;
d_epsilon       = - kappa ./ Plevels .* l2_h.*(epsilon_r - epsilon_h) ./ denom_h;
...%}

%{
% choice 2
d_k2_l2_m       = - (k2_m_r - k2_m_h) .* sigma_m_h .* omega_QG_h_rec_moist ./ denom_h;
d_J_res       = - kappa ./ Plevels .* (l2_r .* J_res_r - l2_h .* J_res_h) ./ denom_h;
d_excessive_1   = - kappa ./ Plevels .* cpd .* (dtheta_dp_ma_r - dtheta_dp_ma_h) .* (l2_h.*omega_h - k2_h.*omega_QG_h_rec) ./ denom_h;
d_excessive_2   = - kappa ./ Plevels .* cpd .* dtheta_dp_ma_h .* ((l2_r.*omega_r - k2_r.*omega_QG_r_rec) - ...
                                                                  (l2_h.*omega_h - k2_h.*omega_QG_h_rec)) ./ denom_h;
d_epsilon       = - kappa ./ Plevels .* (l2_r.*epsilon_r - l2_h.*epsilon_h) ./ denom_h;
%}

d_excessive     = d_excessive_1 + d_excessive_2;
d_sigma_m       = - (sigma_m_r - sigma_m_h) .* k2_m_h .* omega_QG_h_rec_moist ./ denom_h;
d_m2_m          = - (m2_r - m2_h) .* F0.^2 .* omega_QG_h_rec_moist ./ denom_h;
d_denom         = d_sigma_m + d_k2_l2_m + d_m2_m;
d_Adv_m         = - (Adv_r - Adv_h) ./ denom_h;
d_omega         = omega_r - omega_h;
d_rec_m_full = d_sigma_m + d_k2_l2_m + d_m2_m + d_Adv_m + d_J_res;
d_rec_m      = d_sigma_m + d_k2_l2_m + d_m2_m + d_Adv_m;


% code backup from ensemble_analysis.m
%{
if strcmp(v_str, '0.0') && J_OVER_OMEGA

    % version A, subtract J from both sides
    sigma_m_h   = sigma_h + kappa ./ Plevels .* J_h ./ omega_QG_h_rec .* l2_h ./ k2_h;
    sigma_m_r   = sigma_r + kappa ./ Plevels .* J_r ./ omega_QG_r_rec .* l2_r ./ k2_r;
    k2_m_h      = (sigma_h .* k2_h + kappa ./ Plevels .* J_h ./ omega_QG_h_rec .* l2_h) ./ sigma_m_h;
    k2_m_r      = (sigma_r .* k2_r + kappa ./ Plevels .* J_r ./ omega_QG_r_rec .* l2_r) ./ sigma_m_r;
    % version A.1, define sigma_m as sigma + kappa/p*J/omega_QG, problem: this method is not stable in k2
    %sigma_m_h   = sigma_h + kappa ./ Plevels .* J_h ./ omega_QG_h_rec;
    %sigma_m_r   = sigma_r + kappa ./ Plevels .* J_r ./ omega_QG_r_rec;
    % version B
    %sigma_m_h   = max(sigma_h + kappa ./ Plevels .* J_h ./ omega_QG_h_rec .* l2_h ./ k2_h, 0);
    %sigma_m_r   = max(sigma_r + kappa ./ Plevels .* J_r ./ omega_QG_r_rec .* l2_r ./ k2_r, 0);
    %k2_m_h = k2_h;
    %k2_m_r = k2_r;
    latent_h = 0;
    latent_r = 0;
    denom_h     = sigma_m_h .* k2_m_h + F0.^2.*m2_h;
    denom_r     = sigma_m_r .* k2_m_r + F0.^2.*m2_r;
elseif strcmp(v_str, '1.0') || strcmp(v_str, '0.0')
    sigma_m_h   = sigma_h + kappa .* cpd ./ Plevels .* dtheta_dp_ma_h .* l2_h ./ k2_h;
    sigma_m_r   = sigma_r + kappa .* cpd ./ Plevels .* dtheta_dp_ma_r .* l2_r ./ k2_r;
    latent_h    = cpd .* dtheta_dp_ma_h .* omega_QG_h_rec;
    latent_r    = cpd .* dtheta_dp_ma_r .* omega_QG_r_rec;
    %sigma_m_r   = sigma_r + kappa .* cpd ./ Plevels .* dtheta_dp_ma_r;
    %sigma_m_h   = sigma_h + kappa .* cpd ./ Plevels .* dtheta_dp_ma_h;
    %latent_h    = cpd .* dtheta_dp_ma_h .* omega_QG_h_rec .* k2_h./l2_h;
    %latent_r    = cpd .* dtheta_dp_ma_r .* omega_QG_r_rec .* k2_r./l2_r;
    denom_h     = sigma_m_h .* k2_h + F0.^2.*m2_h;
    denom_r     = sigma_m_r .* k2_r + F0.^2.*m2_r;
    k2_m_h      = (sigma_h .* k2_h + kappa ./ Plevels .* J_h ./ omega_QG_h_rec .* l2_h) ./ sigma_m_h;
    k2_m_r      = (sigma_r .* k2_r + kappa ./ Plevels .* J_r ./ omega_QG_r_rec .* l2_r) ./ sigma_m_r;
    k2_m_h = k2_h;
    k2_m_r = k2_r;
elseif strcmp(v_str, '1.1') || strcmp(v_str, '0.1')
    %sigma_m_h   = sigma_h - sigma_star;
    sigma_m_h   = sigma_h + kappa .* cpd ./ Plevels .* dtheta_dp_ma_h;
    sigma_m_r   = sigma_r + kappa .* cpd ./ Plevels .* dtheta_dp_ma_r;
    latent_h    = cpd .* dtheta_dp_ma_h .* omega_QG_h_rec .* k2_h./l2_h;
    latent_r    = cpd .* dtheta_dp_ma_r .* omega_QG_r_rec .* k2_r./l2_r;
    denom_h     = sigma_m_h .* k2_h + F0.^2.*m2_h;
    denom_r     = sigma_m_r .* k2_r + F0.^2.*m2_r;
    k2_m_h = k2_h;
    k2_m_r = k2_r;
end
%}


%{
% linear expansion

d_sigma_m       = (sigma_m_r - sigma_m_h) .* k2_m_h .* ...
                    (Adv_h + kappa ./ plot_level .* l2_h .* (epsilon_h + excessive_h)) ./ denom_h.^2;
d_k2_l2_m       = (k2_m_r - k2_m_h) .* sigma_m_h .* ...
                    (Adv_h + kappa ./ plot_level .* l2_h .* (epsilon_h + excessive_h)) ./ denom_h.^2 ...
                    - kappa ./ plot_level .* (l2_r - l2_h) .* excessive_h ./ denom_h;
d_m2_m          = (m2_r - m2_h) .* repmat(F0, [1, 1, length(plevels)]).^2 .* ...
                    (Adv_h + kappa ./ plot_level .* l2_h .* (epsilon_h + excessive_h)) ./ denom_h.^2;
d_denom         = (denom_r - denom_h) .* (Adv_h + kappa ./ plot_level .* l2_h .* (epsilon_h + excessive_h)) ./ denom_h.^2;
d_Adv_m         = - (Adv_r - Adv_h) ./ denom_h;
d_epsilon       = - kappa ./ plot_level .* (l2_r.*epsilon_r - l2_h.*epsilon_h) ./ denom_h;
d_excessive     = - kappa ./ plot_level .* l2_h .* (excessive_r - excessive_h) ./ denom_h;

if strcmp(v_str, '0.1')
    % choice 1: Paul agreed that this is the preferable choice
    d_k2_l2_m       = - ((k2_m_r - k2_m_h) .* sigma_m_h .* omega_QG_h_rec + ...
                        kappa ./ Plevels .* (l2_r - l2_h) .* J_res_h) ./ denom_h;
    d_J_res       = - kappa ./ Plevels .* l2_h .* (J_res_r - J_res_h) ./ denom_h;
    d_excessive_1   = - kappa ./ Plevels .* cpd .* (dtheta_dp_ma_r - dtheta_dp_ma_h) .* (l2_h.*omega_h - k2_h.*omega_QG_h_rec) ./ denom_h;
    d_excessive_2   = - kappa ./ Plevels .* cpd .* dtheta_dp_ma_h .* l2_h.*(omega_r - omega_h + ...
                                                        k2_h./l2_h.*omega_QG_h_rec - k2_r./l2_r.*omega_QG_r_rec) ./ denom_h;

    
    % choice 2
    %d_k2_l2_m       = - (k2_m_r - k2_m_h) .* sigma_m_h .* omega_QG_h_rec ./ denom_h;
    %d_J_res       = - kappa ./ Plevels .* (l2_r .* J_res_r - l2_h .* J_res_h) ./ denom_h;
    %d_excessive_1   = - kappa ./ Plevels .* cpd .* (dtheta_dp_ma_r - dtheta_dp_ma_h) .* (l2_h.*omega_h - k2_h.*omega_QG_h_rec) ./ denom_h;
    %d_excessive_2   = - kappa ./ Plevels .* cpd .* dtheta_dp_ma_h .* ((l2_r.*omega_r - k2_r.*omega_QG_r_rec) - ...
    %                                                                  (l2_h.*omega_h - k2_h.*omega_QG_h_rec)) ./ denom_h;
    

    d_excessive     = d_excessive_1 + d_excessive_2;
    d_sigma_m       = - (sigma_m_r - sigma_m_h) .* k2_m_h .* omega_QG_h_rec ./ denom_h;
    d_m2_m          = - (m2_r - m2_h) .* repmat(F0, [1, 1, length(plevels)]).^2 .* omega_QG_h_rec ./ denom_h;
    %d_denom         = (denom_r - denom_h) .* (Adv_h + kappa ./ Plevels .* l2_h .* J_res_h) ./ denom_h.^2;
    d_denom         = d_sigma_m + d_k2_l2_m + d_m2_m;
    d_Adv_m         = - (Adv_r - Adv_h) ./ denom_h;
    d_epsilon       = - kappa ./ Plevels .* l2_h .* (epsilon_r - epsilon_h) ./ denom_h;

else
    d_excessive_1   = - kappa ./ Plevels .* cpd .* (dtheta_dp_ma_r - dtheta_dp_ma_h) .* l2_h.*(omega_h - omega_QG_h_rec) ./ denom_h;
    d_excessive_2   = - kappa ./ Plevels .* cpd .* l2_h .* ((omega_r - omega_QG_r_rec) - ...
                                                               (omega_h - omega_QG_h_rec)) .* dtheta_dp_ma_h ./ denom_h;
    d_excessive     = d_excessive_1 + d_excessive_2;
    d_J_res       = d_excessive_1 + d_excessive_2 + d_epsilon;
end
d_omega         = omega_r - omega_h;

if J_OVER_OMEGA
    d_sigma_m   = (sigma_m_r - sigma_m_h) .* k2_m_h .* Adv_h ./ denom_h.^2;
    d_k2_l2_m   = sigma_m_h .* (k2_m_r - k2_m_h) .* Adv_h ./ denom_h.^2;
    d_excessive_1 = d_excessive;
    d_excessive_1(:) = 0;
    d_excessive_2(:) = 0;
end
%}

