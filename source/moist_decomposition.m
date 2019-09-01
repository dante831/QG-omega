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

% dtheta_dp_ma is defined using omega, and sigma_star is defined using omega_QG
sigma_star_h_omega = - kappa .* cpd ./ Plevels .* dtheta_dp_ma_h;
sigma_star_r_omega = - kappa .* cpd ./ Plevels .* dtheta_dp_ma_r;

sigma_m_h   = sigma_h - sigma_star_h;
sigma_m_r   = sigma_r - sigma_star_r;
latent_h    = - sigma_star_h .* Plevels ./ kappa .* omega_QG_h_rec .* k2_star_h./l2_h;
latent_r    = - sigma_star_r .* Plevels ./ kappa .* omega_QG_r_rec .* k2_star_r./l2_r;

k2_m_h = (k2_h .* sigma_h - k2_star_h .* sigma_star_h) ./ sigma_m_h;
k2_m_r = (k2_r .* sigma_r - k2_star_r .* sigma_star_r) ./ sigma_m_r;
denom_h     = sigma_m_h .* k2_m_h + F0.^2.*m2_h;
denom_r     = sigma_m_r .* k2_m_r + F0.^2.*m2_r;

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

d_k2_l2_m       = - ((k2_m_r - k2_m_h) .* sigma_m_h .* omega_QG_h_rec + ...
                    kappa ./ Plevels .* (l2_r - l2_h) .* J_res_h) ./ denom_h;
d_J_res       = - kappa ./ Plevels .* l2_h .* (J_res_r - J_res_h) ./ denom_h;
d_excessive_1   = (sigma_star_r - sigma_star_h) .* (l2_h.*omega_h - k2_h.*omega_QG_h_rec) ./ denom_h;
d_excessive_2   = sigma_star_h .* l2_h.*(omega_r - omega_h + ...
                                                    k2_h./l2_h.*omega_QG_h_rec - k2_r./l2_r.*omega_QG_r_rec) ./ denom_h;
d_epsilon       = - kappa ./ Plevels .* l2_h.*(epsilon_r - epsilon_h) ./ denom_h;
d_excessive     = d_excessive_1 + d_excessive_2;
d_sigma_m       = - (sigma_m_r - sigma_m_h) .* k2_m_h .* omega_QG_h_rec_moist ./ denom_h;
d_m2_m          = - (m2_r - m2_h) .* F0.^2 .* omega_QG_h_rec_moist ./ denom_h;
d_denom         = d_sigma_m + d_k2_l2_m + d_m2_m;
d_Adv_m         = - (Adv_r - Adv_h) ./ denom_h;
d_omega         = omega_r - omega_h;
d_rec_m_full = d_sigma_m + d_k2_l2_m + d_m2_m + d_Adv_m + d_J_res;
d_rec_m      = d_sigma_m + d_k2_l2_m + d_m2_m + d_Adv_m;


