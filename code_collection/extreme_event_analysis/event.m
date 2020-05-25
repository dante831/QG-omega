classdef event
     properties
        event_timespan % in days from 1870?
        event_year 
        event_lonspan % longitude indices to calculate QG_omega
        event_latspan % latitude indices to calculate QG_omega
        computation_x
        computation_y
        event_lon
        event_lat
        event_lon_center
        event_lat_center
        event_level
        
        omega_500hPa_full
        omega_QG_500hPa_full % 15 by 15 omega_QG field
        Adv_500hPa_full

        omega_QG_max % maximum (minus) omega 
        omega_QG_max_p % pressure where maximum occurs
        omega
        omega_QG            % omega_QG during event_timespan, and in event_lonspan & event_latspan
        omega_QG_full       % omega_QG inverted with accurate lower and lateral boundaries, and zero upper boundary
        omega_QG_Adv        % omega_QG inverted with only Adv term on the rhs
        omega_QG_interior   % omega_QG induced by rhs only, with zero boundary condition
        omega_QG_full_b     % omega_QG induced by full omega as boundary condition, with zero rhs
        omega_QG_lateral_b  % omega QG induced by lateral boundaries only, with zero lower and rhs
        omega_QG_lateral_clim_b % same as above, except the lateral boundaries is climatological mean

        
        % diagnostics

        k2
        k2_star % the wavenumber of \nabla^2omega_{QG}*sigma_star
        m2
        l2
        A
        A1
        A2
        A3
        B
        C
        %C1
        %C2
        %Qx
        %Qy
        %dug_dlambda
        %dvg_dlambda
        %dug_dphi
        %dvg_dphi
        %dT_dlambda
        %dT_dphi
        sigma
        sigma_eff
        sigma_accu
        sigma_star % \propto dtheta/dp|theta^*
        dtheta_dp_ma
        dtheta_dp_ma_avg
        T
        T_avg
        %J_400hPa
        J_500hPa
        %J_600hPa
        J_center
        %J1
        %J2
        precip
        f0
        q
    
    end
end
