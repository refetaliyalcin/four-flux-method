function [n_eff_m,k_eff_m] = maxwell(f_v,n_pigment,k_pigment,n_medium,k_medium)

    %maxwell garnett method
    % Yalçin, R. A., & Ertürk, H. (2018). Inverse Design of Spectrally Selective Thickness Sensitive Pigmented Coatings for Solar Thermal Applications. Journal of Solar Energy Engineering, 140(3), 031006. doi:10.1115/1.4039276 
    % doi:10.1115/1.4039276 
    ev_pigment_e_a_=n_pigment.^2-k_pigment.^2; % eq. 2.30
    ev_pigment_e_a__=n_pigment.*k_pigment*2;
    ev_pigment_e_b_=n_medium.^2-k_medium.^2;
    ev_pigment_e_b__=n_medium.*k_medium*2;
    e_b=ev_pigment_e_b_-ev_pigment_e_b__*1i;
    e_a=ev_pigment_e_a_-ev_pigment_e_a__*1i;
    
    e_mg=e_b.*(e_a+2*e_b+2*f_v*(e_a-e_b))./(e_a+2*e_b-f_v*(e_a-e_b)); %Effective medium models for the optical properties of inhomogeneous materials
    e_mg_=real(e_mg);
    e_mg__=-imag(e_mg);
    n_eff_m=sqrt(0.5*(e_mg_+sqrt(e_mg_.^2+e_mg__.^2)));
    k_eff_m=sqrt(0.5*(-e_mg_+sqrt(e_mg_.^2+e_mg__.^2)));

end