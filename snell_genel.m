function reflectance=snell_genel(s_z,n_ortam,k_ortam,n_outside,k_outside)

    cos_teta=abs(s_z);
    
    if (cos_teta>0.9999) %dik gelmissen açiyi karistirma
        reflectance=((n_ortam-n_outside)^2+(k_ortam-k_outside)^2)/((n_ortam+n_outside)^2+(k_ortam+k_outside)^2);
    elseif (cos_teta<(0.0001)) %paralel gelmissen yansit
        reflectance=1;
    else
        sin_teta=sqrt(1-cos_teta*cos_teta);
        carpan=(n_ortam-1i*k_ortam)/(n_outside-1i*k_outside);
        sin_x=sin_teta*abs(carpan);
        if (sin_x>1)
            reflectance=1;
        else
            cos_x=sqrt(1-sin_x*sin_x);
            E_parallel=(cos_teta/cos_x-carpan)/(cos_teta/cos_x+carpan);
            R_parallel=E_parallel*conj(E_parallel);
            E_orth=-(cos_x/cos_teta-carpan)/(cos_x/cos_teta+carpan);
            R_orth=E_orth*conj(E_orth);
            reflectance=real(R_parallel+R_orth)*0.5;
        end
    end
    
 end 

    
       
    