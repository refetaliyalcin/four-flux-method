function result = Find_R_T_size_dist(h,lamda,f_v,phase)


ksi=1;
r_avg_nm=42;
r_avg=r_avg_nm*10^-9;

r_vector_nm=(1:200); %size range to be considered
weight_vector = lognpdf(r_vector_nm,log(r_avg_nm),0.45); %the distribution
r_vector=r_vector_nm*10^-9;
weight_vector=weight_vector/trapz(r_vector,weight_vector); 

Z=h;

IcZ=1;
IdZ=10^-5;

n_matrix=PDMS_n(lamda);
k_matrix=PDMS_k(lamda);

n_air=1;
k_air=0;
if phase==1
    ev_pigment_n= vo2_semi_n_currie(lamda);
    ev_pigment_k= vo2_semi_k_currie(lamda);
end
if phase==2
    ev_pigment_n= vo2_metallic_n_currie(lamda);
    ev_pigment_k= vo2_metallic_k_currie(lamda);
end

ev_pigment_e_a_=ev_pigment_n^2-ev_pigment_k^2; % eq. 2.30
ev_pigment_e_a__=ev_pigment_n*ev_pigment_k*2;
ev_pigment_e_b_=n_matrix^2-k_matrix^2;
ev_pigment_e_b__=n_matrix*k_matrix*2;
e_b=ev_pigment_e_b_-ev_pigment_e_b__*1i;
e_a=ev_pigment_e_a_-ev_pigment_e_a__*1i;
e_mg=e_b*(e_a+2*e_b+2*f_v*(e_a-e_b))/(e_a+2*e_b-f_v*(e_a-e_b)); %Effective medium models for the optical properties of inhomogeneous materials
e_mg_=real(e_mg);
e_mg__=-imag(e_mg);
n_coating=sqrt(0.5*(e_mg_+sqrt(e_mg_^2+e_mg__^2)));
k_coating=sqrt(0.5*(-e_mg_+sqrt(e_mg_^2+e_mg__^2)));
% n_coating=n_matrix;
% k_coating=k_matrix;


rc=snell_genel(cosd(0),n_air,k_air,n_coating,k_coating);
rc_b=snell_genel(cosd(0),n_coating,k_coating,n_air,k_air);

acilar_deg=linspace(0,90,10000);
acilar_rad=acilar_deg*pi/180;
r_d_air_to_coating_ang=zeros(1,length(acilar_deg));
r_d_coating_to_air_ang=zeros(1,length(acilar_deg));
for i=1:length(acilar_deg)
    r_d_air_to_coating_ang(i)=snell_genel(cosd(acilar_deg(i)),n_air,k_air,n_coating,k_coating);
    r_d_coating_to_air_ang(i)=snell_genel(cosd(acilar_deg(i)),n_coating,k_coating,n_air,k_air);
end
r_d_air_to_coating=trapz(acilar_rad,2*cos(acilar_rad).*sin(acilar_rad).*r_d_air_to_coating_ang);
r_d_coating_to_air=trapz(acilar_rad,2*cos(acilar_rad).*sin(acilar_rad).*r_d_coating_to_air_ang);

rd_e=r_d_coating_to_air;
rd_i=r_d_air_to_coating;
rd_b=r_d_coating_to_air;



C_sca_r=zeros(1,length(r_vector));
C_abs_r=zeros(1,length(r_vector));
g_arr_r=zeros(1,length(r_vector));
sigma_c_r=zeros(1,length(r_vector));
for z=1:length(r_vector)
    Area=pi*r_vector(z)^2;
    x=2*pi*r_vector(z)*n_matrix/lamda;
    m=(ev_pigment_n+1i*ev_pigment_k)/n_matrix;
    fonksiyon=Mie(m,x);
    Qabs=fonksiyon(3);
    Qsca=fonksiyon(2);
    C_sca_r(z)=Qsca*Area;
    C_abs_r(z)=Qabs*Area;
    g_arr_r(z)=fonksiyon(5);
    
    
    fe=floor(x+4*x^(1/3)+2);
    a_ve_b=Mie_ab(m,x);
    sigma_c_r(z)=0;
    for k=1:fe
    for j=1:fe
        if mod(j,2)==1 && mod(k,2)==1
         sigma_c_r(z)=sigma_c_r(z)+qnm(j,k)*real(a_ve_b(1,k)*conj(a_ve_b(2,j)));
        end
    end
    end
    for k=2:fe
    for j=1:fe
        if mod(j,2)==1 && mod(k,2)==0
        sigma_c_r(z)=sigma_c_r(z)+pnm(j,k)*real(a_ve_b(1,k)*conj(a_ve_b(1,j))+a_ve_b(2,k)*conj(a_ve_b(2,j)));
        end
    end    
    end
    sigma_c_r(z)=0.5-sigma_c_r(z)*2/(x^2*Qsca);  
end

r_3=trapz(r_vector,(weight_vector.*r_vector.^3)');
V_avg=(4/3)*pi*r_3;
C_sca=trapz(r_vector,(weight_vector.*C_sca_r)');
C_abs=trapz(r_vector,(weight_vector.*C_abs_r)');
alfa=(C_sca)*f_v/V_avg;
beta=(C_abs)*f_v/V_avg+4*pi*k_matrix/lamda;
sigma_c=trapz(r_vector,(weight_vector.*sigma_c_r.*C_sca_r)')./trapz(r_vector,(weight_vector.*C_sca_r)');
% g=trapz(r_vector,(weight_vector.*g_arr_sigma.*C_sca_sigma)')./trapz(r_vector,(weight_vector.*C_sca_sigma)');



sigma_d=sigma_c;

A1=ksi^2*beta*(beta+2*(1-sigma_d)*alfa);%mlg'de bu sqrt of A1
A2=alfa*(ksi*beta*sigma_c+ksi*alfa*(1-sigma_d)+(alfa+beta)*sigma_c);
A3=alfa*(1-sigma_d)*(alfa+beta)*(ksi-1);
A4=ksi*(beta+alfa*(1-sigma_d));
A5=ksi*alfa*(1-sigma_d);
A6=(rc*(1-2*rc_b)+rc_b)/(1-rc*rc_b);
A7=(rd_i+rd_b*(1-rd_e-rd_i))/(1-rd_b*rd_e); 



C1_ust=(1-rc)*IcZ*exp(-(alfa+beta)*Z);
C1_alt=1-rc*A6*exp(-2*(alfa+beta)*Z);
C1=C1_ust/C1_alt;
C2_ust=(1-rc)*A6*IcZ*exp(-(alfa+beta)*Z);
C2_alt=1-rc*A6*exp(-2*(alfa+beta)*Z);
C2=C2_ust/C2_alt;

C5=A2*C1/(A1-(alfa+beta)^2);
C6=A3*C2/(A1-(alfa+beta)^2);

C9=A3*C1/(A1-(alfa+beta)^2);
C10=A2*C2/(A1-(alfa+beta)^2);
A8=C9+C10-A7*(C5+C6);
A9=A7-(A4-sqrt(A1))/A5;
A10=A7-(A4+sqrt(A1))/A5;
A11=(1-rd_i*(A4-sqrt(A1))/A5)*exp(sqrt(A1)*Z);
A12=(1-rd_i*(A4+sqrt(A1))/A5)*exp(-sqrt(A1)*Z);
A13=(1-rd_e)*IdZ-(C5-rd_i*C9)*exp((alfa+beta)*Z)-(C6-rd_i*C10)*exp(-Z*(alfa+beta));
C3=(A10*A13-A8*A12)/(A10*A11-A9*A12);
C4=(A8*A11-A9*A13)/(A10*A11-A9*A12);
C7=C3*(A4-sqrt(A1))/A5;
C8=C4*(A4+sqrt(A1))/A5;

term_1= (IcZ)/(IcZ+IdZ);
term_2= (IdZ)/(IcZ+IdZ);

t_c= 1;
t_d= 1;

%Tcc= (t_c*(1-rc_b)*((1-rc)^2)*exp(-(alpha+beta)*Z))/((1-rc*rc_b)-(rc*(rc+rc_b-2*rc*rc_b)*exp(-2*(alpha+beta)*Z)));
Tcc= ((t_c*(1-rc_b)*((1-rc)^2)*exp(-(alfa+beta)*Z))/((1-rc*rc_b)-(rc*(rc+rc_b-2*rc*rc_b)*exp(-2*(alfa+beta)*Z))))*term_1;

Tdt= (((1-rd_b)*(1-rd_i)*t_d)/((1-rd_e*rd_b)*(IcZ+IdZ)))*(((A10*A13+A8*A11-A8*A12-A9*A13)/(A10*A11-A9*A12))+((A2*C1+A3*C2)/(A1-((alfa+beta)^2))));

%Rc=rc+((1-rc)^2*(rc+rc_b-2*rc*rc_b)*exp(-2*(alfa+beta)*Z))/(1-rc*rc_b-rc*(rc+rc_b-2*rc*rc_b)*exp(-2*(alfa+beta)*Z));
Rc=(rc+((1-rc)^2*(rc+rc_b-2*rc*rc_b)*exp(-2*(alfa+beta)*Z))/(1-rc*rc_b-rc*(rc+rc_b-2*rc*rc_b)*exp(-2*(alfa+beta)*Z)))*term_1;

Rd=rd_e*IdZ/(IcZ+IdZ)+(1-rd_i)*(C7*exp(sqrt(A1)*Z)+C8*exp(-sqrt(A1)*Z)+C9*exp(Z*(alfa+beta))+C10*exp(-Z*(alfa+beta)));
%Rd=rd_e*IdZ/(IcZ+IdZ)+((1-rd_i)/(IcZ+IdZ))*((C7*exp(sqrt(A1)*Z)+C8*exp(-sqrt(A1)*Z)+C9*exp(Z*(alfa+beta))+C10*exp(-Z*(alfa+beta))));

R=Rc+Rd;
T=Tcc+Tdt;

result=[R T];
 