clc
clear all
% close all

set(0, 'DefaultLineLineWidth', 2); %set thickness of all the lines = 2

exp_data

lamda=(400:10:2000)'*10^-9;
lamda_nm=lamda*10^9;

f_v=0.00079;
h = 40*10^-6;

eleman = length(lamda);
sonuc = zeros(eleman,1);

RT_semi=zeros(eleman,2);
RT_metal=zeros(eleman,2);
for j=1:eleman
    RT_semi(j,:)=Find_R_T_size_dist(h,lamda(j),f_v,1)*100;
    RT_metal(j,:)=Find_R_T_size_dist(h,lamda(j),f_v,2)*100;
end

figure('Renderer', 'painters', 'Position', [500 300 428 420])
hAx=gca;
plot(tra_d(:,1)*1000,tra_d(:,2)/100,'-k')
hold on;
plot(lamda_nm,RT_semi(:,2)/100,'ok','MarkerIndices',1:10:length(RT_semi(:,2)))
plot(tra_m(:,1)*1000,tra_m(:,2)/100,'--k')
plot(lamda_nm,RT_metal(:,2)/100,'sk','MarkerIndices',1:10:length(RT_metal(:,2)))
hAx.XColor = [0 0 0];
hAx.YColor = [0 0 0];
hAx.LineWidth = 1.5;
axis square
hLg=legend('Semiconducting, Experiment','Semiconducting, This study','Metallic, Experiment','Metallic, This study','Location','southeast');
hLg.LineWidth=1.5;
hLg.EdgeColor = [0 0 0];
xlabel('Wavelength [nm]')
ylh=ylabel('Transmittance');
ylh.VerticalAlignment	= 'bottom';
xlim([400 2000])
ylim([0 1] )
set(gca,'FontSize',13)
saveas(gcf,'validation.png')
saveas(gcf,'validation.pdf')

solar=I_solar(lamda);
T_lamda_metal_exp=interp1(tra_m(:,1)*1000,tra_m(:,2),lamda_nm);
T_lamda_dielec_exp=interp1(tra_d(:,1)*1000,tra_d(:,2),lamda_nm);
T_lamda_metal_calc=RT_metal(:,2);
T_lamda_dielec_calc=RT_semi(:,2);
w_lum=f_lum(lamda);

T_solar_metal_exp=trapz(lamda_nm,T_lamda_metal_exp.*solar)/trapz(lamda_nm,solar);
T_solar_dielec_exp=trapz(lamda_nm,T_lamda_dielec_exp.*solar)/trapz(lamda_nm,solar);
T_solar_metal_calc=trapz(lamda_nm,T_lamda_metal_calc.*solar)/trapz(lamda_nm,solar);
T_solar_dielec_calc=trapz(lamda_nm,T_lamda_dielec_calc.*solar)/trapz(lamda_nm,solar);

T_vis_metal_exp=trapz(lamda_nm,T_lamda_metal_exp.*solar.*w_lum)/trapz(lamda_nm,solar.*w_lum);
T_vis_dielec_exp=trapz(lamda_nm,T_lamda_dielec_exp.*solar.*w_lum)/trapz(lamda_nm,solar.*w_lum);
T_vis_metal_calc=trapz(lamda_nm,T_lamda_metal_calc.*solar.*w_lum)/trapz(lamda_nm,solar.*w_lum);
T_vis_dielec_calc=trapz(lamda_nm,T_lamda_dielec_calc.*solar.*w_lum)/trapz(lamda_nm,solar.*w_lum);

err_sol_metal=abs(T_solar_metal_exp-T_solar_metal_calc)
err_sol_dielec=abs(T_solar_dielec_exp-T_solar_dielec_calc)
err_vis_metal=abs(T_vis_metal_exp-T_vis_metal_calc)
err_vis_dielec=abs(T_vis_dielec_exp-T_vis_dielec_calc)

DT_solar_calc=abs(T_solar_metal_calc-T_solar_dielec_calc)
DT_solar_exp=abs(T_solar_metal_exp-T_solar_dielec_exp)
DT_vis_calc=abs(T_vis_metal_calc-T_vis_dielec_calc)
DT_vis_exp=abs(T_vis_metal_exp-T_vis_dielec_exp)

max_error=max([err_sol_metal,err_sol_dielec,err_vis_metal,err_vis_dielec])  