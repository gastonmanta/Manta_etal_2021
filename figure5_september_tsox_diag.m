clear all;close all;clc
load msm60_matrix_final


water_masses=({'SW','SACW','AAIW','UCDW','NADW','LCDW','AABW'});

% % rr=gamma-roo;
% g_n_levels_sabrina=[19 26.35 27.1 27.6 27.9 28.12 28.22 35];

%following valla et al 2018

g_n_levels=[19 26.35 27.1 27.6 27.9 28.1 28.1 28.27 35];


[SA, in_ocean] = gsw_SA_from_SP(sal,pres,lon,-34.5);

CT = gsw_CT_from_t(SA,temp,pres);



ind=find(isnan(real_values)==1);
gammaa=gamma;gammaa(ind)=NaN;oxx=ox;oxx(ind)=NaN;

CT(ind)=NaN;
SA(ind)=NaN;
sall=sal;sall(ind)=NaN;tempp=temp;tempp(ind)=NaN;thetaa=theta;thetaa(ind)=NaN;
stations=repmat(1:124,length(ox),1);

wm_text=[25.9 26.8 27.3 27.8 28 28.2 28.33];




figure
h = subplot(1,2,2); 
p = get(h, 'pos');
p(3) = p(3) + 0.07; p(1) = p(1) - 0.07;p(4) = p(4) + 0.04;p(2) = p(2) + 0.02;set(h, 'pos', p);box on; hold on


cmap = jet(124);
for k = 1:124
%  plot(oxx(:,k), gamma(:,k), 'Color', cmap(k, :));hold on
  plot(oxx(:,k), gamma(:,k),'.','Markersize',1, 'Color', cmap(k, :));hold on


end
ylim([25 28.39]);axis ij;
hold on
for i=1:length(g_n_levels)
for k=1:length(water_masses);
line ([147 270],[g_n_levels(i) g_n_levels(i)],'Linestyle','--','Color','k')
text(148.5,wm_text(k),water_masses(k),'Fontsize',12);hold on
end
end
xlim([147 265])
colormap(jet(124))

yticks([26.35,27.1,27.6,27.9,28.1,28.27,35])
xlabel('Oxygen Concentration (mol.kg^-^1)')
ylabel('/gamma^n (kg.m^-^3)')
xticks([160:20:260])

text(0.92,0.075,'b','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',18, 'Edgecolor','k')


cbh=colorbar('Location','eastoutside');;
cbh.Ticks = linspace(0, 1, 9) ; %Create 8 ticks from zero to 1
cbh.TickLabels = {'51.8ºW','46.5ºW','37.2ºW','27.3ºW','17.5ºW','7.7ºW','1.5ºE','11.2ºE','17.9ºE'} ;    %Replace the labels of these 8 ticks

h = subplot(1,2,1); 
p = get(h, 'pos');
p(3) = p(3) + 0.07; p(1) = p(1) - 0.07;p(4) = p(4) + 0.04;p(2) = p(2) + 0.02;set(h, 'pos', p);box on; hold on

cmap = jet(124);

theta_sdiag_background(CT,SA);hold on

for k = 1:124
  plot(SA(:,k), CT(:,k),'.','Markersize',1, 'Color', cmap(k, :));hold on
end
ylim([-1 25]);xlim([34.2 36.7]);%axis ij;

g_n_levels=[26.35 27.1 27.6 27.9 28.1 28.1 28.27];

hold on
contour(SA,CT,gamma,[g_n_levels],'k','linewidth',3)


text(0.006,0.075,'a','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',18, 'Edgecolor','k')

set(findall(gcf,'-property','Markersize'),'Markersize',5)
set(findall(gcf,'-property','FontSize'),'FontSize',13)


% savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/ts_ox_merian_june')
% print(gcf, '-dpng','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/ts_ox_merian_june')
% print(gcf, '-dpdf','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/ts_ox_merian_june')


% theta_new=theta.*land_mask;
% 
% figure
% cmap = jet(124);
% 
% theta_sdiag_background(theta_new,sal);hold on
% 
% for k = 1:124
%   plot(sal(:,k), theta_new(:,k),'.','Markersize',4, 'Color', cmap(k, :));hold on
% end
% 
% ylim([-1 25]);xlim([34.0 36.5]);%axis ij;
% 
% 
%  set(findall(gcf,'-property','FontSize'),'FontSize',13)
% 
% text(0.006,0.075,'c','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',20, 'Edgecolor','k')
