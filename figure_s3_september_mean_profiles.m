%LOAD DATA
clear all;close all;clc;
cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian
load MSM60_ctd_october.mat

load cloro_merian

load vu;v=v/100;u=u/100;

load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian/stations.mat')
 
ae_in=nansum(ae_in);ce_in=nansum(ce_in);
ae_inn=find(ae_in==1);ce_inn=find(ce_in==1);

cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian
% tsg=importdata('msm_060_1_no_headers.tsg');
load stations
load merian_contours
cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian

g_n_levels=[19 26.35 27.1 27.6 27.9 28.1 28.1 28.27 35];
distance_between_stations=sw_dist(lon,lat,'km')*1000;
pp=cumsum(distance_between_stations);pp=pp(end);

[SA, in_ocean] = gsw_SA_from_SP(sal,pres,lon,-34.5);

CT = gsw_CT_from_t(SA,temp,pres);

z=pres;x=lon;pos_grid=lon;


CT=CT.*bucket;
SA=SA.*bucket;
ox=ox.*bucket;



tw1=nmean(CT(:,1:700),2);
tw2=nstd(CT(:,1:700),1,2);

te1=nmean(CT(:,700:end),2);
te2=nstd(CT(:,700:end),1,2);

sw1=nmean(SA(:,1:700),2);
sw2=nstd(SA(:,1:700),1,2);

se1=nmean(SA(:,700:end),2);
se2=nstd(SA(:,700:end),1,2);

ow1=nmean(ox(:,1:700),2);
ow2=nstd(ox(:,1:700),1,2);

oe1=nmean(ox(:,700:end),2);
oe2=nstd(ox(:,700:end),1,2);

y1=z;


% figure; hold on
ha = tight_subplot(1,3,[.01 .01],[.1 .01],[.1 .01]);

% ha = tight_subplot(1,3,[.1 .01],[.1 .01],[.1 .01]);

set(findall(gcf,'-property','FontSize'),'FontSize',11)


axes(ha(1));box on; grid on

plot(-100:-10,-100:-10,'r');hold on;plot(-100:-10,-100:-10,'b');

legend({'West','East'},'Autoupdate','off','location','south','box','off')
boundedline([tw1], y1, tw2,'r', 'orientation', 'horiz', 'alpha', 'nan', 'fill')
boundedline([te1], y1, te2,'b', 'orientation', 'horiz', 'alpha', 'nan', 'fill')
axis ij;ylim([0 5500]);yticks([0:250:5500])
yticklabels([0:250:5500]);ylabel('Pressure (dbar)')

yticklabels({'0','','500','','1000','','1500','','2000','','2500','','3000','','3500','','4000','','4500','','5000',''})

xticks([0:2:24]);%xticklabels([0:2:24])

xlim([-0.5 24])

xticklabels({'0','','4','','8','','12','','16','','20'})
grid on
xlabel('\Theta (ºC)')


ylim([0 5300])

text(0.9,0.065,'a','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')


axes(ha(2));box on; grid on
boundedline([sw1], y1, sw2,'r', 'orientation', 'horiz', 'alpha', 'nan', 'fill')
boundedline([se1], y1, se2,'b', 'orientation', 'horiz', 'alpha', 'nan', 'fill')
axis ij;ylim([0 5500]);yticks([0:200:5500])

xlim([34.3 36.2])


ylim([0 5300])

xticks([34.4:0.2:36.6])

%xticklabels([34:0.2:36.6])

xticklabels({'34.4','','34.8','','35.2','','35.6','','36','','36.4'})


xlabel('S_A (kg.m^-^3)')

text(0.9,0.065,'b','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')


axes(ha(3));box on; grid on
boundedline([ow1], y1, ow2,'r', 'orientation', 'horiz', 'alpha', 'nan', 'fill')
boundedline([oe1], y1, oe2,'b', 'orientation', 'horiz', 'alpha', 'nan', 'fill')
axis ij;ylim([0 5500]);yticks([0:200:5500])

xticks([170:10:270])

xticklabels([170:10:270])


xticklabels({'170','','190','','210','','230','','250','','270'})

ylim([0 5300])

xlabel('O_2 (µmol.kg^-^1)')

text(0.9,0.065,'c','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')


set(findall(gcf,'-property','FontSize'),'FontSize',12)
 set(gcf,'color','w');
set(findall(gcf,'-property','Linewidth'),'Linewidth',.7)

set(gcf,'renderer','Painters')

print(gcf, '-depsc','-painters' ,'-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/paper_merian/Figures_paper_merian/fig_s3_mean_profiles')


savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/paper_merian/Figures_paper_merian/fig_s3_mean_profiles')

print(gcf, '-dpng','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/paper_merian/Figures_paper_merian/fig_s3_mean_profiles')
print(gcf, '-djpeg','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/paper_merian/Figures_paper_merian/fig_s3_mean_profiles')


