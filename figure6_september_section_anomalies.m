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

water_masses=({'TW','SACW','AAIW','UCDW','NADW','LCDW','AABW'});

g_n_levels=[19 26.35 27.1 27.6 27.9 28.1 28.1 28.27 35];
% just fot the plotting of the ctd locations
distance_between_stations=sw_dist(lon,lat,'km')*1000;
pp=cumsum(distance_between_stations);pp=pp(end);

[SA, in_ocean] = gsw_SA_from_SP(sal,pres,lon,-34.5);

CT = gsw_CT_from_t(SA,temp,pres);

gamma_GP = gamma_GP_from_SP_pt(SA,CT,pres,lon,-34.5);

gamma=gamma.*bucket;

gamma_GP=gamma_GP.*bucket;

z=pres;x=lon;pos_grid=lon;


%%%%%%%%%%%%%%%%%%%%%% PLOT SECTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLOT SECTIONS

fontsize=13;
gamma_font_color=[0 0 0];
gamma_color=[0.4 0.4 0.4];

gammas=gamma.*bucket;

gamma=gamma.*bucket;
sep=25;
 

%% anoms 
ct_new=CT.*bucket;
sa_new=SA.*bucket;
ox_new=ox.*bucket;
gamma_new=gamma.*bucket;

ct_anom=(ct_new-(nmean(ct_new,2)));
sa_anom=(sa_new-(nmean(sa_new,2)));
ox_anom=(ox_new-(nmean(ox_new,2)));ox_anom=movmean(ox_anom,33,1);ox_anom=movmean(ox_anom,5,2);




gamma_anom=(gamma_new-(nmean(gamma_new,2)));
dens_anom=(gamma_new-(nmean(gamma_new,2)));


rho = sw_pden(SA,CT,pres,1);
rho_new=rho.*bucket;

rho_anom=(rho_new-(nmean(rho_new,2)));


extent=1000;



%% temperature anomaly in 2
%%%THETA PLOT%%%%.

ind=find(ct_anom<-5);ct_anom(ind)=-5;
ind=find(ct_anom>5);ct_anom(ind)=5;

figure;hold on;

h = subplot(3,1,1); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.035;p(4) = p(4) + 0.15;p(2) = p(2) + 0.07;set(h, 'pos', p);box on; hold on

m_proj('mercator','long', [-52 18],'lat', [-36 -33]);hold on
%[CS,CH]=m_etopo2('contourf',[-3500 -200],'edgecolor','none');caxis([-5500 7000]);colormap(flip(gray(7)));  
m_plot(lon_adcp,lat_adcp,'.k','Markersize',1);hold on
ll=2;

for id_time= 1:length(CEs_out60)
    m_plot(CEs_out60{id_time,id_time,1},CEs_out60{id_time,id_time,2},'c','LineWidth',1)
end

for id_time= 1:length(AEs_out60)
    m_plot(AEs_out60{id_time,id_time,1},AEs_out60{id_time,id_time,2},'m','LineWidth',1)
end

m_grid('xtick',[],'ytick',[],'linestyle','none','tickstyle','dd');

h = subplot(3,1,2); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.035;p(4) = p(4) + 0.14;p(2) = p(2) + 0.15;set(h, 'pos', p);box on; hold on

[CS,CH]=contourf(x,z(1:extent),ct_anom(1:extent,:),[-5:1:-1 -1:0.25:1 1:1:5],'edgecolor','none');
hold on

% contour(x,z,ct_anom,[-1 1],'color',gamma_color);
%[CS,CH]=contour(x,z,ct_anom,[-.5 .5],'color',gamma_color);
% [CS,CH]=contour(x,z,ct_anom,[0 0],'w');

caxis([-5 5])
fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');
axis ij

load('MyColormap.mat')
colormap(mymap)

ylim([-sep 1000]);set(gcf,'color','w');  % otherwise 'print' turns lakes black


[CS,CH]=contour(x,z,ct_anom,[-1 1],'color',gamma_color);
clabel(CS,CH,'manual','FontSize',fontsize,'Color',gamma_color);

plot(lon_stations,-sep,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
plot(lon_stations(ce_inn),-sep,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
plot(lon_stations(ae_inn),-sep,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on

axis ij

[Cx,hContourx] = contour(x,z(1:extent),ct_anom(1:extent,:),[0 0],'w','linewidth',2);
[C,hContour] = contour(pos_grid(1:250),z(1:extent),gammas(1:extent,(1:250)),g_n_levels,'--k', 'showtext','on','labelspacing',150);
clabel(C,hContour,'FontSize',fontsize,'Color','k');
hold on
[C,hContour] = contour(pos_grid(250:end),z(1:extent),gammas(1:extent,(250:end)),g_n_levels,'--k');

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');

ylim([-sep extent]);
xticks([-50:5:15]);xticklabels([]);
set(gca,'tickdir','out');box on

h = subplot(3,1,3); p = get(h, 'pos');p(3) = p(3) + 0.12; p(1) = p(1) - 0.035;p(4) = p(4) + 0.25; p(2) = p(2) - 0.02;set(h, 'pos', p);box on; hold on

[CS,CH]=contourf(x,z,ct_anom,[-5:1:-1 -1:0.25:1 1:1:5],'edgecolor','none');
hold on;axis ij
[C,hContour] = contour(x,z,ct_anom,[0 0],'w','linewidth',2);

% [C,hContour] = contour(pos_grid,z,gamma,g_n_levels,'--k', 'showtext','on','labelspacing',550);
% clabel(C,hContour,'FontSize',fontsize,'Color','k');
% hold on
% [C,hContour] = contour(pos_grid(550:end),z(extent:end),gammas(extent:end,(550:end)),g_n_levels,'--k');

[C,hContour] = contour(pos_grid(1:550),z(1000:end),gammas(1000:end,(1:550)),g_n_levels,'--k', 'showtext','on','linewidth',0.5,'labelspacing',550);
clabel(C,hContour,'FontSize',fontsize,'Color','k');
hold on
[C,hContour] = contour(pos_grid(550:end),z(1000:end),gammas(1000:end,(550:end)),g_n_levels,'--k','linewidth',0.5);

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);

% fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);
ylim([extent+1 5500]);
caxis([-5 5])

[CS,CH]=m_contfbar([.37 .72],.145,CS,CH);
title(CS,' ? anomaly (ºC)','Fontweight','normal')


[CS,CH]=contour(x,z,ct_anom,[-1 1],'color',gamma_color);
clabel(CS,CH,'manual','FontSize',fontsize,'Color',gamma_color);


set(gcf,'color','w');  % otherwise 'print' turns lakes black
set(gca,'tickdir','out')
xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
ylim([extent+1 5500])

xlh = ylabel('Pressure (dbar)');xlh.Position(2) = xlh.Position(2) - 2000;clear xlh.Position(1);clear xlh.Position(2)

set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)


text(0.013,0.145,'b','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',20, 'Edgecolor','k')

% 
%  savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/ct_anom_section_june2')
%  print(gcf, '-dpng','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/ct_anom_section_june2')
%  print(gcf, '-dpdf','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/ct_anom_section_junee2')

%% sal anom new


ind=find(sa_anom<-.9);sa_anom(ind)=-.9;
ind=find(sa_anom>.9);sa_anom(ind)=.9;

figure;hold on;

h = subplot(3,1,1); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.035;p(4) = p(4) + 0.15;p(2) = p(2) + 0.07;set(h, 'pos', p);box on; hold on

m_proj('mercator','long', [-52 18],'lat', [-36 -33]);hold on
%[CS,CH]=m_etopo2('contourf',[-3500 -200],'edgecolor','none');caxis([-5500 7000]);colormap(flip(gray(7)));  
m_plot(lon_adcp,lat_adcp,'.k','Markersize',1);hold on
ll=2;

for id_time= 1:length(CEs_out60)
    m_plot(CEs_out60{id_time,id_time,1},CEs_out60{id_time,id_time,2},'c','LineWidth',1)
end

for id_time= 1:length(AEs_out60)
    m_plot(AEs_out60{id_time,id_time,1},AEs_out60{id_time,id_time,2},'m','LineWidth',1)
end

m_grid('xtick',[],'ytick',[],'linestyle','none','tickstyle','dd');

h = subplot(3,1,2); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.035;p(4) = p(4) + 0.14;p(2) = p(2) + 0.15;set(h, 'pos', p);box on; hold on

[CS,CH]=contourf(x,z(1:extent),sa_anom(1:extent,:),[-.9:.1:-.1 -.1:0.025:.1 .1:.1:.9],'edgecolor','none');
hold on

% contour(x,z,ct_anom,[-1 1],'color',gamma_color);
%[CS,CH]=contour(x,z,ct_anom,[-.5 .5],'color',gamma_color);
% [CS,CH]=contour(x,z,ct_anom,[0 0],'w');

caxis([-.9 .9])

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');
axis ij
colormap(mymap)
ylim([-sep 1000]);set(gcf,'color','w');  % otherwise 'print' turns lakes black


[CS,CH]=contour(x,z,sa_anom,[-.1 .1],'color',gamma_color);
clabel(CS,CH,'manual','FontSize',fontsize,'Color',gamma_color);

plot(lon_stations,-sep,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
plot(lon_stations(ce_inn),-sep,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
plot(lon_stations(ae_inn),-sep,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on

axis ij

[Cx,hContourx] = contour(x,z(1:extent),sa_anom(1:extent,:),[0 0],'w','linewidth',2);
[C,hContour] = contour(pos_grid(1:250),z(1:extent),gammas(1:extent,(1:250)),g_n_levels,'--k', 'showtext','on','labelspacing',150);
clabel(C,hContour,'FontSize',fontsize,'Color','k');
hold on
[C,hContour] = contour(pos_grid(250:end),z(1:extent),gammas(1:extent,(250:end)),g_n_levels,'--k');

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');

ylim([-sep extent]);
xticks([-50:5:15]);xticklabels([]);
set(gca,'tickdir','out');box on

h = subplot(3,1,3); p = get(h, 'pos');p(3) = p(3) + 0.12; p(1) = p(1) - 0.035;p(4) = p(4) + 0.25; p(2) = p(2) - 0.02;set(h, 'pos', p);box on; hold on

[CS,CH]=contourf(x,z,sa_anom,[-.9:.1:-.1 -.1:0.02:.1 .1:.1:.9],'edgecolor','none');
hold on


% [CS,CH]=contourf(x,z,ct_anom,[-5:1:-1 -1:0.25:1 1:1:5],'edgecolor','none');
% hold on;axis ij
[C,hContour] = contour(x,z,sa_anom,[0 0],'w','linewidth',2);
axis ij
ylim([extent+1 5500]);

[C,hContour]=contour(x,z,sa_anom,[-.1 .1],'color',gamma_color);
clabel(C,hContour,'manual','FontSize',fontsize,'Color',gamma_color);

% [C,hContour] = contour(pos_grid,z,gamma,g_n_levels,'--k', 'showtext','on','labelspacing',550);
% clabel(C,hContour,'FontSize',fontsize,'Color','k');
% hold on
% [C,hContour] = contour(pos_grid(550:end),z(extent:end),gammas(extent:end,(550:end)),g_n_levels,'--k');

[C,hContour] = contour(pos_grid(1:550),z(1000:end),gammas(1000:end,(1:550)),g_n_levels,'--k', 'showtext','on','linewidth',0.5,'labelspacing',550);
clabel(C,hContour,'FontSize',fontsize,'Color','k');
hold on
[C,hContour] = contour(pos_grid(550:end),z(1000:end),gammas(1000:end,(550:end)),g_n_levels,'--k','linewidth',0.5);

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);

% fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);

caxis([-.9 .9])



[CS,CH]=m_contfbar([.37 .72],.145,CS,CH);
title(CS,' S_A anomaly (g.kg^-^1)','Fontweight','normal')


set(gcf,'color','w');  % otherwise 'print' turns lakes black
set(gca,'tickdir','out')
xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
ylim([extent+1 5500])

xlh = ylabel('Pressure (dbar)');xlh.Position(2) = xlh.Position(2) - 2000;clear xlh.Position(1);clear xlh.Position(2)

set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)


text(0.013,0.145,'a','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',20, 'Edgecolor','k')

% 
%   savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/sa_anom_section_june')
%    print(gcf, '-dpng','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/sa_anom_section_june')
%  print(gcf, '-dpdf','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/sa_anom_section_june')
% 

%% rho anom


ind=find(rho_anom<-.6);rho_anom(ind)=-.6;
ind=find(rho_anom>.6);rho_anom(ind)=.6;

figure;hold on;

h = subplot(3,1,1); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.035;p(4) = p(4) + 0.15;p(2) = p(2) + 0.07;set(h, 'pos', p);box on; hold on

m_proj('mercator','long', [-52 18],'lat', [-36 -33]);hold on
%[CS,CH]=m_etopo2('contourf',[-3500 -200],'edgecolor','none');caxis([-5500 7000]);colormap(flip(gray(7)));  
m_plot(lon_adcp,lat_adcp,'.k','Markersize',1);hold on
ll=2;

for id_time= 1:length(CEs_out60)
    m_plot(CEs_out60{id_time,id_time,1},CEs_out60{id_time,id_time,2},'c','LineWidth',1)
end

for id_time= 1:length(AEs_out60)
    m_plot(AEs_out60{id_time,id_time,1},AEs_out60{id_time,id_time,2},'m','LineWidth',1)
end

m_grid('xtick',[],'ytick',[],'linestyle','none','tickstyle','dd');

h = subplot(3,1,2); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.035;p(4) = p(4) + 0.14;p(2) = p(2) + 0.15;set(h, 'pos', p);box on; hold on

[CS,CH]=contourf(x,z(1:extent),rho_anom(1:extent,:),[-.6:.1:-.1 -.1:0.025:.1 .1:.1:.6],'edgecolor','none');
hold on

% contour(x,z,ct_anom,[-1 1],'color',gamma_color);
%[CS,CH]=contour(x,z,ct_anom,[-.5 .5],'color',gamma_color);
% [CS,CH]=contour(x,z,ct_anom,[0 0],'w');

caxis([-.6 .6])

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');
axis ij
colormap(mymap)
ylim([-sep 1000]);set(gcf,'color','w');  % otherwise 'print' turns lakes black


[CS,CH]=contour(x,z,rho_anom,[-.1 .1],'color',gamma_color);
clabel(CS,CH,'manual','FontSize',fontsize,'Color',gamma_color);

plot(lon_stations,-sep,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
plot(lon_stations(ce_inn),-sep,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
plot(lon_stations(ae_inn),-sep,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on

axis ij

[Cx,hContourx] = contour(x,z(1:extent),rho_anom(1:extent,:),[0 0],'w','linewidth',2);

[C,hContour] = contour(pos_grid(1:250),z(1:extent),gammas(1:extent,(1:250)),g_n_levels,'--k', 'showtext','on','labelspacing',150);
clabel(C,hContour,'FontSize',fontsize,'Color','k');
hold on
[C,hContour] = contour(pos_grid(250:end),z(1:extent),gammas(1:extent,(250:end)),g_n_levels,'--k');

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');

ylim([-sep extent]);
xticks([-50:5:15]);xticklabels([]);
set(gca,'tickdir','out');box on

h = subplot(3,1,3); p = get(h, 'pos');p(3) = p(3) + 0.12; p(1) = p(1) - 0.035;p(4) = p(4) + 0.25; p(2) = p(2) - 0.02;set(h, 'pos', p);box on; hold on

[CS,CH]=contourf(x,z,rho_anom,[-.6:.1:-.1 -.1:0.02:.1 .1:.1:.6],'edgecolor','none');
hold on


% [CS,CH]=contourf(x,z,ct_anom,[-5:1:-1 -1:0.25:1 1:1:5],'edgecolor','none');
% hold on;
axis ij
[C,hContour] = contour(x,z,rho_anom,[0 0],'w','linewidth',2);
ylim([extent+1 5500]);

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);

[C,hContour]=contour(x,z,rho_anom,[-.1 .1],'color',gamma_color);
clabel(C,hContour,'manual','FontSize',fontsize,'Color',gamma_color);

% [C,hContour] = contour(pos_grid,z,gamma,g_n_levels,'--k', 'showtext','on','labelspacing',550);
% clabel(C,hContour,'FontSize',fontsize,'Color','k');
% hold on
% [C,hContour] = contour(pos_grid(550:end),z(extent:end),gammas(extent:end,(550:end)),g_n_levels,'--k');

[C,hContour] = contour(pos_grid(1:550),z(1000:end),gammas(1000:end,(1:550)),g_n_levels,'--k', 'showtext','on','linewidth',0.5,'labelspacing',550);
clabel(C,hContour,'FontSize',fontsize,'Color','k');
hold on
[C,hContour] = contour(pos_grid(550:end),z(1000:end),gammas(1000:end,(550:end)),g_n_levels,'--k','linewidth',0.5);


% fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);

caxis([-.6 .6])



[CS,CH]=m_contfbar([.37 .72],.145,CS,CH);
title(CS,' ?_? anomaly (kg.m^-^3)','Fontweight','normal')



set(gcf,'color','w');  % otherwise 'print' turns lakes black
set(gca,'tickdir','out')
xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
ylim([extent+1 5500])

xlh = ylabel('Pressure (dbar)');xlh.Position(2) = xlh.Position(2) - 2000;clear xlh.Position(1);clear xlh.Position(2)

set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)


text(0.013,0.145,'d','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',20, 'Edgecolor','k')

% % 
  savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/rho_anom_section_june')
   print(gcf, '-dpng','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/rho_anom_section_june')
 print(gcf, '-dpdf','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/rho_anom_section_june')



%%


%% ox anom


ind=find(ox_anom<-55);ox_anom(ind)=-55;
ind=find(ox_anom>55);ox_anom(ind)=55;



figure;hold on;

h = subplot(3,1,1); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.035;p(4) = p(4) + 0.15;p(2) = p(2) + 0.07;set(h, 'pos', p);box on; hold on

m_proj('mercator','long', [-52 18],'lat', [-36 -33]);hold on
%[CS,CH]=m_etopo2('contourf',[-3500 -200],'edgecolor','none');caxis([-5500 7000]);colormap(flip(gray(7)));  
m_plot(lon_adcp,lat_adcp,'.k','Markersize',1);hold on
ll=2;

for id_time= 1:length(CEs_out60)
    m_plot(CEs_out60{id_time,id_time,1},CEs_out60{id_time,id_time,2},'c','LineWidth',1)
end

for id_time= 1:length(AEs_out60)
    m_plot(AEs_out60{id_time,id_time,1},AEs_out60{id_time,id_time,2},'m','LineWidth',1)
end

m_grid('xtick',[],'ytick',[],'linestyle','none','tickstyle','dd');

h = subplot(3,1,2); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.035;p(4) = p(4) + 0.14;p(2) = p(2) + 0.15;set(h, 'pos', p);box on; hold on

[CS,CH]=contourf(x,z(1:extent),ox_anom(1:extent,:),[-55:5:55],'edgecolor','none');
hold on

% contour(x,z,ct_anom,[-1 1],'color',gamma_color);
%[CS,CH]=contour(x,z,ct_anom,[-.5 .5],'color',gamma_color);
% [CS,CH]=contour(x,z,ct_anom,[0 0],'w');

caxis([-55 55])

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');
axis ij
colormap(mymap)
ylim([-sep 1000]);set(gcf,'color','w');  % otherwise 'print' turns lakes black


[CS,CH]=contour(x,z,ox_anom,[-20 20],'color',gamma_color);
clabel(CS,CH,'manual','FontSize',fontsize,'Color',gamma_color);

plot(lon_stations,-sep,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
plot(lon_stations(ce_inn),-sep,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
plot(lon_stations(ae_inn),-sep,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on

axis ij

[Cx,hContourx] = contour(x,z(1:extent),ox_anom(1:extent,:),[0 0],'w','linewidth',2);

[C,hContour] = contour(pos_grid(1:250),z(1:extent),gammas(1:extent,(1:250)),g_n_levels,'--k', 'showtext','on','labelspacing',150);
clabel(C,hContour,'FontSize',fontsize,'Color','k');
hold on
[C,hContour] = contour(pos_grid(250:end),z(1:extent),gammas(1:extent,(250:end)),g_n_levels,'--k');

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');

ylim([-sep extent]);
xticks([-50:5:15]);xticklabels([]);
set(gca,'tickdir','out');box on

h = subplot(3,1,3); p = get(h, 'pos');p(3) = p(3) + 0.12; p(1) = p(1) - 0.035;p(4) = p(4) + 0.25; p(2) = p(2) - 0.02;set(h, 'pos', p);box on; hold on

[CS,CH]=contourf(x,z,ox_anom,[-55:5:55],'edgecolor','none');
hold on


% [CS,CH]=contourf(x,z,ct_anom,[-5:1:-1 -1:0.25:1 1:1:5],'edgecolor','none');
% hold on;
axis ij

[C,hContour] = contour(x,z,ox_anom,[0 0],'w','linewidth',2);
ylim([extent+1 5500]);

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);

[C,hContour]=contour(x,z,ox_anom,[-20 20],'color',gamma_color);
clabel(C,hContour,'manual','FontSize',fontsize,'Color',gamma_color);

% [C,hContour] = contour(pos_grid,z,gamma,g_n_levels,'--k', 'showtext','on','labelspacing',550);
% clabel(C,hContour,'FontSize',fontsize,'Color','k');
% hold on
% [C,hContour] = contour(pos_grid(550:end),z(extent:end),gammas(extent:end,(550:end)),g_n_levels,'--k');

[C,hContour] = contour(pos_grid(1:550),z(1000:end),gammas(1000:end,(1:550)),g_n_levels,'--k', 'showtext','on','linewidth',0.5,'labelspacing',550);
clabel(C,hContour,'FontSize',fontsize,'Color','k');
hold on
[C,hContour] = contour(pos_grid(550:end),z(1000:end),gammas(1000:end,(550:end)),g_n_levels,'--k','linewidth',0.5);


% fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);

caxis([-55 55])



[CS,CH]=m_contfbar([.37 .72],.145,CS,CH);

title(CS,' Oxygen anomaly (µmol.kg^-^1)','Fontweight','normal')


set(gcf,'color','w');  % otherwise 'print' turns lakes black
set(gca,'tickdir','out')
xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
ylim([extent+1 5500])

xlh = ylabel('Pressure (dbar)');xlh.Position(2) = xlh.Position(2) - 2000;clear xlh.Position(1);clear xlh.Position(2)

set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)


text(0.013,0.145,'c','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',20, 'Edgecolor','k')
% 
% 
  savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/ox_anom_section_june')
   print(gcf, '-dpng','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/ox_anom_section_june')
 print(gcf, '-dpdf','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/ox_anom_section_june')

