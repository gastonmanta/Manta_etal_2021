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


% cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian
% tsg=importdata('msm_060_1_no_headers.tsg');
% 
% sst=ncread('sst_merian_MSM60.nc','analysed_sst');sst=sst-273.18;sst=flip(sst,1);%1 january to 1 february included
% lon_sst=ncread('sst_merian_MSM60.nc','lon');lon_sst=flip(lon_sst);
% lat_sst=ncread('sst_merian_MSM60.nc','lat');%topo=sw_pres(topo',-35);


% load /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/4gaston/vmadcp/msm60os75000_000000_39_hc.mat
% 
% v=squeeze(b.vel(:,2,:));
% u=squeeze(b.vel(:,1,:));
% zc=c.depth(:);
% pg=c.pg;%8 meter interval from 17m to 800
% kk=find(pg<=25);
% v(kk)=NaN;
% u(kk)=NaN;

% %mean between 40 and 100m from 38 khz adcp
% uplot=squeeze(nanmean(u(4:11,:),1));
% vplot=squeeze(nanmean(v(4:11,:),1));
% lon_adcp=b.nav.txy1(2,:);lat_adcp=b.nav.txy1(3,:);


% Neutral Density (Î³, kg/m3) Range Assign for Each Water Mass Layera
% Water mass Î³ range Î¸ (Â°C) S O2 (Î¼mol/kg)
% TW <26.35 19.75
% 36.09 225 35.10 211 34.245 237 34.58 188 34.91 241 34.72 215 34.67 226
% SACW 26.35â€“27.10
% AAIW 27.10â€“27.60
% UCDW 27.60â€“27.90
% NADW 27.90â€“28.10
% LCDW 28.10â€“28.27
% AABW >28.27 0.10
% 12.10 4.05 2.90 2.95 0.72
% Note. The zonally averaged potential temperature (Î¸), salinity (S), and dis- solved oxygen (O2) values at the center of each range computed from all cruises occupied at the SAMBA-W line are given as a reference.
% aTW: Tropical Water; SACW: South Atlantic Central Water; AAIW: Antarctic Intermediate Water; UCDW: Upper Circumpolar Deep Water; NADW: North Atlantic Deep Water; LCDW: Lower Circumpolar Deep Water; AABW: Antarctic Bottom Water.

water_masses=({'TW','SACW','AAIW','UCDW','NADW','LCDW','AABW'});

% rr=gamma-roo;
g_n_levels_sabrina=[19 26.35 27.1 27.6 27.9 28.12 28.22 35];

%following valla et al 2018

g_n_levels=[19 26.35 27.1 27.6 27.9 28.1 28.1 28.27 35];
% Waterdepth(121)=1040;

% lonn=lon;
% lonn(lonn<0)=lonn(lonn<0)+360
% lon_ssh(lon_ssh>180)=lon_ssh(lon_ssh>180)-360;
% [gamma_n, gamma_error_lower, gamma_error_upper] = eos80_legacy_gamma_n(sal_new(1:end-10,:),temp_new(1:end-10,:),press(1:end-10,:),lonn',lat')

% just fot the plotting of the ctd locations
distance_between_stations=sw_dist(lon,lat,'km')*1000;
pp=cumsum(distance_between_stations);pp=pp(end);
% xs=distance_between_stations*123/pp;xs=cumsum(xs); xs=[0; xs];
% xs=xs+1;
% ceros=(1:124)*0;
% Waterdepthh=Waterdepth;Waterdepthh(Waterdepthh>1000)=1000;
% Waterdepthhh=Waterdepth;Waterdepthhh(Waterdepthhh<1000)=1000;

[SA, in_ocean] = gsw_SA_from_SP(sal,pres,lon,-34.5);

CT = gsw_CT_from_t(SA,temp,pres);

gamma_GP = gamma_GP_from_SP_pt(SA,CT,pres,lon,-34.5);

gamma=gamma.*bucket;

gamma_GP=gamma_GP.*bucket;

z=pres;x=lon;pos_grid=lon;


%%%%%%%%%%%%%%%%%%%%%% PLOT SECTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT SECTIONS

fontsize=14;

%% temperature
%%%THETA PLOT%%%%


figure;hold on;

h = subplot(3,1,1); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.15;p(2) = p(2) + 0.07;set(h, 'pos', p);box on; hold on

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
m_gshhs_h('patch',[0.81818 0.77647 0.70909]);



h = subplot(3,1,2); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.14;p(2) = p(2) + 0.15;set(h, 'pos', p);box on; hold on
contourf(x,z(1:1000),theta(1:1000,:),[-0.5:0.5:23.5],'linestyle','none');hold on;%colormap(jet);colorbar;%title('temp');hold on

plot(lon_stations,0,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
plot(lon_stations(ce_inn),0,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
plot(lon_stations(ae_inn),0,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on


for l=1:length(stations)
line([stations(l) stations(l)],[0 max_pres(l)],'Color',[.6 .6 .6],'LineStyle','--')
xline(stations(l),'--','color',[.65 .65 .65]);hold on
end

[C,hContour] = contour(x,z(1:1000),CT(1:1000,:),[0 1 2 3 4 5 7 9 11 13 15 17 19 21 23],'k', 'showtext','on','Labelspacing',2500);
clabel(C,hContour,'FontSize',fontsize,'Color','k');

[Cx,hContourx] = contour(x,z(1:1000),CT(1:1000,:),[6 8 10 12 14 16 18  20 22],'k');

[C,hContour] = contour(x,z(1:1000),gamma(1:1000,:),g_n_levels,'--w', 'showtext','on','linewidth',1,'Labelspacing',2200);
set(hContour,'showtext','on');clabel(C,hContour,'FontSize',fontsize,'FontWeight','bold','Color','w');

% [C,hContour] = contour(x,z(1:1000),gamma(1:1000,:),g_n_levels_sabrina,'--w', 'showtext','on','linewidth',1,'Labelspacing',2300);
% set(hContour,'showtext','on');clabel(C,hContour,'FontSize',12,'Color','w');
% [C,hContour] = contour(x,z(1:1000),gamma_GP(1:1000,:),g_n_levels,'w', 'showtext','on','linewidth',1,'Labelspacing',2300);
% set(hContour,'showtext','on');clabel(C,hContour,'FontSize',fontsize,'FontWeight','bold','Color','w');
% [C,hContour] = contour(x,z(1:1000),gamma_GP(1:1000,:),g_n_levels_sabrina,'w', 'showtext','on','linewidth',1,'Labelspacing',2300);
% set(hContour,'showtext','on');clabel(C,hContour,'FontSize',12,'Color','w');


fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[1000 topo(1) topo 1000 ],[.65 .65 .65],'linestyle','none');
colormap(flipud(brewermap([],'Spectral')));caxis([-0.5 23.5]);ylim([0 1000]);xticks;axis ij;
xticks([-50:5:15]);xticklabels([]);
set(gca,'tickdir','out');box on
set(findall(gcf,'-property','FontSize'),'FontSize',12)

h = subplot(3,1,3); p = get(h, 'pos');p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.25; p(2) = p(2) - 0.02;set(h, 'pos', p);box on; hold on
contourf(x,z(1000:end),CT(1000:end,:),[-0.5:0.5:23.5],'linestyle','none');hold on;colormap(jet);%colorbar;title('temp');hold on

for l=1:length(stations)
line([stations(l) stations(l)],[0 max_pres(l)],'Color',[.6 .6 .6],'LineStyle','--')
%xline(stations(l),'--','color',[.65 .65 .65]);hold on
end

[C,hContour] = contour(x,z(1000:end),CT(1000:end,:),[0 1 2 3 4 5 7 9 11 13 15 17 19 21 23],'k', 'showtext','on','Labelspacing',2500);
clabel(C,hContour,'FontSize',fontsize,'Color','k');
[Cx,hContourx] = contour(x,z(1000:end,:),CT(1000:end,:),[6 8 10 12 14 16 18  20 22],'k');

[C,hContour] = contour(x,z(1000:end),gamma(1000:end,:),g_n_levels,'--w', 'showtext','on','linewidth',1,'Labelspacing',2700);

set(hContour,'showtext','on');clabel(C,hContour,'FontSize',fontsize,'FontWeight','bold','Color','w');axis ij
% 
% [C,hContour] = contour(x,z,gamma,g_n_levels_sabrina,'--w', 'showtext','on','linewidth',1,'Labelspacing',2700);
% set(hContour,'showtext','on');clabel(C,hContour,'FontSize',12,'Color','w');axis ij
% 
% [C,hContour] = contour(x,z,gamma_GP,g_n_levels,'w', 'showtext','on','linewidth',1,'Labelspacing',2700);
% set(hContour,'showtext','on');clabel(C,hContour,'FontSize',12,'Color','w');axis ij
% 
% [C,hContour] = contour(x,z,gamma_GP,g_n_levels_sabrina,'--w', 'showtext','on','linewidth',1,'Labelspacing',2700);
% set(hContour,'showtext','on');clabel(C,hContour,'FontSize',12,'Color','w');axis ij

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);

colormap(flipud(brewermap(49,'Spectral')));caxis([-0.5 23.5])
[ax,h]=m_contfbar([.35 .75],.15,[-0.5:0.5: 23.5],[-0.5:0.5: 23.5]);

title(ax,'?(ºC)','Fontweight','normal')
set(gcf,'color','w');  % otherwise 'print' turns lakes black
set(gca,'tickdir','out')
xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});

xlh = ylabel('Pressure (db)');xlh.Position(1) = xlh.Position(1) + 0.05;xlh.Position(2) = xlh.Position(2) - 2000;clear xlh.Position(1);clear xlh.Position(2)

set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)



savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/ct_section_octobers')
print(gcf, '-djpeg','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/ct_section')




% 
% spd=sqrt(v.*v+u.*u);spd=spd.*bucket;
% 
% 
% hold on
% 
% figure;colormap(jet)
% 
% plot(lon_stations,0,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
% plot(lon_stations(ce_inn),0,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
% plot(lon_stations(ae_inn),0,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on
% 
% contourf(lon,z,spd,[0:.05:.9],'linestyle','none');shading interp;hold on;caxis([0 .9]);axis ij
% contour(lon,z,spd,[.1 .1],'w');
% 
% 
% figure;colormap(jet)
% 
% plot(lon_stations,0,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
% plot(lon_stations(ce_inn),0,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
% plot(lon_stations(ae_inn),0,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on
% 
% contourf(lon,z,vort,[-0.6:.05:.6],'linestyle','none');shading interp;hold on;caxis([0 .9]);axis ij
% contour(lon,z,spd,[.1 .1],'w');
% 


%% salinity
%%%SALINITY PLOT%%%%
sall=SA;

ind=find(sall<34.001);sall(ind)=34;


fontsize=12;
%%%THETA PLOT%%%%
figure;hold on;

h = subplot(3,1,1); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.15;p(2) = p(2) + 0.07;set(h, 'pos', p);box on; hold on

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
m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
% cmocean('gray')


h = subplot(3,1,2); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.14;p(2) = p(2) + 0.15;set(h, 'pos', p);box on; hold on

contourf(x,z(1:1000),sall(1:1000,:),[34:0.05:36.4],'linestyle','none');hold on;%colormap(jet);colorbar;%title('temp');hold on

cmocean('delta',length(34:.1:36.4));caxis([34 36.4]);axis ij
 
% cmocean('delta',length(34:.025:36.5));caxis([34 36.5])


% for l=1:length(stations)
% xline(stations(l),'--','color',[.65 .65 .65]);hold on
% end

% [C,hContour] = contourf(x,z(5:1000),sall(5:1000,:),[34.2:0.1:37],'linestyle','none');axis ij

[C,hContour] = contour(x,z(50:1000),sall(50:1000,:),[34:0.2:36.4],'k');axis ij

clabel(C,hContour,'manual','FontSize',fontsize,'Color','k');



[C,hContour] = contour(x,z(1:1000),gamma(1:1000,:),g_n_levels,'--w', 'showtext','on','linewidth',1,'labelspacing',2700);

set(hContour,'showtext','on');
clabel(C,hContour,'FontSize',fontsize,'Fontweight','bold','Color','w');
 
% [C,hContour] = contour(x,z,gamma,g_n_levels_sabrina,'--w', 'showtext','on','linewidth',1,'Labelspacing',2700);
% set(hContour,'showtext','on');clabel(C,hContour,'FontSize',12,'Color','w');


% cmocean('delta',length(33.5:.1:36.5));caxis([33.5 36.5])
fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);

ylim([0 1000]);xticks;axis ij; hold on;
xticks([-50:5:15]);xticklabels([]);
set(gca,'tickdir','out');box on
set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)
plot(lon_stations,0,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
plot(lon_stations(ce_inn),0,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
plot(lon_stations(ae_inn),0,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on
% cmocean('tarn'),length(33.5:.1:36.5));caxis([33.5 36.5])

h = subplot(3,1,3); p = get(h, 'pos');p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.25; p(2) = p(2) - 0.02;set(h, 'pos', p);box on; hold on

contourf(x,z(1000:end),SA(1000:end,:),[33.5:0.025:37],'linestyle','none');hold on;%colormap(jet);%colorbar;title('temp');hold on

% for l=1:length(stations)
% xline(stations(l),'--','color',[.65 .65 .65]);hold on
% end

% [C,hContour] = contour(x,z(1000:end),sal(1000:end,:),[34.2:0.1:36.5],'k', 'showtext','on','Labelspacing',2500);
% clabel(C,hContour,'FontSize',12,'Color','k');

[C,hContour] = contour(x,z(1000:end),SA(1000:end,:),[34:0.1:36.4],'k');axis ij
clabel(C,hContour,'manual','FontSize',fontsize,'Color','k');


[C,hContour] = contour(x,z(1000:end),SA(1000:end,:),[34.85:0.1:37],'--k');axis ij

% [Cx,hContourx] = contour(x,z(1000:end,:),theta(1000:end,:),[6 8 10 12 14 16 18  20 22],'k');
% [C,hContour] = contour(x,z,gamma,g_n_levels,'--g', 'showtext','on','linewidth',1,'Labelspacing',2700);
% set(hContour,'showtext','on');clabel(C,hContour,'FontSize',12,'Color','g');

[C,hContour] = contour(x,z,gamma,g_n_levels,'--w', 'showtext','on','linewidth',1,'Labelspacing',2700);
set(hContour,'showtext','on');clabel(C,hContour,'FontSize',fontsize,'Fontweight','bold','Color','w');


fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);
axis ij

% cmocean((flipud('tarn')),length(33.5:.1:36.5));caxis([33.5 36.5])
cmocean('delta',length(34:.025:36.4));caxis([34 36.4])

[ax,h]=m_contfbar([.35 .675],.15,[34:.05:36.5],[34:.05:36.5]);
title(ax,'Absolute Salinity (g.kg^-^1)','Fontweight','normal')
set(gcf,'color','w');  % otherwise 'print' turns lakes black
set(gca,'tickdir','out')
xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
xlhh = ylabel('Pressure (db)');xlhh.Position(1) = xlhh.Position(1) + 0.35;xlhh.Position(2) = xlhh.Position(2) - 2000;%clear xlh

set(findall(gcf,'-property','FontSize'),'FontSize',11.5)
set(findall(gcf,'-property','FontWeight'),'FontWeight','normal')


savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/salinity_abs_section_dec')
print(gcf, '-djpeg','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/salinity_abs_section_dec')



print(gcf, '-dpng','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/salinity_abs_section_dec')


%% oxygen
%%%OXYGEN PLOT%%%%

oxplot=movmean(ox,33,1);oxplot=movmean(oxplot,5,2);
ox=oxplot;ox=ox.*bucket;

fontsize=12;
%%%THETA PLOT%%%%
figure;hold on;

h = subplot(3,1,1); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.15;p(2) = p(2) + 0.07;set(h, 'pos', p);box on; hold on

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
m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
% cmocean('gray')


h = subplot(3,1,2); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.14;p(2) = p(2) + 0.15;set(h, 'pos', p);box on; hold on
contourf(x,z(1:1000),ox(1:1000,:),[170:10:260],'linestyle','none');hold on;%colormap(jet);colorbar;%title('temp');hold on
% % 
% for l=1:length(stations)
% xline(stations(l),'--','color',[.65 .65 .65]);hold on
% end

% [C,hContour] = contour(x,z(1:1000),ox(1:1000,:),[165:20:265],'k');%clabel(C,hContour,'FontSize',12,'Color','k');
% [C,hContour] = contour(x,z(50:1000),ox(50:1000,:),[160:20:260],'k', 'showtext','on','labelspacing',2500);clabel(C,hContour,'FontSize',12,'Color','k');
%[C,hContour] = contour(x,z(1:1000),ox(1:1000,:),[175:20:265],'k','showtext','on','labelspacing',2700);clabel(C,hContour,'FontSize',12,'Color','k');
% [C,hContour] = contour(x,z(1:1000),ox(1:1000,:),[175:20:265],'k')%,'showtext','on','labelspacing',2700);clabel(C,hContour,'FontSize',12,'Color','k');
%  [C,hContour] = contour(x,z(1:1000),ox(1:1000,:),[170:20:270],'k', 'showtext','on','labelspacing',2500);axis ij
% clabel(C,hContour,'manual','FontSize',12,'Color','k');
 [C,hContour] = contour(x,z(1:1000),ox(1:1000,:),[170:20:270],'k');axis ij

clabel(C,hContour,'manual','FontSize',fontsize,'Color','k');

colormap(brewermap([],'RdPu'));caxis([170 260])

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);

% clabel( C , h ,'manual') 

[C,hContour] = contour(x,z,gamma,g_n_levels,'--w', 'showtext','on','linewidth',1,'Labelspacing',2700);
set(hContour,'showtext','on');clabel(C,hContour,'FontSize',fontsize,'Fontweight','bold','Color','w');
% set(hContour,'showtext','on');
% clabel(C,hContour,'FontSize',12,'Color','w');

% [C,hContour] = contour(x,z,gamma,g_n_levels_sabrina,'--w', 'showtext','on','linewidth',1,'Labelspacing',2700);
% set(hContour,'showtext','on');clabel(C,hContour,'FontSize',12,'Color','w');


% cmocean('delta')
ylim([0 1000]);xticks;axis ij;
xticks([-50:5:15]);xticklabels([]);
set(gca,'tickdir','out');box on
set(findall(gcf,'-property','FontSize'),'FontSize',12)
colormap(brewermap([],'RdPu'));caxis([165 265])
cmap = colormap((brewermap([],'RdPu')));
plot(lon_stations,0,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
plot(lon_stations(ce_inn),0,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
plot(lon_stations(ae_inn),0,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on

 
h = subplot(3,1,3); p = get(h, 'pos');p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.25; p(2) = p(2) - 0.02;set(h, 'pos', p);box on; hold on



contourf(x,z(1000:end),ox(1000:end,:),[170:10:260],'linestyle','none');hold on;%colormap(jet);%colorbar;title('temp');hold on

% for l=1:length(stations)
% xline(stations(l),'--','color',[.65 .65 .65]);hold on
% end


%[C,hContour] = contour(x,z(1000:end),ox(1000:end,:),[165:20:265],'k');%clabel(C,hContour,'FontSize',12,'Color','k');
%  [C,hContour] = contour(x,z(1000:end),ox(1000:end,:),[170:20:270],'k', 'showtext','on','labelspacing',2500);clabel(C,hContour,'FontSize',12,'Color','k');


 [C,hContour] = contour(x,z(1000:end),ox(1000:end,:),[170:20:270],'k');axis ij
clabel(C,hContour,'manual','FontSize',fontsize,'Color','k');


% [C,hContour] = contour(x,z,gamma,g_n_levels,'--g', 'showtext','on','linewidth',1,'Labelspacing',2700);
% set(hContour,'showtext','on');clabel(C,hContour,'FontSize',12,'Color','g');

[C,hContour] = contour(x,z,gamma,g_n_levels,'--w', 'showtext','on','linewidth',1,'Labelspacing',2700);
set(hContour,'showtext','on');clabel(C,hContour,'FontSize',fontsize,'Fontweight','bold','Color','w');


fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);
axis ij
% cmocean((flipud('tarn')),length(33.5:.1:36.5))
colormap(brewermap([],'RdPu'));caxis([170 260])

[ax,h]=m_contfbar([.38 .68],.15,[170:10:270],[170:10:270]);
title(ax,'Dissolved oxygen (?mol.kg^-^1)','Fontweight','normal')
set(gcf,'color','w');  % otherwise 'print' turns lakes black
set(gca,'tickdir','out')
xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
xlhh = ylabel('Pressure (db)');xlhh.Position(1) = xlhh.Position(1) - 0.2;xlhh.Position(2) = xlhh.Position(2) - 2000;%clear xlh
set(findall(gcf,'-property','FontSize'),'FontSize',12)

% 
%  savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/oxygen_section_octobers_smoothed')
%  print(gcf, '-djpeg','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/oxygen_section')

%%%%%%%%%%%%%%%%%%%%     V %%%%%%%%%%%%%%%%%%%%%%%%


%% meridional vel plot
figure
h = subplot(3,1,1); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.15;p(2) = p(2) + 0.07;set(h, 'pos', p);box on; hold on

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
m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
% cmocean('gray')


h = subplot(3,1,2); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.14;p(2) = p(2) + 0.15;set(h, 'pos', p);box on; hold on
pcolor(pos_grid,z(1:1000),v(1:1000,:));shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on

colormap(flipud(brewermap(11,'RdBu')))
hold on;
%contour(pos_grid,z,v,[0 0], 'linecolor', 'w','linewidth',1)
[C,hContour] = contour(pos_grid,z(1:1000),gamma(1:1000,:),g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',1400);
[C,hContour] = contour(pos_grid,z(1:1000),gamma(1:1000,:),[26.4:0.1:27],'k');
[C,hContour] = contour(pos_grid,z(1:1000),gamma(1:1000,:),[27.2 27.3 27.4 27.5 27.7 27.8 28 28.2],'k');

[C,hContour] = contour(pos_grid,z(1:1000),gamma(1:1000,:),[20.35:0.2:26.35],'Color',[.2 .2 .2]);

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[1000 topo(1) topo 1000 ],[.65 .65 .65],'linestyle','none');
%cmocean('balance',10);
caxis([-.55 .55])
ylim([0 999]);xticks;
% set(gca, 'XTickLabel', {[], [], [], [], [], [], [], []});box on
 xticks([-50  -40 -30 -20 -10 0 10]);xticklabels([]);%xticklabels({'[]','[]','[]','[]','[]','[]','[]','[]'});
% yticks([-800 -600 -400 -200 0]);yticklabels({'800','600','400','200','0'});
set(gca,'tickdir','out'); box on
% for l=1:length(stations)
% xline(stations(l),'--','color',[.65 .65 .65]);hold on
% end

plot(lon_stations,0,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
plot(lon_stations(ce_inn),0,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
plot(lon_stations(ae_inn),0,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on


ylim([0 999]);xticks;axis ij;
xticks([-50:5:15]);xticklabels([]);
set(gca,'tickdir','out');box on
set(findall(gcf,'-property','FontSize'),'FontSize',12)
 
h = subplot(3,1,3); p = get(h, 'pos');p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.25; p(2) = p(2) - 0.02;set(h, 'pos', p);box on; hold on
pcolor(pos_grid,z(1000:end),v(1000:end,:));shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on

for l=1:length(stations)
xline(stations(l),'--','color',[.65 .65 .65]);hold on
end
[C,hContour] = contour(pos_grid,z(1000:end),gamma(1000:end,:),g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',1400);
[C,hContour] = contour(pos_grid,z(1000:end),gamma(1000:end,:),[26.4:0.1:27],'k');
[C,hContour] = contour(pos_grid,z(1000:end),gamma(1000:end,:),[27.2 27.3 27.4 27.5 27.7 27.8 28 28.2],'k');

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);

% cmocean('delta')
ylim([1000 5500]);xticks;axis ij;
xticks([-50:5:15]);xticklabels([]);
set(gca,'tickdir','out');box on

[ax,h]=m_contfbar([.38 .68],.13,[-.55:.1:.55],[-.55:.1:.55]);

caxis([-.55 .55])


title(ax,'Meridional Velocity (m.s^-^1)','Fontweight','normal')
set(gcf,'color','w');  % otherwise 'print' turns lakes black
set(gca,'tickdir','out')
xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
xlhh = ylabel('Pressure (db)');xlhh.Position(2) = xlhh.Position(2) - 2000;xlhh.Position(1) = xlhh.Position(1) - 0.2;%clear xlh
set(findall(gcf,'-property','FontSize'),'FontSize',12)

% savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/v_section_dec')
% print(gcf, '-djpeg','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/v_section_dec')


%% zonal vel plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  U %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uplot=u;uplot(uplot>.55)=.55;uplot(uplot<-.55)=-.56;


figure
h = subplot(3,1,1); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.15;p(2) = p(2) + 0.07;set(h, 'pos', p);box on; hold on

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
m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
% cmocean('gray')


h = subplot(3,1,2); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.14;p(2) = p(2) + 0.15;set(h, 'pos', p);box on; hold on

pcolor(pos_grid,z(1:1000),uplot(1:1000,:));shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on

colormap(flipud(brewermap(11,'RdBu')))
hold on;

%contour(pos_grid,z,v,[0 0], 'linecolor', 'w','linewidth',1)
[C,hContour] = contour(pos_grid,z(1:1000),gamma(1:1000,:),g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',2700);
[C,hContour] = contour(pos_grid,z(1:1000),gamma(1:1000,:),[26.4:0.1:27],'k');
[C,hContour] = contour(pos_grid,z(1:1000),gamma(1:1000,:),[27.2 27.3 27.4 27.5 27.7 27.8 28 28.2],'k');

[C,hContour] = contour(pos_grid,z(1:1000),gamma(1:1000,:),[20.35:0.2:26.35],'Color',[.2 .2 .2]);

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[1000 topo(1) topo 1000 ],[.65 .65 .65],'linestyle','none');
%cmocean('balance',10);
caxis([-.55 .55])
ylim([0 999]);xticks;
% set(gca, 'XTickLabel', {[], [], [], [], [], [], [], []});box on
 xticks([-50  -40 -30 -20 -10 0 10]);xticklabels([]);%xticklabels({'[]','[]','[]','[]','[]','[]','[]','[]'});
% yticks([-800 -600 -400 -200 0]);yticklabels({'800','600','400','200','0'});
set(gca,'tickdir','out'); box on
% for l=1:length(stations)
% xline(stations(l),'--','color',[.65 .65 .65]);hold on
% end

plot(lon_stations,0,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
plot(lon_stations(ce_inn),0,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
plot(lon_stations(ae_inn),0,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on


ylim([0 999]);xticks;axis ij;
xticks([-50:5:15]);xticklabels([]);
set(gca,'tickdir','out');box on
set(findall(gcf,'-property','FontSize'),'FontSize',12)

h = subplot(3,1,3); p = get(h, 'pos');p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.25; p(2) = p(2) - 0.02;set(h, 'pos', p);box on; hold on
pcolor(pos_grid,z(1000:end),u(1000:end,:));shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on

for l=1:length(stations)
xline(stations(l),'--','color',[.65 .65 .65]);hold on
end

[C,hContour] = contour(pos_grid,z(1000:end),gamma(1000:end,:),g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',1400);
[C,hContour] = contour(pos_grid,z(1000:end),gamma(1000:end,:),[26.4:0.1:27],'k');
[C,hContour] = contour(pos_grid,z(1000:end),gamma(1000:end,:),[27.2 27.3 27.4 27.5 27.7 27.8 28 28.2],'k');

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);

% cmocean('delta')
ylim([1000 5500]);xticks;axis ij;
xticks([-50:5:15]);xticklabels([]);
set(gca,'tickdir','out');box on

ax=m_contfbar([.38 .68],.13,[-.55:.1:.55],[-.55:.1:.55],'endpiece', 'yes');

ax=m_contfbar([.38 .68],.13,[-.6:.1:.55],[-.5:.1:.5],'endpiece', 'yes');


%  set(ax,'ytick',([-.5 0 .5]),'yticklabel',({'<-0.5','0','.>0.5'}))    

title(ax,'Zonal Velocity (m.s^-^1)','Fontweight','normal')

set(gcf,'color','w');  % otherwise 'print' turns lakes black
set(gca,'tickdir','out')

xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
xlhh = ylabel('Pressure (db)');xlhh.Position(2) = xlhh.Position(2) - 2000;xlhh.Position(1) = xlhh.Position(1) - 0.2;%clear xlh
set(findall(gcf,'-property','FontSize'),'FontSize',12)
caxis([-.55 .55])

% savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/V_section_dec')
% print(gcf, '-djpeg','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/V_section_dec')



%% fluorescence borders and more
%%%%%%%%%%%%%%%%%%%%%% FLUORESCENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure;hold on;
% 
% h = subplot(3,1,1); 
% p = get(h, 'pos');
% p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.15;p(2) = p(2) + 0.07;set(h, 'pos', p);box on; hold on
% 
% m_proj('mercator','long', [-52 18],'lat', [-36 -33]);hold on
% %[CS,CH]=m_etopo2('contourf',[-3500 -200],'edgecolor','none');caxis([-5500 7000]);colormap(flip(gray(7)));  
% m_plot(lon_adcp,lat_adcp,'.k','Markersize',1);hold on
% ll=2;
% 
% for id_time= 1:length(CEs_out60)
%     m_plot(CEs_out60{id_time,id_time,1},CEs_out60{id_time,id_time,2},'c','LineWidth',1)
% end
% 
% for id_time= 1:length(AEs_out60)
%     m_plot(AEs_out60{id_time,id_time,1},AEs_out60{id_time,id_time,2},'m','LineWidth',1)
% end
% 
% m_grid('xtick',[],'ytick',[],'linestyle','none','tickstyle','dd');
% m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
% 
% h = subplot(3,1,2); 
% p = get(h, 'pos');
% p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) - 0.05;p(2) = p(2) + 0.35;set(h, 'pos', p);box on; hold on
% 
% 
% plot(lon_stations,0,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
% plot(lon_stations(ce_inn),0,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
% plot(lon_stations(ae_inn),0,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on
% 
% for l=1:length(stations)
% xline(stations(l),'--','color',[.65 .65 .65]);hold on
% end
% % aa=nsum(cloro_merian,1);hold on
%  plot(lon,maxz_cloro,'k','linewidth',4);hold on
% %scatter(lon,maxz_cloro,3,aa+140);cmocean('algae')
% jj=9;
% 
% scatter(lon(1:jj:end),maxz_cloro(1:jj:end),8,aa(1:jj:end),'filled');cmocean('speed',5)
% 
% 
% axis ij;xlim([-52 18]);caxis([35 85])
% 
% [ax,h]=m_contfbar([.85 .99],.28,[35:10:85],[35:10:85],'endpiece','yes')
% 
% title(ax,'Flourescence','Fontweight','normal')
% xticklabels(ax,{'Low','','High'})
% ylim([0 170])
% 
% xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
% ylabel('Max. Pres. (db)');%xlhh.Position(1) = xlhh.Position(1) - 0.2;xlhh.Position(2) = xlhh.Position(2) - 2000;%clear xlh
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',12)
% set(findall(gcf,'-property','Fontweight'),'Fontweight','normal')
% 
% 
% savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/fluoresence_section_octobers_smoothed')
% print(gcf, '-djpeg','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/fluoresence_section_octobers_smoothed')
% 
%  %% PLOT BORDERS
% % % 
% % % figure('Renderer', 'painters', 'Position', [100 100 550 600]) 
% % % 
% % % h=subplot(2,2,1)
% % % p = get(h, 'pos'); p(1) = p(1)- .07 ; p(2) = p(2) ;p(3) = p(3) + 0.05;p(4) = p(4) + 0.06;set(h, 'pos', p);box on; hold on
% % % ylabel('Pressure')
% % % 
% % % %contourf(x(1:81),z,v(:,1:81),50,'linestyle','none');hold on;%colormap(jet);colorbar;%title('temp');hold on
% % % pcm(x(1:81),z,v(:,1:81));shading interp;hold on;%colormap(jet);colorbar;%title('temp');hold on
% % % 
% % % hold on;
% % % 
% % % 
% % % for l=1:length(stations)
% % % line([stations(l) stations(l)],[0 max_pres(l)],'Color',[.6 .6 .6],'LineStyle','--')
% % % xline(stations(l),'--','color',[.65 .65 .65]);hold on
% % % end
% % % 
% % % fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');
% % % xlim([-52 -48]);ylim([0 4500])
% % % axis ij
% % % cmocean('balance');caxis([-1 1])
% % % 
% % % h=subplot(2,2,2)
% % % p = get(h, 'pos'); p(1) = p(1)- .07 ; p(2) = p(2) ;p(3) = p(3) + 0.05;p(4) = p(4) + 0.06;set(h, 'pos', p);box on; hold on
% % % pcm(x(1:81),z,u(:,1:81));shading interp;hold on;%colormap(jet);colorbar;%title('temp');hold on
% % % hold on;
% % % 
% % % for l=1:length(stations)
% % % line([stations(l) stations(l)],[0 max_pres(l)],'Color',[.6 .6 .6],'LineStyle','--')
% % % xline(stations(l),'--','color',[.65 .65 .65]);hold on
% % % end
% % % 
% % % fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');
% % % xlim([-52 -48]);ylim([0 4500])
% % % axis ij
% % % %colorbar
% % % cmocean('balance');caxis([-1 1])
% % % 
% % % h=subplot(2,2,3)
% % % p = get(h, 'pos'); p(1) = p(1)- .07 ; p(2) = p(2) ;p(3) = p(3) + 0.05;p(4) = p(4) + 0.06;set(h, 'pos', p);box on; hold on
% % % 
% % % pcm(x(end-81:end),z,v(:,end-81:end));shading interp;hold on;%colormap(jet);colorbar;%title('temp');hold on
% % % 
% % % for l=1:length(stations)
% % % line([stations(l) stations(l)],[0 max_pres(l)],'Color',[.6 .6 .6],'LineStyle','--')
% % % xline(stations(l),'--','color',[.65 .65 .65]);hold on
% % % end
% % % 
% % % hold on;
% % % fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');
% % % xlim([14 18]);ylim([0 4500])
% % % axis ij
% % % %colorbar
% % % cmocean('balance');caxis([-1 1])
% % % 
% % % h=subplot(2,2,4)
% % % p = get(h, 'pos'); p(1) = p(1)- .07 ; p(2) = p(2) ;p(3) = p(3) + 0.05;p(4) = p(4) + 0.06;set(h, 'pos', p);box on; hold on
% % % 
% % % pcm(x(end-81:end),z,u(:,end-81:end));shading interp;hold on;%colormap(jet);colorbar;%title('temp');hold on
% % % 
% % % 
% % % for l=1:length(stations)
% % % line([stations(l) stations(l)],[0 max_pres(l)],'Color',[.6 .6 .6],'LineStyle','--')
% % % xline(stations(l),'--','color',[.65 .65 .65]);hold on
% % % end
% % % 
% % % 
% % % hold on;
% % % fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');
% % % xlim([14 18]);ylim([0 4500])
% % % axis ij
% % % %colorbar
% % % cmocean('balance');caxis([-1 1])
% % % 
% % % 
% % % [ax,h]=m_contfbar([.05 .55],.7,[-1:0.1:1],[-1:0.1:1]);
% % % title(ax,'Meridional Velocity(m.s^-^1)','Fontweight','normal')
% % % 
% % % set(gcf,'color','w');  % otherwise 'print' turns lakes black
% % % set(gca,'tickdir','out')
% % % xticks([-52:1:-48]);xticklabels({'52Â°W','51Â°W','50Â°W','49Â°W','48Â°W'})
% % 
% % 
% % 
% % % 
% % % figure
% % % subplot(2,2,1)
% % % %contourf(x(1:81),z,SA(:,1:81),[33.5:0.05:37],'linestyle','none');hold on;%colormap(jet);colorbar;%title('temp');hold on
% % % pcm(x(1:81),z,SA(:,1:81)'
% % % 
% % % for l=1:length(stations)
% % % line([stations(l) stations(l)],[0 max_pres(l)],'Color',[.6 .6 .6],'LineStyle','--')
% % % xline(stations(l),'--','color',[.65 .65 .65]);hold on
% % % end
% % % 
% % % 
% % % hold on;
% % % fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');
% % % xlim([-52 -48]);ylim([0 4500])
% % % axis ij
% % % 
% % % 
% % % subplot(2,2,2)
% % % contourf(x(end-81:end),z,SA(:,end-81:end),50,'linestyle','none');hold on;%colormap(jet);colorbar;%title('temp');hold on
% % % hold on;
% % % fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');
% % % xlim([14 18]);ylim([0 4500])
% % % axis ij
% % % colorbar
% % % cmocean('delta');colorbar
% % % axis ij
% % % colorbar
% % % 
% % % subplot(2,2,3)
% % % contourf(x(1:81),z,CT(:,1:81),[-0.50:0.05:25],'linestyle','none');hold on;%colormap(jet);colorbar;%title('temp');hold on
% % % hold on;
% % % fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');
% % % xlim([-52 -48]);ylim([0 4500]);
% % % axis ij
% % % cmocean('delta');colorbar
% % % 
% % % 
% % % 
% % % subplot(2,2,4)
% % % contourf(x(end-81:end),z,CT(:,end-81:end),[-0.50:0.05:25],'linestyle','none');hold on;%colormap(jet);colorbar;%title('temp');hold on
% % % hold on;
% % % fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');
% % % xlim([14 18]);ylim([0 4500])
% % % axis ij
% % % colorbar
% % % 
% % 
% % 
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%     V %%%%%%%%%%%%%%%%%%%%%%%%
% % %figure
% % % h = subplot(3,1,1); 
% % % p = get(h, 'pos');
% % % p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.15;p(2) = p(2) + 0.07;set(h, 'pos', p);box on; hold on
% % % 
% % % m_proj('mercator','long', [-52 -48],'lat', [-36 -33]);hold on
% % % %[CS,CH]=m_etopo2('contourf',[-3500 -200],'edgecolor','none');caxis([-5500 7000]);colormap(flip(gray(7)));  
% % % m_plot(lon_adcp,lat_adcp,'.k','Markersize',1);hold on
% % % ll=2;
% % % 
% % % for id_time= 1:length(CEs_out60)
% % %     m_plot(CEs_out60{id_time,id_time,1},CEs_out60{id_time,id_time,2},'c','LineWidth',1)
% % % end
% % % 
% % % for id_time= 1:length(AEs_out60)
% % %     m_plot(AEs_out60{id_time,id_time,1},AEs_out60{id_time,id_time,2},'m','LineWidth',1)
% % % end
% % % 
% % % m_grid('xtick',[],'ytick',[],'linestyle','none','tickstyle','dd');
% % % m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
% % % % cmocean('gray')
% % 
% % 
% % 
% % figure
% % pcm(pos_grid,z,SA);shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on
% % 
% % % colormap(flipud(brewermap(21,'RdBu')))
% % 
% % hold on;
% % %contour(pos_grid,z,v,[0 0], 'linecolor', 'w','linewidth',1)
% % gammaa=gamma.*bucket;
% % [C,hContour] = contour(pos_grid,z,gammaa,g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',400);
% % %[C,hContour] = contour(pos_grid,z,gammaa,[26.4:0.1:27],'k');
% % % [C,hContour] = contour(pos_grid,z,gammaa,[27.2 27.3 27.4 27.5 27.7 27.8 28 28.2],'k');
% % % [C,hContour] = contour(pos_grid,z,gammaa,[20.35:0.2:26.35],'Color',[.2 .2 .2]);
% % 
% % fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500 ],[.65 .65 .65],'linestyle','none');
% % %cmocean('balance',10);
% % % caxis([-.85 .85])
% % %ylim([0 999]);
% % xticks;
% % xlim([-52 -48]);ylim([0 4500])
% % % set(gca, 'XTickLabel', {[], [], [], [], [], [], [], []});box on
% %  xticks([-50:1:-48]);xticklabels([]);%xticklabels({'[]','[]','[]','[]','[]','[]','[]','[]'});
% % % yticks([-800 -600 -400 -200 0]);yticklabels({'800','600','400','200','0'});
% % set(gca,'tickdir','out'); box on
% % % for l=1:length(stations)
% % % xline(stations(l),'--','color',[.65 .65 .65]);hold on
% % % end
% % 
% % % plot(lon_stations,0,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
% % % plot(lon_stations(ce_inn),0,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
% % % plot(lon_stations(ae_inn),0,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on
% % 
% % for l=1:length(stations)
% % xline(stations(l),'--','color',[.65 .65 .65]);hold on
% % end
% % 
% % %ylim([0 999]);
% % xticks;axis ij;
% % 
% % 
% % xticks([-52:1:-48]);xticklabels({'52Â°W','51Â°W','50Â°W','49Â°W','48Â°W'});
% % 
% % set(gca,'tickdir','out');box on
% % 
% % [ax,h]=m_contfbar([.38 .68],.13,[-.85:.1:.85],[-.85:.1:.85]);
% % 
% % caxis([-.85 .85])
% % 
% % 
% % title(ax,'Meridional Velocity (m.s^-^1)','Fontweight','normal')
% % set(gcf,'color','w');  % otherwise 'print' turns lakes black
% % set(gca,'tickdir','out')
% % xlhh = ylabel('Pressure (db)');xlhh.Position(2) = xlhh.Position(2) - 2000;xlhh.Position(1) = xlhh.Position(1) - 0.2;%clear xlh
% % set(findall(gcf,'-property','FontSize'),'FontSize',12)
% % 
% % savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/v_section_dec')
% % print(gcf, '-djpeg','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/v_section_dec')
% % 
% % 
% % 
% % 
% % figure
% % pcm(pos_grid,z,v);shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on
% % 
% % colormap(flipud(brewermap(21,'RdBu')))
% % 
% % hold on;
% % %contour(pos_grid,z,v,[0 0], 'linecolor', 'w','linewidth',1)
% % gammaa=gamma.*bucket;
% % [C,hContour] = contour(pos_grid,z,gammaa,g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',400);
% % %[C,hContour] = contour(pos_grid,z,gammaa,[26.4:0.1:27],'k');
% % % [C,hContour] = contour(pos_grid,z,gammaa,[27.2 27.3 27.4 27.5 27.7 27.8 28 28.2],'k');
% % % [C,hContour] = contour(pos_grid,z,gammaa,[20.35:0.2:26.35],'Color',[.2 .2 .2]);
% % 
% % fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500 ],[.65 .65 .65],'linestyle','none');
% % %cmocean('balance',10);
% % caxis([-.85 .85])
% % %ylim([0 999]);
% % xticks;
% % xlim([-52 -48]);ylim([0 4500])
% % % set(gca, 'XTickLabel', {[], [], [], [], [], [], [], []});box on
% %  xticks([-50:1:-48]);xticklabels([]);%xticklabels({'[]','[]','[]','[]','[]','[]','[]','[]'});
% % % yticks([-800 -600 -400 -200 0]);yticklabels({'800','600','400','200','0'});
% % set(gca,'tickdir','out'); box on
% % % for l=1:length(stations)
% % % xline(stations(l),'--','color',[.65 .65 .65]);hold on
% % % end
% % 
% % % plot(lon_stations,0,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
% % % plot(lon_stations(ce_inn),0,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
% % % plot(lon_stations(ae_inn),0,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on
% % 
% % for l=1:length(stations)
% % xline(stations(l),'--','color',[.65 .65 .65]);hold on
% % end
% % 
% % %ylim([0 999]);
% % xticks;axis ij;
% % 
% % 
% % xticks([-52:1:-48]);xticklabels({'52Â°W','51Â°W','50Â°W','49Â°W','48Â°W'});
% % 
% % set(gca,'tickdir','out');box on
% % 
% % [ax,h]=m_contfbar([.08 .48],.13,[-.85:.1:.85],[-.85:.1:.85]);
% % 
% % caxis([-.85 .85])
% % 
% % 
% % title(ax,'Meridional Velocity (m.s^-^1)','Fontweight','normal')
% % set(gcf,'color','w');  % otherwise 'print' turns lakes black
% % set(gca,'tickdir','out')
% % xlhh = ylabel('Pressure (db)');xlhh.Position(2) = xlhh.Position(2) - 2000;xlhh.Position(1) = xlhh.Position(1) - 0.2;%clear xlh
% % set(findall(gcf,'-property','FontSize'),'FontSize',12)
% 
% %%
% 
% oxplot=movmean(ox,40,1,'omitnan');oxplot=movmean(oxplot,20,2,'omitnan');
% 
% figure
% % h=subplot(2,2,1)
% % p = get(h, 'pos'); p(1) = p(1)- .07 ; p(2) = p(2) ;p(3) = p(3) + 0.05;p(4) = p(4) + 0.06;set(h, 'pos', p);box on; hold on
% contourf(pos_grid,z,oxplot,[170:10:260],'linestyle','none');axis ij%colormap(jet);colorbar;%title('temp');hold on
% colormap(brewermap([],'RdPu'));caxis([170 260])
% % clabel( C , h ,'manual') 
% hold on
% [C,hContour] = contour(x,z,gamma,g_n_levels,'--w', 'showtext','on','linewidth',1,'Labelspacing',200);
% set(hContour,'showtext','on');clabel(C,hContour,'FontSize',fontsize,'Fontweight','normal','Color','w');
% fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500 ],[.65 .65 .65],'linestyle','none');
% xticks;
% xlim([-52 -48]);ylim([0 4500])
% for l=1:length(stations)
% xline(stations(l),'--','color',[.65 .65 .65]);hold on
% end
% xticks;axis ij;
% xticks([-52:1:-48]);xticklabels({'52Â°W','51Â°W','50Â°W','49Â°W','48Â°W'});
% set(gca,'tickdir','in');box on
% 
% [C,hContour] = contour(x,z,oxplot,[170:10:270],'k');axis ij
% clabel(C,hContour,'manual','FontSize',fontsize,'Color','k');
% 
% [C,hContour] = contour(x,z,gamma,g_n_levels,'--w', 'showtext','on','linewidth',1,'Labelspacing',2700);
% set(hContour,'showtext','on');clabel(C,hContour,'FontSize',fontsize,'Fontweight','bold','Color','w');
% 
% fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');
% ylim([0 4500]);
% axis ij
% 
% colormap(brewermap([],'RdPu'));caxis([170 260])
% 
% [ax,h]=m_contfbar([.1 .45],.1,[170:10:260],[170:10:260]);
% title(ax,'Dissolved oxygen (Î¼mol.kg^-^1)','Fontweight','normal')
% set(gcf,'color','w');  % otherwise 'print' turns lakes black
% set(gca,'tickdir','out')
% ylabel('Pressure (db)');
% set(findall(gcf,'-property','FontSize'),'FontSize',12)
% hold on;
% 
% savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/western_boundry_oxygen')
% 
% 
% 
% %% more borders
% figure
% % h=subplot(2,2,1)
% % p = get(h, 'pos'); p(1) = p(1)- .07 ; p(2) = p(2) ;p(3) = p(3) + 0.05;p(4) = p(4) + 0.06;set(h, 'pos', p);box on; hold on
% contourf(pos_grid,z,SA,[34:0.025:36.4],'linestyle','none');axis ij%colormap(jet);colorbar;%title('temp');hold on
% cmocean('delta');hold on
% [C,hContour] = contour(x,z,gamma,g_n_levels,'--w', 'showtext','on','linewidth',1,'Labelspacing',200);
% set(hContour,'showtext','on');clabel(C,hContour,'FontSize',fontsize,'Fontweight','normal','Color','w');
% 
% xticks;
% xlim([-52 -48]);ylim([0 4500])
% for l=1:length(stations)
% xline(stations(l),'--','color',[.65 .65 .65]);hold on
% end
% 
% xticks;axis ij;
% caxis([34 36.4])
% 
% xticks([-52:1:-48]);xticklabels({'52Â°W','51Â°W','50Â°W','49Â°W','48Â°W'});
% set(gca,'tickdir','in');box on
% 
% [C,hContour] = contour(x,z,SA,[34:0.1:36.4],'k');axis ij
% clabel(C,hContour,'manual','FontSize',fontsize,'Color','k');
% 
% % [C,hContour] = contour(x,z,gamma,g_n_levels,'--w', 'showtext','on','linewidth',1,'Labelspacing',2700);
% % set(hContour,'showtext','on');clabel(C,hContour,'FontSize',fontsize,'Fontweight','bold','Color','w');
% % 
% fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500 ],[.65 .65 .65],'linestyle','none');
% 
% ylim([0 4500]);
% axis ij
% 
% 
% 
% [ax,h]=m_contfbar([.1 .45],.1,[34:0.1:36.4],[34:0.1:36.4]);
% title(ax,'Absolute Salinity','Fontweight','normal')
% set(gcf,'color','w');  % otherwise 'print' turns lakes black
% set(gca,'tickdir','out')
% ylabel('Pressure (db)');
% set(findall(gcf,'-property','FontSize'),'FontSize',12)
% hold on;
% 
% savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/estern_boundry_sal')
% 
% 
% 
% 
% 
% figure
% % h=subplot(2,2,1)
% % p = get(h, 'pos'); p(1) = p(1)- .07 ; p(2) = p(2) ;p(3) = p(3) + 0.05;p(4) = p(4) + 0.06;set(h, 'pos', p);box on; hold on
% contourf(pos_grid,z,oxplot,[170:10:260],'linestyle','none');axis ij%colormap(jet);colorbar;%title('temp');hold on
% colormap(brewermap([],'RdPu'));caxis([170 260])
% % clabel( C , h ,'manual') 
% hold on
% [C,hContour] = contour(x,z,gamma,g_n_levels,'--w', 'showtext','on','linewidth',1,'Labelspacing',200);
% set(hContour,'showtext','on');clabel(C,hContour,'FontSize',fontsize,'Fontweight','normal','Color','w');
% fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500 ],[.65 .65 .65],'linestyle','none');
% xticks;
% xlim([-52 -48]);ylim([0 4500])
% for l=1:length(stations)
% xline(stations(l),'--','color',[.65 .65 .65]);hold on
% end
% xticks;axis ij;
% xticks([-52:1:-48]);xticklabels({'52Â°W','51Â°W','50Â°W','49Â°W','48Â°W'});
% set(gca,'tickdir','in');box on
% 
% [C,hContour] = contour(x,z,oxplot,[170:10:270],'k');axis ij
% clabel(C,hContour,'manual','FontSize',fontsize,'Color','k');
% 
% [C,hContour] = contour(x,z,gamma,g_n_levels,'--w', 'showtext','on','linewidth',1,'Labelspacing',2700);
% set(hContour,'showtext','on');clabel(C,hContour,'FontSize',fontsize,'Fontweight','bold','Color','w');
% 
% fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');
% ylim([0 4500]);
% axis ij
% 
% colormap(brewermap([],'RdPu'));caxis([170 260])
% 
% [ax,h]=m_contfbar([.1 .45],.1,[170:10:260],[170:10:260]);
% title(ax,'Dissolved oxygen (Î¼mol.kg^-^1)','Fontweight','normal')
% set(gcf,'color','w');  % otherwise 'print' turns lakes black
% set(gca,'tickdir','out')
% ylabel('Pressure (db)');
% set(findall(gcf,'-property','FontSize'),'FontSize',12)
% hold on;
% 
% savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/western_boundry_oxygen')
% % 
% % 
% 
% %%
% figure
% % h=subplot(2,2,1)
% % p = get(h, 'pos'); p(1) = p(1)- .07 ; p(2) = p(2) ;p(3) = p(3) + 0.05;p(4) = p(4) + 0.06;set(h, 'pos', p);box on; hold on
% contourf(pos_grid,z,CT,[-0.5:0.25:22],'linestyle','none');axis ij%colormap(jet);colorbar;%title('temp');hold on
% hold on
% [C,hContour] = contour(x,z,gamma,g_n_levels,'--w', 'showtext','on','linewidth',1,'Labelspacing',200);
% set(hContour,'showtext','on');clabel(C,hContour,'FontSize',fontsize,'Fontweight','normal','Color','w');
% colormap(flipud(brewermap([],'Spectral')));caxis([-0.5 23.5]);ylim([0 1000]);xticks;axis ij;
% 
% xticks;
% xlim([-52 -48]);ylim([0 4500])
% for l=1:length(stations)
% xline(stations(l),'--','color',[.65 .65 .65]);hold on
% end
% 
% xticks;axis ij;
% 
% xticks([-52:1:-48]);xticklabels({'52Â°W','51Â°W','50Â°W','49Â°W','48Â°W'});
% set(gca,'tickdir','in');box on
% 
% [C,hContour] = contour(x,z,CT,[-0.5:0.5:22],'k');axis ij
% clabel(C,hContour,'manual','FontSize',fontsize,'Color','k');
% 
% % [C,hContour] = contour(x,z,gamma,g_n_levels,'--w', 'showtext','on','linewidth',1,'Labelspacing',2700);
% % set(hContour,'showtext','on');clabel(C,hContour,'FontSize',fontsize,'Fontweight','bold','Color','w');
% % 
% fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500 ],[.65 .65 .65],'linestyle','none');
% 
% ylim([0 4500]);
% axis ij
% 
% 
% 
% [ax,h]=m_contfbar([.1 .45],.1,[-0.5:0.5:22],[-0.5:0.5:22]);
% title(ax,'Î˜ (ÂºC)','Fontweight','normal')
% set(gcf,'color','w');  % otherwise 'print' turns lakes black
% set(gca,'tickdir','out')
% ylabel('Pressure (db)');
% set(findall(gcf,'-property','FontSize'),'FontSize',12)
% hold on;
% 
% savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/estern_boundry_temp')
% 
% 
% 
% xticks([14:1:18]);xticklabels({'14Â°E','15Â°E','16Â°E','17Â°E','18Â°E'});
% 

% 
