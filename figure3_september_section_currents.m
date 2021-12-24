%% LOAD DATA
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

gammas=gamma;

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
% Note. The zonally averaged potential temperature (Î¸), salinity (S), and dis- solved ygen (O2) values at the center of each range computed from all cruises occupied at the SAMBA-W line are given as a reference.
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

fontsize=13;
gamma_font_color=[0 0 0];
gamma_color=[0.1 0.1 0.1];

gammas=gamma.*bucket;

gamma=gamma.*bucket;
sep=25;
%%%%%%%%%%%%%%%%%%%%%% PLOT SECTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CONSTRAIN DATA

 u(u>.5)=.5;u(u<-.5)=-.5;
v(v>.5)=.5;v(v<-.5)=-.5;
load corrientes_colormap_final.mat

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

%contourf(pos_grid,z(1:1000),v(1:1000,:),[-8:0.05:8],'linestyle','none');%shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on

%contourf(pos_grid,z(1:1000),v(1:1000,:),[-8:.1:-.5 -.5:0.05:.5 .5:0.1:8],'linestyle','none');%shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on
axis ij
pcolor(pos_grid,z(1:1000),v(1:1000,:));shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on

%colormap(flipud(brewermap(11,'RdBu')))

% % ax = gca;
% % corrientes_colormap_final = colormap(ax);
% % save('corrientes_colormap_final','corrientes_colormap_final')



load corrientes_colormap_final.mat
 colormap(corrientes_colormap_final)
% cmocean('balance',21)
%colormap(mymap)
hold on;
% contour(pos_grid,z,v,[0 0], 'linecolor', 'w')
% hold on
% [CSs,CHs]=contour(pos_grid,z(1:1000),v(1:1000,:),[-.5 .5],'color',gamma_color);%shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on
% clabel(CSs,CHs,'manual','FontSize',fontsize,'Color',gamma_color);
% 
% 

% [C,hContour] = contour(pos_grid,z(1:1000),gamma(1:1000,:),g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',1400);
% clabel(C,hContour,'FontSize',fontsize,'Color','k');


[C,hContour] = contour(pos_grid(1:250),z(1:1000),gammas(1:1000,(1:250)),g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',150);
clabel(C,hContour,'FontSize',fontsize,'Color','k');
hold on
[C,hContour] = contour(pos_grid(250:end),z(1:1000),gammas(1:1000,(250:end)),g_n_levels,'k','linewidth',1);


% [C,hContour] = contour(pos_grid,z(1:1000),gamma(1:1000,:),[26.4:0.1:27],'k');
% [C,hContour] = contour(pos_grid,z(1:1000),gamma(1:1000,:),[27.2 27.3 27.4 27.5 27.7 27.8 28 28.2],'k');
% 
% [C,hContour] = contour(pos_grid,z(1:1000),gamma(1:1000,:),[20.35:0.2:26.35],'Color',[.2 .2 .2]);

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[1000 topo(1) topo 1000 ],[.65 .65 .65],'linestyle','none');
%cmocean('balance',10);
caxis([-.5 .5])
ylim([-sep 999]);xticks;
% set(gca, 'XTickLabel', {[], [], [], [], [], [], [], []});box on
 xticks([-50  -40 -30 -20 -10 0 10]);xticklabels([]);%xticklabels({'[]','[]','[]','[]','[]','[]','[]','[]'});
% yticks([-800 -600 -400 -200 0]);yticklabels({'800','600','400','200','0'});
set(gca,'tickdir','out'); box on
% for l=1:length(stations)
% xline(stations(l),'--','color',[.65 .65 .65]);hold on
% end

plot(lon_stations,-sep,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
plot(lon_stations(ce_inn),-sep,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
plot(lon_stations(ae_inn),-sep,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on


ylim([-sep 1000]);xticks;axis ij;
xticks([-50:5:15]);xticklabels([]);
set(gca,'tickdir','out');box on
set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)

 
h = subplot(3,1,3); p = get(h, 'pos');p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.25; p(2) = p(2) - 0.02;set(h, 'pos', p);box on; hold on

%contourf(pos_grid,z(1000:end),v(1000:end,:),[-.8:0.05:.8],'linestyle','none');%shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on
%[CS,CH]=contourf(pos_grid,z(1000:end),v(1000:end,:),[-.8:.1:-.5 -.45:0.05:.45 .5:0.1:.81],'linestyle','none');%shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on


% [CS,CH]=contour(pos_grid,z(1000:end),v(1000:end,:),[-0.2 -0.1 0.1 0.2],'linestyle','none');
% clabel(CS,CH,'manual','FontSize',fontsize,'Color',gamma_color);



pcolor(pos_grid,z(1000:end),v(1000:end,:));shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on
axis ij
% for l=1:length(stations)
% xline(stations(l),'--','color',[.65 .65 .65]);hold on
% end

%  contour(pos_grid,z(1000:end),v(1000:end,:),[0 0],'color',[.9 .9 .9]);%shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);


[C,hContour] = contour(pos_grid(1:250),z(1000:4500),gammas(1000:4500,(1:250)),g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',150);
clabel(C,hContour,'FontSize',fontsize,'Color','k');
hold on
[C,hContour] = contour(pos_grid(250:end),z(1000:4500),gammas(1000:4500,(250:end)),g_n_levels,'k','linewidth',1);


%  [C,hContour] = contour(pos_grid,z(1000:end),gammas(1000:end,:),g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',4000);
%  clabel(C,hContour,'FontSize',fontsize,'Color','k');
% 

% [C,hContour] = contour(pos_grid,z(1000:end),gamma(1000:end,:),g_n_levels,'k','linewidth',1);
% clabel(C,hContour,'manual','FontSize',fontsize,'Color','k');


% [C,hContour] = contour(pos_grid,z(1000:end),gamma(1000:end,:),[26.4:0.1:27],'k');
% [C,hContour] = contour(pos_grid,z(1000:end),gamma(1000:end,:),[27.2 27.3 27.4 27.5 27.7 27.8 28 28.2],'k');


% cmocean('delta')
ylim([1001 5500]);xticks;axis ij;
xticks([-50:5:15]);xticklabels([]);
set(gca,'tickdir','out');box on



% colormap(mymap2)



caxis([-.5 .5])

%[ax,h]=m_contfbar([.38 .68],.15,CS,CH,'endpiece','none');


% [ax,h]=m_contfbar([.38 .68],.15,[-1:.05:1],[-1:.05:1]);

box on


%title(ax,'Meridional Velocity (m.s^-^1)','Fontweight','normal')
set(gcf,'color','w');  % otherwise 'print' turns lakes black
set(gca,'tickdir','out')

xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
xlhh = ylabel('Pressure (dbar)');xlhh.Position(2) = xlhh.Position(2) - 2000;xlhh.Position(1) = xlhh.Position(1) - 0.2;%clear xlh


set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)


text(0.013,0.145,'a','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',20, 'Edgecolor','k')

% savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september/figure3_september_section_currentsv')
% print(gcf, '-dpng','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september/figure3_september_section_currentsv')
% print(gcf, '-djpeg','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september/figure3_september_section_currentsv')


%% u 


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

%contourf(pos_grid,z(1:1000),v(1:1000,:),[-8:0.05:8],'linestyle','none');%shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on

%contourf(pos_grid,z(1:1000),v(1:1000,:),[-8:.1:-.5 -.5:0.05:.5 .5:0.1:8],'linestyle','none');%shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on
axis ij
pcolor(pos_grid,z(1:1000),u(1:1000,:));shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on

%colormap(flipud(brewermap(11,'RdBu')))
load corrientes_colormap_final.mat
 colormap(corrientes_colormap_final)
% cmocean('balance',21)
%colormap(mymap)
hold on;
% contour(pos_grid,z,v,[0 0], 'linecolor', 'w')
% hold on
% [CSs,CHs]=contour(pos_grid,z(1:1000),v(1:1000,:),[-.5 .5],'color',gamma_color);%shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on
% clabel(CSs,CHs,'manual','FontSize',fontsize,'Color',gamma_color);
% 
% 

% [C,hContour] = contour(pos_grid,z(1:1000),gamma(1:1000,:),g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',1400);
% clabel(C,hContour,'FontSize',fontsize,'Color','k');


[C,hContour] = contour(pos_grid(1:250),z(1:1000),gammas(1:1000,(1:250)),g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',150);
clabel(C,hContour,'FontSize',fontsize,'Color','k');
hold on
[C,hContour] = contour(pos_grid(250:end),z(1:1000),gammas(1:1000,(250:end)),g_n_levels,'k','linewidth',1);


% [C,hContour] = contour(pos_grid,z(1:1000),gamma(1:1000,:),[26.4:0.1:27],'k');
% [C,hContour] = contour(pos_grid,z(1:1000),gamma(1:1000,:),[27.2 27.3 27.4 27.5 27.7 27.8 28 28.2],'k');
% 
% [C,hContour] = contour(pos_grid,z(1:1000),gamma(1:1000,:),[20.35:0.2:26.35],'Color',[.2 .2 .2]);

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[1000 topo(1) topo 1000 ],[.65 .65 .65],'linestyle','none');
%cmocean('balance',10);
caxis([-.5 .5])
ylim([-sep 999]);xticks;
% set(gca, 'XTickLabel', {[], [], [], [], [], [], [], []});box on
 xticks([-50  -40 -30 -20 -10 0 10]);xticklabels([]);%xticklabels({'[]','[]','[]','[]','[]','[]','[]','[]'});
% yticks([-800 -600 -400 -200 0]);yticklabels({'800','600','400','200','0'});
set(gca,'tickdir','out'); box on
% for l=1:length(stations)
% xline(stations(l),'--','color',[.65 .65 .65]);hold on
% end

plot(lon_stations,-sep,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
plot(lon_stations(ce_inn),-sep,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
plot(lon_stations(ae_inn),-sep,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on


ylim([-sep 1000]);xticks;axis ij;
xticks([-50:5:15]);xticklabels([]);
set(gca,'tickdir','out');box on
set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)

 
h = subplot(3,1,3); p = get(h, 'pos');p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.25; p(2) = p(2) - 0.02;set(h, 'pos', p);box on; hold on

%contourf(pos_grid,z(1000:end),v(1000:end,:),[-.8:0.05:.8],'linestyle','none');%shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on
%[CS,CH]=contourf(pos_grid,z(1000:end),v(1000:end,:),[-.8:.1:-.5 -.45:0.05:.45 .5:0.1:.81],'linestyle','none');%shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on


% [CS,CH]=contour(pos_grid,z(1000:end),v(1000:end,:),[-0.2 -0.1 0.1 0.2],'linestyle','none');
% clabel(CS,CH,'manual','FontSize',fontsize,'Color',gamma_color);



pcolor(pos_grid,z(1000:end),u(1000:end,:));shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on
axis ij
% for l=1:length(stations)
% xline(stations(l),'--','color',[.65 .65 .65]);hold on
% end

%  contour(pos_grid,z(1000:end),v(1000:end,:),[0 0],'color',[.9 .9 .9]);%shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);


[C,hContour] = contour(pos_grid(1:250),z(1000:4500),gammas(1000:4500,(1:250)),g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',150);
clabel(C,hContour,'FontSize',fontsize,'Color','k');
hold on
[C,hContour] = contour(pos_grid(250:end),z(1000:4500),gammas(1000:4500,(250:end)),g_n_levels,'k','linewidth',1);


%  [C,hContour] = contour(pos_grid,z(1000:end),gammas(1000:end,:),g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',4000);
%  clabel(C,hContour,'FontSize',fontsize,'Color','k');
% 

% [C,hContour] = contour(pos_grid,z(1000:end),gamma(1000:end,:),g_n_levels,'k','linewidth',1);
% clabel(C,hContour,'manual','FontSize',fontsize,'Color','k');


% [C,hContour] = contour(pos_grid,z(1000:end),gamma(1000:end,:),[26.4:0.1:27],'k');
% [C,hContour] = contour(pos_grid,z(1000:end),gamma(1000:end,:),[27.2 27.3 27.4 27.5 27.7 27.8 28 28.2],'k');


% cmocean('delta')
ylim([1001 5500]);xticks;axis ij;
xticks([-50:5:15]);xticklabels([]);
set(gca,'tickdir','out');box on



% colormap(mymap2)



caxis([-.5 .5])

%[ax,h]=m_contfbar([.38 .68],.15,CS,CH,'endpiece','none');


% [ax,h]=m_contfbar([.38 .68],.15,[-1:.05:1],[-1:.05:1]);

box on


%title(ax,'Meridional Velocity (m.s^-^1)','Fontweight','normal')
set(gcf,'color','w');  % otherwise 'print' turns lakes black
set(gca,'tickdir','out')

xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
xlhh = ylabel('Pressure (dbar)');xlhh.Position(2) = xlhh.Position(2) - 2000;xlhh.Position(1) = xlhh.Position(1) - 0.2;%clear xlh


set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)


text(0.013,0.145,'b','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',20, 'Edgecolor','k')
% 
% savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september/figure3_september_section_currentsu')
% print(gcf, '-dpng','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september/figure3_september_section_currentsu')
% print(gcf, '-djpeg','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september/figure3_september_section_currentsu')
% 

