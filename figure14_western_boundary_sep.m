cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian
% tsg=importdata('msm_060_1_no_headers.tsg');
load stations
load merian_contours
cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian

cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian
load MSM60_ctd_october.mat

load argo_merian_gridded

water_masses=({'TW','SACW','AAIW','UCDW','NADW','LCDW','AABW'});

pos_grid=lon;
pos_gridd=pos_grid(:,1:end-1)+0.05/2;

cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian
load figure_7a


cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/

load new_AEs_max60
load geostrophic_vel_281


v_new=v.*bucket;

u_new=u.*bucket;


ind=find(dens_anom>0.6);dens_anom(ind)=0.6;
ind=find(dens_anom<-0.6);dens_anom(ind)=-0.6;



smooth=20;
oxx=movmean(ox,smooth,1);oxx=movmean(oxx,smooth,2);oxx=movmean(oxx,smooth,1);

saa=movmean(SA,smooth*3,1);saa=movmean(saa,smooth,2);saa=movmean(SA,smooth*3,1);

% ind=find(saa(:,1:40)>35.05 & saa(:,1:40)<35.2);saa(ind)=35.1;
% ind=find(saa(:,1:40)<35.05 & saa(:,1:40)>35);saa(ind)=35;

g_n_levels=[ 26.35 27.1 27.6 27.9 28.1 28.27];


gamma=movmedian(gamma,5,1);



%% figure

% figure
% 
% subplot(2,2,1)
figure('Renderer', 'painters', 'Position', [200 200 700 550]) 

h = subplot(2,2,1); 
p = get(h, 'pos');
p(3) = p(3) + 0.06; p(1) = p(1) - 0.06;p(4) = p(4) + 0.06;p(2) = p(2) + 0.08;set(h, 'pos', p);box on; hold on

m_proj('mercator','long', [-52 -46.4],'lat', [-35.499 -33.501]);hold on
text(0.02,0.13,'a','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

m_plot(lon_adcp,lat_adcp,'k');hold on
ll=2;
for id_time= 1:length(CEs_out60)
    m_plot(CEs_out60{id_time,id_time,1},CEs_out60{id_time,id_time,2},'c','LineWidth',1)
end
% 
% for id_time= 1:length(CEs_max60)
%     m_plot(CEs_max60{id_time,id_time,1},CEs_max60{id_time,id_time,2},'b','LineWidth',1)
% end

for id_time= 15:length(AEs_out60)
    m_plot(AEs_out60{id_time,id_time,1},AEs_out60{id_time,id_time,2},'m','LineWidth',1)
end
jj=120;

 for uhj=15:length(dates_ces_traj60_crossed)

ind=find(dates_ces_traj60_crossed(uhj)==cy_traj_msm60{uhj,7});
if ind>jj
%   m_plot(cy_traj_msm60{uhj,5}(ind-jj:ind),cy_traj_msm60{uhj,6}(ind-jj:ind),'linewidth',1,'color','k'); hold on
m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.','MarkerSize',18,'color','c');
else 
%  m_plot(cy_traj_msm60{uhj,5}(1:ind),cy_traj_msm60{uhj,6}(1:ind),'linewidth',1,'color','k'); hold on
m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.','MarkerSize',18,'color','c');

end
end


hold on
m_contour(lon_ssh,lat_ssh,ssh_composite',[-.2:0.05:1],'color',[.55 .55 .55]);%colormap(jet(ww));

h2=m_quiver(lonn38,latt38,uu38,vv38,1,'LineWidth',1,'Color','k','ShowArrowhead','on','AutoScale','on','AutoScaleFactor',0.8)
h3=m_quiver(lonn38,latt38,uvgos,vvgos,1,'LineWidth',1,'Color',[.5 .5 .5],'ShowArrowhead','on','AutoScale','on','AutoScaleFactor',0.8)

m_plot(lon_stations,lat_stations,'dk','Markersize',4);hold on

set(h2,'AutoScale','on', 'AutoScaleFactor',1.2)
m_grid('xticklabels',[],'xtick',[-60:1:20], 'ytick',[-35.5:0.5:-33.5],'linestyle','none');


axes('position',[0.0657,0.5041,0.3986,0.2414]);

contourf(pos_grid(1:113),pres(1:2000),dens_anom(1:2000,1:113),[-.6:.025:.6],'linestyle','none');shading interp;axis ij;hold on

ylim([0 1300]);xticks;

% [C,hContour] = contour(pos_grid,z(1:1500),gamma(1:1500,:),[25:0.2:29.5],'k');

[C,hContour] = contour(pos_grid,z(1:1500),gamma(1:1500,:),[22:0.15:26.2 26.5:0.15:26.95 27.25:0.15:28],'color',[.6 .6 .6],'linewidth',.4);
[C,hContour] = contour(pos_grid,z(1:5000),gamma(1:5000,:),g_n_levels,'color',[.1 .1 .1], 'showtext','on','linewidth',1,'labelspacing',350);


caxis([-.6 .6]);cmocean('curl',13);

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(113)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');
xlim([-52 -46.4]);ylim([0 1200])
xlim([-52 -46.6])
set(gcf,'color','w');  % otherwise 'print' turns lakes black
set(gca,'tickdir','out')

% 
% cb=colorbar;
%  cb.Position = cb.Position + 1e-10;
ylabel('Pressure (dbar)')

cmocean('curl',13);

[ax,h]=m_contfbar([0.07 .475],.2,[-.6:.1:.6],[-.6:.1:.6]);%cmocean('curl',14);

title(ax,'\sigma_0 anomaly (kg.m^-^3)','Fontweight','normal')
%\rho. ?.

text(0.02,0.13,'b','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')
xticklabels([])


h =subplot(2,2,3)
p = get(h, 'pos');
p(3) = p(3) + 0.06; p(1) = p(1) - 0.06;p(4) = p(4) + 0.1;p(2) = p(2) - 0.06;

set(h, 'pos', p);box on; hold on
pcolor(pos_gridd(1:113),pres(1:4800),u_new(1:4800,1:113));shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on

ylabel('Pressure (dbar)')
%pcolor(pos_gridd(1:200),pres(1:1500),geostrophic_velocity_281_ek_balanced(1:1500,1:200));shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on

hold on

%cmocean('balance',21);caxis([-.8 .8]);%colorbar

for l=1:length(stations)
xline(stations(l),'--','color',[.65 .65 .65]);hold on
end
% [C,hContour] = contour(pos_grid,z(1:1500),gamma(1:1500,:),g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',400);


fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(113)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');
xlim([-52 -46.4]);ylim([0 4800])

xticklabels({'52ºW','51ºW','50ºW','49ºW','48ºW','47ºW'})

text(0.02,0.07,'d','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

%text(0.13,0.08,'Zonal-LADCP','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','none')


cmocean('balance',21);
[ax,h]=m_contfbar([0.07 .475],.1,[-0.8:0.1:0.8],[-0.8:0.1:0.8]);
title(ax,'U_A_D_C_P_s(m.s^-^1)','Fontweight','normal')
caxis([-.75 .75]);%colorbar

%colormap(flipud(brewermap(21,'RdBu')))
hold on;



h = subplot(2,2,2); 
p = get(h, 'pos');
p(3) = p(3) + 0.04; p(1) = p(1) - 0.04;p(4) = p(4) + 0.08;p(2) = p(2) - 0.05;set(h, 'pos', p);box on; hold on

%pcolor(pos_gridd(1:200),pres(1:1500),geostrophic_velocity_281_ek_balanced(1:1500,1:200));shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on
pcolor(pos_gridd(1:113),pres(1:4800),geostrophic_velocity_281_ek_balanced(1:4800,1:113));shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on


hold on

cmocean('balance',21);caxis([-.8 .8]);%colorbar
% for l=1:length(stations)
% xline(stations(l),'--','color',[.65 .65 .65]);hold on
% end
% [C,hContour] = contour(pos_grid,z(1:1500),gamma(1:1500,:),g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',400);
% hold on
% [C,hContour] = contour(pos_gridd(1:113),pres(1:4800),gamma(1:4800,1:113),[26:0.2:29],'k');
hold on
% [C,hContour] = contour(pos_gridd(1:113),pres(1:4800),gamma(1:4800,1:113),[26.35 27.1 27.6 27.9 28.1 28.27],'k','labelspacing',300);
% 
% th = clabel(C,hContour);


[C,hContour] = contour(pos_grid,z,gamma,[22:0.15:26.2 26.5:0.15:26.95 27.25:0.15:28],'color',[.6 .6 .6],'linewidth',.4);
[C,hContour] = contour(pos_grid,z(1:3500),gamma(1:3500,:),g_n_levels(1:end-1),'color',[.1 .1 .1], 'showtext','on','linewidth',1,'labelspacing',300); hold on
[C,hContour] = contour(pos_grid,z(3950:end),gamma(3950:end,:),[g_n_levels(end) g_n_levels(end)],'color',[.1 .1 .1], 'showtext','on','linewidth',1,'labelspacing',140);


for l=1:length(stations)
xline(stations(l),'--','color',[.65 .65 .65]);hold on
end
% [C,hContour] = contour(pos_grid(1:113),z(1000:4800),saa(1000:4800,1:113),[34.9 34.9],'k', 'showtext','on','labelspacing',1400);
% 
% [C,hContour] = contour(pos_grid(1:113),z(3950:4800),CT(3950:4800,1:113),[-2 0],'r', 'showtext','on','labelspacing',100);



fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(113)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');
xlim([-52 -46.4]);ylim([0 4800])

ylabel('Pressure (dbar)')



cmocean('balance',21);
[ax,h]=m_contfbar([0.07 .475],.1,[-0.8:0.1:0.8],[-0.8:0.1:0.8]);
title(ax,'V_g_e_o_s (m.s^-^1)','Fontweight','normal')
caxis([-.75 .75]);%colorbar


% 
% cb=colorbar;
% pos=get(cb,'Position');
% set(cb,'Position',pos+[0.06,0,0,0]);
% ylabel(cb, 'Velocity (m.s^-^1)')

text(0.03,0.08,'c','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

%text(0.13,0.08,'Meridional-Geostrophy','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','none')

xticklabels({'52ºW','51ºW','50ºW','49ºW','48ºW','47ºW'})

h = subplot(2,2,4); 
p = get(h, 'pos');
p(3) = p(3) + 0.04; p(1) = p(1) - 0.04;p(4) = p(4) + 0.09;p(2) = p(2) - 0.05;set(h, 'pos', p);box on; hold on


% p(3) = p(3) + 0.06; p(1) = p(1) - 0.06;p(4) = p(4) + 0.06;p(2) = p(2) - 0.05;set(h, 'pos', p);box on; hold on

%pcolor(pos_gridd(1:200),pres(1:1500),geostrophic_velocity_281_ek_balanced(1:1500,1:200));shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on

v_cut=v_new(1001:end,:);v_cutt=horzcat(v_cut(:,1:24),v_cut(:,24),v_cut(:,25),v_cut(:,25:28),v_cut(:,28),v_cut(:,29),v_cut(:,29:end));

v_neww=vertcat(v_new(1:1000,1:113),v_cutt(:,1:113));



pcolor(pos_gridd(1:113),pres(1:4800),v_neww(1:4800,:));shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on

for l=1:length(stations)
xline(stations(l),'--','color',[.65 .65 .65]);hold on
end

[C,hContour] = contour(pos_grid,z,gamma,[22:0.15:26.2 26.5:0.15:26.95 27.25:0.15:28],'color',[.6 .6 .6],'linewidth',.4);
[C,hContour] = contour(pos_grid,z(1:3500),gamma(1:3500,:),g_n_levels(1:end-1),'color',[.1 .1 .1], 'showtext','on','linewidth',1,'labelspacing',300); hold on
[C,hContour] = contour(pos_grid,z(3950:end),gamma(3950:end,:),[g_n_levels(end) g_n_levels(end)],'color',[.1 .1 .1], 'showtext','on','linewidth',1,'labelspacing',140);




hold on
cmocean('balance',21);caxis([-.8 .8]);%colorbar
% cb=colorbar;
% pos=get(cb,'Position');
% set(cb,'Position',pos+[0.06,0,0,0]);
% ylabel(cb, 'Velocity (m.s^-^1)')

% [C,hContour] = contour(pos_grid(1:113),z(1:4800),SA(1:4800,1:113),[34.5 34.9 35.1],'k', 'showtext','on','linewidth',1,'labelspacing',400);

% 
% [C,hContour] = contour(pos_grid(1:113),z(800:4800),oxx(800:4800,1:113),[220 240],'k', 'showtext','on','labelspacing',400);
% 
% 
% [C,hContour] = contour(pos_grid(1:113),z(3950:4800),CT(3950:4800,1:113),[-2 0],'r', 'showtext','on','labelspacing',100);

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(113)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');
xlim([-52 -46.4]);ylim([0 4800])

xticklabels({'52ºW','51ºW','50ºW','49ºW','48ºW','47ºW'})

% for l=1:length(stations)
% xline(stations(l),'--','color',[.65 .65 .65]);hold on
% end
% [C,hContour] = contour(pos_grid,z(1:1500),gamma(1:1500,:),g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',400);
% hold on
% [C,hContour] = contour(pos_gridd(1:113),pres(1:4800),gamma(1:4800,1:113),[26:0.2:29],'k');

text(0.03,0.08,'e','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

% text(0.13,0.08,'Meridional-LADCP','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','none')

cmocean('balance',21);
[ax,h]=m_contfbar([0.07 .475],.1,[-0.8:0.1:0.8],[-0.8:0.1:0.8]);
title(ax,'V_A_D_C_P_s (m.s^-^1)','Fontweight','normal')
caxis([-.75 .75]);%colorbar

ylabel('Pressure (dbar)')

set(findall(gcf,'-property','FontSize'),'FontSize',11)
set(findall(gcf,'-property','Linewidth'),'Linewidth',.6)

set(gcf,'color','w');

print(gcf, '-dpng','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september/figure14_september_western_boundary')

print(gcf, '-djpeg','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september/figure14_september_western_boundary')

savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september/figure_western_boundary_sep')

%close all


