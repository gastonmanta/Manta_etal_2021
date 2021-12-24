
% 
% load('eddies_merian_contours_traj_and_argo_end.mat', 'Anticyclonic_Trajectories')
% load('eddies_merian_contours_traj_and_argo_end.mat', 'Cyclonic_Trajectories')

cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian
load MSM60_ctd_october.mat
 load('eddies_merian_contours_traj_and_argo.mat');

cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian
load dates_crossed
dates_ces_traj60_crossed=[dates_ces_traj60_crossed datenum('09-jan-2017')]
cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian

load figure8_merian


% load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian/lost_traj.mat')
% cd /Users/gaston/Documents/Phd/Database_South_Atl/New_Remi_database_new
%% %FIRST PART OK


AEss={'A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12'}

CEss={'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13'}

CEs_lon=[13.5 8 -1 -5 -7.5 -11 -14.5 -17.5 -27.5 -31 -35 -47.5 -51.5];

AEs_lon=[15.6 13 3 -3.9 -10.2 -13 -19.5  -28.9 -33.9 -39.6 -43 -48.6]


figure('rend','painters','pos',[100 100 800 350]); 

h=subplot(2,1,1)

p = get(h, 'pos');

p(3) = p(3) +0.14 ; p(1) = p(1) - 0.085;p(4) = p(4) + 0.065;p(2) = p(2) + 0.04;set(h, 'pos', p);box on; hold on

text(0.015,0.15,'a','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

% m_proj('mercator','long', [-60 20],'lat', [-41.5 -32]);hold on
m_proj('mercator','long', [-60 20],'lat', [-41.5 -32]);hold on

% m_proj('mercator','long', [-60 20],'lat', [-41.5 -23.5]);hold on

 %[CS,CH]=m_etopo2('contourf',[-3500 -200],'edgecolor','none');caxis([-5500 7000]);colormap(flip(gray(7)));  
%  [CS,CH]=m_etopo2('contour',[-3500 -200],'color',[.5 .5 .5]);%caxis([-5500 7000]);%colormap(flip(gray(7)));  
% [CS,CH]=m_etopo2('contourf',[-3500 -200],'edgecolor','none');caxis([-5500 7000]);colormap(flip(gray(7)));  

m_plot(lon_adcp,lat_adcp,'--k','LineWidth',1); hold on


for id_time= 1:length(CEs_out60)
    m_plot(CEs_out60{id_time,id_time,1},CEs_out60{id_time,id_time,2},'c','LineWidth',1.5)
end



% end
for uhj=1:size(cy_traj_msm60,1)
    m_plot(cy_traj_msm60{uhj,5},cy_traj_msm60{uhj,6},'k','linewidth',1);hold on
%     m_plot(cy_traj_msm60{uhj,5},cy_traj_msm60{uhj,6},'.c','markersize',13);hold on
%  m_plot(cy_traj_msm60{uhj,5},cy_traj_msm60{uhj,6},5,cy_traj_msm60{uhj,7});
end


jj=9000;
 for uhj=1:length(dates_ces_traj60_crossed)

ind=find(dates_ces_traj60_crossed(uhj)==cy_traj_msm60{uhj,7});
% if ind>jj
%   m_plot(cy_traj_msm60{uhj,5}(ind-jj:ind),cy_traj_msm60{uhj,6}(ind-jj:ind),'linewidth',1,'color','b'); hold on
% m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.b','MarkerSize',8);
% else 
 %m_plot(cy_traj_msm60{uhj,5}(1:ind),cy_traj_msm60{uhj,6}(1:ind),'linewidth',1,'color','b'); hold on
m_scatter(cy_traj_msm60{uhj,5}(1:ind),cy_traj_msm60{uhj,6}(1:ind),5,cy_traj_msm60{uhj,7}(1:ind)); hold on

m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.c','MarkerSize',8);

end



% 
% for uhj=1:size(lost_traj,1)
%     m_plot(lost_traj{uhj,5},lost_traj{uhj,6},'c','linewidth',2);hold on
% %     m_plot(cy_traj_msm60{uhj,5},cy_traj_msm60{uhj,6},'.c','markersize',13);hold on
%  m_scatter(lost_traj{uhj,5},lost_traj{uhj,6},5,lost_traj{uhj,7});
% end



% 
% for uhj=1:size(ay_traj_msm60,1)
%     m_plot(ay_traj_msm60{uhj,5},ay_traj_msm60{uhj,6},'m','linewidth',2);hold on
% %     m_plot(ay_traj_msm60{uhj,5},ay_traj_msm60{uhj,6},'.m','markersize',13);hold on
%  m_scatter(ay_traj_msm60{uhj,5},ay_traj_msm60{uhj,6},5,ay_traj_msm60{uhj,7});
% end
% 
% for id_time= 1:length(AEs_out60)
%     m_plot(AEs_out60{id_time,id_time,1},AEs_out60{id_time,id_time,2},'m','LineWidth',1)
% end
% 

for i=1:length(CEs_lon)
m_text(CEs_lon(i),-36.2,CEss{i},'color',[0 1 1],'fontsize',16,'fontweight','bold'); hold on
m_text(CEs_lon(i),-36.2,CEss{i},'color',[.4 .4 .4],'fontsize',16,'linewidth',.25)

end



m_gshhs_h('patch', [0.81818 0.77647 0.70909]);colorbar

m_grid('linestyle','none')
colormap(flip(parula(20)));caxis([736699-150 736699+20])

colorbar;cbdate('mmm/yy');


%colormap(jet(3));caxis([736702-1 736702+1]);cbdate


h=subplot(2,1,2)

p = get(h, 'pos');

p(3) = p(3) +0.14 ; p(1) = p(1) - 0.085;p(4) = p(4) + 0.15;p(2) = p(2) + 0.0;set(h, 'pos', p);box on; hold on

m_proj('mercator','long', [-60 20],'lat', [-40 -23.5]);hold on

 %[CS,CH]=m_etopo2('contourf',[-3500 -200],'edgecolor','none');caxis([-5500 7000]);colormap(flip(gray(7)));  
% [CS,CH]=m_etopo2('contour',[-3500 -200],'color',[.5 .5 .5]);%caxis([-5500 7000]);%colormap(flip(gray(7)));  


% for uhj=1:size(cy_traj_msm60,1)
%     m_plot(cy_traj_msm60{uhj,5},cy_traj_msm60{uhj,6},'c','linewidth',2);hold on
% %     m_plot(cy_traj_msm60{uhj,5},cy_traj_msm60{uhj,6},'.c','markersize',13);hold on
%  m_scatter(cy_traj_msm60{uhj,5},cy_traj_msm60{uhj,6},5,cy_traj_msm60{uhj,7});
% end
% 
% 
% for id_time= 1:length(CEs_out60)
%     m_plot(CEs_out60{id_time,id_time,1},CEs_out60{id_time,id_time,2},'c','LineWidth',1)
% end

m_plot(lon_adcp,lat_adcp,'--k','LineWidth',1); hold on

for id_time= 1:length(AEs_out60)
    m_plot(AEs_out60{id_time,id_time,1},AEs_out60{id_time,id_time,2},'m','LineWidth',1.5)
end


% end
for uhj=1:size(ay_traj_msm60,1)
    m_plot(ay_traj_msm60{uhj,5},ay_traj_msm60{uhj,6},'k','linewidth',1);hold on
%     m_plot(cy_traj_msm60{uhj,5},cy_traj_msm60{uhj,6},'.c','markersize',13);hold on
%  m_plot(cy_traj_msm60{uhj,5},cy_traj_msm60{uhj,6},5,cy_traj_msm60{uhj,7});
end


jj=9000;
 for uhj=1:length(dates_aes_traj60_crossed)

ind=find(dates_aes_traj60_crossed(uhj)==ay_traj_msm60{uhj,7});
% if ind>jj
%   m_plot(cy_traj_msm60{uhj,5}(ind-jj:ind),cy_traj_msm60{uhj,6}(ind-jj:ind),'linewidth',1,'color','b'); hold on
% m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.b','MarkerSize',8);
% else 
 %m_plot(cy_traj_msm60{uhj,5}(1:ind),cy_traj_msm60{uhj,6}(1:ind),'linewidth',1,'color','b'); hold on
m_scatter(ay_traj_msm60{uhj,5}(1:ind),ay_traj_msm60{uhj,6}(1:ind),5,ay_traj_msm60{uhj,7}(1:ind)); hold on

% m_plot(ay_traj_msm60{uhj,5}(ind),ay_traj_msm60{uhj,6}(ind),'.r','MarkerSize',8);

end


% for uhj=1:size(ay_traj_msm60,1)
%     m_plot(ay_traj_msm60{uhj,5},ay_traj_msm60{uhj,6},'k','linewidth',2);hold on
% %     m_plot(ay_traj_msm60{uhj,5},ay_traj_msm60{uhj,6},'.m','markersize',13);hold on
%  m_scatter(ay_traj_msm60{uhj,5},ay_traj_msm60{uhj,6},5,ay_traj_msm60{uhj,7});
% end

for i=1:length(AEs_lon)
m_text(AEs_lon(i),-36.7,AEss{i},'color',[1 0 1],'fontsize',16,'fontweight','bold'); hold on
m_text(AEs_lon(i),-36.7,AEss{i},'color',[.4 .4 .4],'fontsize',16)

end

m_gshhs_h('patch', [0.81818 0.77647 0.70909]);colorbar


m_grid('linestyle','none');
% colormap(flip(parula(20)));caxis([736699-150 736699+20])

cmocean('solar');caxis([734875 736699-150])






cbdate('mmm/yy');
grid off
%colormap(jet(3));caxis([736702-1 736702+1]);cbdate

%m_plot(8.70,-35.875,'or','markerfacecolor','r','markersize',6)

text(0.008,0.12,'b','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

set(findall(gcf,'-property','FontSize'),'FontSize',12)


% cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september



% print(gcf, '-depsc','-painters' ,'-r600','figure11_september_eddy_trajs3')
% 
% 
% print(gcf, '-dpng','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/rho_anom')

% % savefig('Fig_8merian')