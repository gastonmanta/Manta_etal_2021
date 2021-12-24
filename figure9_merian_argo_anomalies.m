% load('/Users/gaston/Documents/Phd/Database_South_Atl/New_Remi_database/Argo_proxy.mat')
% 
% load /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/msm60_matrix_final.mat
%  %%%PARTICULARES%%%%%
%  
%  argodates=[Argo_proxy{:,6}]';argox=[Argo_proxy{:,4}]';argoy=[Argo_proxy{:,5}]';
%  
%  %CHOOSE DATES AND AREA TO FIND EDDIES AND ARGOS
%  initdate=datenum(2000,1,31); finaldate=datenum(2019,1,31);
%  Wlim=-53;Elim=18; Nlim=-33.5;Slim=-35.5;
%  
% 
% %applies the restrictions
% ind=find([Argo_proxy{:,6}]>=initdate & [Argo_proxy{:,6}]<=finaldate & [Argo_proxy{:,4}]>= Wlim ...
%     & [Argo_proxy{:,4}]<= Elim & [Argo_proxy{:,5}]>= Slim & [Argo_proxy{:,5}]<= Nlim);
% 
% argo_proxy_cut=Argo_proxy(ind,:);%generate a new restricted data base, hay X perfiles argos en ese dominio
% 
% argox=[argo_proxy_cut{:,4}];argoy=[argo_proxy_cut{:,5}];
% 
% 
% ind=find([argo_proxy_cut{:,16}]==0);
% argo_out_eddie=argo_proxy_cut(ind,:);
% 
% 
% 
% % tt=NaN(198,length(argo_out_eddie));
% % 
% % for i=1:length(argo_out_eddie)
% %         tt(:,i)=argo_out_eddie{i,8};
% % end
% % 
% % maxEl=198;%max number of elements of a cell element in your cell
% % 
% % C=cell2mat(cellfun(@(x) [cell2mat(x) zeros(1,maxEl-numel(x))],c,'uni',0))
% % 
% % % rather awkward cell array of cell arrays of numeric arrays ...
% % % engine
% % c2 = cellfun(@(x) [x{:}], c, 'un',0) % convert to cell array of numeric array
% % 
% % % [m, tf] = padcat(c{:}); % concatenate, pad rows with NaNs
% % % m(~tf) = 0 % replace NaNs by zeros
% % 
% 
% dates=datevec([argo_out_eddie{:,6}])
% ind=find(dates(:,2)<3 | dates(:,2)>11 );
% argo_out_eddie_sum=argo_out_eddie(ind,:);
% 
% 
% c=argo_out_eddie_sum(:,8);tempa=padcat(c{:});
% c=argo_out_eddie_sum(:,9);sala=padcat(c{:});
% 
% c=argo_out_eddie_sum(:,7);deptha=padcat(c{:});
% 
% 
% lona=[argo_out_eddie_sum{:,4}];
% lata=[argo_out_eddie_sum{:,5}];
% 
% % vq = griddata(x,y,z,v,xq,yq,zq)
% % vq = griddata(x,y,z,v,xq,yq,zq)
% % vq = griddata(x,y,v,xq,yq)
% % lonna=repmat(lona,201,1)
% % tempa(isnan(tempa))=0
% 
% 
% 
% 
% temp_argo_merian = griddata(lona,(-2000:10:0),tempa,lon,z(1:2000)');
% sal_argo_merian = griddata(lona,(-2000:10:0),sala,lon,z(1:2000)');
% 
% 
% figure;
% subplot(1,2,1)
% pcolor(temp_argo_merian);shading interp
% colorbar
% subplot(1,2,2)
% pcolor(sal_argo_merian);shading interp
% colorbar
% colormap(jet)
% title('argo climatology along 34.5S without eddies')
% 
% 
% figure;
%  m_proj('mercator','lon',[-60 20], 'lat',[-37.5 -31.5]);hold on
% % m_plot(argox,argoy,'.','Markersize',5,'linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k')
% % m_plot([argo_out_eddie{:,4}],[argo_out_eddie{:,5}],'.','Markersize',5,'linewidth',2,'MarkerEdgeColor','r','MarkerFaceColor','k')
% m_plot([argo_out_eddie_sum{:,4}],[argo_out_eddie_sum{:,5}],'.','Markersize',5,'linewidth',2,'MarkerEdgeColor','m','MarkerFaceColor','k')
% m_grid
% m_coast
% 
% depth_argo=z(1:2000);
% 
% temp_argo_merian=flip(temp_argo_merian);sal_argo_merian=flip(sal_argo_merian);
% 
% 
% 
% tempp=temp(1:1000,:)-temp_argo_merian(1:1000,:);
% 
% sall=sal(1:1000,:)-sal_argo_merian(1:1000,:);
% 
% tempp=temp(1:2000,:)-temp_argo_merian(1:2000,:);
% 
% sall=sal(1:2000,:)-sal_argo_merian(1:2000,:);
% 
% figure;
% subplot(2,3,3)
% pcolor(tempp);shading interp;axis ij;
% title('anomaly');colorbar;caxis([-6 6]);colormap(jet)
% subplot(2,3,1)
% pcolor(temp(1:1000,:));shading interp;axis ij
% title('temperature MSM60');colorbar;caxis([0 25])
% subplot(2,3,2)
% pcolor(temp_argo_merian(1:1000,:));shading interp;axis ij
% title('argo clim without eddies');colorbar;caxis([0 25])
% 
% subplot(2,3,6)
% pcolor(sall);shading interp;axis ij;
% title('anomaly');colorbar;caxis([-1 1]);colormap(jet)
% subplot(2,3,4)
% pcolor(sal(1:1000,:));shading interp;axis ij
% title('salinity MSM60');colorbar;caxis([33.8 36.8])
% subplot(2,3,5)
% pcolor(sal_argo_merian(1:1000,:));shading interp;axis ij
% title('argo clim without eddies');colorbar;caxis([33.8 36.8])
% 
% 
% 
% savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/argo_clim_anom_merian')
% print(gcf, '-djpeg','-r300','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/argo_clim_anom_merian')
% 
% 
% clearvars -except temp_argo_merian sal_argo_merian lon depth_argo
% %save /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/argo_clim_summer_345S
% 
% 
% 
% depth_argo=-depth_argo;
% 
% pres_argo=sw_pres(depth_argo,-35);pres_argo=pres_argo;
% 
% 
% sal_argo=griddata(lon',depth_argo,sal_argo_merian,lon',pres_argo);
% 
% temp_argo=griddata(lon',depth_argo,temp_argo_merian,lon',pres_argo);
% 
% 
% ptemp_argo_merian=sw_ptmp(sal_argo,temp_argo,pres_argo,-35);
% 
% 
% ptemp_argo=fillmissing(ptemp_argo_merian,'linear',2);
% 
% sal_argo=fillmissing(sal_argo,'linear',2);
% 
% 
% %%%% NOW INTERACT WITH THE SECTIONS%%%
% lonn=lon;
% cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian
% load MSM60_ctd_october.mat
% 
% 
% ptemp_ar=griddata(lonn,pres_argo,ptemp_argo,lon,pres(1:2000));
% 
% sal_ar=griddata(lonn,pres_argo,sal_argo,lon,pres(1:2000));
% 
% dens_ar=sw_pden(sal_ar,ptemp_ar,pres(1:2000),1);   
% %sw_pden(S,T,P,PR)
% 
% 
% clearvars -except ptemp_ar sal_ar dens_ar pres_argo
% 
% save argo_merian_gridded ptemp_ar sal_ar dens_ar pres_argo

cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian
% tsg=importdata('msm_060_1_no_headers.tsg');
load stations
load merian_contours
cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian

cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian
load MSM60_ctd_october.mat



load argo_merian_gridded


% AEss={'A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12'}
% 
% CEss={'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13'}

AEss={'1','2','3','4','5','6','7','8','9','10','11','12'}

CEss={'1','2','3','4','5','6','7','8','9','10','11','12','13'}



CEs_lon=[13.5 8 -1 -5 -7.5 -11 -14.5 -17.5 -27.5 -31 -37 -46.5 -51.5];

AEs_lon=[15.6 12 3 -3.9 -9.2 -13 -19.5  -28.9 -33.9 -39.6 -43 -48.6]



water_masses=({'TW','SACW','AAIW','UCDW','NADW','LCDW','AABW'});

% % rr=gamma-roo;
% g_n_levels_sabrina=[19 26.35 27.1 27.6 27.9 28.12 28.22 35];

%following valla et al 2018

g_n_levels=[19 26.35 27.1 27.6 27.9 28.1 28.27 35];
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

dens=sw_pden(SA,CT,pres,1);   
% dens=sw_pden(sal,theta,pres,1);   

dens_anom=dens(1:2000,:)-dens_ar;


% dens_anom_corrected=dens_anom-nmedian(dens_anom(:));

x=lon;pos_grid=x;

% figure;hold on;
% pcolor(lon,pres(1:2000),dens_anom_corrected);shading interp;axis ij
% fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);
% axis ij
% axis tight;ylim([50 1800])
%  caxis([-.55 .55]);cmocean('balance',11)
% colorbar
% 
% print(gcf, '-djpeg','-r300','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/density_anomaly_from_argo_climatology')
% savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/density_anomaly_from_argo_climatology')

% dens_anom_cor=movmean(dens_anom_corrected,10,1);
% dens_anom_cor=movmean(dens_anom_cor,5,10);dens_anom_cor=movmean(dens_anom_cor,10,1);


dens_anom_cor=movmean(dens_anom,5,10);dens_anom_cor=movmean(dens_anom_cor,10,1);


% figure
% h = subplot(2,1,1); 
% p = get(h, 'pos');
% p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.15;p(2) = p(2) - 0.12;set(h, 'pos', p);box on; hold on
% % pcolor(lon,pres(1:2000),dens_anom_corrected);caxis([-.65 .65]);shading interp;axis ij;cmocean('curl');hold on
% 
% contourf(lon,pres(1:2000),dens_anom_cor,[-.6:.05:.6],'linestyle','none');shading interp;axis ij;cmocean('curl');hold on
% 
% 
% C=contour(lon,pres(1:2000),dens_anom_cor,[-.6:.1:.6],'k');shading interp;axis ij
% % 
% % clabel(C,'manual')
% 
% 
% % contourf(lon,pres(1:2000),dens_anom_cor,[-.65:.1:.65],'linestyle','none');shading interp;axis ij
% 
% % 
% fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);
% %colorbar
% axis tight;ylim([0 1800])
% %cmocean('curl',13);%cmocean('diff',13);cmocean('tarn',13)
% 
% clabel(C,'manual')
% %colorbar
% xticks([-50:5:15]);xticklabels({'50Â°W','45Â°W','40Â°W','35Â°W','30Â°W','25Â°W','20Â°W','15Â°W','10Â°W','5Â°W','0Â°E','5Â°E','10Â°E','15Â°E'});
% ylabel('Pressure (db)')
% plot(lon_stations,0,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
% plot(lon_stations(ce_inn),0,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
% plot(lon_stations(ae_inn),0,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on
% 
% [ax,h]=m_contfbar([.15 .375],.13,[-.65:0.1: .65],[-.65:0.1: .65]);
% title(ax,'Ï_Î¸ anomaly (kg.m^-^3)','Fontweight','normal')
% set(gcf,'color','w');  % otherwise 'print' turns lakes black
% set(gca,'tickdir','out')
% 
% 
% print(gcf, '-djpeg','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/density_anomaly_from_argo_climatology')
% savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/density_anomaly_from_argo_climatology')



% tt=theta(1:1500,:);
% 
% 
% ctt=CT(1:1500,:);
% 
% ss=SA(1:1500,:);
% 
% plot(ptemp_ar(:),tt(:),'.')
% 
% plot(sal_ar(:),ss(:),'.')
% 
% % plot(ptemp_ar(:),ctt(:),'.')
% % 
% % plot(sal_ar(:),ss(:),'.')
% 
% offset_argo_sal=nmean(nmean(sal_ar(1:1500,:)-ss))
% 
% offset_argo_ptemp=nmean(nmean(ptemp_ar(1:1500,:)-tt))
% 
% 
% offset_argo_sal=nmean(nmean(sal_ar(1:1500,:)-as(1:1500,:)))
% 
% offset_argo_ptemp=nmean(nmean(ptemp_ar(1:1500,:)-ct(1:1500,:)))
% 
% 
% sal_corrected=sal_ar-offset_argo_sal;
% 
% temp_corrected=ptemp_ar-offset_argo_ptemp;
% 
% theta_anom=theta(1:2000,:)-temp_corrected;
% 
% sal_anom=sal(1:2000,:)-sal_corrected;
% 
% 
% 
theta_anom=CT(1:2000,:)-ptemp_ar;

sal_anom=SA(1:2000,:)-sal_ar;

ae_in=nansum(ae_in);ce_in=nansum(ce_in);
ae_inn=find(ae_in==1);ce_inn=find(ce_in==1);




figure
h = subplot(4,1,1); 
p = get(h, 'pos');
%p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.05;p(2) = p(2) + 0.02;set(h, 'pos', p);box on; hold on
p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.05;p(2) = p(2) + 0.045;set(h, 'pos', p);box on; hold on

m_proj('mercator','long', [-52 18],'lat', [-36 -32]);hold on
[CS,CH]=m_etopo2('contourf',[-3500 -200],'edgecolor','none');caxis([-5500 7000]);colormap(flip(gray(7)));  
m_plot(lon_adcp,lat_adcp,'.k','Markersize',1);hold on

m_plot(lona,lata,'.','color',[.3 .3 .3],'Markersize',3);hold on

ll=2;

for id_time= 1:length(CEs_out60)
    m_plot(CEs_out60{id_time,id_time,1},CEs_out60{id_time,id_time,2},'c','LineWidth',1)
end

for id_time= 1:length(AEs_out60)
    m_plot(AEs_out60{id_time,id_time,1},AEs_out60{id_time,id_time,2},'m','LineWidth',1)
end



for i=1:length(CEs_lon)
m_text(CEs_lon(i),-32.5,CEss{i},'color',[0 .8 .8],'fontsize',11,'fontweight','bold'); hold on

end


for i=1:length(AEs_lon)
m_text(AEs_lon(i),-32.5,AEss{i},'color',[.8 0 .8],'fontsize',11,'fontweight','bold'); hold on
end

m_grid('xtick',[],'ytick',[],'linestyle','none','tickstyle','dd');
m_gshhs_h('patch',[0.81818 0.77647 0.70909]);


h = subplot(4,1,2); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.1;p(2) = p(2) + 0.05;set(h, 'pos', p);box on; hold on

contourf(lon,pres(1:2000),dens_anom,[-.65:.05:.65],'linestyle','none');shading interp;axis ij;hold on

% contour(lon,pres(1:2000),gamma(1:2000,:),[g_n_levels,'k','showtext','on');
contour(lon,pres(1:2000),gamma(1:2000,:),[26.3500   27.1000   27.6000   27.9000   28.1000   28.2700  ],'k','showtext','on','Labelspacing',1400);

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);
axis tight;ylim([0 1800]);cmocean('curl',13);
[ax,h]=m_contfbar([.15 .375],.2,[-.65:0.05: .65],[-.65:0.05: .65]);cmocean('curl',27);
title(ax,' \sigma_0¸ anomaly (kg.m^-^3)','Fontweight','normal')
%clabel(C,'manual')
%colorbar
xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
ylabel('Pressure (dbar)')

plot(lon_stations,0,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
plot(lon_stations(ce_inn),0,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
plot(lon_stations(ae_inn),0,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on

set(gcf,'color','w');  % otherwise 'print' turns lakes black
set(gca,'tickdir','out')
cmocean('curl',27);
text(0.98,0.15,'a','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')


h = subplot(4,1,3); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.095;p(2) = p(2) + 0.00;set(h, 'pos', p);box on; hold on


contourf(lon,pres(1:2000),theta_anom,[-6:0.5:6],'linestyle','none');shading interp;axis ij
hold on;
caxis([-6 6])
contour(lon,pres(1:2000),gamma(1:2000,:),[26.3500   27.1000   27.6000   27.9000   28.1000   28.2700  ],'k','showtext','on','Labelspacing',1400);

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);

axis tight;ylim([0 1800]);cmocean('balance',25);
[ax,h]=m_contfbar([.15 .375],.2,[-6:0.5:6],[-6:0.5:6]);cmocean('balance',13);
title(ax,'\Theta anomaly (°C)','Fontweight','normal')

%clabel(C,'manual')
%colorbar
xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
ylabel('Pressure (dbar)')
%colormap(brewermap(20,'RdBu'));
plot(lon_stations,0,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on

plot(lon_stations(ce_inn),0,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
plot(lon_stations(ae_inn),0,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on

set(gcf,'color','w');  % otherwise 'print' turns lakes black
set(gca,'tickdir','out')

text(0.98,0.15,'b','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')



h = subplot(4,1,4); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.095;p(2) = p(2) - 0.05;set(h, 'pos', p);box on; hold on

contourf(lon,pres(1:2000),sal_anom,[-1:0.1:1],'linestyle','none');shading interp;axis ij
caxis([-1 1]); hold on
contour(lon,pres(1:2000),gamma(1:2000,:),[26.3500   27.1000   27.6000   27.9000   28.1000   28.2700  ],'k','showtext','on','Labelspacing',1400);

fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);

axis tight;ylim([0 1800]);cmocean('tarn',20);

[ax,h]=m_contfbar([.15 .375],.2,[-1:0.1:1],[-1:0.1:1]);cmocean('tarn',20);
title(ax,'S_A anomaly (g.kg^-^1)','Fontweight','normal')
%clabel(C,'manual')
%colorbar
xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
ylabel('Pressure (dbar)')

plot(lon_stations,0,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
plot(lon_stations(ce_inn),0,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
plot(lon_stations(ae_inn),0,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on

set(gcf,'color','w');  % otherwise 'print' turns lakes black
set(gca,'tickdir','out')

text(0.98,0.15,'c','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')
% 
set(findall(gcf,'-property','FontSize'),'FontSize',11.5)
% set(findall(gcf,'-property','Linewidth'),'Linewidth',.7)



% 
   print(gcf, '-dpng','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september/anomaly_from_argo_climatologys')
print(gcf, '-djpeg','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september/anomaly_from_argo_climatologys')

 savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september/anomaly_from_argo_climatologys')

print(gcf, '-djpeg','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september/anomaly_from_argo_climatologys')


% 
