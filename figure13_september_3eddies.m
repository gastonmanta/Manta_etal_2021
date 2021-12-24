 %% ALtimery data
% cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian
% tsg=importdata('msm_060_1_no_headers.tsg');
% 
% load stations
% load merian_contours
% 
% 
% %ncdisp('sst_merian_MSM602.nc')
% sst=ncread('sst_merian_MSM602.nc','analysed_sst');sst=sst-273.18;sst=flip(sst,1);%1 january to 1 february included
% lon_sst=ncread('sst_merian_MSM602.nc','lon');lon_sst=flip(lon_sst);
% lat_sst=ncread('sst_merian_MSM602.nc','lat');
% time_sst2=ncread('sst_merian_MSM602.nc','time');
% 
% % cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian
% % load eddies_merian_contours_traj_and_argo
% 
% 
% msst=nanmean(sst,3);msubsst=sst(:,:,end-1)-sst(:,:,1);
% 
% 
% dates=tsg(:,1:6);daten=datenum(dates);
% lon_tsg=tsg(:,8);lat_tsg=tsg(:,7);sst_tsg=nanmean(tsg(:,9:10),2);
% 
% 
% figure;%hold on
% m_proj('ortho','lat',-34.5','long',-15');
% m_plot(lon_tsg,lat_tsg,'.r','Markersize',2);
% m_coast('patch', [.7 .7 .7]);
% % m_grid('xtick',[-34.5 -34.5])xticklabels',[],'yticklabels',[])
%  m_grid('xticklabels',[],'xticks',[],'yticklabels',[],'yticks',[])
% 
% 
% m_grid('linest','-','xticklabels',[],'yticklabels',[]);
% 
% 
% 
% 
% cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/ssh_merian
% aa=dir('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/ssh_merian')
% ncfiles = dir('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/ssh_merian/*.nc');   % get all nc files in the folder 
% nfiles = length(ncfiles) ;   % total number of files 
% % ncdisp('dt_global_allsat_phy_l4_20170101_20170530.nc')
% 
% % loop for each file 
% 
% for K = 1 : nfiles
%   filename = ncfiles(K).name;  %just the name
%   %%get the vairable you want from nc file 
%   adt(:,:,K)= ncread(filename,'adt') ; 
%   u(:,:,K)= ncread(filename,'ugos') ; 
%     v(:,:,K)= ncread(filename,'vgos') ; 
%    lon_ssh= ncread(filename,'longitude') ;
%    lat_ssh= ncread(filename,'latitude') ;% doc ncread 
%   %%Append the data where you want 
% end
% 
% 
% 
% cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian
% 
% 
% 
% 
% lon_ssh(lon_ssh>180)=lon_ssh(lon_ssh>180)-360;
% [lon_ssh,id_sort]=sort(lon_ssh);
% adt=adt(id_sort,:,:);
% 
% u=u(id_sort,:,:);
% v=v(id_sort,:,:);
% 
% 
% 
% adt=flip(adt,1);lon_ssh=flip(lon_ssh);
% u=flip(u,1);v=flip(v,1);
% 
% % lon_ssh=lon_ssh(480:800);lat_ssh=lat_ssh(120:320);adt=adt(480:800,120:320,:);
% 
% 
% %create a composite of 1km sst of sst days changing along cruise track
% day_lon=dates(:,3);
% ind=find(day_lon(2:end)-day_lon(1:end-1)>0);
% lon_date=lon_tsg(ind);%lon day are the longitudes when a day change
% lon_date=round(lon_date);%round it to 2 decimal places to match sst 1km longiudes
% date_lon=day_lon(ind);
% 
% yw=find(ismember(lon_sst,lon_date));
% yy=vertcat(1,yw);yy(end)=length(lon_sst);
% 
% ywh=find(ismember(lon_ssh,lon_date-.125));
% 
% yyh=vertcat(1,ywh);yyh(end)=length(lon_ssh);
% 
% 
% % for i=1:27
% %  sst_composite{i} = squeeze(sst(yy(i):yy(i+1)-1,:,date_lon(i)));
% % end
% %%%
% for i=1:27
%  sst_composite{i} = squeeze(sst(yy(i):yy(i+1)-1,:,date_lon(i)));
%   ssh_composite{i} = squeeze(adt(yyh(i):yyh(i+1)-1,:,date_lon(i)));
%     u_composite{i} = squeeze(u(yyh(i):yyh(i+1)-1,:,date_lon(i)));
%         v_composite{i} = squeeze(v(yyh(i):yyh(i+1)-1,:,date_lon(i)));
% end
% 
%  sst_composite= vertcat(sst_composite{:}); sst_composite=vertcat(sst_composite,sst_composite(end,:));
%  ssh_composite= vertcat(ssh_composite{:}); ssh_composite=vertcat(ssh_composite,ssh_composite(end,:));
%  u_composite= vertcat(u_composite{:}); u_composite=vertcat(u_composite,u_composite(end,:));
%  v_composite= vertcat(v_composite{:}); v_composite=vertcat(v_composite,v_composite(end,:));
% 
% clear sst
% 
% b200=load('200.dat');b1000=load('1000.dat');
% 
% 
% load /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/4gaston/vmadcp/msm60os75000_000000_39_hc.mat
% 
% v=squeeze(b.vel(:,2,:));
% u=squeeze(b.vel(:,1,:));
% zc=c.depth(:);
% pg=c.pg;%8 meter interval from 17m to 800
% kk=find(pg<=25);
% v(kk)=NaN;
% u(kk)=NaN;
% 
% %mean between 40 and 100m from 38 khz adcp
% uplot=squeeze(nanmean(u(4:11,:),1));
% vplot=squeeze(nanmean(v(4:11,:),1));
% lon_adcp=b.nav.txy1(2,:);lat_adcp=b.nav.txy1(3,:);
% 
% 
% 
% spacing = 36/2;
% smoothing=36;%1 is 10 minutes of data
% uuplot=movmean(uplot,smoothing);vvplot=movmean(vplot,smoothing);
% lon38=lon_adcp(1:spacing:end);
% lat38=lat_adcp(1:spacing:end);
% 
% u38=uuplot(1:spacing:end);
% v38=vvplot(1:spacing:end);
% 
% 
% 
% [lon_sshx,lat_sshy]=meshgrid(lon_ssh,lat_ssh);
% lx=lon_sshx';ly=lat_sshy';
% 
% lat_d=[-37.5 -31.5];
% lon_d=[lon_date lon_date];
% 
% 
% ae_in=nansum(ae_in);ce_in=nansum(ce_in);
% ae_inn=find(ae_in==1);ce_inn=find(ce_in==1);
% 
% 
% load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian/eddies_merian.mat')
% 
% 
% lon_ssh=double(lon_ssh);lat_ssh=double(lat_ssh);u_composite=double(u_composite);
% 
% ugos=griddata(lon_ssh,lat_ssh',u_composite',lon_adcp,lat_adcp);
% 
% vgos=griddata(lon_ssh,lat_ssh',v_composite',lon_adcp,lat_adcp);
% 
% % figure;plot(lon_adcp,vgos);
% % hold on
% % plot(lon_adcp,vvplot);axis tight; ylim([-.7 .7])
% 
% 
% uvgos=ugos(1:spacing:end);
% vvgos=vgos(1:spacing:end);
% 
% 
% % hold on
% % plot
% % ugos=squeeze(nmean(u_composite(:,223:224),2))
% % 
% %  ugos=griddata(lon_sshx',lat_sshy,u_composite,lon_adcp,lat_adcp);
% 
% 
% 
% %%%%%%%%%%FIGURES%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% % lonn38=[lon38 -59.9];latt38=[lat38 -31.7];uu38=[u38 0.5];vv38=[v38 0];
% % ww=30;
% % 
% 
% %adds a refernce arrow
% %lonn38=[lon38 -60.2];latt38=[lat38 -31.3];uu38=[u38 0.5];vv38=[v38 0];
% %
% 
% ww=30;
% 
% 
% 
% %yi; comes from the colocalization script, and it's the index in with the
% %CEs have trajectories. The other's don't
% cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian
% load('eddies_merian_contours_traj_and_argo.mat')
%  cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian
% % load('eddies_merian_with_lost_traj.mat')
% load('lost_traj.mat')
% 
% 
% yi=[1;2;5;7;8;10;11;12;13;14;15;16;17;18;19;20];
% 
% date_num_CEs_out60_yi=date_num_CEs_out60(yi);
% 
% for i=1:length(yi)
% ind(i)=find(date_num_CEs_out60_yi(i)==cy_traj_msm60{i,7});
% dates_ces_traj60_crossed(i)=cy_traj_msm60{i,7}(ind(i));
% end
% 
% 
% yi=[1:16];
% date_num_AEs_out60_yi=date_num_AEs_out60(yi);
% dates_aes_traj60_crossed=date_num_AEs_out60;
% 
% % 
% % for i=1:length(yi)
% % ind(i)=find(date_num_AEs_out60_yi(i)==ay_traj_msm60{i,7});
% % dates_ces_traj60_crossed(i)=cy_traj_msm60{i,7}(ind(i));
% % end
% 
% 
% 
% %yi; comes from the colocalization script, and it's the index in with the
% %CEs have trajectories. The other's don't
% % yi=[1;2;5;7;8;10;11;12;13;14;15;16;17;18;19;20];
% % date_num_CEs_out60_yi=date_num_CEs_out60(yi);
% % date_num_CEs_out60_yi=date_num_CEs_out60(yi);
% % 
% % date_num_CEs_out60_yi=cy_traj_msm60{i,7}(yi);
% 
% % 
% % for i=1:length(yi)
% % ind(i)=find(date_num_CEs_out60_yi(i)==cy_traj_msm60{i,7});
% % dates_ces_traj60_crossed(i)=cy_traj_msm60{i,7}(ind(i));
% % end
% 
% 
% % lonn38=[lon38 -59.9];latt38=[lat38 -31.7];uu38=[u38 0.5];vv38=[v38 0];
% % ww=30;
% % 
% 
% spacing = 36/2;
% smoothing=36;%1 is 10 minutes of data
% uuplot=movmean(uplot,smoothing);vvplot=movmean(vplot,smoothing);
% lon38=lon_adcp(1:spacing:end);
% lat38=lat_adcp(1:spacing:end);
% u38=uuplot(1:spacing:end);
% v38=vvplot(1:spacing:end);
% % 
% % lonn38=[lon38 -59.9];latt38=[lat38 -31.7];uu38=[u38 0.5];vv38=[v38 0];uvgos=[uvgos 0.5];vvgos=[vvgos 0];
% % ww=30;
% 
% %adds a refernce arrow
% %lonn38=[lon38 -60.2];latt38=[lat38 -31.3];uu38=[u38 0.5];vv38=[v38 0];
% lonn38=[lon38 -58.9];latt38=[lat38 -30.7];uu38=[u38 0.5];vv38=[v38 0];uvgos=[uvgos 0.5];vvgos=[vvgos 0];
% 
% 
% ww=30;
% 
% 
% %dates_aes_traj60_crossed=ay_traj_msm60{i,7};
% 
% 
% latmin=-39;latmax=-30;
% 
% 
% 
% clear msubsst msst kk dates 
% 
% %% section data
% cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian
% load MSM60_ctd_october.mat
% 
% load cloro_merian
% 
% load vu;v=v/100;u=u/100;
% 
% load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian/stations.mat')
%  
% ae_in=nansum(ae_in);ce_in=nansum(ce_in);
% ae_inn=find(ae_in==1);ce_inn=find(ce_in==1);
% 
% cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian
% % tsg=importdata('msm_060_1_no_headers.tsg');
% load stations
% load merian_contours
% cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian
% 
% 
% % cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian
% % tsg=importdata('msm_060_1_no_headers.tsg');
% % 
% % sst=ncread('sst_merian_MSM60.nc','analysed_sst');sst=sst-273.18;sst=flip(sst,1);%1 january to 1 february included
% % lon_sst=ncread('sst_merian_MSM60.nc','lon');lon_sst=flip(lon_sst);
% % lat_sst=ncread('sst_merian_MSM60.nc','lat');%topo=sw_pres(topo',-35);
% 
% 
% water_masses=({'TW','SACW','AAIW','UCDW','NADW','LCDW','AABW'});
% 
% % rr=gamma-roo;
% g_n_levels_sabrina=[19 26.35 27.1 27.6 27.9 28.12 28.22 35];
% 
% %following valla et al 2018
% 
% g_n_levels=[19 26.35 27.1 27.6 27.9 28.1 28.1 28.27 35];
% % Waterdepth(121)=1040;
% 
% % lonn=lon;
% % lonn(lonn<0)=lonn(lonn<0)+360
% % lon_ssh(lon_ssh>180)=lon_ssh(lon_ssh>180)-360;
% % [gamma_n, gamma_error_lower, gamma_error_upper] = eos80_legacy_gamma_n(sal_new(1:end-10,:),temp_new(1:end-10,:),press(1:end-10,:),lonn',lat')
% 
% % just fot the plotting of the ctd locations
% distance_between_stations=sw_dist(lon,lat,'km')*1000;
% pp=cumsum(distance_between_stations);pp=pp(end);
% % xs=distance_between_stations*123/pp;xs=cumsum(xs); xs=[0; xs];
% % xs=xs+1;
% % ceros=(1:124)*0;
% % Waterdepthh=Waterdepth;Waterdepthh(Waterdepthh>1000)=1000;
% % Waterdepthhh=Waterdepth;Waterdepthhh(Waterdepthhh<1000)=1000;
% 
% [SA, in_ocean] = gsw_SA_from_SP(sal,pres,lon,-34.5);
% 
% CT = gsw_CT_from_t(SA,temp,pres);
% 
% gamma_GP = gamma_GP_from_SP_pt(SA,CT,pres,lon,-34.5);
% 
% gamma=gamma.*bucket;
% 
% gamma_GP=gamma_GP.*bucket;
% 
% z=pres;x=lon;pos_grid=lon;
% 
% fontsize=14;
% 
% %% argo data
% 
% 
% load argo_merian_gridded
% 
% 
% % AEss={'A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12'}
% % 
% % CEss={'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13'}
% 
% AEss={'1','2','3','4','5','6','7','8','9','10','11','12'};
% 
% CEss={'1','2','3','4','5','6','7','8','9','10','11','12','13'};
% 
% 
% 
% CEs_lon=[13.5 8 -1 -5 -7.5 -11 -14.5 -17.5 -27.5 -31 -37 -46.5 -51.5];
% 
% AEs_lon=[15.6 12 3 -3.9 -9.2 -13 -19.5  -28.9 -33.9 -39.6 -43 -48.6];
% 
% 
% 
% water_masses=({'TW','SACW','AAIW','UCDW','NADW','LCDW','AABW'});
% 
% % % rr=gamma-roo;
% % g_n_levels_sabrina=[19 26.35 27.1 27.6 27.9 28.12 28.22 35];
% 
% %following valla et al 2018
% 
% g_n_levels=[19 26.35 27.1 27.6 27.9 28.1 28.27 35];
% % Waterdepth(121)=1040;
% 
% % lonn=lon;
% % lonn(lonn<0)=lonn(lonn<0)+360
% % lon_ssh(lon_ssh>180)=lon_ssh(lon_ssh>180)-360;
% % [gamma_n, gamma_error_lower, gamma_error_upper] = eos80_legacy_gamma_n(sal_new(1:end-10,:),temp_new(1:end-10,:),press(1:end-10,:),lonn',lat')
% 
% % just fot the plotting of the ctd locations
% distance_between_stations=sw_dist(lon,lat,'km')*1000;
% pp=cumsum(distance_between_stations);pp=pp(end);
% % xs=distance_between_stations*123/pp;xs=cumsum(xs); xs=[0; xs];
% % xs=xs+1;
% % ceros=(1:124)*0;
% % Waterdepthh=Waterdepth;Waterdepthh(Waterdepthh>1000)=1000;
% % Waterdepthhh=Waterdepth;Waterdepthhh(Waterdepthhh<1000)=1000;
% 
% % [SA, in_ocean] = gsw_SA_from_SP(sal,pres,lon,-34.5);
% % 
% % CT = gsw_CT_from_t(SA,temp,pres);
% % 
% % gamma_GP = gamma_GP_from_SP_pt(SA,CT,pres,lon,-34.5);
% % 
% % 
% % gamma=gamma.*bucket;
% % 
% % gamma_GP=gamma_GP.*bucket;
% 
% dens=sw_pden(SA,CT,pres,1);   
% % dens=sw_pden(sal,theta,pres,1);   
% 
% dens_anom=dens(1:2000,:)-dens_ar;
% 
% 
% % dens_anom_corrected=dens_anom-nmedian(dens_anom(:));
% 
% x=lon;pos_grid=x;
% 
% 
% %dens_anom_cor=movmean(dens_anom,5,10);dens_anom_cor=movmean(dens_anom_cor,10,1);
% 
% 
% theta_anom=CT(1:2000,:)-ptemp_ar;
% 
% sal_anom=SA(1:2000,:)-sal_ar;
% 
% ae_in=nansum(ae_in);ce_in=nansum(ce_in);
% ae_inn=find(ae_in==1);ce_inn=find(ce_in==1);
% 


%%

cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian
load figure_7a


cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/

load new_AEs_max60
load geostrophic_vel_281

ind=find(dens_anom>0.6);dens_anom(ind)=0.6;ind=find(dens_anom<-0.6);dens_anom(ind)=-0.6;


gamma=movmedian(gamma,5,1);
cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian


%% the figure, finally


lon1w=-32.8;
lon1e=-28.8;

lon2w=-14.5;
lon2e=-10.5;

lon3w=11.75 ;
lon3e=15.75;

cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian
figure('Renderer', 'painters', 'Position', [200 200 800 750]) ; ha = tight_subplot(4,3,[.01 .01],[.05 .01],[.06 .08])

%% 1st part
qq=-0.5;
axes(ha(1)); 

m_proj('mercator','long', [lon1w lon1e],'lat', [-35.5 -33.5]);hold on

text(0.95,0.145,'a','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')


m_plot(lon_stations,lat_stations,'dk','Markersize',4);hold on


m_plot(lon_adcp,lat_adcp,'k');hold on

ll=2;

for id_time= 15%1:length(CEs_out60)
    m_plot(CEs_out60{id_time,id_time,1},CEs_out60{id_time,id_time,2},'c','LineWidth',1)
end


for id_time=15% :length(CEs_max60)
    m_plot(CEs_max60{id_time,id_time,1},CEs_max60{id_time,id_time,2},'b','LineWidth',1)
end

for id_time= 1:length(AEs_out60)
    m_plot(AEs_out60{id_time,id_time,1},AEs_out60{id_time,id_time,2},'m','LineWidth',1)
end


% m_text(lon1e-0.7,-33.7,'C10','color','b','fontsize',16,'fontweight','bold'); hold on

title('C10')

jj=120;


 for uhj=1:length(dates_ces_traj60_crossed)

ind=find(dates_ces_traj60_crossed(uhj)==cy_traj_msm60{uhj,7});
if ind>jj
%   m_plot(cy_traj_msm60{uhj,5}(ind-jj:ind),cy_traj_msm60{uhj,6}(ind-jj:ind),'linewidth',1,'color','k'); hold on
m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.','MarkerSize',18,'color','b');
else 
%  m_plot(cy_traj_msm60{uhj,5}(1:ind),cy_traj_msm60{uhj,6}(1:ind),'linewidth',1,'color','k'); hold on
m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.','MarkerSize',18,'color','b');

end
end


hold on
m_contour(lon_ssh,lat_ssh,ssh_composite',[-.2:0.05:1],'color',[.55 .55 .55]);%colormap(jet(ww));


h2=m_quiver(lonn38,latt38,uu38,vv38,1,'LineWidth',1.2,'Color','k','ShowArrowhead','on','AutoScale','on','AutoScaleFactor',0.8)
h3=m_quiver(lonn38,latt38,uvgos,vvgos,1,'LineWidth',1.2,'Color',[.5 .5 .5],'ShowArrowhead','on','AutoScale','on','AutoScaleFactor',0.8)
set(h2,'AutoScale','on', 'AutoScaleFactor',1.2)

%m_grid('xticklabels',[],'xtick',[-60:1:20], 'ytick',[-35.5:0.5:-33.5],'linestyle','none');
m_grid('xtick',[-60:1:20], 'ytick',[-35.5:0.5:-33.5],'linestyle','none');


axes(ha(2)); 


%m_proj('mercator','long', [ -14 -11],'lat', [-35.5 -33.5]);hold on

m_proj('mercator','long', [lon2w lon2e],'lat', [-35.5 -33.5]);hold on

%[CS,CH]=m_etopo2('contourf',[-3500 -200],'edgecolor','none');caxis([-5500 7000]);colormap(flip(gray(7)));  
% m_plot(lon_adcp,lat_adcp,'.k','Markersize',1);hold on

text(0.95,0.145,'b','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')


m_plot(lon_stations,lat_stations,'dk','Markersize',4);hold on

m_plot(lon_adcp,lat_adcp,'k');hold on

ll=2;

for id_time= 2:length(CEs_out60)
    m_plot(CEs_out60{id_time,id_time,1},CEs_out60{id_time,id_time,2},'c','LineWidth',1)
end


for id_time= 2:length(CEs_max60)
    m_plot(CEs_max60{id_time,id_time,1},CEs_max60{id_time,id_time,2},'b','LineWidth',1)
end


for id_time= 1:length(AEs_out60)
    m_plot(AEs_out60{id_time,id_time,1},AEs_out60{id_time,id_time,2},'m','LineWidth',1)
end


for id_time= 1:length(AEs_max60n)
    m_plot(AEs_max60n{id_time,id_time,1},AEs_max60n{id_time,id_time,2},'r','LineWidth',1)
end



jj=120;

 for uhj=1:length(dates_aes_traj60_crossed)

ind=find(dates_aes_traj60_crossed(uhj)==ay_traj_msm60{uhj,7});
if ind>jj
%   m_plot(ay_traj_msm60{uhj,5}(ind-jj:ind),ay_traj_msm60{uhj,6}(ind-jj:ind),'color','k','linewidth',1); hold on
m_plot(ay_traj_msm60{uhj,5}(ind),ay_traj_msm60{uhj,6}(ind),'.','MarkerSize',18,'color','r');
else 
%  m_plot(ay_traj_msm60{uhj,5}(1:ind),ay_traj_msm60{uhj,6}(1:ind),'m','linewidth',1,'color','k'); hold on
m_plot(ay_traj_msm60{uhj,5}(ind),ay_traj_msm60{uhj,6}(ind),'.','MarkerSize',18,'color','r');
end
end

 for uhj=1:length(dates_ces_traj60_crossed)

ind=find(dates_ces_traj60_crossed(uhj)==cy_traj_msm60{uhj,7});
if ind>jj
%   m_plot(cy_traj_msm60{uhj,5}(ind-jj:ind),cy_traj_msm60{uhj,6}(ind-jj:ind),'linewidth',1,'color','k'); hold on
m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.','MarkerSize',18,'color','b');
else 
%  m_plot(cy_traj_msm60{uhj,5}(1:ind),cy_traj_msm60{uhj,6}(1:ind),'linewidth',1,'color','k'); hold on
m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.','MarkerSize',18,'color','b');

end
end


title('A6')


hold on
m_contour(lon_ssh,lat_ssh,ssh_composite',[-.2:0.05:1],'color',[.55 .55 .55]);%colormap(jet(ww));

h2=m_quiver(lonn38,latt38,uu38,vv38,1,'LineWidth',1.2,'Color','k','ShowArrowhead','on','AutoScale','on','AutoScaleFactor',0.8)
h3=m_quiver(lonn38,latt38,uvgos,vvgos,1,'LineWidth',1.2,'Color',[.5 .5 .5],'ShowArrowhead','on','AutoScale','on','AutoScaleFactor',0.8)

%  m_plot(lon_stations,lat_stations,'dk','Markersize',4,'MarkerFacecolor','k');hold on
% 
%  m_plot(lon_stations(ae_inn),lat_stations(ae_inn),'dm','Markersize',5,'MarkerFacecolor','m');hold on
%  m_plot(lon_stations(ce_inn),lat_stations(ce_inn),'dc','Markersize',5,'MarkerFacecolor','c');hold on

% h2 = quiver(x, y, dxda, dyda)
set(h2,'AutoScale','on', 'AutoScaleFactor',1)

%m_grid('xticklabels',[],'xtick',[-60:1:20], 'ytick',[-35.5:0.5:-33.5],'yticklabels',[],'linestyle','none');

%
m_grid('xtick',[-60:1:20], 'ytick',[-35.5:0.5:-33.5],'yticklabels',[],'linestyle','none');
%m_grid('xtick',[-60:1:20], 'ytick',[-35.5:0.5:-33.5],'yaxislocation','right','linestyle','none');

% 0.070172811059908,0.420238095643569,0.435483870967742,0.158888401581088
% 0.1268,0.4591,0.3375,0.1081

% 
% axes('position',[0.089285714285714,0.460829493087558,0.375014285714285,0.08294930875576]);
% 
% 
% plot(lon_ssh,ssh_composite(:,222));
% xlim([-14 -11]);xticklabels({})
% 
% ylim([0.5 0.85])
% 
% ylabel('ADT(m)')
% 


axes(ha(3));
 
%m_proj('mercator','long', [12 15.5],'lat', [-35.5 -33.5]);hold on

m_proj('mercator','long', [lon3w lon3e],'lat', [-35.5 -33.5]);hold on

text(0.95,0.145,'c','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')


m_plot(lon_adcp,lat_adcp,'k');hold on
ll=2;
for id_time= 2%:length(CEs_out60)
    m_plot(CEs_out60{id_time,id_time,1},CEs_out60{id_time,id_time,2},'c','LineWidth',1)
end

for id_time= 2%:length(CEs_max60)
    m_plot(CEs_max60{id_time,id_time,1},CEs_max60{id_time,id_time,2},'b','LineWidth',1)
end

for id_time= 1:length(AEs_out60)
    m_plot(AEs_out60{id_time,id_time,1},AEs_out60{id_time,id_time,2},'m','LineWidth',1)
end
jj=120;

 for uhj=1%:length(dates_ces_traj60_crossed)

ind=find(dates_ces_traj60_crossed(uhj)==cy_traj_msm60{uhj,7});
if ind>jj
%   m_plot(cy_traj_msm60{uhj,5}(ind-jj:ind),cy_traj_msm60{uhj,6}(ind-jj:ind),'linewidth',1,'color','k'); hold on
m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.','MarkerSize',18,'color','b');
else 
%  m_plot(cy_traj_msm60{uhj,5}(1:ind),cy_traj_msm60{uhj,6}(1:ind),'linewidth',1,'color','k'); hold on
m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.','MarkerSize',18,'color','b');

end
end

hold on
m_contour(lon_ssh,lat_ssh,ssh_composite',[-.2:0.05:1],'color',[.55 .55 .55]);%colormap(jet(ww));

h2=m_quiver(lonn38,latt38,uu38,vv38,1,'LineWidth',1.2,'Color','k','ShowArrowhead','on','AutoScale','on','AutoScaleFactor',0.8)
h3=m_quiver(lonn38,latt38,uvgos,vvgos,1,'LineWidth',1.2,'Color',[.5 .5 .5],'ShowArrowhead','on','AutoScale','on','AutoScaleFactor',0.8)

set(h2,'AutoScale','on', 'AutoScaleFactor',1.2)

m_plot(lon_stations,lat_stations,'dk','Markersize',4);hold on

%m_grid('xticklabels',[],'xtick',[-60:1:20], 'ytick',[-35.5:0.5:-33.5],'yaxislocation','right','linestyle','none');

%m_grid('xtick',[-60:1:20], 'ytick',[-35.5:0.5:-33.5],'yaxislocation','right','linestyle','none');


m_grid('xtick',[-60:1:20], 'ytick',[-35.5:0.5:-33.5],'yticklabels',[],'linestyle','none');

title('C1')



%% 2nd part


axes(ha(4));hold on
text(0.92,0.125,'d','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

pcolor(pos_grid,z(1:1500),v(1:1500,:));shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on

colormap(flipud(brewermap(11,'RdBu')))
hold on;

caxis([-.75 .75])

xlim([lon1w lon1e])

ylim([0 1300]);xticks;

[C,hContour] = contour(pos_grid,z(1:1500),gamma(1:1500,:),[22:0.15:26.2 26.5:0.15:26.95 27.25:0.15:28],'color',[.6 .6 .6],'linewidth',.4);
[C,hContour] = contour(pos_grid,z(1:5000),gamma(1:5000,:),g_n_levels,'color',[.1 .1 .1], 'showtext','on','linewidth',1,'labelspacing',350);

ylabel('Pressure (dbar)')

xticks([-34 -33  -32 -31 -30 -29 -28]);set(gca,'Xticklabel',[]) ; box on
% xticklabels({'33¬∫W','32¬∫W','31¬∫W','30¬∫W','29¬∫W'})

yticks([0 200 400 600 800 1000 1200]);set(gca,'yticklabels',[0:200:1200]); box on


box on
grid on



axes(ha(5)); hold on
text(0.92,0.125,'e','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')


pcolor(pos_grid,z(1:1500),v(1:1500,:));shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on

colormap(flipud(brewermap(11,'RdBu')))
hold on;
caxis([-.8 .8])
xlim([lon2w lon2e]);

ylim([0 1300]);


[C,hContour] = contour(pos_grid,z(1:1500),gamma(1:1500,:),[22:0.15:26.2 26.5:0.15:26.95 27.25:0.15:28],'color',[.6 .6 .6],'linewidth',.4);
[C,hContour] = contour(pos_grid,z(1:5000),gamma(1:5000,:),g_n_levels,'color',[.1 .1 .1], 'showtext','on','linewidth',1,'labelspacing',240);


yticks([0 200 400 600 800 1000 1200]);set(gca,'yticklabels',[]); box on

xticks([-14 -13 -12 -11]);set(gca,'Xticklabel',[]) ; box on


axes(ha(6)); hold on

text(0.92,0.125,'f','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

pcolor(pos_grid,z(1:1500),v(1:1500,:));shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on

colormap(flipud(brewermap(11,'RdBu')))
hold on;

caxis([-.75 .75])

xlim([lon3w lon3e])

ylim([0 1300]);xticks;

[C,hContour] = contour(pos_grid,z(1:1500),gamma(1:1500,:),[22:0.15:26.2 26.5:0.15:26.95 27.25:0.15:28],'color',[.6 .6 .6],'linewidth',.4);
[C,hContour] = contour(pos_grid,z(1:5000),gamma(1:5000,:),g_n_levels,'color',[.1 .1 .1], 'showtext','on','linewidth',1,'labelspacing',1400);

xticks([12 13 14 15]);set(gca,'Xticklabel',[]) ; box on

yticks([0 200 400 600 800 1000 1200]);set(gca,'Yticklabel',[]) ; box on

cb=colorbar;

pos=get(cb,'Position');set(cb,'Position',pos+[0.049,0,0,0]);
xlhg = ylabel(cb, 'Vadcp (m.s^-^1)');xlhg.Position(1) = xlhg.Position(1) - 0.2;%xlh.Position(2) = xlh.Position(2) - 2000;clear xlh.Position(1);clear xlh.Position(2)



axes(ha(7));hold on

text(0.92,0.125,'g','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

pos_gridd=pos_grid(:,1:end-1)+0.05/2;

pcolor(pos_gridd,z(1:1500),geostrophic_velocity_281_ek_balanced(1:1500,:));shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on

hold on


for l=1:length(stations)
xline(stations(l),'--','color',[.65 .65 .65]);hold on
end
%contour(pos_grid,z,v,[0 0], 'linecolor', 'w','linewidth',1)

% [C,hContour] = contour(pos_grid,z(1:1500),gamma(1:1500,:),g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',400);

[C,hContour] = contour(pos_grid,z(1:1500),gamma(1:1500,:),[22:0.15:26.2 26.5:0.15:26.95 27.25:0.15:28],'color',[.6 .6 .6],'linewidth',.4);
[C,hContour] = contour(pos_grid,z(1:5000),gamma(1:5000,:),g_n_levels,'color',[.1 .1 .1], 'showtext','on','linewidth',1,'labelspacing',350);


colormap(flipud(brewermap(11,'RdBu')))
hold on;
caxis([-.75 .75])

xlim([lon1w lon1e])

ylim([0 1300]);xticks;


xticks([-33  -32 -31 -30 -29]);set(gca,'xticklabels',[]); box on
%xticklabels({'33∫W','32∫W','31¬∫W','30¬∫W','29¬∫W'})
yticks([0 200 400 600 800 1000 1200]);set(gca,'yticklabels',[0:200:1200]); box on
ylabel('Pressure (dbar)')

% cb=colorbar;
%  cb.Position = cb.Position + 1e-10;
% ylabel(cb, 'Velocity (m.s^-^1)')

%yticklabels({})


% cmocean('delta')
% %contour(pos_grid,z,v,[0 0], 'linecolor', 'w','linewidth',1)
% [C,hContour] = contour(pos_grid,z(1:1500),gamma(1:1500,:),g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',400);
% % 
%  [C,hContour] = contour(pos_grid,z(1:1500),gamma(1:1500,:),[26.4:0.1:27],'k');
% 
% 


% p = get(h, 'pos');
% p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.1;p(2) = p(2) + 0.05;set(h, 'pos', p);box on; hold on

grid on

box on

axes(ha(8));hold on

text(0.92,0.125,'h','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

pos_gridd=pos_grid(:,1:end-1)+0.05/2;

pcolor(pos_gridd,z(1:1500),geostrophic_velocity_281_ek_balanced(1:1500,:));shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on

hold on


for l=1:length(stations)
xline(stations(l),'--','color',[.65 .65 .65]);hold on
end
%contour(pos_grid,z,v,[0 0], 'linecolor', 'w','linewidth',1)

% [C,hContour] = contour(pos_grid,z(1:1500),gamma(1:1500,:),g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',400);

[C,hContour] = contour(pos_grid,z(1:1500),gamma(1:1500,:),[22:0.15:26.2 26.5:0.15:26.95 27.25:0.15:28],'color',[.6 .6 .6],'linewidth',.4);
[C,hContour] = contour(pos_grid,z(1:5000),gamma(1:5000,:),g_n_levels,'color',[.1 .1 .1], 'showtext','on','linewidth',1,'labelspacing',350);


colormap(flipud(brewermap(11,'RdBu')))
hold on;
caxis([-.75 .75])

xlim([lon2w lon2e])

ylim([0 1300]);xticks;

xticks([-14 -13 -12 -11]);set(gca,'Xticklabel',[]) ; box on

% xticks([-33  -32 -31 -30 -29]);set(gca,'xticklabels',[]); box on
%xticklabels({'33∫W','32∫W','31¬∫W','30¬∫W','29¬∫W'})
yticks([0 200 400 600 800 1000 1200]);set(gca,'yticklabels',[]); box on









axes(ha(9));hold on 

text(0.92,0.125,'i','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

pos_gridd=pos_grid(:,1:end-1)+0.05/2;

pcolor(pos_gridd,z(1:1500),geostrophic_velocity_281_ek_balanced(1:1500,:));shading interp;axis ij%colormap(jet);colorbar;%title('temp');hold on

hold on

for l=1:length(stations)
xline(stations(l),'--','color',[.65 .65 .65]);hold on
end
%contour(pos_grid,z,v,[0 0], 'linecolor', 'w','linewidth',1)

% [C,hContour] = contour(pos_grid,z(1:1500),gamma(1:1500,:),g_n_levels,'k', 'showtext','on','linewidth',1,'labelspacing',400);
hold on

[C,hContour] = contour(pos_grid,z(1:1500),gamma(1:1500,:),[22:0.15:26.2 26.5:0.15:26.95 27.25:0.15:28],'color',[.6 .6 .6],'linewidth',.4);
[C,hContour] = contour(pos_grid,z(1:5000),gamma(1:5000,:),g_n_levels,'color',[.1 .1 .1], 'showtext','on','linewidth',1,'labelspacing',1400);

colormap(flipud(brewermap(11,'RdBu')))
hold on;
caxis([-.75 .75])

xlim([lon3w lon3e])

ylim([0 1300]);xticks;

yticks([0 200 400 600 800 1000 1200]);set(gca,'yticklabels',[]); box on

xticks([12 13 14 15]);xticklabels([]); box on

cb=colorbar;

pos=get(cb,'Position');set(cb,'Position',pos+[0.049,0,0,0]);
xlhgg = ylabel(cb, 'Vgeos. (m.s^-^1)');xlhgg.Position(1) = xlhgg.Position(1) - 0.2;





%% density part
%10
axes(ha(10)); hold on
text(0.92,0.125,'j','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

contourf(pos_grid,pres(1:2000),dens_anom,[-.6:.025:.6],'linestyle','none');shading interp;axis ij;hold on

xlim([lon1w lon1e])

ylim([0 1300]);xticks;


for l=1:length(stations)
xline(stations(l),'--','color',[.65 .65 .65]);hold on
end

[C,hContour] = contour(pos_grid,z(1:1500),gamma(1:1500,:),[22:0.15:26.2 26.5:0.15:26.95 27.25:0.15:28],'color',[.6 .6 .6],'linewidth',.4);
[C,hContour] = contour(pos_grid,z(1:5000),gamma(1:5000,:),g_n_levels,'color',[.1 .1 .1], 'showtext','on','linewidth',1,'labelspacing',350);
caxis([-.6 .6])


set(gcf,'color','w');  % otherwise 'print' turns lakes black
set(gca,'tickdir','out')

cmocean('curl',13);
% 
% cb=colorbar;
%  cb.Position = cb.Position + 1e-10;
ylabel('Pressure (dbar)')

% ylabel(cb,'œÅ_Œ∏ anomaly (kg.m^-^3)')

% xticks([-33  -32 -31 -30 -29]);xticklabels({'33¬∫W','32¬∫W','31¬∫W','30¬∫W','29¬∫W'})

xticks([-33  -32 -31 -30 -29]); 
xticklabels({'33∫W','32∫W','31∫W','30∫W','29∫W'});box on
yticks([0 200 400 600 800 1000 1200]);set(gca,'yticklabels',[0:200:1200]); box on



axes(ha(11)); 

 contourf(pos_grid-0.5,pres(1:2000),dens_anom,[-.65:.025:.65],'linestyle','none');shading interp;axis ij;hold on

% contour(lon,pres(1:2000),gamma(1:2000,:),25,'k');


ylim([0 1300]);xticks;

for l=1:length(stations)
xline(stations(l),'--','color',[.65 .65 .65]);hold on
end


[C,hContour] = contour(pos_grid,z(1:1500),gamma(1:1500,:),[22:0.15:26.2 26.5:0.15:26.95 27.25:0.15:28],'color',[.6 .6 .6],'linewidth',.4);
[C,hContour] = contour(pos_grid,z(1:5000),gamma(1:5000,:),g_n_levels,'color',[.1 .1 .1], 'showtext','on','linewidth',1,'labelspacing',350);

caxis([-.6 .6])

cmocean('curl',13);


text(0.92,0.125,'k','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')


xlim([lon2w lon2e]);

% xticks([-14 -13 -12 -11]);xticklabels({'14¬∫W','13¬∫W','12¬∫W','11¬∫W'})

yticks([0 200 400 600 800 1000 1200]);set(gca,'yticklabels',[]); box on

xticks([-14 -13 -12 -11]);xticklabels({'14∫W','13∫W','12∫W','11∫W'}); box on


axes(ha(12)); hold on

contourf(pos_grid,pres(1:2000),dens_anom,[-.6:.025:.6],'linestyle','none');shading interp;axis ij;hold on

xlim([lon3w lon3e])

ylim([0 1300]);xticks;


for l=1:length(stations)
xline(stations(l),'--','color',[.65 .65 .65]);hold on
end 


[C,hContour] = contour(pos_grid,z(1:1500),gamma(1:1500,:),[22:0.15:26.2 26.5:0.15:26.95 27.25:0.15:28],'color',[.6 .6 .6],'linewidth',.4);
[C,hContour] = contour(pos_grid,z(1:5000),gamma(1:5000,:),g_n_levels,'color',[.1 .1 .1], 'showtext','on','linewidth',1,'labelspacing',1400);


caxis([-.6 .6])


cmocean('curl',13);

cb=colorbar;

pos=get(cb,'Position');set(cb,'Position',pos+[0.049,0,0,0]);

xlh = ylabel(cb,'/sigma anomaly (kg.m^-^3)');xlh.Position(1) = xlh.Position(1) - 0.2;%xlh.Position(2) = xlh.Position(2) - 2000;clear xlh.Position(1);clear xlh.Position(2)

text(0.92,0.125,'l','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')


xticks([12 13 14 15]);xticklabels({'12∫E','13∫E','14∫E','15∫E'})


grid on
 box on




%% saving and tunning
set(gcf,'color','w');  % otherwise 'print' turns lakes black
set(gca,'tickdir','in')

set(findall(gcf,'-property','FontSize'),'FontSize',12)

set(findall(gcf,'-property','FontType'),'Fonttype','normal')
set(findall(gcf,'-property','Fontname'),'Fontname','Helvetica')

%set(findall(gcf,'-property','Linewidth'),'Linewidth',.7)

print(gcf, '-dpng','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/paper_merian/Figures_paper_merian/Figure_eddies_sep')

close all