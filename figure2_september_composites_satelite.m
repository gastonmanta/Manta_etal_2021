%% PLOT SST DATA, EDDIE CONTOURS, AND SADCP
cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian
tsg=importdata('msm_060_1_no_headers.tsg');

load stations
load merian_contours


%ncdisp('sst_merian_MSM602.nc')
sst=ncread('sst_merian_MSM602.nc','analysed_sst');sst=sst-273.18;sst=flip(sst,1);%1 january to 1 february included
lon_sst=ncread('sst_merian_MSM602.nc','lon');lon_sst=flip(lon_sst);
lat_sst=ncread('sst_merian_MSM602.nc','lat');
time_sst2=ncread('sst_merian_MSM602.nc','time');

% cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian
% load eddies_merian_contours_traj_and_argo


msst=nanmean(sst,3);msubsst=sst(:,:,end-1)-sst(:,:,1);


dates=tsg(:,1:6);daten=datenum(dates);
lon_tsg=tsg(:,8);lat_tsg=tsg(:,7);sst_tsg=nanmean(tsg(:,9:10),2);


figure;%hold on
m_proj('ortho','lat',-34.5','long',-15');
m_plot(lon_tsg,lat_tsg,'.r','Markersize',2);
m_coast('patch', [.7 .7 .7]);
% m_grid('xtick',[-34.5 -34.5])xticklabels',[],'yticklabels',[])
 m_grid('xticklabels',[],'xticks',[],'yticklabels',[],'yticks',[])


m_grid('linest','-','xticklabels',[],'yticklabels',[]);




cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/ssh_merian
aa=dir('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/ssh_merian')
ncfiles = dir('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/ssh_merian/*.nc');   % get all nc files in the folder 
nfiles = length(ncfiles) ;   % total number of files 
% ncdisp('dt_global_allsat_phy_l4_20170101_20170530.nc')

% loop for each file 

for K = 1 : nfiles
  filename = ncfiles(K).name;  %just the name
  %%get the vairable you want from nc file 
  adt(:,:,K)= ncread(filename,'adt') ; 
  u(:,:,K)= ncread(filename,'ugos') ; 
    v(:,:,K)= ncread(filename,'vgos') ; 
   lon_ssh= ncread(filename,'longitude') ;
   lat_ssh= ncread(filename,'latitude') ;% doc ncread 
  %%Append the data where you want 
end



cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian




lon_ssh(lon_ssh>180)=lon_ssh(lon_ssh>180)-360;
[lon_ssh,id_sort]=sort(lon_ssh);
adt=adt(id_sort,:,:);

u=u(id_sort,:,:);
v=v(id_sort,:,:);



adt=flip(adt,1);lon_ssh=flip(lon_ssh);
u=flip(u,1);v=flip(v,1);

% lon_ssh=lon_ssh(480:800);lat_ssh=lat_ssh(120:320);adt=adt(480:800,120:320,:);


%create a composite of 1km sst of sst days changing along cruise track
day_lon=dates(:,3);
ind=find(day_lon(2:end)-day_lon(1:end-1)>0);
lon_date=lon_tsg(ind);%lon day are the longitudes when a day change
lon_date=round(lon_date);%round it to 2 decimal places to match sst 1km longiudes
date_lon=day_lon(ind);

yw=find(ismember(lon_sst,lon_date));
yy=vertcat(1,yw);yy(end)=length(lon_sst);

ywh=find(ismember(lon_ssh,lon_date-.125));

yyh=vertcat(1,ywh);yyh(end)=length(lon_ssh);


% for i=1:27
%  sst_composite{i} = squeeze(sst(yy(i):yy(i+1)-1,:,date_lon(i)));
% end
%%%
for i=1:27
 sst_composite{i} = squeeze(sst(yy(i):yy(i+1)-1,:,date_lon(i)));
  ssh_composite{i} = squeeze(adt(yyh(i):yyh(i+1)-1,:,date_lon(i)));
    u_composite{i} = squeeze(u(yyh(i):yyh(i+1)-1,:,date_lon(i)));
        v_composite{i} = squeeze(v(yyh(i):yyh(i+1)-1,:,date_lon(i)));
end

 sst_composite= vertcat(sst_composite{:}); sst_composite=vertcat(sst_composite,sst_composite(end,:));
 ssh_composite= vertcat(ssh_composite{:}); ssh_composite=vertcat(ssh_composite,ssh_composite(end,:));
 u_composite= vertcat(u_composite{:}); u_composite=vertcat(u_composite,u_composite(end,:));
 v_composite= vertcat(v_composite{:}); v_composite=vertcat(v_composite,v_composite(end,:));

clear sst

b200=load('200.dat');b1000=load('1000.dat');


load /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/4gaston/vmadcp/msm60os75000_000000_39_hc.mat

v=squeeze(b.vel(:,2,:));
u=squeeze(b.vel(:,1,:));
zc=c.depth(:);
pg=c.pg;%8 meter interval from 17m to 800
kk=find(pg<=25);
v(kk)=NaN;
u(kk)=NaN;

%mean between 40 and 100m from 38 khz adcp
uplot=squeeze(nanmean(u(4:11,:),1));
vplot=squeeze(nanmean(v(4:11,:),1));
lon_adcp=b.nav.txy1(2,:);lat_adcp=b.nav.txy1(3,:);



spacing = 36/2;
smoothing=36;%1 is 10 minutes of data
uuplot=movmean(uplot,smoothing);vvplot=movmean(vplot,smoothing);
lon38=lon_adcp(1:spacing:end);
lat38=lat_adcp(1:spacing:end);

u38=uuplot(1:spacing:end);
v38=vvplot(1:spacing:end);



[lon_sshx,lat_sshy]=meshgrid(lon_ssh,lat_ssh);
lx=lon_sshx';ly=lat_sshy';

lat_d=[-37.5 -31.5];
lon_d=[lon_date lon_date];


ae_in=nansum(ae_in);ce_in=nansum(ce_in);
ae_inn=find(ae_in==1);ce_inn=find(ce_in==1);


load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian/eddies_merian.mat')


lon_ssh=double(lon_ssh);lat_ssh=double(lat_ssh);u_composite=double(u_composite);

ugos=griddata(lon_ssh,lat_ssh',u_composite',lon_adcp,lat_adcp);

vgos=griddata(lon_ssh,lat_ssh',v_composite',lon_adcp,lat_adcp);

figure;plot(lon_adcp,vgos);
hold on
plot(lon_adcp,vvplot);axis tight; ylim([-.7 .7])


uvgos=ugos(1:spacing:end);
vvgos=vgos(1:spacing:end);


% hold on
% plot
% ugos=squeeze(nmean(u_composite(:,223:224),2))
% 
%  ugos=griddata(lon_sshx',lat_sshy,u_composite,lon_adcp,lat_adcp);



%%%%%%%%%%FIGURES%%%%%%%%%%%%%%%%%%%%%%%%%%




% lonn38=[lon38 -59.9];latt38=[lat38 -31.7];uu38=[u38 0.5];vv38=[v38 0];
% ww=30;
% 

%adds a refernce arrow
%lonn38=[lon38 -60.2];latt38=[lat38 -31.3];uu38=[u38 0.5];vv38=[v38 0];
%

ww=30;



%yi; comes from the colocalization script, and it's the index in with the
%CEs have trajectories. The other's don't
cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian
load('eddies_merian_contours_traj_and_argo.mat')
 cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian
% load('eddies_merian_with_lost_traj.mat')
load('lost_traj.mat')


yi=[1;2;5;7;8;10;11;12;13;14;15;16;17;18;19;20];

date_num_CEs_out60_yi=date_num_CEs_out60(yi);

for i=1:length(yi)
ind(i)=find(date_num_CEs_out60_yi(i)==cy_traj_msm60{i,7});
dates_ces_traj60_crossed(i)=cy_traj_msm60{i,7}(ind(i));
end


yi=[1:16];
date_num_AEs_out60_yi=date_num_AEs_out60(yi);
dates_aes_traj60_crossed=date_num_AEs_out60;

% 
% for i=1:length(yi)
% ind(i)=find(date_num_AEs_out60_yi(i)==ay_traj_msm60{i,7});
% dates_ces_traj60_crossed(i)=cy_traj_msm60{i,7}(ind(i));
% end



%yi; comes from the colocalization script, and it's the index in with the
%CEs have trajectories. The other's don't
% yi=[1;2;5;7;8;10;11;12;13;14;15;16;17;18;19;20];
% date_num_CEs_out60_yi=date_num_CEs_out60(yi);
% date_num_CEs_out60_yi=date_num_CEs_out60(yi);
% 
% date_num_CEs_out60_yi=cy_traj_msm60{i,7}(yi);

% 
% for i=1:length(yi)
% ind(i)=find(date_num_CEs_out60_yi(i)==cy_traj_msm60{i,7});
% dates_ces_traj60_crossed(i)=cy_traj_msm60{i,7}(ind(i));
% end


% lonn38=[lon38 -59.9];latt38=[lat38 -31.7];uu38=[u38 0.5];vv38=[v38 0];
% ww=30;
% 

spacing = 36/2;
smoothing=36;%1 is 10 minutes of data
uuplot=movmean(uplot,smoothing);vvplot=movmean(vplot,smoothing);
lon38=lon_adcp(1:spacing:end);
lat38=lat_adcp(1:spacing:end);
u38=uuplot(1:spacing:end);
v38=vvplot(1:spacing:end);
% 
% lonn38=[lon38 -59.9];latt38=[lat38 -31.7];uu38=[u38 0.5];vv38=[v38 0];uvgos=[uvgos 0.5];vvgos=[vvgos 0];
% ww=30;

%adds a refernce arrow
%lonn38=[lon38 -60.2];latt38=[lat38 -31.3];uu38=[u38 0.5];vv38=[v38 0];
lonn38=[lon38 -58.9];latt38=[lat38 -30.7];uu38=[u38 0.5];vv38=[v38 0];uvgos=[uvgos 0.5];vvgos=[vvgos 0];


ww=30;


%dates_aes_traj60_crossed=ay_traj_msm60{i,7};


latmin=-39;latmax=-30;

%m_proj('mercator','long', [-60 20],'lat', [-39.5 -30.5]);hold on
%m_proj('mercator','long', [-60 20],'lat', [-44.5 -23.5]);hold on

%m_proj('mercator','long', [-60 20],'lat', [latmin latmax]);hold on
%figure('rend','painters','pos',[100 100 1150 600]); 

% figure; 
% 
% m_proj('mercator','long', [-60 20],'lat', [latmin latmax]);hold on
% 
% %m_proj('mercator','long', [-60 20],'lat', [latmin latmax]);hold on
% 
% 
%  [CS,CH]=m_etopo2('contourf',[-3500 -200],'edgecolor','none');caxis([-5500 7000]);colormap(flip(gray(7)));  
% 
% m_contour(lon_ssh,lat_ssh,ssh_composite',[-1:0.05:1],'color',[.55 .55 .55]);%colormap(jet(ww));
% 
% %  m_plot(lon_adcp,lat_adcp,'.k','Markersize',1);hold on
% 
%  %m_plot(lon_stations(CEs_in),lat_stations(CEs_in),'dg','Markersize',4);hold on
% 
% % for qq=1:length(lon_date)
% %     m_line(lon_d(qq,:),lat_d,'color',[.7 .7 .7],'LineWidth',1.5,'linestyle','--')
% % end
% ll=2;
% % h2=m_vec(1.3,lon38,lat38,u38,v38,'shaftwidth',0.5,'headlength',3,'headwidth',3);
% 
% for id_time= 1:length(CEs_out60)
%     m_plot(CEs_out60{id_time,id_time,1},CEs_out60{id_time,id_time,2},'c','LineWidth',1)
% end
% 
% for id_time= 1:length(AEs_out60)
%     m_plot(AEs_out60{id_time,id_time,1},AEs_out60{id_time,id_time,2},'m','LineWidth',1)
% end
% 
% 
% % m_grid('xtick',lon_date(1:17,'ytick',[-39 -37 -34.5 -32 -30],'linestyle','none','tickstyle','dd');
% 
% m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
% m_text(-59.9,-32.8,'0.5 m.sˆ-ˆ1','verticalalignment','bottom')
% 
% % h2=m_quiver(lonn38,latt38,uu38,vv38,1,'LineWidth',1,'Color','k','AutoScale','on','AutoScaleFactor',0.8)
% % h2=m_quiver(lonn38,latt38,uu38,vv38,1,'LineWidth',1,'Color','r','AutoScale','on','AutoScaleFactor',0.8)
% 
% 
% h3=m_quiver(lonn38,latt38,uvgos,vvgos,1,'LineWidth',1,'Color',[.5 .5 .5],'ShowArrowhead','off','AutoScale','on','AutoScaleFactor',0.8)
% 
% h2=m_quiver(lonn38,latt38,uu38,vv38,1,'LineWidth',1,'Color','k','ShowArrowhead','off','AutoScale','on','AutoScaleFactor',0.8)
% 
%  m_plot(lon_stations,lat_stations,'dk','Markersize',3,'MarkerFacecolor','k');hold on
% 
%  m_plot(lon_stations(ae_inn),lat_stations(ae_inn),'dm','Markersize',4,'MarkerFacecolor','m');hold on
%  m_plot(lon_stations(ce_inn),lat_stations(ce_inn),'dc','Markersize',4,'MarkerFacecolor','c');hold on
% 
% % h2 = quiver(x, y, dxda, dyda)
% set(h2,'AutoScale','on', 'AutoScaleFactor',1)
% 
% 
% m_grid('xtick',lon_date,'ytick',[-42 -38 -34.5 -30 -26],'tickstyle', 'dd' );
% 
% 
% 
% jj=120;
% 
%  for uhj=1:length(dates_ces_traj60_crossed)
% 
% ind=find(dates_ces_traj60_crossed(uhj)==cy_traj_msm60{uhj,7});
% if ind>jj
%   m_plot(cy_traj_msm60{uhj,5}(ind-jj:ind),cy_traj_msm60{uhj,6}(ind-jj:ind),'linewidth',1,'color',[.7 1 1]); hold on
% m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.','MarkerSize',8,'color',[.7 1 1]);
% else 
%  m_plot(cy_traj_msm60{uhj,5}(1:ind),cy_traj_msm60{uhj,6}(1:ind),'linewidth',1,'color',[.7 1 1]); hold on
% m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.','MarkerSize',8,'color',[.7 1 1]);
% 
% end
% end
% 
% %datevec(lost_traj{:,7}(1:7))
% 
% m_plot(lost_traj{:,5}(1:6),lost_traj{:,6}(1:6),'linewidth',1,'color',[.7 1 1]);
% 
% 
% jj=120;
% 
%  for uhj=1:length(dates_aes_traj60_crossed)
% 
% ind=find(dates_aes_traj60_crossed(uhj)==ay_traj_msm60{uhj,7});
% if ind>jj
%   m_plot(ay_traj_msm60{uhj,5}(ind-jj:ind),ay_traj_msm60{uhj,6}(ind-jj:ind),'color',[1 .7 1],'linewidth',1); hold on
% m_plot(ay_traj_msm60{uhj,5}(ind),ay_traj_msm60{uhj,6}(ind),'.','MarkerSize',8,'color',[1 .7 1]);
% else 
%  m_plot(ay_traj_msm60{uhj,5}(1:ind),ay_traj_msm60{uhj,6}(1:ind),'m','linewidth',1,'color',[1 .7 1]); hold on
% m_plot(ay_traj_msm60{uhj,5}(ind),ay_traj_msm60{uhj,6}(ind),'.','MarkerSize',8,'color',[1 .7 1]);
% 
% end
% end


clear msubsst msst kk dates 

%% fig 1 a adt

figure('Renderer', 'painters', 'Position', [10 10 1000 150])
m_proj('mercator','long', [-60 20],'lat', [latmin latmax]);hold on
text(0.008,0.25,'a','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

[CS,CH]=m_etopo2('contourf',[-3500 -200],'edgecolor','none');caxis([-5500 7000]);colormap(flip(gray(7)));  

% m_contour(lon_ssh,lat_ssh,ssh_composite',25,'color',[.55 .55 .55]);%colormap(jet(ww));
% m_grid('xtick',lon_date,'tickstyle', 'dd' );

m_grid('xtick',lon_date,'xticklabels',[],'linestyle','none');

[C,hContour] = m_contour(lon_ssh,lat_ssh,ssh_composite',[-.2:0.1:1],'color',[.55 .55 .55]);
% clabel(C,hContour,'manual','FontSize',12,'Color',[.55 .55 .55]);


% m_contour(lon_ssh,lat_ssh,ssh_composite',[-2.05:0.1:1.05],'color',[.55 .55 .55]);%


% m_contour(lon_ssh,lat_ssh,ssh_composite',[0.4 0.4],'color','r');%
% m_contour(lon_ssh,lat_ssh,ssh_composite',[0.5 0.5],'color','r');%
% m_contour(lon_ssh,lat_ssh,ssh_composite',[0.7 0.7],'color','b');%

%m_contour(lon_ssh,lat_ssh,ssh_composite',[0.4 .4],'color',[.5 .5 .5],'linewidth',1);%colormap(jet(ww));
% m_contour(lon_ssh,lat_ssh,ssh_composite',[0.5 .5],'color',[.5 .5 .5],'linewidth',1);%colormap(jet(ww));

m_plot(lon_adcp,lat_adcp,'.k','Markersize',1);hold on

 %m_plot(lon_stations(CEs_in),lat_stations(CEs_in),'dg','Markersize',4);hold on

% for qq=1:length(lon_date)
%     m_line(lon_d(qq,:),lat_d,'color',[.7 .7 .7],'LineWidth',1.5,'linestyle','--')
% end

ll=2;
% h2=m_vec(1.3,lon38,lat38,u38,v38,'shaftwidth',0.5,'headlength',3,'headwidth',3);


for id_time= 1:length(CEs_out60)
    m_plot(CEs_out60{id_time,id_time,1},CEs_out60{id_time,id_time,2},'c','LineWidth',1)

end

for id_time= 1:length(AEs_out60)
    m_plot(AEs_out60{id_time,id_time,1},AEs_out60{id_time,id_time,2},'m','LineWidth',1)
end


AEss={'A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12'}

CEss={'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13'}

CEs_lon=[13.5 8 -1 -5 -7.5 -11 -14.5 -17.5 -27.5 -31 -35 -47.5 -51.5];

AEs_lon=[15.6 13 3 -3.9 -10.2 -13 -19.5  -28.9 -33.9 -39.6 -43 -48.6]


for i=1:length(CEs_lon)
m_text(CEs_lon(i),-31.2,CEss{i},'color',[0 1 1],'fontsize',16,'fontweight','bold'); hold on
m_text(CEs_lon(i),-31.2,CEss{i},'color',[.4 .4 .4],'fontsize',16,'linewidth',.25)

end


for i=1:length(AEs_lon)
m_text(AEs_lon(i),-36.7,AEss{i},'color',[1 0 1],'fontsize',16,'fontweight','bold'); hold on
m_text(AEs_lon(i),-36.7,AEss{i},'color',[.4 .4 .4],'fontsize',16)

end

% m_grid('xtick',lon_date(1:17,'ytick',[-39 -37 -34.5 -32 -30],'linestyle','none','tickstyle','dd');

m_gshhs_h('patch',[0.81818 0.77647 0.70909]);

m_text(-58.65,-32.1,'0.5 m.s^-^1','verticalalignment','bottom')

% h2=m_quiver(lonn38,latt38,uu38,vv38,1,'LineWidth',1,'Color','k','AutoScale','on','AutoScaleFactor',0.8)
% h2=m_quiver(lonn38,latt38,uu38,vv38,1,'LineWidth',1,'Color','r','AutoScale','on','AutoScaleFactor',0.8)
%h3=m_quiver(lonn38,latt38,uvgos,vvgos,1,'LineWidth',1,'Color',[.5 .5 .5],'ShowArrowhead','off','AutoScale','on','AutoScaleFactor',0.8)

h2=m_quiver(lonn38,latt38,uu38,vv38,1,'LineWidth',1,'Color','k','ShowArrowhead','off','AutoScale','on','AutoScaleFactor',0.8)

 m_plot(lon_stations,lat_stations,'dk','Markersize',4,'MarkerFacecolor','k');hold on

 m_plot(lon_stations(ae_inn),lat_stations(ae_inn),'dm','Markersize',5,'MarkerFacecolor','m');hold on
 m_plot(lon_stations(ce_inn),lat_stations(ce_inn),'dc','Markersize',5,'MarkerFacecolor','c');hold on

% h2 = quiver(x, y, dxda, dyda)
set(h2,'AutoScale','on', 'AutoScaleFactor',1.2)



jj=120;

 for uhj=1:length(dates_ces_traj60_crossed)

ind=find(dates_ces_traj60_crossed(uhj)==cy_traj_msm60{uhj,7});
if ind>jj
  m_plot(cy_traj_msm60{uhj,5}(ind-jj:ind),cy_traj_msm60{uhj,6}(ind-jj:ind),'linewidth',1,'color',[.55 1 1]); hold on
m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.','MarkerSize',8,'color',[.55 1 1]);
else 
 m_plot(cy_traj_msm60{uhj,5}(1:ind),cy_traj_msm60{uhj,6}(1:ind),'linewidth',1,'color',[.55 1 1]); hold on
m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.','MarkerSize',8,'color',[.55 1 1]);

end
end
m_plot(lost_traj{:,5}(1:6),lost_traj{:,6}(1:6),'linewidth',1,'color',[.55 1 1]);

%datevec(lost_traj{:,7}(1:7))


jj=120;

 for uhj=1:length(dates_aes_traj60_crossed)

ind=find(dates_aes_traj60_crossed(uhj)==ay_traj_msm60{uhj,7});
if ind>jj
  m_plot(ay_traj_msm60{uhj,5}(ind-jj:ind),ay_traj_msm60{uhj,6}(ind-jj:ind),'color',[1 .55 1],'linewidth',1); hold on
m_plot(ay_traj_msm60{uhj,5}(ind),ay_traj_msm60{uhj,6}(ind),'.','MarkerSize',8,'color',[1 .55 1]);
else 
 m_plot(ay_traj_msm60{uhj,5}(1:ind),ay_traj_msm60{uhj,6}(1:ind),'m','linewidth',1,'color',[1 .55 1]); hold on
m_plot(ay_traj_msm60{uhj,5}(ind),ay_traj_msm60{uhj,6}(ind),'.','MarkerSize',8,'color',[1 .55 1]);

end
end


 m_plot(lost_traj{:,5}(1:6),lost_traj{:,6}(1:6),'linewidth',1,'color',[.7 1 1]);


set(findall(gcf,'-property','FontSize'),'FontSize',12)

%savefig('Figure_1a_merian')
% 
print(gcf, '-dpng','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september/figure2_set_composites_adt')%  
print(gcf, '-djpeg','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september/figure2_set_composite_adt')%  
% 
% savefig('fig1a_sadcp')
% 
% print(gcf, '-djpeg','-r600','/Users/gaston/Documents/eureca_sadcp_adt')
% jj=120;
% 
%  for uhj=1:length(dates_ces_traj60_crossed)
% 
% ind=find(dates_ces_traj60_crossed(uhj)==cy_traj_msm60{uhj,7});
% if ind>jj
%   m_plot(cy_traj_msm60{uhj,5}(ind-jj:ind),cy_traj_msm60{uhj,6}(ind-jj:ind),'linewidth',1,'color',[.7 1 1]); hold on
% m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.','MarkerSize',8,'color',[.7 1 1]);
% else 
%  m_plot(cy_traj_msm60{uhj,5}(1:ind),cy_traj_msm60{uhj,6}(1:ind),'linewidth',1,'color',[.7 1 1]); hold on
% m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.','MarkerSize',8,'color',[.7 1 1]);
% 
% end
% end
% 
% %datevec(lost_traj{:,7}(1:7))
% 
% 
% 
% jj=120;
% 
%  for uhj=1:length(dates_aes_traj60_crossed)
% 
% ind=find(dates_aes_traj60_crossed(uhj)==ay_traj_msm60{uhj,7});
% if ind>jj
%   m_plot(ay_traj_msm60{uhj,5}(ind-jj:ind),ay_traj_msm60{uhj,6}(ind-jj:ind),'color',[1 .7 1],'linewidth',1); hold on
% m_plot(ay_traj_msm60{uhj,5}(ind),ay_traj_msm60{uhj,6}(ind),'.','MarkerSize',8,'color',[1 .7 1]);
% else 
%  m_plot(ay_traj_msm60{uhj,5}(1:ind),ay_traj_msm60{uhj,6}(1:ind),'m','linewidth',1,'color',[1 .7 1]); hold on
% m_plot(ay_traj_msm60{uhj,5}(ind),ay_traj_msm60{uhj,6}(ind),'.','MarkerSize',8,'color',[1 .7 1]);
% 
% end
% end
% 
% 

%% SST fig

%%%SST%%%%%%%%
b200=load('200.dat');b1000=load('1000.dat');

gy=35;

figure('Renderer', 'painters', 'Position', [10 10 1000 150])
m_proj('mercator','long', [-60 20],'lat', [latmin latmax]);hold on
text(0.008,0.25,'b','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')


m_proj('mercator','long', [-60 20],'lat', [latmin latmax]);hold on
m_pcolor(lon_sst,lat_sst,sst_composite');shading interp;colormap(flipud(brewermap(49,'Spectral')));%caxis([-0.5 23.5])

% m_contourf(lon_sst,lat_sst,sst_composite',20,'linestyle','none');colormap(jet(20));hold on

m_line(b200(:,1),b200(:,2),'color',[0.81818 0.77647 0.70909],'LineWidth',1.5); hold on    % Area outline
%m_line(b1000(:,1),b1000(:,2),'color',[.7 .7 .7],'LineWidth',1.5); hold on    % Area outline


% m_gshhs_h('patch', [.7 .7 .7]);%m_grid('xtick',lon_date,'ytick',[-37.5 -34.5 -31.5],'tickstyle', 'dd' );
%title('sst 31/jan/17 - sst 1/jan/17');
tsg_space=20;
m_plot(lon_tsg(250:tsg_space:end),lat_tsg(250:tsg_space:end),'.k','Markersize',gy-5);hold on
m_scatter(lon_tsg(250:tsg_space:end),lat_tsg(250:tsg_space:end),gy,sst_tsg(250:tsg_space:end),'filled','MarkerEdgeColor','none');%colormap(jet(20));

caxis([15.5 25.5]);hold on; %colorbar

% cc = colorbar;
% cc.Label
% cc.Label.String = 'SST (°C)';
% figure

% m_proj('mercator','long', [-60 20],'lat', [latmin latmax]);hold on

%m_proj('mercator','long', [-60 20],'lat', [latmin latmax]);hold on


% [CS,CH]=m_etopo2('contourf',[-3500 -200],'edgecolor','none');caxis([-5500 7000]);colormap(flip(gray(7)));  


hold on
% m_contour(lon_ssh,lat_ssh,ssh_composite',[-.2:0.05:1],'color',[.55 .55 .55]);%colormap(jet(ww));


% [C,hContour] = m_contour(lon_ssh,lat_ssh,ssh_composite',[-.2:0.1:1],'color',[.55 .55 .55]);
% clabel(C,hContour,'manual','FontSize',12,'Color',[.55 .55 .55]);

m_contour(lon_ssh,lat_ssh,ssh_composite',[-2.05:0.1:1.05],'color',[.55 .55 .55]);%

hold on
m_grid('xtick',lon_date(1:2:end),'tickstyle', 'dd' );

m_gshhs_h('patch',[0.81818 0.77647 0.70909]);

% m_contour(lon_ssh,lat_ssh,ssh_composite',[0.4 0.4],'color','r');%
% 

[C,hContour] = m_contour(lon_ssh,lat_ssh,ssh_composite',[0.5 0.5],'color',[.3 .3 .3],'linewidth',1);

clabel(C,hContour,'manual','FontSize',12,'Color','k');




% 
% m_contour(lon_ssh,lat_ssh,ssh_composite',[0.7 0.7],'color','b');%





% %m_contour(lon_ssh,lat_ssh,ssh_composite',[0.4 .4],'color',[.5 .5 .5],'linewidth',1);%colormap(jet(ww));
% m_contour(lon_ssh,lat_ssh,ssh_composite',[0.5 .5],'color',[.5 .5 .5],'linewidth',1);%colormap(jet(ww));
% 
% m_contour(lon_ssh,lat_ssh,ssh_composite',[0.7 .7],'color',[.5 .5 .5],'linewidth',3);%colormap(jet(ww));
% 
% m_contour(lon_ssh,lat_ssh,ssh_composite',[0.55 .55],'color',[.5 .5 .5],'linewidth',1);%colormap(jet(ww));

% 
% 
% m_contour(lon_ssh,lat_ssh,ssh_composite',[0.5 .5],'color','r','linewidth',2);%colormap(jet(ww));
% 
% m_contour(lon_ssh,lat_ssh,ssh_composite',[0.6 .6],'color','k','linewidth',2);%colormap(jet(ww));
% 
% m_contour(lon_ssh,lat_ssh,ssh_composite',[0.6 .6],'color',[.5 .5 .5],'linewidth',);%colormap(jet(ww));

%m_contour(lon_ssh,lat_ssh,ssh_composite',[0.6 .6],'color',[.45 .45 .45],'linewidth',1);%colormap(jet(ww));

%  m_plot(lon_adcp,lat_adcp,'.k','Markersize',1);hold on

 %m_plot(lon_stations(CEs_in),lat_stations(CEs_in),'dg','Markersize',4);hold on

% for qq=1:length(lon_date)
%     m_line(lon_d(qq,:),lat_d,'color',[.7 .7 .7],'LineWidth',1.5,'linestyle','--')
% end
ll=2;
% h2=m_vec(1.3,lon38,lat38,u38,v38,'shaftwidth',0.5,'headlength',3,'headwidth',3);

for id_time= 1:length(CEs_out60)
    m_plot(CEs_out60{id_time,id_time,1},CEs_out60{id_time,id_time,2},'c','LineWidth',1)
end

for id_time= 1:length(AEs_out60)
    m_plot(AEs_out60{id_time,id_time,1},AEs_out60{id_time,id_time,2},'m','LineWidth',1)
end


% m_grid('xtick',lon_date(1:17,'ytick',[-39 -37 -34.5 -32 -30],'linestyle','none','tickstyle','dd');

% m_text(-59.9,-32.8,'0.5 m.sˆ-ˆ1','verticalalignment','bottom')

% h2=m_quiver(lonn38,latt38,uu38,vv38,1,'LineWidth',1,'Color','k','AutoScale','on','AutoScaleFactor',0.8)
% h2=m_quiver(lonn38,latt38,uu38,vv38,1,'LineWidth',1,'Color','r','AutoScale','on','AutoScaleFactor',0.8)

% 
% h3=m_quiver(lonn38,latt38,uvgos,vvgos,1,'LineWidth',1,'Color',[.2 .2 .2],'ShowArrowhead','off','AutoScale','on','AutoScaleFactor',0.8)
% 
% % h2=m_quiver(lonn38,latt38,uu38,vv38,1,'LineWidth',1,'Color','k','ShowArrowhead','off','AutoScale','on','AutoScaleFactor',0.8)
% 
%  m_plot(lon_stations,lat_stations,'dk','Markersize',3,'MarkerFacecolor','k');hold on
% 
%  m_plot(lon_stations(ae_inn),lat_stations(ae_inn),'dm','Markersize',4,'MarkerFacecolor','m');hold on
%  m_plot(lon_stations(ce_inn),lat_stations(ce_inn),'dc','Markersize',4,'MarkerFacecolor','c');hold on
% 
% % h2 = quiver(x, y, dxda, dyda)
% set(h2,'AutoScale','on', 'AutoScaleFactor',1)



% m_grid('xtick',lon_date,'tickstyle', 'none' );

% jj=12000;
% 
%  for uhj=1:length(dates_ces_traj60_crossed)
% 
% ind=find(dates_ces_traj60_crossed(uhj)==cy_traj_msm60{uhj,7});
% if ind>jj
%   m_plot(cy_traj_msm60{uhj,5}(ind-jj:ind),cy_traj_msm60{uhj,6}(ind-jj:ind),'linewidth',1,'color',[.7 1 1]); hold on
% m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.','MarkerSize',8,'color',[.7 1 1]);
% else 
%  m_plot(cy_traj_msm60{uhj,5}(1:ind),cy_traj_msm60{uhj,6}(1:ind),'linewidth',1,'color',[.7 1 1]); hold on
% m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.','MarkerSize',8,'color',[.7 1 1]);
% 
% end
% end

%datevec(lost_traj{:,7}(1:7))

m_plot(lost_traj{:,5}(1:6),lost_traj{:,6}(1:6),'linewidth',1,'color','k');


% jj=120;
% 
%  for uhj=1:length(dates_aes_traj60_crossed)
% 
% ind=find(dates_aes_traj60_crossed(uhj)==ay_traj_msm60{uhj,7});
% if ind>jj
%   m_plot(ay_traj_msm60{uhj,5}(ind-jj:ind),ay_traj_msm60{uhj,6}(ind-jj:ind),'color',[1 .7 1],'linewidth',1); hold on
% m_plot(ay_traj_msm60{uhj,5}(ind),ay_traj_msm60{uhj,6}(ind),'.','MarkerSize',8,'color',[1 .7 1]);
% else 
%  m_plot(ay_traj_msm60{uhj,5}(1:ind),ay_traj_msm60{uhj,6}(1:ind),'m','linewidth',1,'color','m'); hold on
% m_plot(ay_traj_msm60{uhj,5}(ind),ay_traj_msm60{uhj,6}(ind),'.','MarkerSize',8,'color','m');
% 
% end
% end


jj=120;

 for uhj=1:length(dates_ces_traj60_crossed)

ind=find(dates_ces_traj60_crossed(uhj)==cy_traj_msm60{uhj,7});
if ind>jj
  m_plot(cy_traj_msm60{uhj,5}(ind-jj:ind),cy_traj_msm60{uhj,6}(ind-jj:ind),'linewidth',1,'color','k'); hold on
m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.','MarkerSize',8,'color','k');
else 
 m_plot(cy_traj_msm60{uhj,5}(1:ind),cy_traj_msm60{uhj,6}(1:ind),'linewidth',1,'color','k'); hold on
m_plot(cy_traj_msm60{uhj,5}(ind),cy_traj_msm60{uhj,6}(ind),'.','MarkerSize',8,'color','k');

end
end

%datevec(lost_traj{:,7}(1:7))

m_plot(lost_traj{:,5}(1:6),lost_traj{:,6}(1:6),'linewidth',1,'color','k');


jj=120;

anit_color=[.2 .2 .2]
 for uhj=1:length(dates_aes_traj60_crossed)

ind=find(dates_aes_traj60_crossed(uhj)==ay_traj_msm60{uhj,7});
if ind>jj
  m_plot(ay_traj_msm60{uhj,5}(ind-jj:ind),ay_traj_msm60{uhj,6}(ind-jj:ind),'color',anit_color,'linewidth',1); hold on
m_plot(ay_traj_msm60{uhj,5}(ind),ay_traj_msm60{uhj,6}(ind),'.','MarkerSize',8,'color',anit_color);
else 
 m_plot(ay_traj_msm60{uhj,5}(1:ind),ay_traj_msm60{uhj,6}(1:ind),'m','linewidth',1,'color',anit_color); hold on
m_plot(ay_traj_msm60{uhj,5}(ind),ay_traj_msm60{uhj,6}(ind),'.','MarkerSize',8,'color',anit_color);

end
end





[ax,h]=m_contfbar([.0078 .085],.598,[15:1:26],[15:1:26]);
title(ax,'SST (�C)','Fontweight','normal')

set(findall(gcf,'-property','FontSize'),'FontSize',12)


print(gcf, '-dpng','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september/figure2_set_composite_sst')%  
print(gcf, '-djpeg','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september/figure2_set_composite_sst')%  

% savefig('Figure_1b_merian')
% 
%  print(gcf, '-djpeg','-r600','/Users/gaston/Documents/Figure_1b_merian')

