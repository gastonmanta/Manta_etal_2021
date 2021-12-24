clear all; close all;clc

cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian
load MSM60_ctd_november.mat

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

% geo_strf_dyn_height = gsw_geo_strf_dyn_height(SA,CT,pres,0);
% [geostrophic_velocity0, mid_lat, mid_long] =gsw_geostrophic_velocity(geo_strf_dyn_height,lon,lat,0);

% Neutral Density (Î³, kg/m3) Range Assign for Each Water Mass Layers
% Water mass Î³ range Î¸ (Â°C) S O2 (Î¼mol/kg)
% TW <26.35 19.75
% 36.09 225 35.10 211 34.245 237 34.58 188 34.91 241 34.72 215 34.67 226
% SACW 26.35â€“27.10
% AAIW 27.10â€“27.60
% UCDW 27.60â€“27.90
% NADW 27.90â€“28.10
% LCDW 28.10â€“28.27
% AABW >28.27 0.10
% Note. The zonally averaged potential temperature (Î¸), salinity (S), and dis- solved oxygen (O2) values at the center of each range computed from all cruises occupied at the SAMBA-W line are given as a reference.
% aTW: Tropical Water; SACW: South Atlantic Central Water; AAIW: Antarctic Intermediate Water; UCDW: Upper Circumpolar Deep Water; NADW: North Atlantic Deep Water; LCDW: Lower Circumpolar Deep Water; AABW: Antarctic Bottom Water.

load('argo_merian_gridded.mat')

ptemp_comp=cat(1,ptemp_ar(1:1200,:),CT(1201:end,:));
sal_comp=cat(1,sal_ar(1:1200,:),SA(1201:end,:));
% % SA=sal_comp;CT=ptemp_comp;

spd=sqrt(u.*u+v.*v);

spd_z=nanmean(spd,2);[val_min,pres_min]=min(spd_z); %2204 minimun spd 0.0630

bucket(end,:)=NaN;buckett(:,1040)=buckett(:,1039);

%4579 m2 is the area of each pixel (1m.0.05degrees), bucket land mask and
%1e6 to get it into sverdrups, so it's not just the area what it means
length_grid_point=sw_dist([-34.5 -34.5],[-30 -30.05],'km')*1000


area=length_grid_point*buckett/1e6;
 
%BUCKET IS THE LAND MASK
temp_new=temp.*bucket;sal_new=sal.*bucket;ox_new=ox.*bucket;u_new=u.*bucket;v_new=v.*bucket;
gamma_new=gamma.*bucket;ro_new=ro.*bucket;theta_new=theta.*bucket;spd_new=spd.*bucket;

water_masses=({'TW','SACW','AAIW','UCDW','NADW','LCDW','AABW'})

% g_n_levels_sabrina=[19 26.35 27.1 27.6 27.9 28.12 28.22 35];

%following valla et al 2018
g_n_levels=[19 26.35 27.1 27.6 27.9 28.1 28.27 35];
g_n_levels=[26.35 27.1 27.6 27.9 28.1 28.27];

%hernandez guerra
 g_n_levels_hg=[21 26.14 26.45 27.0 27.23 27.58 27.84 28.04 28.1 28.15 28.23 29];

xx=1:123;tt=31; 
pres=1:5501;%pres=repmat(pres,1401,1);pres=pres';


%GEOSTROPHIC VELOCITY WITH RESPECT TO 3400db AND 0db
geo_strf_dyn_height_3400 = gsw_geo_strf_dyn_height(SA,CT,pres',3400);
geostrophic_velocity_3400=gsw_geostrophic_velocity(geo_strf_dyn_height_3400,lon,lat,0);

geo_strf_dyn_height_0 = gsw_geo_strf_dyn_height(SA,CT,pres',0);
geo_strf_dyn_height_0 = gsw_geo_strf_dyn_height(sal_new,theta_new,pres',0);

[geostrophic_velocity_0, mid_lat, mid_long]=gsw_geostrophic_velocity(geo_strf_dyn_height_0,lon,lat,0);


%MEAN AND MEDIAN DEPTH OF GAMMA 28.1
for i=1:length(lon)
ind=find(gamma(:,i)>28.098 & gamma(:,i)<28.102);
pres281(i)=nmedian(pres(ind));
ppres281(i)=nmean(pres(ind));
v281(i)=nmedian(v_new(ind,i));
vv281(i)=nmean(v_new(ind,i));
end

%CALCLUATES GEOSTROPHIC VEL AT GAMMA 28.1
for i=1:length(lon)-1
ind=find(gamma(:,i)>28.098 & gamma(:,i)<28.102);
ladcp281(i)=nmean(v_new(ind,i));
geostrophic_velocity_281(i)=nmean(geostrophic_velocity0(ind,i));
end



%GEOSTROPHIC Velocity with respect to gamma 28.1

geostrophic_velocity281r=repmat(geostrophic_velocity_281,5501,1);
geostrophic_velocity_281=geostrophic_velocity0-geostrophic_velocity281r; clear geostrophic_velocity281r;


%REMOVING A COUPLE OF SPIKES

geostrophic_velocity_281(:,441)=((geostrophic_velocity_281(:,440)+geostrophic_velocity_281(:,442))./2);
geostrophic_velocity_281(:,1392)=((geostrophic_velocity_281(:,1391)+geostrophic_velocity_281(:,1393))./2);
% 

% geostrophic_velocity_281(:,1380:end)=movmean(geostrophic_velocity_281(:,1380:end),5,1);
% geostrophic_velocity_281(:,1:20)=movmean(geostrophic_velocity_281(:,1:20),5,1);
% geostrophic_velocity_2811=geostrophic_velocity_281.*buckett;


%REMOVING A COUPLE OF SPIKES
geostrophic_velocity_3400(:,441)=((geostrophic_velocity_3400(:,440)+geostrophic_velocity_3400(:,442))./2);

%geostrophic_velocity_3400(:,1380:end)=movmean(geostrophic_velocity_3400(:,1380:end),5,1);
geostrophic_velocity_3400(:,1392)=((geostrophic_velocity_3400(:,1391)+geostrophic_velocity_3400(:,1393))./2);


geostrophic_velocity_3400(:,1:4)=0;geostrophic_velocity_3400(:,end-4:end)=0;
geostrophic_velocity_3400(geostrophic_velocity_3400>1.5)=1.5;
geostrophic_velocity_3400(geostrophic_velocity_3400<-1.5)=-1.5;



geostrophic_velocity_281(:,1:4)=0;geostrophic_velocity_281(:,end-4:end)=0;
geostrophic_velocity_281(geostrophic_velocity_281>1.5)=1.5;
geostrophic_velocity_281(geostrophic_velocity_281<-1.5)=-1.5;

% geostrophic_velocity_3400(:,1:20)=movmean(geostrophic_velocity_3400(:,1:20),5,1);
% geostrophic_velocity_3400=geostrophic_velocity_3400.*buckett;


% geostrophic_velocity_3400(:,441)=((geostrophic_velocity_281(:,440)+geostrophic_velocity_281(:,442))./2);
% 
% geostrophic_velocity_281(:,1380:end)=movmean(geostrophic_velocity_281(:,1380:end),5,1);
% geostrophic_velocity_281(:,1:20)=movmean(geostrophic_velocity_281(:,1:20),5,1);
% geostrophic_velocity_2811=geostrophic_velocity_281.*buckett;

%TRANSPORT, BASICALLY MULTYPLYING PER area and taking to Sverdrups
area=4579*buckett/1e6;

%the volume contribution of each pixel
0.004579;
%the amount of pixels to repart the transport
1400*50;
%the ekman transport
-0.42;

% 
%%

%ADDS THE EKMAN COMPONENT, AND ALSO CREATES A BALANCED BY ADDING BAROTROPIC
%EQUAL VELOCITY
geostrophic_velocity_281_ek=geostrophic_velocity_281;
geostrophic_velocity_281_ek(1:100,:)=geostrophic_velocity_281_ek(1:100,:)-0.0011/2;
geostrophic_velocity_281_ek=geostrophic_velocity_281_ek.*buckett;
geostrophic_velocity_281_ek_balanced=geostrophic_velocity_281_ek-nsum(nsum(geostrophic_velocity_281_ek))./nsum(nsum((~isnan(geostrophic_velocity_281_ek))));


geostrophic_velocity_3400_ek=geostrophic_velocity_3400;
geostrophic_velocity_3400_ek(1:100,:)=geostrophic_velocity_3400_ek(1:100,:)-0.0011/2;
geostrophic_velocity_3400_ek=geostrophic_velocity_3400_ek.*buckett;
geostrophic_velocity_3400_ek_balanced=geostrophic_velocity_3400_ek-nsum(nsum(geostrophic_velocity_3400_ek))./nsum(nsum((~isnan(geostrophic_velocity_3400_ek))));



v_new_ek=v_new;
v_new_ek(1:100,:)=v_new_ek(1:100,:)-0.0041/2;
v_new_ek=v_new_ek.*bucket;
v_new_ek_balanced=v_new_ek-nsum(nsum(v_new_ek))./nsum(nsum((~isnan(v_new_ek))));


transportvgos28_balanced=geostrophic_velocity_281_ek_balanced.*area;%transport per area
transportvgos281_balanced=nansum(transportvgos28_balanced,2);transportvgos28_balanced=movmean(transportvgos28_balanced,10);%averaged in z


transportvgos28_ek=geostrophic_velocity_281.*area;%transport per area
transportvgos281_ek=nansum(transportvgos28_ek,2);transportvgos28_ek=movmean(transportvgos28_ek,10);%averaged in z


% % 

transportvgos3400a_balanced=geostrophic_velocity_3400_ek_balanced.*area;%transport per area
transportvgos3400_balanced=nansum(transportvgos3400a_balanced,2);transportvgos3400_balanced=movmean(transportvgos3400_balanced,10);%averaged in z

transportvgos3400a=geostrophic_velocity_3400.*area;%transport per area
transportvgos3400=nansum(transportvgos3400a,2);transportvgos3400=movmean(transportvgos3400,10);%averaged in z


transportva=v_new_ek(:,2:end).*area;%transport per area
transportv=nansum(transportva,2);transportvv=movmean(transportv,100);%averaged in z


transportva_balanced=v_new_ek_balanced(:,2:end).*area;%transport per area
transportv_balanced=nansum(transportva_balanced,2);transportvv_balanced=movmean(transportv_balanced,100);%averaged in z




%THE SAME BUT REPLACING THE UPPER 1500m with argo clim data
load('argo_merian_gridded.mat')

ct_ar=cat(1,ptemp_ar(1:1500,:),CT(1501:end,:));
sal_ar=cat(1,sal_ar(1:1500,:),SA(1501:end,:));
ct_argo=cat(2,CT(:,1:20),ct_ar(:,21:1401-20),CT(:,end-19:end));
sal_argo=cat(2,SA(:,1:20),sal_ar(:,21:1401-20),SA(:,end-19:end));



%GEOSTROPHIC VELOCITY WITH RESPECT TO 3400db
geo_strf_dyn_height_3400argo = gsw_geo_strf_dyn_height(sal_argo,ct_argo,pres',3400);
geostrophic_velocity_3400argo=gsw_geostrophic_velocity(geo_strf_dyn_height_3400argo,lon,lat,0);

geostrophic_velocity_3400argo_ek=geostrophic_velocity_3400argo;
geostrophic_velocity_3400argo_ek(1:100,:)=geostrophic_velocity_3400argo_ek(1:100,:)-0.0011/2;
geostrophic_velocity_3400argo_ek=geostrophic_velocity_3400argo_ek.*buckett;
geostrophic_velocity_3400argo_ek_balanced=geostrophic_velocity_3400argo_ek-nsum(nsum(geostrophic_velocity_3400argo_ek))./nsum(nsum((~isnan(geostrophic_velocity_3400argo_ek))));

transportvgos3400argo=geostrophic_velocity_3400argo_ek.*area;
transportvgos3400argo_balanced=geostrophic_velocity_3400argo_ek_balanced.*area;





% %following valla et al 2018
g_n_levels=[19 26.35 27.1 27.6 27.9 28.1 28.27 38 38];



for i=2:length(g_n_levels)
ind=find(gammaa(:,1:end-1)<g_n_levels(i) & gammaa(:,1:end-1)>g_n_levels(i-1));

transport_by_layer_hg(i)=nansum(transportva(ind));
transport_by_layer_hg_geos3400(i)=nansum(transportvgos3400a(ind));
transport_by_layer_hg_geos3400argo(i)=nansum(transportvgos3400argo(ind));
transport_by_layer_hg_geos281(i)=nansum(transportvgos28_ek(ind));


transport_by_layer_hgb(i)=nansum(transportva_balanced(ind));
transport_by_layer_hg_geos3400b(i)=nansum(transportvgos3400a_balanced(ind));
transport_by_layer_hg_geos3400argob(i)=nansum(transportvgos3400argo_balanced(ind));
transport_by_layer_hg_geos281b(i)=nansum(transportvgos28_balanced(ind));

end



%just the cumulative sum of the transports
tt=nancumsum(transport_by_layer_hg,'omitnan');
tts=nancumsum(transport_by_layer_hg_geos281);
tthg=nancumsum(transport_by_layer_hg_geos3400);
tthga=nancumsum(transport_by_layer_hg_geos3400argo);


%b stands for balanced at the end
ttb=nancumsum(transport_by_layer_hgb);
ttsb=nancumsum(transport_by_layer_hg_geos281b);
tthgb=nancumsum(transport_by_layer_hg_geos3400b);
tthgab=nancumsum(transport_by_layer_hg_geos3400argob);
%%

%% AND NOW THE PLOTS BEGINS
%PLOT SST DATA, EDDIE CONTOURS, AND SADCP
cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian
tsg=importdata('msm_060_1_no_headers.tsg');

load stations
load merian_contours


%ncdisp('sst_merian_MSM602.nc')
sst=ncread('sst_merian_MSM602.nc','analysed_sst');sst=sst-273.18;sst=flip(sst,1);%1 january to 1 february included
lon_sst=ncread('sst_merian_MSM602.nc','lon');lon_sst=flip(lon_sst);
lat_sst=ncread('sst_merian_MSM602.nc','lat');
time_sst2=ncread('sst_merian_MSM602.nc','time');
msst=nanmean(sst,3);msubsst=sst(:,:,end-1)-sst(:,:,1);
dates=tsg(:,1:6);daten=datenum(dates);
lon_tsg=tsg(:,8);lat_tsg=tsg(:,7);sst_tsg=nanmean(tsg(:,9:10),2);

%WORLD MAP WITH CRUISE TRACK, NEVER KNOW WHEN NEED IT
% figure;%hold on
% m_proj('ortho','lat',-34.5','long',-15');
% m_plot(lon_tsg,lat_tsg,'.r','Markersize',2);
% m_coast('patch', [.7 .7 .7]);
% % m_grid('xtick',[-34.5 -34.5])xticklabels',[],'yticklabels',[])
%  m_grid('xticklabels',[],'xticks',[],'yticklabels',[],'yticks',[])
% m_grid('linest','-','xticklabels',[],'yticklabels',[]);




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
  ugos(:,:,K)= ncread(filename,'ugos') ; 
    vgos(:,:,K)= ncread(filename,'vgos') ; 
   lon_ssh= ncread(filename,'longitude') ;
   lat_ssh= ncread(filename,'latitude') ;% doc ncread 
  %%Append the data where you want 
end



cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian




lon_ssh(lon_ssh>180)=lon_ssh(lon_ssh>180)-360;
[lon_ssh,id_sort]=sort(lon_ssh);
adt=adt(id_sort,:,:);

ugos=ugos(id_sort,:,:);
vgos=vgos(id_sort,:,:);
adt=flip(adt,1);lon_ssh=flip(lon_ssh);
ugos=flip(ugos,1);vgos=flip(vgos,1);

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
 %sst_composite{i} = squeeze(sst(yy(i):yy(i+1)-1,:,date_lon(i)));
  ssh_composite{i} = squeeze(adt(yyh(i):yyh(i+1)-1,:,date_lon(i)));
    u_composite{i} = squeeze(ugos(yyh(i):yyh(i+1)-1,:,date_lon(i)));
        v_composite{i} = squeeze(vgos(yyh(i):yyh(i+1)-1,:,date_lon(i)));
end

% sst_composite= vertcat(sst_composite{:}); sst_composite=vertcat(sst_composite,sst_composite(end,:));
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


% uvgos=ugos(1:spacing:end);
% vvgos=vgos(1:spacing:end);

v_new100=nmean(v_new(1:100,:),1);



% save ('figure3a_surface_currents','lon_adcp','lon','vvplot','vgos')

v_new100(442:452)=movmean(v_new100(442:452),10);

vvplot(442:452)=movmean(v_new100(442:452),10);

cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian
load MSM60_ctd_october.mat

load vu;v=v/100;u=u/100;

load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian/stations.mat')
 
ae_in=nansum(ae_in);ce_in=nansum(ce_in);
ae_inn=find(ae_in==1);ce_inn=find(ce_in==1);


water_masses=({'TW','SACW','AAIW','UCDW','NADW','LCDW','AABW'});

% rr=gamma-roo;
g_n_levels_sabrina=[19 26.35 27.1 27.6 27.9 28.12 28.22 35];

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


gamma=gamma.*bucket;

z=pres;x=lon;pos_grid=lon;




%%

% FIGURE 3

figure('Renderer', 'painters', 'Position', [200 200 800 650]) 

h = subplot(2,2,1); 
p = get(h, 'pos');
p(3) = p(3) + 0.06; p(1) = p(1) - 0.05;p(4) = p(4) + 0.06;p(2) = p(2) - 0.02;set(h, 'pos', p);box on; hold on
text(0.02,0.98,'a','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

%plot(transportv,1:length(transportv),'color',[.5 .5 .5]);hold on;
plot(transportvv,1:length(transportv),'k');hold on;
plot(transportvgos281_ek,1:length(transportvgos281_ek),'--r');hold on;
plot(transportvgos3400,1:length(transportvgos3400),'--b');hold on;
%plot(nsum(transportvgos3400argo,1:length(transportvgos3400),'--g');hold on;
plot(nsum(transportvgos3400argo,2),pres,'--g');hold on

%plot(transportv_balanced,1:length(transportv),'color',[.5 .5 .5],'linestyle','--');hold on;
plot(transportvv_balanced,1:length(transportv),'--k');hold on;
% plot(transportvv,1:length(transportv),'k','linewidth',1);hold on;
plot(transportvgos281_balanced,1:length(transportvgos281_ek),'r','linewidth',1);hold on;
plot(transportvgos3400_balanced,1:length(transportvgos3400),'b','linewidth',1);hold on;
% plot(transportvgos3400_balanced,1:length(transportvgos3400),'b','linewidth',1);hold on;
plot(nsum(transportvgos3400argo_balanced,2),pres,'g','linewidth',1);hold on
axis ij;xline(0);
%title('transport per 1db (Sv)')
ylim([0 5500]);xlim([-0.05 0.05]);ylabel('Pressure (db)');xlabel('Transport (Sv)')
xticks([-0.04:0.02:0.04])
%  savefig('figure_3a')
%  print(gcf, '-dpdf','-r600','/Users/gaston/Documents/figure_3a')
% 
legend('LADCP','Î³=28.1','3400db','3400db ARGO','LADCP adj','Î³=28.1 adj','3400db adj','3400db ARGOadj','location','southwest','AutoUpdate', 'off','box','off')

grid on
h = subplot(2,2,2); 
p = get(h, 'pos');
p(3) = p(3) + 0.06; p(1) = p(1) +0.01;p(4) = p(4) + 0.06;p(2) = p(2) - 0.02;set(h, 'pos', p);box on; hold on
text(0.02,0.98,'b','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

tt=nancumsum(transport_by_layer_hg,'omitnan');
tts=nancumsum(transport_by_layer_hg_geos281);
tthg=nancumsum(transport_by_layer_hg_geos3400);

%plot(transportv_balanced,1:length(transportv),'color',[.5 .5 .5],'linestyle','--');hold on;
plot(nancumsum(transportvv_balanced),1:length(transportv),'--k');hold on;
% plot(transportvv,1:length(transportv),'k','linewidth',1);hold on;
plot(nancumsum(transportvgos281_balanced),1:length(transportvgos3400),'r','linewidth',1);hold on;
plot(nancumsum(transportvgos3400_balanced),1:length(transportvgos3400),'b','linewidth',1);hold on;
% plot(transportvgos3400_balanced,1:length(transportvgos3400),'b','linewidth',1);hold on;
plot(nancumsum(nsum(transportvgos3400argo_balanced,2)),1:length(transportvgos3400),'g','linewidth',1);hold on
% legend('LADCP adj','Î³=28.1 adj','3400db adj','3400db ARGOadj','location','south','AutoUpdate', 'off')
axis ij;xline(0);
%title('transport per 1db (Sv)')


ylim([0 5500]);
%xlim([-0.05 0.05]);
ylabel('Pressure (db)');xlabel('Cumulative transport (Sv)')

grid on

h = subplot(2,2,3); 
p = get(h, 'pos');
p(3) = p(3) + 0.06; p(1) = p(1) - 0.05;p(4) = p(4) + 0.05;p(2) = p(2) - 0.02;set(h, 'pos', p);box on; hold on
text(0.02,0.98,'c','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

plot(transport_by_layer_hg,1:length(transport_by_layer_hg),'k--');hold on;
plot(transport_by_layer_hg_geos281,1:length(transport_by_layer_hg),'r--');hold on;
plot(transport_by_layer_hg_geos3400,1:length(transport_by_layer_hg),'b--');hold on;
plot(transport_by_layer_hg_geos3400argo,1:length(transport_by_layer_hg),'g--');hold on;

plot(transport_by_layer_hgb,1:length(transport_by_layer_hg),'k-');hold on;
plot(transport_by_layer_hg_geos281b,1:length(transport_by_layer_hg),'r-');hold on;
plot(transport_by_layer_hg_geos3400b,1:length(transport_by_layer_hg),'b-');hold on;
plot(transport_by_layer_hg_geos3400argob,1:length(transport_by_layer_hg),'g-');hold on;


% plot(transport_by_layer_hg_geos281bar,1:length(transport_by_layer_hg_geos281bar));hold on;
% plot(transport_by_layer_hg_geos3400bar,1:length(transport_by_layer_hg_geos3400bar));hold on;

axis ij;xline(0);

%yticks('manual')
yticklabels({'','SW','SACW','AAIW','UCDW','NADW','LCDW','AABW',''})
grid on
%

transport_msm_by_layer=transport_by_layer_hg_geos281b;
transport_msm_by_laye_imb=transport_by_layer_hg_geos281;

% title('transport by layer (Sv)')
% 
% %legend('ladcp','28.1','3400','281bar','3400bar')
xlabel('Transport (Sv)')
h = subplot(2,2,4); 
p = get(h, 'pos');
%p(3) = p(3) + 0.06; p(1) = p(1) - 0.05;p(4) = p(4) + 0.02;p(2) = p(2) + 0.03;set(h, 'pos', p);box on; hold on
p(3) = p(3) + 0.06; p(1) = p(1) +0.01;p(4) = p(4) + 0.05;p(2) = p(2) -0.02;set(h, 'pos', p);box on; hold on
text(0.02,0.98,'d','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

% plot(tt,1:length(g_n_levels),'k+--');hold on;
% plot(tts,1:length(g_n_levels),'r+--');hold on;
% plot(tthg,1:length(g_n_levels),'b+--');hold on;
% plot(tthga,1:length(g_n_levels),'g+--');hold on;
plot(ttb,1:length(ttb),'k+-');hold on;
plot(ttsb,1:length(ttb),'r+-');hold on;
plot(tthgb,1:length(ttb),'b+-');hold on;
plot(tthgab,1:length(ttb),'g+-');hold on;

axis ij;xline(0);

yticklabels({'','SW','SACW','AAIW','UCDW','NADW','LCDW','AABW',''})
grid on

xlabel('Cumulative transport (Sv)');

set(findall(gcf,'-property','FontSize'),'FontSize',14)
set(findall(gcf,'-property','Linewidth'),'Linewidth',.7)

% savefig('Fig3')
% print(gcf, '-dpdf','-r600','/Users/gaston/Documents/figure_3')
% 
% print(gcf, '-djpeg','-r600','/Users/gaston/Documents/figure_3')
% 
% print(gcf, '-dpng','-r600','/Users/gaston/Documents/figure_3')







for i=2:length(g_n_levels)
for j=1:1400
ind=find(gammaa(:,j)<g_n_levels(i) & gammaa(:,j)>g_n_levels(i-1));

transport_by_layer_hg(i,j)=nansum(transportva(ind,j));
transport_by_layer_hg_geos3400(i,j)=nansum(transportvgos3400a(ind,j));
transport_by_layer_hg_geos3400argo(i,j)=nansum(transportvgos3400argo(ind,j));
transport_by_layer_hg_geos281(i,j)=nansum(transportvgos28_ek(ind,j));


transport_by_layer_hgb(i,j)=nansum(transportva_balanced(ind,j));
transport_by_layer_hg_geos3400b(i,j)=nansum(transportvgos3400a_balanced(ind,j));
transport_by_layer_hg_geos3400argob(i,j)=nansum(transportvgos3400argo_balanced(ind,j));
transport_by_layer_hg_geos281b(i,j)=nansum(transportvgos28_balanced(ind,j));

end
end

% PLOT SECTIONS


%% Fig 3.

%%%THETA PLOT%%%%
figure;hold on;

h = subplot(4,1,1); 
p = get(h, 'pos');

p(3) = p(3) +.12 ; p(1) = p(1) - 0.04;p(4) = p(4) ;p(2) = p(2) + 0.07;set(h, 'pos', p);box on; hold on

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

h = subplot(4,1,2); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.09;p(2) = p(2) + 0.08;set(h, 'pos', p);box on; hold on
text(0.975,0.2,'a','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

plot(lon,v_new100,'k');
hold on
plot(lon_adcp,vvplot,'b');
hold on
plot(lon_adcp,vgos,'r');
axis tight; ylim([-.8 .8])
ylabel('V (m.s^-^1)')
grid on

legend('Hidrography','SADCP','Altimetry','location','north','AutoUpdate', 'off','orientation','horizontal','box','off')
yline(0,'color',[.5 .5 .5])

cmap = jet(7);ind=find(cmap==1);cmap(ind)=0.7;

h = subplot(4,1,3); 
p = get(h, 'pos');

p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.075;p(2) = p(2) + 0.05;set(h, 'pos', p);box on; hold on

%p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.155;p(2) = p(2) - 0.03;set(h, 'pos', p);box on; hold on

text(0.975,0.2,'b','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

% for i=2:8
% plot(lon(1:end-1),transport_by_layer_hg_geos3400b(i,:)*4,'Color', cmap(i-1, :));
% hold on
% 
% end
% %ylim([-2 2])
% 
% plot(lon(1:end-1),nsum(transportvgos3400a_balanced),'k');

% for i=2:8
% plot(lon(1:end-1),transport_by_layer_hg_geos281b(i,:)*4,'Color', cmap(i-1, :));
% hold on
% 
% end
%ylim([-2 2])



plot(lon(1:end-1),nsum(transportvgos28_balanced),'k');

ylabel('Transport (Sv)')

% plot(lon(1:end-1),nsum(transportva_balanced),'k');
yline(0,'color',[.5 .5 .5])


% legend('SW','SACW','AAIW','UCDW', newline,'NADW','LCDW','AABW','Total','location','south','AutoUpdate', 'off','orientation','horizontal','box','off');


%legend({['blue' newline 'line'],'red line'})
hold on;axis tight
axis tight;

% ylim([-5.1 3])

grid on
xticklabels({'50ºW','40ºW','30ºW','20ºW','10ºW','0ºW','10ºE'})


h = subplot(4,1,4); 
p = get(h, 'pos');
p(3) = p(3) + 0.12; p(1) = p(1) - 0.04;p(4) = p(4) + 0.155;p(2) = p(2) - 0.06;set(h, 'pos', p);box on; hold on

%text(0.01,0.96,'c','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

text(0.975,0.15,'c','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')



aa=ncumsum(transport_by_layer_hg_geos3400b,2)*2;
aa=ncumsum(transport_by_layer_hg_geos281b,2)*2;


for i=2:8
plot(lon(1:end-1),aa(i,:),'Color', cmap(i-1, :));
hold on

end

%plot(lon(1:end-1),cumsum(nsum(transportvgos3400a_balanced)),'k');
plot(lon(1:end-1),cumsum(nsum(transport_by_layer_hg_geos281b)),'k');

% plot(lon(1:end-1),nsum(transportva_balanced),'k');
yline(0,'color',[.5 .5 .5])
hold on;axis tight
axis tight;



ylim([-50 30])



xticklabels({'50ºW','40ºW','30ºW','20ºW','10ºW','0ºW','10ºE'})
 
ylabel('Cum. transport (Sv)');

legend('SW(.2)','SACW(.2)','AAIW(.2)','UCDW(.2)','NADW(.2)','LCDW(.2)','AABW(.2)','Total','location','north','AutoUpdate', 'off','orientation','horizontal','box','off','NumColumns',4);


grid on
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(findall(gcf,'-property','Linewidth'),'Linewidth',.7)



% savefig('figure3c_horizontal_transport')
% 
% print(gcf, '-dpdf','-r600','/Users/gaston/Documents/figure3c_horizontal_transport')
% print(gcf, '-djpeg','-r600','/Users/gaston/Documents/figure3c_horizontal_transport')
% print(gcf, '-dpng','-r600','/Users/gaston/Documents/figure3c_horizontal_transport')
% 



%%



% ctt=CT+273.15;


%the area of each pixel
%4579;
%heat capacity of seawater
%3850;
%mean reference density
%1027.5;


% heat_transport_balanced = nsum(nsum((geostrophic_velocity_3400_ek_balanced.*CT(:,1:end - 1))*1027.5*3850*4579))

heat_transport_balanced = nsum(nsum((geostrophic_velocity_281_ek_balanced.*CT(:,1:end - 1))*1027.5*3850*4579))



heat_transport_balanced = nsum(nsum((geostrophic_velocity_281_ek_balanced.*CT(:,1:end - 1))*1027.5*3850*4579))

heat_transport_281 = (geostrophic_velocity_281_ek_balanced.*CT(:,1:end - 1)*1027.5*3850*4579);

% heat_transport_balanced = nsum(nsum((geostrophic_velocity_281_ek_balanced.*theta_new(:,1:end - 1)).*area*1027.5*3850.1e6))

heat_transport=nsum(nsum((geostrophic_velocity_3400_ek.*CT(:,1:end-1))*1027.5*3850*4579))

heat_transportv=nsum(nsum(v_new_ek_balanced(:,1:end-1).*CT(:,1:end-1).*area*1027.5*3850.1e6))

ss=1-(sal_new./(nmean(nmean(sal_new))));

ss=1-(SA./(nmean(nmean(SA))));



fresshwater_transport_28=ss(:,1:end-1).*transportvgos28_balanced;

press=repmat(pres,1,1401);

for i=2:length(g_n_levels)

ind=find(gammaa<g_n_levels(i) & gammaa>g_n_levels(i-1));

temp_by_layer_hg(i)=nanmean(CT(ind));
sal_by_layer_hg(i)=nanmean(SA(ind));

ox_by_layer_hg(i)=nanmean(ox(ind));

pres_by_layer_hg(i)=nanmean(press(ind));

transport_by_layer_hgg(i)=nansum(transportvgos3400a_balanced(ind));


heat_transport_by_layer_hgg(i)=nansum(heat_transport_281(ind));

fresh_transport_by_layer_hgg(i)=nansum(fresshwater_transport_28(ind));

end


area=4579*buckett;

% heat_vgos3400=3900.*CT(:,1:end-1).*ro_new(:,1:end-1).*geostrophic_velocity_3400_ek_balanced.*area;
% 
% %heat_vgos3400=3900.*theta(:,1:end-1).*ro_new(:,1:end-1).*geostrophic_velocity_3400_ek_balanced.*area;
% 
% 
% heat_vgos3400=3900.*CT(:,1:end-1).*ro_new(:,1:end-1).*geostrophic_velocity_3400_ek_balanced.*area;
% 
% heat_vgos3400=3900.*CT(:,1:end-1).*ro_new(:,1:end-1).*geostrophic_velocity_3400_ek.*area;
% 
% nsum(nsum(heat_vgos3400))


heat_vgos281=3900.*CT(:,1:end-1).*ro_new(:,1:end-1).*geostrophic_velocity_281.*area;

3900*3.26*1025*0.1*4579

nsum(nsum(heat_vgos281))



% heat_transport281balance=geostrophic_velocity_281_ek_balanced.*area.*1027.5.*3900.*tt(:,1:end-1);nsum(nsum(heat_transport281balance))


nsum(transportvgos281_ek)


nsum(transport_by_layer_hg_geos281b,2)



nsum(nsum(fresshwater_transport_28))

nmean(nmean(theta_new.*bucket))

nmean(nmean(sal_new.*bucket))

nmean(nmean(SA.*bucket))

nmean(nmean(CT.*bucket))

nmean(nmean(ox_new.*bucket))

%% transports 
cd  /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/go_ship
%load('transport_a11.mat');
load('transport_a10.mat');%transport_by_layer_hg_geos281b_a10=[ transport_by_layer_hg_geos281b_a10 0];

load('transport_a9.mat');

 load('goship_a10.mat', 'sal_a10');load('goship_a9.mat', 'theta_a9')
load('transport_a10.mat', 'transportvgos28_balanced_a10')


%% load a11 data from a11
cd /Users/gaston/Documents/Phd/Datos_CTD_phd/GO_SHIP/A10_Molly/0_a10_33RO20110926_nc_ctd
projectdir2 = '/Users/gaston/Documents/Phd/Datos_CTD_phd/GO_SHIP/A11_Saunders/0_a11_nc_ctd';

info2 = dir( fullfile(projectdir2, '*.nc') );
n_files2 = length(info2);
filenames2 = fullfile( projectdir2, {info2.name} );
ctda11= struct;

for i= 1 : n_files2;
f = filenames2{i};
ctda11.temp{i}=ncread(f,'temperature');
ctda11.pres{i}=ncread(f,'pressure');
ctda11.sal{i}=ncread(f,'salinity');
ctda11.oxy{i}=ncread(f,'oxygen');
ctda11.lat{i}=ncread(f,'latitude');
ctda11.lon{i}=ncread(f,'longitude');

end

tempa11 = cell2mat(ctda11.temp(:));
sala11 = cell2mat(ctda11.sal(:));
presa11 = cell2mat(ctda11.pres(:));
oxya11 = cell2mat(ctda11.oxy(:));
lona11 = cell2mat(ctda11.lon(:));
lata11 = cell2mat(ctda11.lat(:));

[SA_a11, in_ocean] = gsw_SA_from_SP(sala11,presa11,-30,-45);
CT_a11 = gsw_CT_from_t(SA_a11,tempa11,presa11);
ptempa11 = gsw_pt0_from_t(SA_a11,tempa11,presa11);

%% data from a10
cd /Users/gaston/Documents/Phd/Datos_CTD_phd/GO_SHIP/A10_Molly/0_a10_33RO20110926_nc_ctd
projectdir2 = '/Users/gaston/Documents/Phd/Datos_CTD_phd/GO_SHIP/A10_Molly/0_a10_33RO20110926_nc_ctd';

info2 = dir( fullfile(projectdir2, '*.nc') );
n_files2 = length(info2);
filenames2 = fullfile( projectdir2, {info2.name} );
ctda10= struct;

for i= 1 : n_files2;
f = filenames2{i};
ctda10.temp{i}=ncread(f,'temperature');
ctda10.pres{i}=ncread(f,'pressure');
ctda10.sal{i}=ncread(f,'salinity');
ctda10.oxy{i}=ncread(f,'oxygen');
ctda10.lat{i}=ncread(f,'latitude');
ctda10.lon{i}=ncread(f,'longitude');

end

% tempa10 = cell2mat(ctda10.temp(:));
% sala10= cell2mat(ctda10.sal(:));
% presa10 = cell2mat(ctda10.pres(:));
% oxya10 = cell2mat(ctda10.oxy(:));
% lona10 = cell2mat(ctda10.lon(:));
% lata10 = cell2mat(ctda10.lat(:));

tempa10 = cell2mat(ctda10.temp(:));
sala10= cell2mat(ctda10.sal(:));
presa10 = cell2mat(ctda10.pres(:));
oxya10 = cell2mat(ctda10.oxy(:));
lona10 = cell2mat(ctda10.lon(:));
lata10 = cell2mat(ctda10.lat(:));

[SA_a10, in_ocean] = gsw_SA_from_SP(sala10,presa10,-30,-30);

CT_a10 = gsw_CT_from_t(SA_a10,tempa10,presa10);
ptempa10 = gsw_pt0_from_t(SA_a10,tempa10,presa10);


% ss_a10=1-(SA_a10./(nmean(nmean(SA_a10))));
% 
% fresshwater_transport_28_a10=ss_a10(:,1:end-1).*transportvgos28_balanced_a10;

%% data from a9
cd /Users/gaston/Documents/Phd/Datos_CTD_phd/GO_SHIP/A9_King/cchdo_740H20090307/only_a9_profiles
projectdir9 = '/Users/gaston/Documents/Phd/Datos_CTD_phd/GO_SHIP/A9_King/cchdo_740H20090307/only_a9_profiles';

info9 = dir( fullfile(projectdir9, '*.nc') );
n_files9 = length(info9);
filenames9 = fullfile( projectdir9, {info9.name} );
ctda9= struct;

for i= 1 : n_files9;
f = filenames9{i};
ctda9.temp{i}=ncread(f,'temperature');
ctda9.pres{i}=ncread(f,'pressure');
ctda9.sal{i}=ncread(f,'salinity');
ctda9.oxy{i}=ncread(f,'oxygen');
ctda9.lat{i}=ncread(f,'latitude');
ctda9.lon{i}=ncread(f,'longitude');

end

tempa9 = cell2mat(ctda9.temp(:));
sala9= cell2mat(ctda9.sal(:));
presa9 = cell2mat(ctda9.pres(:));
oxya9 = cell2mat(ctda9.oxy(:));
lona9 = cell2mat(ctda9.lon(:));
lata9 = cell2mat(ctda9.lat(:));

[SA_a9, in_ocean] = gsw_SA_from_SP(sala9,presa9,-30,-24);

CT_a9 = gsw_CT_from_t(SA_a9,tempa9,presa9);
ptempa9 = gsw_pt0_from_t(SA_a9,tempa9,presa9);

%% data from merian
%creation_date = ncreadatt('msm_060_1_ctd_126.nc','/','rodb_hdr')

projectdir = '/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/msm60ctd_nc';
info = dir( fullfile(projectdir, '*.nc') );
n_files = length(info);
filenames = fullfile( projectdir, {info.name} );

ctdamerian= struct;


for i= 1 : n_files;
f = filenames{i};
ctdmerian.temp{i}=ncread(f,'TEMP');
ctdmerian.pres{i}=ncread(f,'PRES');
ctdmerian.sal{i}=ncread(f,'PSAL');
ctdmerian.oxy{i}=ncread(f,'DOX2');
ctdmerian.lat{i}=ncread(f,'LATITUDE');
ctdmerian.lon{i}=ncread(f,'LONGITUDE');

end


tempmerian = cell2mat(ctdmerian.temp(:));
salmerian= cell2mat(ctdmerian.sal(:));
presmerian = cell2mat(ctdmerian.pres(:));
oxymerian = cell2mat(ctdmerian.oxy(:));
lonmerian= cell2mat(ctdmerian.lon(:));
latmerian= cell2mat(ctdmerian.lat(:));



[SA_merian, in_ocean] = gsw_SA_from_SP(salmerian,presmerian,-35,-35);

CTmerian = gsw_CT_from_t(SA_merian,tempmerian,presmerian);

ptempmerian = gsw_pt0_from_t(SA_merian,tempmerian,presmerian);




clear gamma gamma_GP geo_str* geos* sal_argo sal_comp sal ss
%% TS diag
% 
% figure
% tsdiagw(ptempmerian,SA_merian); hold on
% % xlim([34.2 36.6]);ylim([-.5 22])
% xlim([34.2 36.4]);ylim([-.5 20])
% % xlim([34.2 35.2]);ylim([-.5 8])
% scatter2 = scatter(SA_a9,ptempa9,.5,[0    1    1.0000],'filled'); 
% scatter2.MarkerFaceAlpha = 0.2;scatter1.MarkerEdgeAlpha = 0.2;
% 
% scatter2 = scatter(SA_a10,ptempa10,.5,[0    0.6    1.0000],'filled'); 
% scatter2.MarkerFaceAlpha = 0.2;scatter1.MarkerEdgeAlpha = 0.2;
% 
% scatter2 = scatter(SA_merian,ptempmerian,.5,[0 0 1],'filled'); hold on
% scatter2.MarkerFaceAlpha = 0.2;scatter2.MarkerEdgeAlpha = 0.2;
% 
% scatter3 = scatter(SA_a11,ptempa11,.5,[0 0 .5],'filled'); hold on
% scatter3.MarkerFaceAlpha = 0.2;scatter3.MarkerEdgeAlpha = 0.2;
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',12)
% set(findall(gcf,'-property','Linewidth'),'Linewidth',.7)


%%

g_n_levels=[19 26.35 27.1 27.6 27.9 28.1 28.27 38 38];


for i=2:length(g_n_levels)

ind=find(gammaa<g_n_levels(i) & gammaa>g_n_levels(i-1));

temp_by_layer_hg(i)=nanmean(CT(ind));
sal_by_layer_hg(i)=nanmean(SA(ind));

ox_by_layer_hg(i)=nanmean(ox(ind));

pres_by_layer_hg(i)=nanmean(press(ind));

transport_by_layer_28(i)=nansum(transportvgos3400a_balanced(ind));
end








%% the figure 

figure('Renderer', 'painters', 'Position', [200 200 800 650]) 
h = subplot(2,2,3); 
p = get(h, 'pos');
p(3) = p(3) + 0.05; p(1) = p(1) - 0.05;p(4) = p(4) + 0.05;p(2) = p(2) - 0.02;set(h, 'pos', p);box on; hold on
text(0.02,0.98,'c','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')


transport_by_layer_hg_geos281b_a10(1)=NaN;transport_by_layer_hg_geos281b_a10(end)=NaN;

transport_by_layer_hg_geos281b_a9(1)=NaN;transport_by_layer_hg_geos281b_a9(end)=NaN;
transport_msm_by_layer(1)=NaN;transport_msm_by_layer(end)=NaN;



hold on;
plot(transport_by_layer_hg_geos281b_a9,1:length(transport_msm_by_layer),'+-','color', [0    .8    1.0000]);hold on;

plot(transport_by_layer_hg_geos281b_a10,1:length(transport_msm_by_layer),'b+-');hold on;

plot(transport_msm_by_layer,1:length(transport_msm_by_layer),'+-','color',[0         0    0.5]);hold on;
% plot(transport_by_layer_hg_geos281b_a11,1:length(transport_msm_by_layer),'+-','color',[0         0    0.5])
hold on;

% plot(transport_by_layer_hg_geos281bar,1:length(transport_by_layer_hg_geos281bar));hold on;
% plot(transport_by_layer_hg_geos3400bar,1:length(transport_by_layer_hg_geos3400bar));hold on;
ylim([1.5 8.5])
axis ij;xline(0);

%yticks('manual')
yticklabels({'SW','SACW','AAIW','UCDW','NADW','LCDW','AABW',''})
grid on


legend('24.5ºS','30ºS','34.5ºS','location','southwest','AutoUpdate', 'off','box','off')



% title('transport by layer (Sv)')
% 
% %legend('ladcp','28.1','3400','281bar','3400bar')
xlabel('Transport (Sv)')
xlim([-20 15])


h = subplot(2,2,4); 
p = get(h, 'pos');
p(3) = p(3) + 0.05; p(1) = p(1) - 0.05;p(4) = p(4) + 0.05;p(2) = p(2) - 0.02;set(h, 'pos', p);box on; hold on
text(0.02,0.98,'d','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

tsdiagw(ptempmerian,SA_merian); hold on
% xlim([34.2 36.6]);ylim([-.5 22])
xlim([34.3 36.1]);ylim([-.5 18])
% xlim([34.2 35.2]);ylim([-.5 8])

scatter2 = scatter(SA_a9,ptempa9,.5,[0    .8    1.0000],'filled'); 
scatter2.MarkerFaceAlpha = 0.2;scatter1.MarkerEdgeAlpha = 0.2;


scatter2 = scatter(SA_a10,ptempa10,.5,'b','filled'); 
scatter2.MarkerFaceAlpha = 0.2;scatter1.MarkerEdgeAlpha = 0.2;

scatter2 = scatter(SA_merian,ptempmerian,.5,[0 0 .5],'filled'); hold on
scatter2.MarkerFaceAlpha = 0.2;scatter2.MarkerEdgeAlpha = 0.2;

grid on

 set(gca,'Xtick',34.3:0.5:36.3)
 set(gca,'XtickLabel',34.5:2:36.5)


xticks([34.4:0.2:36.4])
yticks([0:1:25])




% scatter3 = scatter(SA_a11,ptempa11,.5,[0 0 .5],'filled'); hold on
% scatter3.MarkerFaceAlpha = 0.2;scatter3.MarkerEdgeAlpha = 0.2;


% axes('position',[0.705,0.1,0.2,0.2]);
% 
% tsdiagw(ptempmerian,SA_merian); hold on
% % xlim([34.2 36.6]);ylim([-.5 22])
% % xlim([34.2 36.4]);ylim([-.5 20])
%  xlim([34.2 35.2]);ylim([-.5 6])
% scatter2 = scatter(SA_a10,ptempa10,.5,[0    0.7461    1.0000],'filled'); 
% scatter2.MarkerFaceAlpha = 0.2;scatter1.MarkerEdgeAlpha = 0.2;
% 
% scatter2 = scatter(SA_merian,ptempmerian,.5,[0 0 1],'filled'); hold on
% scatter2.MarkerFaceAlpha = 0.2;scatter2.MarkerEdgeAlpha = 0.2;
% 
% scatter3 = scatter(SA_a11,ptempa11,.5,[0 0 .5],'filled'); hold on
% scatter3.MarkerFaceAlpha = 0.2;scatter3.MarkerEdgeAlpha = 0.2;
% 
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',12)
% set(findall(gcf,'-property','Linewidth'),'Linewidth',.7)

hs = subplot(2,2,2); 

text(0.02,0.98,'b','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')
hold on
plot(lon(1:end-1),cumsum(nsum(transportvgos28_balanced_a9)),'color',[0 .8 1]);

hold on
plot(lon(1:end-1),cumsum(nsum(transportvgos28_balanced)),'color',[0 0 .5]);
hold on
plot(lon(1:end-1),cumsum(nsum(transportvgos28_balanced_a10)),'b');
pos_grida11=[-60:0.05:18];
% hold on
% plot(pos_grida11(1:end-1),cumsum(nsum(transportvgos28_balanced_a11)),'color',[0 0 .5]);


% plot(lon(1:end-1),nsum(transportva_balanced),'k');
yline(0,'color',[.5 .5 .5])
hold on;axis tight
axis tight;box on

xlim([-60 20])

xticks([-60:10:20])
xticklabels({'60ºW','50ºW','40ºW','30ºW','20ºW','10ºW','0ºE','10ºE','20ºE'})
 
ylabel('Cum. transport (Sv)');
grid on


h = subplot(2,2,1); 




p = get(h, 'pos');
%[[x0,y0,width,height]

p(3) = p(3) + 0.5; p(1) = p(1) - 0.06;p(4) = p(4) - 0.1;p(2) = p(2) + 0.15;set(h, 'pos', p);box on; hold on

m_proj('mercator','long',[-60 20],'lat',[-36.5 -21])

 [CS,CH]=m_etopo2('contourf',[-3500 -200],'edgecolor','none');caxis([-5500 7000]);colormap(flip(gray(7)));  

%  [CS,CH]=m_etopo2('contour',[-3500 -3500],'color',[.5 .5 .5]);%caxis([-5500 7000]);%colormap(flip(gray(7)));  
% hold on
%  [CS,CH]=m_etopo2('contour',[-200 -200],'color',[.2 .2 .2]);%caxis([-5500 7000]);%colormap(flip(gray(7)));  


m_gshhs_h('patch',[0.81818 0.77647 0.70909]);

m_plot(cell2mat(ctda9.lon),cell2mat(ctda9.lat),'*','color',[0 .8 1])

m_plot(cell2mat(ctdmerian.lon),cell2mat(ctdmerian.lat),'*','color',[0 0 .5])

m_plot(cell2mat(ctda10.lon),cell2mat(ctda10.lat),'b*')

% m_plot(cell2mat(ctda11.lon),cell2mat(ctda11.lat),'*','color',[0 0 .5])

text(0.02,0.98,'a','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')


th=m_annotation('textarrow',[-27 -27],[-34.8 -32.5],'String','0.23 PW','color', [0.6 0 0]);
th=m_annotation('textarrow',[-23.5 -23.5],[-30.7 -27.5],'String','0.65 PW','color', [0.6 0 0]);
th=m_annotation('textarrow',[-19 -19],[-24.9 -21.6],'String','1.08 PW','color', [0.6 0 0]);


% th=m_annotation('textarrow',[-27 -27],[-34.7 -32.5],'String','0.18 PW','color', [0.6 0 0]);
% th=m_annotation('textarrow',[-23.5 -23.5],[-30.6 -27],'String','0.80 PW','color', [0.6 0 0]);
% th=m_annotation('textarrow',[-19 -19],[-24.8 -21.5],'String','1.20 PW','color', [0.6 0 0]);


th=m_annotation('textarrow',[-2 -2],[-34.8 -32],'String','0.21 Sv','color',[0 0.4 0]);
th=m_annotation('textarrow',[-4.5 -4.5],[ -29.5 -31.2],'String','-0.02 Sv','color',[0 0.4 0]);
th=m_annotation('textarrow',[-6 -6], [-23 -25.9],'String','-0.23 Sv','color',[0 0.4 0]);




m_grid('xtick',[-60:10:20],'xticklabels',[],'linestyle','none')


% m_gshhs_h('patch', [0.81818 0.77647 0.70909]);colorbar


%[[x0,y0,width,height]

ps = get(hs, 'pos');

ps(3) = ps(3) + 0.5; ps(1) = ps(1) - 0.5;ps(4) = ps(4) - 0.15;ps(2) = ps(2) - 0.045;set(hs, 'pos', ps);box on; hold on



set(findall(gcf,'-property','FontSize'),'FontSize',13)
set(findall(gcf,'-property','Linewidth'),'Linewidth',.7)


set(gcf,'color','w');



cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september

print(gcf, '-djpeg','-r600','figure_15_september_go_ship')

print(gcf, '-depsc','-painters' ,'-r600','figure_15_september_go_ship')


% 

% savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/figure_9_merian')
% print(gcf, '-dpng','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/go_ship/figure_9_merian')



lonn=lon(1:end-1)+0.05/2;

pp=1:5501;pp=pp';

Ctt=griddata(lon,pp,CT,lonn,pp);

% Ctt=griddata(lon,pres_a9,CT_a9,lonn,pres_a9);

 area=4579*buckett/1e6;

heat_transport=geostrophic_velocity_281_ek_balanced.*Ctt*1027.5*3900.*area*1e6;
nsum(nsum(heat_transport))


heat_transport=geostrophic_velocity_281_ek_balanced.*Ctt*1027.5*3900.*area*1e6;
nsum(nsum(heat_transport))

heat_transport=geostrophic_velocity_281_ek_balanced.*CT(:,1:end-1)*1027.5*3900.*area(:,1:end-1)*1e6;

nsum(nsum(heat_transport))


heat_transport=geostrophic_velocity_281_ek_balanced.*Ctt*1027.5*3900.*length_grid_point;
nsum(nsum(heat_transport))


heat_transportv=v_new_ek_balanced(:,1:end-1).*Ctt*1027.5*3900.*length_grid_point;
nsum(nsum(heat_transportv))


heat_transportv=v_new_ek_balanced(:,1:end-1).*Ctt*1027.5*3900.*length_grid_point;
nsum(nsum(heat_transportv))



heat_transportv=v_new_ek_balanced(:,1:end-1).*Ctt*1027.5*3900.*length_grid_point;
nsum(nsum(heat_transportv))

% sal_anom=1-(SA_a9-nmean(SA_a9.*bucket,2));

salm=nmean(nmean(SA.*bucket));

sal_anom=1-(SA/salm);

sal_anom=sal_anom.*bucket;

pp=1:5501;pp=pp';

sal_anomm=griddata(lon,pp,sal_anom,lonn,pp);

freshwater_transport=geostrophic_velocity_281_ek_balanced.*sal_anomm.*length_grid_point;

nsum(nsum(freshwater_transport))



nsum(nsum(transportvgos28_balanced(1:2000,41:50))





