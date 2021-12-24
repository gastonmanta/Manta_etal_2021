clear all; close all;clc

cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian
load MSM60_ctd_november.mat
load SA_CT_sined

%fixing a bug
 SA_sined(1000:end,900:980)=SA(1000:end,900:980); SA_sined(:,1:50)=SA(:,1:50);SA_sined(:,end-50:end)=SA(:,end-50:end);
 CT_sined(1000:end,900:980)=CT(1000:end,900:980); CT_sined(:,1:50)=CT(:,1:50);CT_sined(:,end-50:end)=CT(:,end-50:end);


% %activate calculations without eddies
% SA=SA_sined;CT=CT_sined;

load vu;v=v/100;u=u/100;

load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian/stations.mat')
 
ae_in=nansum(ae_in);ce_in=nansum(ce_in);
ae_inn=find(ae_in==1);ce_inn=find(ce_in==1);

water_masses=({'TW','SACW','AAIW','UCDW','NADW','LCDW','AABW'});

% geo_strf_dyn_height = gsw_geo_strf_dyn_height(SA,CT,pres,0);
% [geostrophic_velocity0, mid_lat, mid_long] =gsw_geostrophic_velocity(geo_strf_dyn_height,lon,lat,0);

load('argo_merian_gridded.mat')


argo_replace=1700;
ptemp_comp=cat(1,ptemp_ar(1:argo_replace,:),CT(argo_replace+1:end,:));
sal_comp=cat(1,sal_ar(1:argo_replace,:),SA(argo_replace+1:end,:));


% ptemp_comp=cat(1,ptemp_ar(1:argo_replace,:),CT(argo_replace+1:end,:));
% sal_comp=cat(1,sal_ar(1:argo_replace,:),SA(argo_replace+1:end,:));

%THE SHALLOWER THAN 2000M I keep the original
sal_comp(:,1:20)=SA(:,1:20);sal_comp(:,end-20:end)=SA(:,end-20:end);

ptemp_comp(:,1:20)=CT(:,1:20);ptemp_comp(:,end-20:end)=CT(:,end-20:end);


% % SA=sal_comp;CT=ptemp_comp;
% spd=sqrt(u.*u+v.*v);
% spd_z=nanmean(spd,2);[val_min,pres_min]=min(spd_z); %2204 minimun spd 0.0630

bucket(end,:)=NaN;buckett(:,1040)=buckett(:,1039);

%4579 m2 is the area of each pixel (1m.0.05degrees), bucket land mask and
%1e6 to get it into sverdrups, so it's not just the area what it means
length_grid_point=sw_dist([-34.5 -34.5],[-30 -30.05],'km')*1000
areaa=length_grid_point*buckett/1e6;
area=length_grid_point*bucket/1e6;

%BUCKET IS THE LAND MASK
temp_new=temp.*bucket;sal_new=sal.*bucket;ox_new=ox.*bucket;u_new=u.*bucket;v_new=v.*bucket;
gamma_new=gamma.*bucket;ro_new=ro.*bucket;theta_new=theta.*bucket;%spd_new=spd.*bucket;

water_masses=({'TW','SACW','AAIW','UCDW','NADW','LCDW','AABW'});

%following valla et al 2018
g_n_levels=[26.35 27.1 27.6 27.9 28.1 28.27];


xx=1:123;tt=31; 
pres=1:5501;%pres=repmat(pres,1401,1);pres=pres';


%GEOSTROPHIC VELOCITY WITH RESPECT TO 3400dbar AND 0dbar
geo_strf_dyn_height_3400 = gsw_geo_strf_dyn_height(SA,CT,pres',3400);
geostrophic_velocity_3400=gsw_geostrophic_velocity(geo_strf_dyn_height_3400,lon,lat,0);


%GEOSTROPHIC VELOCITY WITH RESPECT TO 3400dbar AND 0dbar
geo_strf_dyn_height_3400_sined = gsw_geo_strf_dyn_height(SA_sined,CT_sined,pres',3400);
geostrophic_velocity_3400_sined=gsw_geostrophic_velocity(geo_strf_dyn_height_3400_sined,lon,lat,0);

geo_strf_dyn_height_0 = gsw_geo_strf_dyn_height(SA,CT,pres',0);
[geostrophic_velocity_0, mid_lat, mid_long]=gsw_geostrophic_velocity(geo_strf_dyn_height_0,lon,lat,0);

geo_strf_dyn_height_0_sined = gsw_geo_strf_dyn_height(SA_sined,CT_sined,pres',0);
[geostrophic_velocity_0_sined, mid_lat, mid_long]=gsw_geostrophic_velocity(geo_strf_dyn_height_0_sined,lon,lat,0);


%SAME WITH ARGO
geo_strf_dyn_height_0argo = gsw_geo_strf_dyn_height(sal_comp,ptemp_comp,pres',0);
[geostrophic_velocity_0argo, mid_lat, mid_long]=gsw_geostrophic_velocity(geo_strf_dyn_height_0argo,lon,lat,0);

geostrophic_velocity_3400argo = gsw_geo_strf_dyn_height(sal_comp,ptemp_comp,pres',3400);
geostrophic_velocity_3400argo=gsw_geostrophic_velocity(geo_strf_dyn_height_3400,lon,lat,0);

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
geostrophic_velocity_281argo(i)=nmean(geostrophic_velocity_0argo(ind,i));
geostrophic_velocity_281_sined(i)=nmean(geostrophic_velocity_0_sined(ind,i));
end

ladcp281=movmean(ladcp281,2);
%GEOSTROPHIC Velocity with respect to gamma 28.1

geostrophic_velocity281r=repmat(geostrophic_velocity_281,5501,1);
geostrophic_velocity_281=geostrophic_velocity0-geostrophic_velocity281r; clear geostrophic_velocity281r;

%without eddies
geostrophic_velocity281r_sined=repmat(geostrophic_velocity_281_sined,5501,1);
geostrophic_velocity_281_sined=geostrophic_velocity_0_sined-geostrophic_velocity281r_sined; clear geostrophic_velocity281r_sined;

%GEOSTROPHIC Velocity from argo with respect to gamma 28.1
geostrophic_velocity281rargo=repmat(geostrophic_velocity_281argo,5501,1);
geostrophic_velocity_281argo=geostrophic_velocity_0argo-geostrophic_velocity281rargo; clear geostrophic_velocity281r;


%REMOVING SPIKES

% geostrophic_velocity_281(:,441)=((geostrophic_velocity_281(:,440)+geostrophic_velocity_281(:,442))./2);
% geostrophic_velocity_281(:,1392)=((geostrophic_velocity_281(:,1391)+geostrophic_velocity_281(:,1393))./2);
% 
% geostrophic_velocity_281argo(:,441)=((geostrophic_velocity_281argo(:,440)+geostrophic_velocity_281argo(:,442))./2);
% geostrophic_velocity_281argo(:,1392)=((geostrophic_velocity_281argo(:,1391)+geostrophic_velocity_281argo(:,1393))./2);
% 
% geostrophic_velocity_3400(:,441)=((geostrophic_velocity_3400(:,440)+geostrophic_velocity_3400(:,442))./2);
% 
% %geostrophic_velocity_3400(:,1380:end)=movmean(geostrophic_velocity_3400(:,1380:end),5,1);
% geostrophic_velocity_3400(:,1392)=((geostrophic_velocity_3400(:,1391)+geostrophic_velocity_3400(:,1393))./2);
% 
% geostrophic_velocity_3400(:,1:4)=0;geostrophic_velocity_3400(:,end-4:end)=0;
% geostrophic_velocity_3400(geostrophic_velocity_3400>1.5)=1.5;
% geostrophic_velocity_3400(geostrophic_velocity_3400<-1.5)=-1.5;
% 
% geostrophic_velocity_281(:,1:4)=0;geostrophic_velocity_281(:,end-4:end)=0;
% geostrophic_velocity_281(geostrophic_velocity_281>1.5)=1.5;
% geostrophic_velocity_281(geostrophic_velocity_281<-1.5)=-1.5;
% 
% geostrophic_velocity_281argo(:,1:4)=0;geostrophic_velocity_281argo(:,end-4:end)=0;
% geostrophic_velocity_281argo(geostrophic_velocity_281argo>1.5)=1.5;
% geostrophic_velocity_281argo(geostrophic_velocity_281argo<-1.5)=-1.5;

% geostrophic_velocity_3400(:,1:20)=movmean(geostrophic_velocity_3400(:,1:20),5,1);
% geostrophic_velocity_3400=geostrophic_velocity_3400.*buckett;
% geostrophic_velocity_3400(:,441)=((geostrophic_velocity_281(:,440)+geostrophic_velocity_281(:,442))./2);
% geostrophic_velocity_281(:,1380:end)=movmean(geostrophic_velocity_281(:,1380:end),5,1);
% geostrophic_velocity_281(:,1:20)=movmean(geostrophic_velocity_281(:,1:20),5,1);
% geostrophic_velocity_2811=geostrophic_velocity_281.*buckett;

% %TRANSPORT, BASICALLY MULTYPLYING PER area and taking to Sverdrups
 area=4579*buckett/1e6;
% 
% %the volume contribution of each pixel
% 0.004579;
% %the amount of pixels to repart the transport
% 1400*50;
% %the ekman transport
% -0.42;
% -0.0013 m/s is the velocity to add in the upper 50 dbar to obtain -0.42
% Ekman transport
% 
%%

%ADDS THE EKMAN COMPONENT, AND ALSO CREATES A BALANCED BY ADDING BAROTROPIC
%EQUAL VELOCITY
geostrophic_velocity_281_ek=geostrophic_velocity_281;
geostrophic_velocity_281_ek(1:50,:)=geostrophic_velocity_281_ek(1:50,:)-0.0013;
geostrophic_velocity_281_ek=geostrophic_velocity_281_ek.*buckett;
geostrophic_velocity_281_ek_balanced=geostrophic_velocity_281_ek-nsum(nsum(geostrophic_velocity_281_ek))./nsum(nsum((~isnan(geostrophic_velocity_281_ek))));
% geostrophic_velocity_281_ek_g=griddata(mid_long(1,:),pres,geostrophic_velocity_281_ek,lon',pres);
% 

geostrophic_velocity_281_eka=geostrophic_velocity_281argo;
geostrophic_velocity_281_eka(1:50,:)=geostrophic_velocity_281_eka(1:50,:)-0.0013;
geostrophic_velocity_281_eka=geostrophic_velocity_281_eka.*buckett;
geostrophic_velocity_281_ek_balanced_argo=geostrophic_velocity_281_eka-nsum(nsum(geostrophic_velocity_281_eka))./nsum(nsum((~isnan(geostrophic_velocity_281_eka))));


geostrophic_velocity_281_ek_sined=geostrophic_velocity_281_sined;
geostrophic_velocity_281_ek_sined(1:50,:)=geostrophic_velocity_281_sined(1:50,:)-0.0013;
geostrophic_velocity_281_ek_sined=geostrophic_velocity_281_ek_sined.*buckett;
geostrophic_velocity_281_ek_balanced_sined=geostrophic_velocity_281_ek_sined-nsum(nsum(geostrophic_velocity_281_ek_sined))./nsum(nsum((~isnan(geostrophic_velocity_281_ek_sined))));


geostrophic_velocity_281_ek_ladcp=geostrophic_velocity_281;
geostrophic_velocity_281_ek_ladcp(1:50,:)=geostrophic_velocity_281_ek(1:50,:)-0.0013;
geostrophic_velocity_281_ek_ladcp=geostrophic_velocity_281_ek.*buckett;
ind=find(isnan(ladcp281));ladcp281(ind)=0;
geostrophic_velocity_281_ek_ladcp=(geostrophic_velocity_281_ek+ladcp281).*buckett;
geostrophic_velocity_281_ek_balanced_ladcp=geostrophic_velocity_281_ek_ladcp-nsum(nsum(geostrophic_velocity_281_ek_ladcp))./nsum(nsum((~isnan(geostrophic_velocity_281_ek_ladcp))));

vb=fillmissing(v_new,'nearest',1);vb=nmean(vb(end-20:end,:));vb(1390:end)=0;vb(1:15)=0;
geostrophic_velocity_281_ek_balanced_bottom_track=(geostrophic_velocity_281_ek+vb(1:end-1).*buckett);

geostrophic_velocity_281_ek_balanced_bottom_track=geostrophic_velocity_281_ek_balanced_bottom_track-nsum(nsum(geostrophic_velocity_281_ek_balanced_bottom_track))./nsum(nsum((~isnan(geostrophic_velocity_281_ek_balanced_bottom_track))));


geostrophic_velocity_3400_ek_sined=geostrophic_velocity_3400_sined;
geostrophic_velocity_3400_ek_sined(1:50,:)=geostrophic_velocity_3400_ek_sined(1:50,:)-0.0013;
geostrophic_velocity_3400_ek_sined=geostrophic_velocity_3400_ek_sined.*buckett;
geostrophic_velocity_3400_ek_balanced_sined=geostrophic_velocity_3400_ek_sined-nsum(nsum(geostrophic_velocity_3400_ek_sined))./nsum(nsum((~isnan(geostrophic_velocity_3400_ek_sined))));

geostrophic_velocity_3400_ek=geostrophic_velocity_3400;
geostrophic_velocity_3400_ek(1:50,:)=geostrophic_velocity_3400_ek(1:50,:)-0.0013;
geostrophic_velocity_3400_ek=geostrophic_velocity_3400_ek.*buckett;
geostrophic_velocity_3400_ek_balanced=geostrophic_velocity_3400_ek-nsum(nsum(geostrophic_velocity_3400_ek))./nsum(nsum((~isnan(geostrophic_velocity_3400_ek))));


geostrophic_velocity_3400argo_ek=geostrophic_velocity_3400argo;
geostrophic_velocity_3400argo_ek(1:50,:)=geostrophic_velocity_3400argo_ek(1:50,:)-0.0013;
geostrophic_velocity_3400argo_ek=geostrophic_velocity_3400argo_ek.*buckett;
geostrophic_velocity_3400_ek_balanced_argo=geostrophic_velocity_3400argo_ek-nsum(nsum(geostrophic_velocity_3400argo_ek))./nsum(nsum((~isnan(geostrophic_velocity_3400argo_ek))));


% v_new_ek=v_new;
% v_new_ek(1:50,:)=v_new_ek(1:50,:)-0.0013;
% v_new_ek=v_new_ek.*bucket;
% v_new_ek_balanced=v_new_ek-nsum(nsum(v_new_ek))./nsum(nsum((~isnan(v_new_ek))));


transportvgos28ek_balanced=geostrophic_velocity_281_ek_balanced.*areaa;%transport per area

transportvgos28ek_balanced_sined=geostrophic_velocity_281_ek_balanced_sined.*areaa;%transport per area

transportvgos28ek_balanced_argo=geostrophic_velocity_281_ek_balanced_argo.*areaa;%transport per area

transportvgos28_ek=geostrophic_velocity_281.*areaa;%transport per area

transportvgos3400=geostrophic_velocity_3400_ek.*areaa;%transport per area

transportvgos3400_balanced=geostrophic_velocity_3400_ek_balanced.*areaa;%transport per area

transportvgos3400_balanced_sined=geostrophic_velocity_3400_ek_balanced_sined.*areaa;%transport per area

transportvgos3400_balanced_argo=geostrophic_velocity_3400_ek_balanced_argo.*areaa;%transport per area

areav=length_grid_point*bucket/1e6;
transportv=v_new.*areav;%transport per area

transportvgos28ek_balanced_ladcp=geostrophic_velocity_281_ek_balanced_ladcp.*areaa;%transport per area

transportvgos28ek_balanced_bottom_track=geostrophic_velocity_281_ek_balanced_bottom_track.*areaa;

%% AMOC MHT and FWT


%the area of each pixel
%4579;
%heat capacity of seawater
%3850;
%mean reference density
%1027.5;

amoc28b=nmax(nmax(ncumsum(nsum(transportvgos28ek_balanced,2))))

amoc28b_sined=nmax(nmax(ncumsum(nsum(transportvgos28ek_balanced_sined,2))))

amoc28a=nmax(nmax(ncumsum(nsum(transportvgos28ek_balanced_argo,2))))

% amoc28=nmax(nmax(ncumsum(nsum(transportvgos28_ek,2))))

amoc3400b=nmax(nmax(ncumsum(nsum(transportvgos3400_balanced,2))))

amoc3400_sined=nmax(nmax(ncumsum(nsum(transportvgos3400_balanced_sined,2))))

amoc3400a=nmax(nmax(ncumsum(nsum(transportvgos3400_balanced_argo,2))))

amocv=nmax(nmax(ncumsum(nsum(transportv,2))))

amoc_ladcp=nmax(nmax(ncumsum(nsum(transportvgos28ek_balanced_ladcp,2))))

amoc_bottom_track=(nmax(ncumsum(nsum(transportvgos28ek_balanced_bottom_track,2))))



% amoc3400a=nmax(nmax(ncumsum(nsum(transportvgos3400a_balanced,2))))


% amoc=nmax(nmax(ncumsum(nsum(transportvgos28ek_balanced_bottom_track,2))))
% 

% amoc=nmax(nmax(ncumsum(nsum(geostrophic_velocity_3400argo_ek,2))))

% pp=1:5501;
% Ctt=griddata(lon',pp,CT,lonn',pp);
% heat_transport_balanced = nsum(nsum((geostrophic_velocity_281_ek_balancedd.*CT*1027.5*3850*4579)))
%MHT vgos

heat_transport_balanced28 = nsum(nsum((geostrophic_velocity_281_ek_balanced.*CT(:,1:end - 1))*1027.5*3850*4579))


heat_transport_balanced28a = nsum(nsum((geostrophic_velocity_281_ek_balanced_argo.*ptemp_comp(:,1:end - 1))*1027.5*3850*4579))

heat_transport_balanced28sined = nsum(nsum((geostrophic_velocity_281_ek_balanced_sined.*CT_sined(:,1:end - 1))*1027.5*3850*4579))

heat_transport_balanced3400 = nsum(nsum((geostrophic_velocity_3400_ek_balanced.*CT(:,1:end - 1))*1027.5*3850*4579))

heat_transport_balanced3400a = nsum(nsum((geostrophic_velocity_3400_ek_balanced_argo.*ptemp_comp(:,1:end - 1))*1027.5*3850*4579))

heat_transport_balanced3400sined = nsum(nsum((geostrophic_velocity_3400_ek_balanced_sined.*CT_sined(:,1:end - 1))*1027.5*3850*4579))

heat_transport_balanced_ladcp = nsum(nsum((geostrophic_velocity_281_ek_balanced_ladcp.*CT(:,1:end - 1))*1027.5*3850*4579))

%MHT ladcp
heat_transport_v = nsum(nsum((v_new(:,1:end - 1).*CT(:,1:end - 1))*1027.5*3850*4579))

%MHT argo
heat_transport=nsum(nsum((geostrophic_velocity_3400_ek_balanced_argo.*CT(:,1:end-1))*1027.5*3850*4579))



%ss=1-(sal_new./(nmean(nmean(sal_new))));

ss=1-(SA./(nmean(nmean(SA))));

ss_sined=1-(SA_sined./(nmean(nmean(SA_sined))));

ssa=1-(sal_comp./(nmean(nmean(sal_comp))));

fresshwater_transport_28=nsum(nsum(ss(:,1:end-1).*transportvgos28ek_balanced))

fresshwater_transport_28a=nsum(nsum(ssa(:,1:end-1).*transportvgos28ek_balanced_argo))

fresshwater_transport_28_sined=nsum(nsum(ss_sined(:,1:end-1).*transportvgos28ek_balanced_sined))

fresshwater_transport_3400=nsum(nsum(ss(:,1:end-1).*transportvgos3400_balanced))

fresshwater_transport_3400a=nsum(nsum(ssa(:,1:end-1).*transportvgos3400_balanced_argo))

fresshwater_transport_3400_sined=nsum(nsum(ss_sined(:,1:end-1).*transportvgos3400_balanced_sined))


fresshwater_transport_v=nsum(nsum(ss.*transportv))

fresshwater_transport_28_ladcp=nsum(nsum(ss(:,1:end-1).*transportvgos28ek_balanced_ladcp))


 % BOTTOM TRACK bottom track

SAA=fillmissing(SA,'nearest',1);CTT=fillmissing(CT,'nearest',1);
press=1:5501';

geo_strf_dyn_height_5500 = gsw_geo_strf_dyn_height(SAA,CTT,press',5500);

% geo_strf_dyn_height_0 = gsw_geo_strf_dyn_height(sal_new,theta_new,pres',0);
[geostrophic_velocity_5500, mid_lat, mid_long]=gsw_geostrophic_velocity(geo_strf_dyn_height_5500,lon,lat,0);

% geostrophic_velocity_00= fillmissing(geostrophic_velocity_0,'nearest',1);

v_neww=fillmissing(v_new,'nearest',1);

vbottom_track=nmean(v_neww(end-50:end,:),1);vbottom_track(1:10)=0;vbottom_track(end-10:end)=0;


vel_bottom_track=geostrophic_velocity_5500+vbottom_track(1:end-1);


geostrophic_velocity_5500_ek=vel_bottom_track;
geostrophic_velocity_5500_ek(1:50,:)=vel_bottom_track(1:50,:)-0.0013;

transport_bottom_track=geostrophic_velocity_5500_ek.*areaa;

transport_bottom_track_balanced=transport_bottom_track-nsum(nsum(transport_bottom_track))./nsum(nsum((~isnan(transport_bottom_track))));


heat_transport_bottom_tracK = nsum(nsum((transport_bottom_track_balanced.*CT(:,1:end - 1))*1027.5*3850*1e6))

fresshwater_transport__bottom_tracK=nsum(nsum(ss(:,1:end-1).*transport_bottom_track_balanced))

amoc_bottom_track=(nmax(ncumsum(nsum(transportvgos28ek_balanced_bottom_track,2))))
% transportv=nansum(transportva,2);transportvv=movmean(transportv,100);%averaged in z

% 
% transportva_balanced=v_new_ek_balanced(:,2:end).*area;%transport per area
% transportv_balanced=nansum(transportva_balanced,2);transportvv_balanced=movmean(transportv_balanced,100);%averaged in z



g_n_levels=[19 26.35 27.1 27.6 27.9 28.1 28.27 38 38];
transportvv=transportv(:,1:end-1);
heat_transport_balanced28 = geostrophic_velocity_281_ek_balanced.*CT(:,1:end - 1)*1027.5*3850*4579;

fresshwater_transport_28=ss(:,1:end-1).*transportvgos28ek_balanced;

transportvgos28ek_balanced_west=transportvgos28ek_balanced;transportvgos28ek_balanced_west(:,700:end)=NaN;

transportvgos28ek_balanced_east=transportvgos28ek_balanced;transportvgos28ek_balanced_east(:,1:700)=NaN;



for i=2:length(g_n_levels)
ind=find(gammaa(:,1:end-1)<g_n_levels(i) & gammaa(:,1:end-1)>g_n_levels(i-1));

transport_by_layer_hg(i)=nansum(transportvv(ind));
transport_by_layer_hg_geos3400(i)=nansum(transportvgos3400_balanced_argo(ind));
transport_by_layer_hg_geos3400argo(i)=nansum(transportvgos28ek_balanced_argo(ind));
transport_by_layer_hg_geos281(i)=nansum(transportvgos28ek_balanced(ind));
transport_by_layer_hg_geos281imb(i)=nansum(transportvgos28_ek(ind));

transport_by_layer_hg_geos281a(i)=nansum(transportvgos28ek_balanced_argo(ind));
transport_by_layer_hg_geos281ladcp(i)=nansum(transportvgos28ek_balanced_ladcp(ind));

transport_by_layer_bottom_track(i)=nansum(transport_bottom_track_balanced(ind));


transport_by_layer_hg_geos281_west(i)=nansum(transportvgos28ek_balanced_west(ind));

transport_by_layer_hg_geos281_east(i)=nansum(transportvgos28ek_balanced_east(ind));


heat_transport_by_layer_hg_geos281(i)=nansum(heat_transport_balanced28(ind));
fwt_transport_by_layer_hg_geos281(i)=nansum(fresshwater_transport_28(ind));
% ttww(i)=nansum(transport_bottom_track_balanced(ind));
% transport_by_layer_hgb(i)=nansum(transportva_balanced(ind));
% transport_by_layer_hg_geos3400b(i)=nansum(transportvgos3400a_balanced(ind));
% transport_by_layer_hg_geos3400argob(i)=nansum(transportvgos3400argo_balanced(ind));
% transport_by_layer_hg_geos281b(i)=nansum(transportvgos28_balanced(ind));

end





heat_transport_balanced28f = geostrophic_velocity_281_ek_balanced.*CT(:,1:end - 1)*1027.5*3850*4579;

heat_transport_balanced28fsined = geostrophic_velocity_281_ek_balanced_sined.*CT_sined(:,1:end - 1)*1027.5*3850*4579;

% figure;pcolor(heat_transport_balanced28f);shading interp;axis ij
% 
% ylim([0 1000])
% caxis([-3e11 3e11])
% colormap(rednblue)




% figure;
% plot(nsum(heat_transport_balanced28f(1:1700,:),1))
% hold on
% plot(nsum(heat_transport_balanced28fsined(1:1700,:),1))



%% new fig 3

% FIGURE 3

figure('Renderer', 'painters', 'Position', [200 200 800 650]) 

h = subplot(2,2,1); 
p = get(h, 'pos');
p(3) = p(3) + 0.06; p(1) = p(1) - 0.05;p(4) = p(4) + 0.06;p(2) = p(2) - 0.02;set(h, 'pos', p);box on; hold on
text(0.02,0.1,'a','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

transport_by_layer_hg([1 end])=NaN;
transport_by_layer_hg_geos3400([1 end])=NaN;
transport_by_layer_hg_geos281a([1 end])=NaN;
transport_by_layer_hg_geos281([1 end])=NaN;
transport_by_layer_hg_geos281imb([1 end])=NaN;
transport_by_layer_bottom_track([1 end])=NaN;
% plot(transport_by_layer_hgb,1:length(transport_by_layer_hg),'k-');hold on;

plot(transport_by_layer_hg,1:length(transport_by_layer_hg),'k+-');hold on;
plot(transport_by_layer_hg_geos3400,1:length(transport_by_layer_hg),'b+-');hold on;
plot(transport_by_layer_hg_geos281a,1:length(transport_by_layer_hg),'o-','color',[1 0.7 0]);hold on;
plot(transport_by_layer_hg_geos281,1:length(transport_by_layer_hg),'r+-');hold on;

plot(transport_by_layer_hg_geos281imb,1:length(transport_by_layer_hg),'r+--');hold on;
plot(transport_by_layer_bottom_track,1:length(transport_by_layer_hg),'+--','color',[0 0.4 0]);hold on;

% plot(transport_by_layer_hg_geos281,1:length(transport_by_layer_hg),'r-');hold on;
%plot(transport_by_layer_hg_geos281ladcp,1:length(transport_by_layer_hg),'-','color',[1.0000 0.8398 0]);hold on;

legend('ADCPs','3400dbar','\gamma28.1 kg.m^-^3 Argo','\gamma28.1 kg.m^-^3','\gamma28.1 kg.m^-^3 imb.','Bottom Track','location','northwest','AutoUpdate', 'off','box','off')


% plot(transport_by_layer_hg_geos3400,1:length(transport_by_layer_hg),'b--');hold on;
% plot(transport_by_layer_hg_geos3400argo,1:length(transport_by_layer_hg),'--','color',[0 .5 0]);hold on;

% plot(transport_by_layer_hg_geos281a,1:length(transport_by_layer_hg),'y-','linewidth',2);hold on;

% plot(transport_by_layer_hg_geos281bar,1:length(transport_by_layer_hg_geos281bar));hold on;
% plot(transport_by_layer_hg_geos3400bar,1:length(transport_by_layer_hg_geos3400bar));hold on;

axis ij;xline(0);

%yticks('manual')
yticklabels({'SW','SACW','AAIW','UCDW','NADW','LCDW','AABW',''})
grid on
ylim([1.5 8.5])


%ylim([1.7 8.2])

% transport_msm_by_layer=transport_by_layer_hg_geos281b;
% transport_msm_by_laye_imb=transport_by_layer_hg_geos281;

xlim([-12 12])
% title('transport by layer (Sv)')

% %legend('ladcp','28.1','3400','281bar','3400bar')
xlabel('Transport (Sv)')
% xticks([-10:2:10])


grid on
h = subplot(2,2,2); 
p = get(h, 'pos');
p(3) = p(3) + 0.06; p(1) = p(1) +0.01;p(4) = p(4) + 0.06;p(2) = p(2) - 0.02;set(h, 'pos', p);box on; hold on


text(0.02,0.1,'b','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

%load transport_west_east

transport_by_layer_hg_geos281([1 end])=NaN;
transport_by_layer_hg_geos281_west([1 end])=NaN;
transport_by_layer_hg_geos281_east([1 end])=NaN;
transport_by_layer_hg_geos281_east([1 end-1])=NaN;


plot(transport_by_layer_hg_geos281,1:length(transport_by_layer_hg),'+k-');hold on;

plot(transport_by_layer_hg_geos281_west,1:length(transport_by_layer_hg),'or-');hold on;

plot(transport_by_layer_hg_geos281_east,1:length(transport_by_layer_hg),'ob-');hold on;


axis ij;xline(0);

%yticks('manual')
yticklabels({'SW','SACW','AAIW','UCDW','NADW','LCDW','AABW',''})
grid on
%

ylim([1.5 8.5])

xlim([-13.5 13.5])
% xticks([-12:2:12])

% title('transport by layer (Sv)')
% 
legend('Total','West','East','location','southeast','AutoUpdate', 'off','box','off')
xlabel('Transport (Sv)')



% 
% %plot(transportv,1:length(transportv),'color',[.5 .5 .5]);hold on;
% plot(transportvv,1:length(transportv),'k');hold on;
% plot(transportvgos281_ek,1:length(transportvgos281_ek),'--r');hold on;
% plot(transportvgos3400,1:length(transportvgos3400),'--b');hold on;
% %plot(nsum(transportvgos3400argo,1:length(transportvgos3400),'--g');hold on;
% plot(nsum(transportvgos3400argo,2),pres,'--g');hold on
% 
% %plot(transportv_balanced,1:length(transportv),'color',[.5 .5 .5],'linestyle','--');hold on;
% plot(transportvv_balanced,1:length(transportv),'--k');hold on;
% % plot(transportvv,1:length(transportv),'k','linewidth',1);hold on;
% plot(transportvgos281_balanced,1:length(transportvgos281_ek),'r','linewidth',1);hold on;
% plot(transportvgos3400_balanced,1:length(transportvgos3400),'b','linewidth',1);hold on;
% % plot(transportvgos3400_balanced,1:length(transportvgos3400),'b','linewidth',1);hold on;
% plot(nsum(transportvgos3400argo_balanced,2),pres,'g','linewidth',1);hold on
% 
% % plot(transportvgos281_eka,pres,'y','linewidth',2);hold on
% 
% axis ij;xline(0);
% %title('transport per 1db (Sv)')
% ylim([0 5500]);xlim([-0.05 0.05]);ylabel('Pressure (dbar)');xlabel('Transport (Sv)')
% xticks([-0.04:0.02:0.04])
%  savefig('figure_3a')
%  print(gcf, '-dpdf','-r600','/Users/gaston/Documents/figure_3a')
% 
% legend('LADCP','\gamma28.1','3400db','3400db ARGO','LADCP adj','\gamma28.1 adj','3400db adj','3400db ARGOadj','location','southwest','AutoUpdate', 'off','box','off')





grid on



h = subplot(2,2,3); 
p = get(h, 'pos');
p(3) = p(3) + 0.06; p(1) = p(1) - 0.05;p(4) = p(4) + 0.05;p(2) = p(2) - 0.02;set(h, 'pos', p);box on; hold on
text(0.02,0.1,'c','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

% % plot(tt,1:length(g_n_levels),'k+--');hold on;
% % plot(tts,1:length(g_n_levels),'r+--');hold on;
% % plot(tthg,1:length(g_n_levels),'b+--');hold on;
% % plot(tthga,1:length(g_n_levels),'g+--');hold on;
% plot(ttb,1:length(ttb),'k+-');hold on;
% plot(ttsb,1:length(ttb),'r+-');hold on;
% plot(tthgb,1:length(ttb),'b+-');hold on;
% plot(tthgab,1:length(ttb),'g+-');hold on;
% 
% 
% xlabel('Cumulative transport (Sv)');
% 

min_max=cat(1,transport_by_layer_hg_geos281,transport_by_layer_hg_geos3400...
,transport_by_layer_hg_geos281a,transport_by_layer_bottom_track);


minn=min(min_max);
maxx=max(min_max);

meann=mean(min_max);

meanns=mean([minn ;maxx]);

x = 1:9;

barh(x,transport_by_layer_hg_geos281,'facecolor',[.6 .6 .6])                
hold on

er=errorbar(meanns,x,((maxx-minn)/2),'horizontal')

% bar(x,transport_by_layer_hg_geos281b)                
% er = errorbar(transport_by_layer_hg_geos281b,transport_by_layer_hg_geos281b-minn,transport_by_layer_hg_geos281b+maxx,'.','horizontal');

er.Color = [0 0 0];                            
er.LineStyle = 'none';  

% flip axis
% hold off
% barh(transport_by_layer_hg_geos281b)
ylim([1 9])
 axis ij;xline(0);
% 
yticklabels({'','SW','SACW','AAIW','UCDW','NADW','LCDW','AABW',''})
grid on

% xlim([-10.9 10.9])

xlim([-11.6 11.6])

axis ij

xticks([-10:2:10])
xlabel('Transport (Sv)')



h = subplot(2,2,4); 
p = get(h, 'pos');
%p(3) = p(3) + 0.06; p(1) = p(1) - 0.05;p(4) = p(4) + 0.02;p(2) = p(2) + 0.03;set(h, 'pos', p);box on; hold on
p(3) = p(3) + 0.06; p(1) = p(1) +0.01;p(4) = p(4) + 0.05;p(2) = p(2) -0.02;set(h, 'pos', p);box on; hold on


text(0.02,0.1,'d','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

tt=nancumsum(transport_by_layer_hg,'omitnan');
tts=nancumsum(transport_by_layer_hg_geos281);
tthg=nancumsum(transport_by_layer_hg_geos3400);

% tthh=nancumsum(transport_by_layer_hg_geos_);

%plot(transportv_balanced,1:length(transportv),'color',[.5 .5 .5],'linestyle','--');hold on;
% plot(nancumsum(transportvv_balanced),1:length(transportv),'--k');hold on;
% plot(transportvv,1:length(transportv),'k','linewidth',1);hold on;
% plot(nancumsum(transportvgos3400_balanced),1:length(transportvgos3400),'b','linewidth',1);hold on;
% plot(transportvgos3400_balanced,1:length(transportvgos3400),'b','linewidth',1);hold on;

% plot(nancumsum(nsum(transportvgos28_eka,2)),1:length(transportvgos3400),'g','linewidth',1);hold on
% 
% plot(nancumsum(nsum(transportvgos3400argo_balanced,2)),1:length(transportvgos3400),'g','linewidth',1);hold on
% 
 aa=nancumsum(nsum(transportvgos28ek_balanced,2));
 ab=nancumsum(nsum(transportvgos3400_balanced,2));
 ac=nancumsum(nsum(transport_bottom_track_balanced,2));
 ad=nancumsum(nsum(transportvgos28ek_balanced_argo,2));

 aladcp=nancumsum(nsum(transportv,2));

% hold on
% 
% fill(aa,ab,1:length(transportvgos3400))
% 
% fill(aa,ab,'k')
% 
% 
% inBetween = [aa, fliplr(ab)];
% fill(x2, inBetween, 'g');


x = 1:length(transportvgos3400);
curve1 = aa;
curve2 = ab;
curve3 = ac;

x2 = [x, fliplr(x)];

x2 = [x; x];

inBetween = [curve1, fliplr(curve2)];
inBetween2 = [curve1, fliplr(curve3)];
 hold on
%  aa=nancumsum(nsum(transportvgos28ek_balanced,2));
%  ab=nancumsum(nsum(transportvgos3400_balanced,2));
%  ac=nancumsum(nsum(transport_bottom_track_balanced,2));
%  ad=nancumsum(nsum(transportvgos28ek_balanced_argo,2));
plot(aladcp,1:length(transportvgos3400),'k');hold on;

plot(aa,1:length(transportvgos3400),'r');hold on;
plot(ab,1:length(transportvgos3400),'b');hold on;
plot(ac,1:length(transportvgos3400),'color',[0 0.4 0]);hold on;
plot(ad,1:length(transportvgos3400),'color',[1 0.7 0]);hold on;

legend('ADCPs','\gamma28.1 kg.m^-^3','3400dbar','Bottom Track','\gamma28.1 kg.m^-^3  Argo','location','southeast','AutoUpdate', 'off','box','off')

% 
% fill(inBetween',x2,[.7 .7 .7],'edgecolor',[.6 .6 .6]);hold on
% fill(inBetween2',x2,[.7 .7 .7],'edgecolor',[.6 .6 .6]);

plot(aa,1:length(transportvgos3400),'r');hold on;
plot(ab,1:length(transportvgos3400),'b');hold on;
plot(ac,1:length(transportvgos3400),'color',[0 0.4 0]);hold on;
plot(ad,1:length(transportvgos3400),'color',[1 0.7 0]);hold on;
% legend('LADCP adj','γ=28.1 adj','3400db adj','3400db ARGOadj','location','south','AutoUpdate', 'off')
axis ij;xline(0);
%title('transport per 1db (Sv)')
% hold on
%  plot(nancumsum(nsum(transportvgos281_eka,2)),1:length(transportvgos3400),'k','linewidth',1);hold on

ylim([0 5500]);
xlim([-5 18]);
ylabel('Pressure (dbar)');xlabel('Cumulative transport (Sv)')

grid on



set(gcf,'color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',13)
set(findall(gcf,'-property','Linewidth'),'Linewidth',.8)

 print(gcf, '-depsc','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_november/figure7_transport_vertical')

print(gcf, '-djpeg','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_november/figure7_transport_vertical')

print(gcf, '-dpng','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_november/figure7_transport_vertical')

savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_november/figure7_transport_vertical')


%% ARGO AND HEAT TRANSPORT DIFFERENCES PART


% nsum(nsum(heat_transport_balanced28a(1:1700,:)))
% 
% nsum(nsum(heat_transport_balanced28(1:1700,:)))
% 
% 
% nsum(nsum(heat_transport_balanced28a(1700:end,:)))
% 
% nsum(nsum(heat_transport_balanced28(1700:end,:)))
% 
% 
% nmean(nmean(geostrophic_velocity_281_ek_balanced_argo(1:1700,:)))
% 
% nmean(nmean(geostrophic_velocity_281_ek_balanced(1:1700,:)))
% 
% 
% nmean(nmean(geostrophic_velocity_281_ek_balanced_argo(1:1700,:)))-nmean(nmean(geostrophic_velocity_281_ek_balanced(1:1700,:)))
% 
% nmean(nmean(ptemp_comp(1:1700,:)))



%% EDDY TRANSPORT C1
% 
% eddy_transport_vgos=transportvgos28ek_balanced(1:1500,1290:1335);
% eddy_transport_v=transportv(1:1500,1290:1335);
% 
% % eddy_center 13.83 lon(1318)
% 
% eddy_transport_vgos=transportvgos28ek_balanced(1:1000,1290:1335);
% 
% eddy_transport_vgos_argo=transportvgos28ek_balanced_argo(1:1000,1290:1335);
% 
% 
% eddy_transport_v=transportv(1:1000,1290:1335);
% 
% 
% eddy_transport_sined=transportvgos28ek_balanced_sined(1:1000,1290:1335);
% 
% 
% nmean(nmean(CT(1:1000,1290:1335)))
% 
% nmean(nmean(CT_sined(1:1000,1290:1335)))
% 
% nmean(nmean(ptemp_ar(1:1000,1290:1335)))
% 
% anom=CT(1:1000,1290:1335)-ptemp_ar(1:1000,1290:1335);
% 
% 
% nmean(nmean(sqrt(v(1:1000,1290:1335).^2)))
% 
% nmean(nmean(sqrt(geostrophic_velocity_281_ek_balanced(1:1000,1290:1335).^2)))
% nmean(nmean((geostrophic_velocity_281_ek_balanced(1:1000,1290:1335))))
% 
% nmean(nmean(geostrophic_velocity_281_ek_balanced_sined(1:1000,1290:1335)))
% 
% nmean(nmean(sqrt(geostrophic_velocity_281_ek_balanced_argo(1:1000,1290:1335).^2)))
% 
% 
% nmean(nmean(sqrt(geostrophic_velocity_281_ek_balanced_argo(1:1000,1290:1335).^2)))
% nmean(nmean((geostrophic_velocity_281_ek_balanced_argo(1:1000,1290:1335))))
% 
% 
% nmean(nmean(sqrt(geostrophic_velocity_281_ek_balanced_sined(1:1000,1290:1335).^2)))
% 
% 
% brasil_transport=nsum(nsum(transportvgos28ek_balanced(:,1:50)))
% 
% brasil_transportv=nsum(nsum(transportv(:,1:50)))
% 
% 
% 
% nmean(nmean(u(1:1000,1290:1335)))
% 
% 
% nmean(nmean(v(1:1000,1290:1335)))
% 
% 
% aa=transportv(1:1000,1290:1335).*anom(1:1000,:);
% nsum(nsum((aa*1027.5*3850*1e6)))
% 
% mht_v= ((transportv.*CT.*1027.5*3850*1e6));
% mht_v_eddy=nsum(nsum(mht_v(1:1000,1290:1335)))
% 
% fwt_eddy=nsum(nsum(fresshwater_transport_28(1:1000,1290:1335)))
% 
% fwt_eddy_argo=nsum(nsum(fresshwater_transport_28a(1:1000,1290:1335)))
% 
% fwt_eddy=nsum(nsum(ss(1:1000,1290:1335).*transportvgos28ek_balanced(1:1000,1290:1335)))
% 
% fwt_eddy_a=nsum(nsum(ssa(1:1000,1290:1335).*transportvgos28ek_balanced_argo(1:1000,1290:1335)))
% 
% fwt_eddyva=nsum(nsum(ssa(1:1000,1290:1335).*transportvgos28ek_balanced(1:1000,1290:1335)))
% 
% fwt_eddyva=nsum(nsum(ssa(1:1000,1290:1335).*transportvgos28ek_balanced(1:1000,1290:1335)))
% 
% fwt_eddyav=nsum(nsum(ss(1:1000,1290:1335).*transportvgos28ek_balanced_argo(1:1000,1290:1335)))
% 
% 
% fwt_eddy=nsum(nsum(fresshwater_transport_28(1:1000,1290:1335)))
% 
% mht_v=nsum(nsum(transportv(1:1000,1290:1335)))
% 
% 
% 
% 
% nsum(nsum((transportvgos28ek_balanced.*CT(:,1:end - 1))*1027.5*3850*1e6));
% 
% mht_vgos_eddy=mht_vgos(1:1000,1318-30:1318+30)
% 
% 
% 
% (210634*1000*0.5)/1e6
% 
% 
% 
% mht_v= (transportv.*CT.*1027.5*3850*1e6);
% mht_v_eddy=nsum(nsum(mht_v(1:1000,1290:1335)))
% 
% 
% mht_v= (transportv.*CT.*1027.5*3850*1e6);
% mht_v_eddy=nsum(nsum(mht_v(1:1000,1290:1335)))
% 
% 
% 
% mht_v= (transportvgos28ek_balanced.*CT(:,1:end-1).*1027.5*3850*1e6);
% 
% mht_v_eddy=nsum(nsum(mht_v(1:1000,1290:1335)))
% 
% 
% mht_v= (transportvgos28ek_balanced_argo.*ptemp_comp(:,1:end-1).*1027.5*3850*1e6);
% 
% mht_v_eddy=nsum(nsum(mht_v(1:1000,1290:1335)))
% 
% 
% eddy_v=nmean(nmean(v(1:1000,1290:1335)))
% 
% 
% eddy_u=nmean(nmean(u(1:1000,1290:1335)))
% 
% 
% load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/eddies_merian_contours_traj_and_argo.mat', 'cy_traj_msm60')
% % 2017           1           6   C1 was crossed 123 in traje
% 
% 
% in_lat=cy_traj_msm60{1,4}(123);
% in_lon=cy_traj_msm60{1,3}(123);
% 
% endd=10;
% end_lat=cy_traj_msm60{1,4}(123-endd)
% end_lon=cy_traj_msm60{1,3}(123-endd)
% 
% 
% 
% in_lat=cy_traj_msm60{1,6}(123);
% in_lon=cy_traj_msm60{1,5}(123);
% 
% endd=10;
% end_lat=cy_traj_msm60{1,6}(123-endd)
% end_lon=cy_traj_msm60{1,5}(123-endd)
% 
% 
% [dist, angle]=sw_dist([in_lat end_lat],[in_lon end_lon],'km')
% 
% speed=dist/endd
% spd_ms=speed/(60*60*24)*1000
% 
% anglee=270+(-angle)
% 
% [uu,vv]=uvposta(spd_ms,anglee)
% 
% 
% 
% %heading 281 at 10 km per day
% 
