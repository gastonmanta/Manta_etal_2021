close all;clear all;

addpath('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/4gaston/obana')
datedir = dir('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/msm60ladcp_lad');
filenames = {datedir.name};files=filenames(3:end);%files=cell2str(files)
files=files';files(end-1)=[]; 

stastr='[2:5 7:29 31:126 128]';%why do you take out the 127?
sta_l=str2num(stastr);


%% lADCP %%
% path('/Users/gaston/Documents/MATLAB/Johannes_Merian/');
% stastr='[2:5 7:29 31:126 128]';
% sta_l=str2num(stastr);
cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/msm60ladcp_lad
[lat_l,lon_l,z_l,u_l,v_l]=rodbload(...
files,['Latitude:Longitude:z:u:v']);
cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/
%todos los valores de ladcp
[lon_l i]=sort(lon_l);% sorts by longitude from w to e
sta_l=sta_l(i);
lat_l=lat_l(i);
z_l=z_l(:,i);
u_l=u_l(:,i).*100;% asume that takes from m to cm per s
v_l=v_l(:,i).*100;

i=find(z_l<40); %eliminates the first 30 m, why not just the first 10 or 20? 
u_l(i)=NaN;
v_l(i)=NaN;


%535(10 meter intervals)x124 stations and first 30 m removed)

%% vmadcp %%
 
load /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/4gaston/vmadcp/msm60os38000_000000_32_hc.mat

v=squeeze(b.vel(:,2,:));v=squeeze(b.vel(:,2,:))*100;
u=squeeze(b.vel(:,1,:));u=squeeze(b.vel(:,1,:))*100;

zc=c.depth(:);%32 meter interval from 53m to 1621
pg=c.pg;%what is pg and why is erased at less than 25? ping? and something with land?
kk=find(pg<=25);
v(kk)=NaN;
u(kk)=NaN;


v_adcp1=[];u_adcp1=[];
for i=1:length(v);
  v_adcp1=[v_adcp1 interp1(zc,v(:,i),[10:10:1600]')];%interpolates in 10, interval from 10 to 1600
  u_adcp1=[u_adcp1 interp1(zc,u(:,i),[10:10:1600]')];
end
zc_adcp1=[10:10:1600]';

xyt_adcp1=[b.nav.txy1(2,:); b.nav.txy1(3,:); b.nav.txy1(1,:)];%lon lat time measurements


[lon_adcp1 i]=sort(b.nav.txy1(2,:));% sorts by longitude from w to e
u_adcp1=u_adcp1(:,i);% asume that takes from m to cm per s
v_adcp1=v_adcp1(:,i);

clear u v b c zc

load /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/4gaston/vmadcp/msm60os75000_000000_39_hc.mat

v=squeeze(b.vel(:,2,:))*100;
u=squeeze(b.vel(:,1,:))*100;
zc=c.depth(:);
pg=c.pg;%8 meter interval from 17m to 800
kk=find(pg<=25);
v(kk)=NaN;
u(kk)=NaN;

v_adcp2=[];u_adcp2=[];
for i=1:length(v);
  v_adcp2=[v_adcp2 interp1(zc,v(:,i),[10:10:800]')];
  u_adcp2=[u_adcp2 interp1(zc,u(:,i),[10:10:800]')];
end
zc_adcp2=10:10:800;
xyt_adcp2=[b.nav.txy1(2,:); b.nav.txy1(3,:); b.nav.txy1(1,:)];


[lon_adcp2 i]=sort(b.nav.txy1(2,:));% sorts by longitude from w to e
u_adcp2=u_adcp2(:,i);% asume that takes from m to cm per s
v_adcp2=v_adcp2(:,i);



pos_grid=[-52:0.02:18];
depth_grid=[0:10:1300];


%REMOVES THE NOISY DATA IN THE SHELF
u_adcp1(:,1:59)=NaN;
v_adcp1(:,1:59)=NaN;



%JUST A PLOT OF THE 3 SOURCES OF MESUREMENTS
figure(1);
subplot(3,1,1)
pcolor(lon_l,1:100,v_l(1:100,:));shading interp
title('V ladcp (cm/s)');ylabel('Depth (10m intervals)')
caxis([-60 60]);axis ij;colorbar
cmocean('balance',11);ylim([0 60]);xlim([lon_adcp2(1) lon_adcp2(end)])
subplot(3,1,2)
pcolor(lon_adcp1,1:100,v_adcp1(1:100,:));shading interp
caxis([-60 60]);axis ij;colorbar;ylabel('Depth (10m intervals)')
cmocean('balance',11)
title('V adcp 38 (cm/s)');ylim([0 60]);xlim([lon_adcp2(1) lon_adcp2(end)])
subplot(3,1,3)
pcolor(lon_adcp2,1:80,v_adcp2);shading interp
title('V adcp 75 (cm/s)');ylim([0 60]);xlim([lon_adcp2(1) lon_adcp2(end)])
shading interp;ylabel('Depth (10m intervals)')
cmocean('balance',11)
caxis([-60 60]);axis ij;colorbar






% t1=nanmean(v_adcp1(1:80,:),2);
% t2=nanmean(v_adcp2(1:80,:),2);
% t3=nanmean(v_l(1:80,:),2);
t1=nanmean(v_adcp1(:,:),2);
t2=nanmean(v_adcp2(:,:),2);
t3=nanmean(v_l(:,:),2);

figure(2);plot(t1,1:length(t1));hold on
plot(t2,1:length(t2)); 
plot(t3,1:length(t3)); axis ij
legend('adcp38','adcp75','ladcp')
ylim([0 535])
xline(0,'--')
ylabel('depth (10m intervals)')
xlabel('Mean V (cm/s)');title('First 600m measurents averaged by depth')



lat=[xyt_adcp1(2,:) xyt_adcp2(2,:)];%lat from both sadcp
lon=[xyt_adcp1(1,:) xyt_adcp2(1,:)];%lon from both sadcp
tim=[xyt_adcp1(3,:)+julian(2017,1,1) xyt_adcp2(3,:)+julian(2017,1,1)];%times
u=[u_adcp1 [u_adcp2;NaN.*ones(80,length(u_adcp2))]];%horzcat concatenates both sacdps, to the 2nd 80 nans value are added to each column in order to get them to the same size
v=[v_adcp1 [v_adcp2;NaN.*ones(80,length(v_adcp2))]];
z=[zc_adcp1];%10m interval from 10 to 1600m


%NOW I WILL GRID DATA 

pos_grid=[-52:0.02:18];

% [pos2,zc2]=meshgrid(pos_grid2,depth_grid);

depth_grid=[0:10:500];

u2=griddata(lon_adcp2,zc_adcp2,u_adcp2,pos_grid',zc_adcp2);
v2=griddata(lon_adcp2,zc_adcp2,v_adcp2,pos_grid',zc_adcp2);


u1=griddata(lon_adcp1(1:end-2),zc_adcp1,u_adcp1(:,1:end-2),pos_grid,zc_adcp1);
v1=griddata(lon_adcp1(1:end-2),zc_adcp1,v_adcp1(:,1:end-2),pos_grid,zc_adcp1);


figure
pcolor(pos_grid,1:80,(v1(1:80,:)-v2));shading interp
title('V adcp 38 -V adcp 75 (cm/s)');ylim([0 60]);xlim([lon_adcp2(1) lon_adcp2(end)])
shading interp;ylabel('Depth (10m intervals)')
cmocean('balance',11)
caxis([-15 15]);axis ij;colorbar

% fills the gaps of the ladcp due to differences in depth with the average
% of the 2 near, I''l doi only for the deep ocean
  for ll=2:123
        for kk=1:535
    v_int(kk,ll)= nanmean(v_l(kk,ll-1:ll+1),2);
    u_int(kk,ll)= nanmean(u_l(kk,ll-1:ll+1),2);

       end
  end

%take only the part of the deep ocean
v_intt=horzcat([v_l(:,1:19)  v_int(:,20:end-20) v_l(:,end-20:end)]) ;
u_intt=horzcat([u_l(:,1:19)  u_int(:,20:end-20) u_l(:,end-20:end)]) ;

u_new=u_l;
v_new=v_l;

ind=find(isnan(v_l));v_new(ind)=v_intt(ind);
ind=find(isnan(u_l));u_new(ind)=u_intt(ind);




ulad=griddata(lon_l,10:10:5350,u_new,pos_grid',10:10:5350);
vlad=griddata(lon_l,10:10:5350,v_new,pos_grid',10:10:5350);

%fills adcp2 missing values with adcp1
vv=v1(1:80,:);ind=find(isnan(v2));v2(ind)=vv(ind);
uu=u1(1:80,:);ind=find(isnan(u2));u2(ind)=uu(ind);

u_end= vertcat(u2(1:34,:), u1(35:119,:), ulad(120:end,:));
v_end= vertcat(v2(1:34,:), v1(35:119,:), vlad(120:end,:));

%fills adcp1 misssing with ctd
ind=find(isnan(v_end));v_end(ind)=vlad(ind);
ind=find(isnan(u_end));u_end(ind)=ulad(ind);



load /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/msm60_dvs.mat

kk=find(depth>6000);
depth(kk)=NaN;clear kk

lon=longitude;
for ip=1:length(pos_grid)
    kk=find(lon>=pos_grid(ip)-0.025 & lon<=pos_grid(ip)+0.025);
    topo(ip)=nmean(depth(kk));
end

% lon_topo=pos_grid;

topoo=fillmissing(topo,'linear')


v_endd=movmean(v_end,5,2);
v_endd=movmean(v_endd,2,1);
v_endd=movmean(v_endd,5,2);
v_endd=movmean(v_endd,2,1);


u_endd=movmean(u_end,5,2);
u_endd=movmean(u_endd,2,1);
u_endd=movmean(u_endd,5,2);
u_endd=movmean(u_endd,2,1);



t4=nanmean(v_end(:,:),2);
t5=nanmean(v_endd(:,:),2);


figure;
pcolor(pos_grid,1:535,v_endd);shading interp; axis ij;
cmocean('balance',23);caxis([-55 55]);colorbar
hold on;
contour(pos_grid,1:535,v_endd,[0 0], 'linecolor', [.7 .7 .7])
plot(pos_grid,topoo/10,'k');axis ij;hold on;

v_endw=v_endd-nanmean(nanmean(v_endd));
t6=nanmean(v_endw(:,:),2);



figure;plot(t1,1:length(t1));hold on
plot(t2,1:length(t2)); 
plot(t3,1:length(t3)); 
plot(t4,1:length(t4));axis ij
plot(t5,1:length(t5));hold on
plot(t6,1:length(t6));hold on
legend('adcp38','adcp75','ladcp','v_merged','v_merged_smoothed','v_merged_smoothed_corrected')
ylim([0 535]);xlim([-3 3])
xline(0,'--')
yline(34,'.-');yline(119,'.-')
ylabel('depth (10m intervals)')
xlabel('Mean V (cm/s)');title('First 600m measurents averaged by depth')




% v_endd=fillmissing(v_end,'linear',1,'EndValues','nearest');




v_end_new=griddata(pos_grid,10:10:5350,v_endd,pos_grid',1:5500);

u_end_new=griddata(pos_grid,10:10:5350,u_endd,pos_grid',1:5500);

% %IF NECESARY fills the gaps, not now...
u_end_new=fillmissing(u_end_new,'linear',1,'EndValues','nearest');
v_end_new=fillmissing(v_end_new,'linear',1,'EndValues','nearest');


aa=1:5500;aa=aa';
pres=sw_pres(aa,-35);clear aa;


vv=griddata(pos_grid',pres,v_end_new,pos_grid',0:5500);

uu=griddata(pos_grid',pres,u_end_new,pos_grid',0:5500);

tt=topo'
topo=sw_pres(tt,-35);topo=fillmissing(topo,'linear');

bucket=ones(5501,length(topo));


for i=1:length(topo);
    bucket(topo(i):end,i)=NaN;
end

vv=vv.*bucket;
uu=uu.*bucket;

save v_end_new vv uu pos_grid topo