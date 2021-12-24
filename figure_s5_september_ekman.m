% load msm60_dvs.mat
% 
% 
% % wind_spd=wind_spd/(3^0.1);
% 
% 
% [u,v]=uvposta(wind_spd,wind_dir);
% 
% 
% [taux,tauy]=windstress(u,v);
% 
% nanmean(taux)
% 
% %coriolis and density
% f=2*2*pi/(24*3600)*sind(-34.5);ro=1030;
% 
% %
% %averaged by longitude
% taux_lon=nanmean(taux,2);
% tauy_lon=nanmean(tauy,2);
% 
% %averaged by time
% taux_t=squeeze(nanmean(taux,1));taux_t=squeeze(nmean(taux_t,1));
% tauy_t=nanmean(tauy,1);
% 
% 
% distance_transect=sw_dist([longitude(1) longitude(end)],[-34.5 -34.5],'km')*1000;
% 
% % for i=1:length(lon_tau)-1
% % distance_transects(i)=sw_dist([lon_tau(i) lon_tau(i+1)],[-34.5 -34.5],'km')*1000;
% % end
% distance_transects=(ones(1,length(taux_lon))+27780); %checlk if this should be n-1
% 
% 
% 
% ekman_transport_ship_wind=((-taux./(ro*f)).*27781)./1e6;
% ekman_transport_temporal_ship_wind=((-taux_t./(ro*f)).*distance_transect)./1e6;


cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/ekman_merian
ncdisp('ekman_temporal_series_merian.nc')

taux=ncread('ekman_temporal_series_merian.nc','surface_downward_eastward_stress');

tauy=ncread('ekman_temporal_series_merian.nc','surface_downward_northward_stress');
lon_tau=ncread('ekman_temporal_series_merian.nc','lon');
lat_tau=ncread('ekman_temporal_series_merian.nc','lat');

tim=ncread('ekman_temporal_series_merian.nc','time');


% timm=tim (1)-datenum(1992,1,1,6,0,0)
% 
% timm (1)-datenum(1992,1,1,6,0,0)

timm=tim/24+693962;timm=double(timm);
tim=datevec(timm);







% 

taux345=squeeze(taux(:,3,:));
tauy345=squeeze(tauy(:,3,:));


% t1 = datetime(1992,1,1,6,0,0);
% t2 = datetime(2018,12,31,18,0,0);
% t3 = t1:t2;
% t3 = t1:hours(6):t2;
% time=datevec(t3);
% t=datenum(t3);

% time= 
%coriolis and density
f=2*2*pi/(24*3600)*sind(-34.5);ro=1030;
%
tau_lon=nanmean(taux345,2);
tauy_lon=nanmean(tauy345,2);
distance_transect=sw_dist([lon_tau(1) lon_tau(end)],[-34.5 -34.5],'km')*1000;

%figure;plot(tau_lon);

tau = squeeze(nanmean(squeeze(nanmean(taux345, 1)), 2));

ekman_transport=((-taux345./(ro*f)).*distance_transect)./1e6;


ekman_transport_temporal=nanmean(ekman_transport,1);

rr=4;
ekt=movmean(ekman_transport_temporal,rr);ekt=ekt(1:rr:end);

rrm=4*30;
ektm=movmean(ekman_transport_temporal,rrm);ektm=ektm(1:rr:end);

figure;plot(ekt);hold on; plot(ektm,'k');axis tight;

% figure;plot(ektm,'k');axis tight;


tx=movmean(tau_lon,rr);tx=tx(1:rr:end);
ty=movmean(tauy_lon,rr);ty=ty(1:rr:end);

figure;
quiver(lon_tau(1:rr:end),lon_tau(1:rr:end)*0,tx,ty,'k'); hold on
 yyaxis right; hold on
 plot(lon_tau(1:rr:end),tx,'r');hold on;xlabel('Longitude')
 plot(lon_tau(1:rr:end),ty,'b');ylabel('taux (r) and tauy (b)');hold on
% quiver(lon_tau(1:rr:end),lon_tau(1:rr:end)*0,tx,ty,'k'); hold on


legend('wind vectors', 'taux', 'tauy','location','northwest')
title(' IFREMER CERSAT Global Blended Mean Wind Fields')

ekman=nmean(ekman_transport_temporal(35076:35187))


% t1 = datetime(2007,5,1);
% t2 = datetime(2018,11,1);
% t3 = t1:t2;
% t3 = t1:months(1):t2;
% time=datevec(t3);
% t=datenum(t3);



% 
% % dates=datevec(tim);
% 
% ekman_transport=((-taux345./(ro*f)).*distance_transect)./1e6;
% 
for i=1:12;
ind=find(tim(:,2)==i);
ekman_monthly(:,i)=nanmean(ekman_transport_temporal(:,ind),2);
end

% 
% % ekman_transport_temporal=((-taux_t./(ro*f)).*distance_transect)./1e6;;
% % 
% % ekman_transport_temporal=((-taux345./(ro*f)).*distance_transect)./1e6;;
% 
% % figure;plot(time,ekman_transport_temporal,'k');axis tight;datetick('x')




ekm=nanmean(ekman_monthly,1);ekm=repmat(ekm,1,27);%ekm=ekm(5:138+4);


rrm=4*30;

ektm=movmean(ekman_transport_temporal(24081:end),rrm);

ektm=ektm(1:rrm:end);

% tt=timm(24081:rrm:end)


y=tim(:,1); m=tim(:,2);
k=0;


for j = 1992:2018
  for i= 1:12
      k=k+1;
      idx=find(y==j&m==i);
      ekman_mensual(k)=nmean(ekman_transport_temporal(idx));
  end
end


%23123 arranca 2008

% ektmm=ektm(end-length(time)+1:end);



%anom=ekman_transport_temporal-ekm;
dn=datenum(1992,[1:((2018-1992)+1)*12].',1);
dates=datevec(dn);dates=dates(:,1:3);%dates=dates(5:end-2,:);
time=datenum(dates);



h =figure('Renderer', 'painters', 'Position', [200 200 550 200])

p = get(h, 'pos');
p(3) = p(3) + 0.4; p(1) = p(1) - 0.0;p(4) = p(4) + 0.075;p(2) = p(2) + 0.03;set(h, 'pos', p);
box on; hold on

%subplot(2,1,1)

plot(time(204:end),ekm(204:end),':k');

hold on;
ylabel('Meridional Ekman Transport (Sv)')
plot(time(204:end)',ekman_mensual(204:end),'k');hold on

plot(time(end-23),ekman_mensual(end-23),'ob');
plot(time(end-23),-1.8,'or');
legend( 'Climatology','Times series','Satelite','Ship wind','Orientation','horizontal','location','northwest','AutoUpdate', 'off','orientation','horizontal','box','off')

xlim([733743 737395])

xticks([time(204):30:time(end)])

datetick('x','mm/yy');
axis tight;


yline(0,'color', [.5 .5 .5]);grid on

grid on
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(findall(gcf,'-property','Linewidth'),'Linewidth',.7)
set(gcf,'color','w'); 
ylim([-2.8 12.8]);yticks([-2:2:12])






savefig ('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september/figure_s5_september_ekman')

print(gcf, '-depsc','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september/figure_s5_september_ekman')
print(gcf, '-djpeg','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/figure_s5_september_ekman')

%subplot(2,1,2)

% 
% figure
% shade_anomaly(time,anom');%legend('Anomalies (Sv)','jan/17')
% hold on; plot(time(117),anom(117),'or')
% 
% 
% datetick('x','mm/yy');axis tight;%ylabel('Anomalies (Sv)')
% 
% 
% figure
% imagesc(lon_tau,time,ekman_transport');ax.CLim = [-15 15];colormap(rednblue);colorbar
% datetick('y', 'mm/yy');xlabel('Longitude ∫'); ylabel('Time'); title('Monthly Ekman transport (Sv) at 34.5∫S every 0.25∫') ; axis tight
% 
% 
% 
% figure;
% %subplot(2,1,1)
% plot(time,ekman_transport_temporal,'k');hold on;
% ylabel('Ekman Transport (Sv)')
% plot(time,ekm,':k')
% plot(time(117),ekman_transport_temporal(117),'or');legend('Times series', 'Climatology','Jan/17(-.24)','Orientation','horizontal')
% axis tight;datetick('x','mm/yy');axis tight;
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %coriolis and density
% f=2*2*pi/(24*3600)*sind(-34.5);ro=1030;
% %
% %averaged by longitude
% taux_lon=nanmean(taux,2);
% tauy_lon=nanmean(tauy,2);
% 
% %averaged by time
% taux_t=nanmean(taux,1);
% tauy_t=nanmean(tauy,1);
% 
% 
% distance_transect=sw_dist([lon_tau(1) lon_tau(end)],[-34.5 -34.5],'km')*1000;
% 
% % for i=1:length(lon_tau)-1
% % distance_transects(i)=sw_dist([lon_tau(i) lon_tau(i+1)],[-34.5 -34.5],'km')*1000;
% % end
% distance_transects=(ones(1,length(lon_tau))+27780); %checlk if this should be n-1
% 
% 
% 
% ekman_transport=((-taux./(ro*f)).*27781)./1e6;
% ekman_transport_temporal=((-taux_t./(ro*f)).*distance_transect)./1e6;;
% 
% 
% 
% for i=1:12;
% ind=find(dates(:,2)==i);
% ekman_monthly(:,i)=nanmean(ekman_transport(:,ind),2)
% end
% 
% 
% % for i=1:12
% %     ekman_seasonal=
% % end
% % y=dates(:,1); m=dates(:,2);
% % k=0;
% 
% 
% 
% 
% 
% 
% % ekman_transport_temporal=((-taux_t./(ro*f)).*distance_transect)./1e6;
% % 
% % figure;plot(time,ekman_transport_temporal,'k');axis tight;datetick('x')
% 
% 
% ekm=nansum(ekman_monthly,1);ekm=repmat(ekm,1,12);ekm=ekm(5:138+4);
% 
% anom=ekman_transport_temporal-ekm;
% 
% 
% dn=datenum(2007,[1:((2018-2007)+1)*12].',1);
% dates=datevec(dn);dates=dates(:,1:3);dates=dates(5:end-2,:);
% 
% time=datenum(dates);
% 
% 
% figure;
% %subplot(2,1,1)
% plot(time,ekman_transport_temporal,'k');hold on;
% ylabel('Ekman Transport (Sv)')
% plot(time,ekm,':k')
% plot(time(117),ekman_transport_temporal(117),'ob');
% plot(time(117),-1.8,'or');
% legend('Times series', 'Climatology','Satelite','Ship Wind','Orientation','horizontal')
% axis tight;datetick('x','mm/yy');axis tight;
% 
% 
% 
% %subplot(2,1,2)
% figure
% shade_anomaly(time,anom');%legend('Anomalies (Sv)','jan/17')
% hold on; plot(time(117),anom(117),'or')
% datetick('x','mm/yy');axis tight;ylabel('Anomalies (Sv)')
% 
% 
% figure
% imagesc(lon_tau,time,ekman_transport');ax.CLim = [-15 15];colormap(rednblue);colorbar
% datetick('y', 'mm/yy');xlabel('Longitude ∫'); ylabel('Time'); title('Monthly Ekman transport (Sv) at 34.5∫S every 0.25∫') ; axis tight
% 
% 
% 
% load('C:\Users\gaston\Desktop\Phd\Datos_CTD_phd\Johannes_Merian\msm60_dvs.mat')
% [u_insitu,v_insitu]=uvposta(wind_spd,wind_dir);
% 
% [Tx,Ty]=windstress(u_insitu,v_insitu);
% 
% 
% t1 = datetime(2007,5,1);
% t2 = datetime(2018,11,1);
% t3 = t1:t2;
% t3 = t1:months(1):t2;
% time=datevec(t3);
% t=datenum(t3);
% 
% dn=datenum(2007,[1:((2018-2007)+1)*12].',1);
% dates=datevec(dn);dates=dates(:,1:3);dates=dates(5:end-2,:);
% time=datenum(dates);
% 
% 
% 
% % t1 = datetime(2007,5,1);
% % t2 = datetime(2018,11,1);
% % t3 = t1:t2;
% % t3 = t1:months(1):t2;
% % time=datevec(t3);
% % t=datenum(t3);
% 
% dn=datenum(2007,[1:((2018-2007)+1)*12].',1);
% dates=datevec(dn);dates=dates(:,1:3);dates=dates(5:end-2,:);
% time=datenum(dates);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% load msm60_dvs.mat
% 
% %Hellmann height correction with alfa=0.12 wind measured at 30m
% wind_spd=wind_spd/(3^0.1);
% 
% [u,v]=uvposta(wind_spd,wind_dir);
% [taux,tauy]=windstress(u,v);
% 
% nanmean(taux)
% 
% %coriolis and density
% f=2*2*pi/(24*3600)*sind(-34.5);ro=1030;
% 
% %
% %averaged by longitude
% taux_lon=nanmean(taux,2);
% tauy_lon=nanmean(tauy,2);
% 
% %averaged by time
% taux_t=nanmean(taux,1);
% tauy_t=nanmean(tauy,1);
% 
% 
% distance_transect=sw_dist([longitude(1) longitude(end)],[-34.5 -34.5],'km')*1000;
% 
% % for i=1:length(lon_tau)-1
% % distance_transects(i)=sw_dist([lon_tau(i) lon_tau(i+1)],[-34.5 -34.5],'km')*1000;
% % end
% distance_transects=(ones(1,length(lon_tau))+27780); %checlk if this should be n-1
% 
% 
% 
% ekman_transport=((-taux./(ro*f)).*27781)./1e6;
% ekman_transport_temporal=((-taux_t./(ro*f)).*distance_transect)./1e6;
% 
% 
% 
% for i=1:12;
% ind=find(dates(:,2)==i);
% ekman_monthly(:,i)=nanmean(ekman_transport(:,ind),2)
% end
% 
% 
% ekman_transport_temporal=((-taux_t./(ro*f)).*distance_transect)./1e6;;
% figure;plot(time,ekman_transport_temporal,'k');axis tight;datetick('x')
% 
% ekm=nansum(ekman_monthly,1);ekm=repmat(ekm,1,12);ekm=ekm(5:138+4);
% anom=ekman_transport_temporal-ekm;
% figure;
% %subplot(2,1,1)
% plot(time,ekman_transport_temporal,'k');hold on;
% ylabel('Ekman Transport (Sv)')
% plot(time,ekm,':k')
% plot(time(117),ekman_transport_temporal(117),'or');legend('Times series', 'Climatology','Jan/17(-.24)','Orientation','horizontal')
% axis tight;datetick('x','mm/yy');axis tight;
% %subplot(2,1,2)
% figure
% shade_anomaly(time,anom');%legend('Anomalies (Sv)','jan/17')
% hold on; plot(time(117),anom(117),'or')
% datetick('x','mm/yy');axis tight;ylabel('Anomalies (Sv)')
% 
% 
% figure
% imagesc(lon_tau,time,ekman_transport');ax.CLim = [-15 15];colormap(rednblue);colorbar
% datetick('y', 'mm/yy');xlabel('Longitude ∫'); ylabel('Time'); title('Monthly Ekman transport (Sv) at 34.5∫S every 0.25∫') ; axis tight
% 
% 
% 