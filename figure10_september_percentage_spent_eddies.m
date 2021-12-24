% 
% main_dir='/Users/gaston/Desktop/adts_eddies'
% 
% list=dir('adt_*.mat');
% 
% apo=zeros(282,1,length(list));
% 
% dates=zeros(1,length(list))';
% 
% % for i=1:length(list)
% % load(list(i).name)
% % Anti_pres_out=Anti_pres_out./Anti_pres_out;
% % apoo=Anti_pres_out(72:353,82:83);apoo=squeeze(nmean(apoo,2));apoo=apoo./apoo;
% % apo(:,:,i)=apoo;
% % %dates(i)=dates;
% % end
% % 
% % 
% % for i=1:length(list)
% % load(list(i).name);
% % dates(i)=date_num;
% % end
% 
% 
% cpo=zeros(282,1,length(list));
% 
% % for i=1:length(list)
% % load(list(i).name)
% % Cyclo_pres_out=Cyclo_pres_out./Cyclo_pres_out;
% % cpoo=Cyclo_pres_out(72:353,82:83);cpoo=squeeze(nmean(cpoo,2));cpoo=cpoo./cpoo;
% % cpo(:,:,i)=cpoo;
% % end
% % 
% 
% for i=1:length(list)
% load(list(i).name)
% Anti_pres_out=Anti_pres_out./Anti_pres_out;
% apoo=Anti_pres_out(72:353,82:83);apoo=squeeze(nmean(apoo,2));apoo=apoo./apoo;
% apo(:,:,i)=apoo;
% Cyclo_pres_out=Cyclo_pres_out./Cyclo_pres_out;
% cpoo=Cyclo_pres_out(72:353,82:83);cpoo=squeeze(nmean(cpoo,2));cpoo=cpoo./cpoo;
% cpo(:,:,i)=cpoo;
% dates(i)=date_num;
% 
% %dates(i)=dates;
% end
% 

% cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian
% load percentage_eddies apo cpo; 



cd /Users/gaston/Documents/adts_eddies

main_dir='/Users/gaston/Documents/adts_eddies'

list=dir('adt_*.mat');

apo=zeros(290,1,length(list));

dates=zeros(1,length(list))';

for i=1:length(list)
load(list(i).name)
Anti_pres_out=Anti_pres_out./Anti_pres_out;
apoo=Anti_pres_out(68:357,83);apoo=squeeze(nmean(apoo,2));apoo=apoo./apoo;
apo(:,:,i)=apoo;
dates(i)=date_num;
Cyclo_pres_out=Cyclo_pres_out./Cyclo_pres_out;
cpoo=Cyclo_pres_out(68:357,83);cpoo=squeeze(nmean(cpoo,2));cpoo=cpoo./cpoo;
cpo(:,:,i)=cpoo;
end

% 
% for i=1:length(list)
% load(list(i).name);
% dates(i)=date_num;
% end

% 
% cpo=zeros(282,1,length(list));
% 
% for i=1:length(list)
% load(list(i).name)
% Cyclo_pres_out=Cyclo_pres_out./Cyclo_pres_out;
% cpoo=Cyclo_pres_out(72:353,83);cpoo=squeeze(nmean(cpoo,2));cpoo=cpoo./cpoo;
% cpo(:,:,i)=cpoo;
% end
% 
% 
% for i=1:length(list)
% load(list(i).name)
% Anti_pres_out=Anti_pres_out./Anti_pres_out;
% apoo=Anti_pres_out(72:353,82:83);apoo=squeeze(nmean(apoo,2));apoo=apoo./apoo;
% apo(:,:,i)=apoo;
% Cyclo_pres_out=Cyclo_pres_out./Cyclo_pres_out;
% cpoo=Cyclo_pres_out(72:353,82:83);cpoo=squeeze(nmean(cpoo,2));cpoo=cpoo./cpoo;
% cpo(:,:,i)=cpoo;
% dates(i)=date_num;
% 
% %dates(i)=dates;
% end




cpo=squeeze(cpo);apo=squeeze(apo);


percentage_apo=nsum(apo,2)/length(apo)*100;

percentage_cpo=nsum(cpo,2)/length(cpo)*100;

datess=datevec(dates);

ind=find(datess(:,2)>11 | datess(:,2)<3);
summer_dates=datess(ind,:);
apo_summer=apo(:,ind);
cpo_summer=cpo(:,ind);
percentage_apo_summer=nsum(apo_summer,2)/length(apo_summer)*100;
percentage_cpo_summer=nsum(cpo_summer,2)/length(cpo_summer)*100;


ind=find(datess(:,2)>5 & datess(:,2)<9);
winter_dates=datess(ind,:)

apo_winter=apo(:,ind);
cpo_winter=cpo(:,ind);
percentage_apo_winter=nsum(apo_winter,2)/length(apo_winter)*100;
percentage_cpo_winter=nsum(cpo_winter,2)/length(cpo_winter)*100;





figure;
subplot(3,1,1)
plot(X(68:357),percentage_apo,'r')
hold on
plot(X(68:357),percentage_cpo,'b')
axis tight
ylabel({'% of time spent with ';'an eddy at 34.5ºS'})
legend('Anticyclonic eddies (19%)','Cyclonic eddies (14%)','orientation','horizontal','location','best','AutoUpdate', 'off','box','off')
xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
title('All days (n= 9292)','fontweight','normal')
ylim([0 50]);grid on

subplot(3,1,2)

plot(X(68:357),percentage_apo_summer,'r')
hold on
plot(X(68:357),percentage_cpo_summer,'b')
axis tight
ylabel({'% of time spent with ';'an eddy at 34.5ºS'})
% legend('Anticyclonic eddies (24%)','Cyclonic eddies (18%)','orientation','horizontal','location','best','AutoUpdate', 'off','box','off')
xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
title('Summer','fontweight','normal')
ylim([0 50]);grid on

subplot(3,1,3)

plot(X(68:357),percentage_apo_winter,'r')
hold on
plot(X(68:357),percentage_cpo_winter,'b')
axis tight
ylabel({'% of time spent with ';'an eddy at 34.5ºS'})
% legend('Anticyclonic eddies (24%)','Cyclonic eddies (18%)','orientation','horizontal','location','best','AutoUpdate', 'off','box','off')
xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
title('Winter','fontweight','normal')
ylim([0 50]);grid on




pos_grid=lon;

figure;

% yyaxis right
% fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');ylim([1001 5500]);
% axis ij
% yyaxis left


plot(X(72:353),percentage_apo,'r')
hold on
plot(X(72:353),percentage_cpo,'b')
axis tight


ylabel({'% of time spent with ';'an eddy at 34.5ºS'})
legend('Anticyclonic eddies (24%)','Cyclonic eddies (18%)','orientation','horizontal','location','best','AutoUpdate', 'off','box','off')
xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
title('All days (n= 9292)','fontweight','normal')
ylim([-10 50]);grid on
yyaxis right
plot(lon,-topo,'k')

ylim([-5500 4000])




figure


% savefig('eddies_presence_clim')
% 
% 
% 
% cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian
% % 
%  save percentage_eddies34375 apo cpo percentage_cpo percentage_apo




load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian/stations.mat')
 
ae_in=nansum(ae_in);ce_in=nansum(ce_in);
ae_inn=find(ae_in==1);ce_inn=find(ce_in==1);

cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian
% tsg=importdata('msm_060_1_no_headers.tsg');
load stations
load merian_contours
cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian


load vu;v=v/100;u=u/100;

cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian
load MSM60_ctd_october.mat



pcoo=cat(2,percentage_apo_summer,-percentage_cpo_summer)


pa=repmat(percentage_apo_summer,1,8);


pc=repmat(-percentage_cpo_summer,1,8);

pcc=round(pc);

pcc = ceil(pc/10)*10;

paa = ceil(pa/10)*10;
% [aa,bb,cc]=meshgrid(X(72:353),-34.3750,percentage_cpo_summer)
% 
% cx=cc(1,1,:)





% figure
% m_proj('mercator','long', [-52 18],'lat', [-37 -32]);hold on
% [CS,CH]=m_etopo2('contour',[-3500 -200],'color',[.2 .2 .2]);%caxis([-5500 7000]);%colormap(flip(gray(7)));  
% hold on
% 
% %ii=3
% %m_pcolor(X(72:353),Y(82-ii:83+ii),pcoo');shading interp
% 
% ii=7
% m_image(X(72:353),Y(83:83+ii),pcc');shading interp
% hold on
% 
% m_image(X(72:353),Y(82-ii:82),paa');shading interp
% 
% caxis([-50 50])
% 
% colormap(rednblue)
% colorbar
% m_grid
% 
% 
% 
% 
% 
% 
% %[CS,CH]=m_etopo2('contourf',[-3500 -200],'edgecolor','none');caxis([-5500 7000]);colormap(flip(gray(7)));  
% 
% 
% m_plot(lon_adcp,lat_adcp,'k');hold on
% ll=2;
% 
% for id_time= 1:length(CEs_out60)
%     m_plot(CEs_out60{id_time,id_time,1},CEs_out60{id_time,id_time,2},'c','LineWidth',1)
% end
% 
% for id_time= 1:length(AEs_out60)
%     m_plot(AEs_out60{id_time,id_time,1},AEs_out60{id_time,id_time,2},'m','LineWidth',1)
% end
% 
% m_grid('xtick',[],'ytick',[],'linestyle','none','tickstyle','dd');
% m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
% 



figure('Renderer', 'painters', 'Position', [200 200 800 750]) 
ha = tight_subplot(3,1,[.01 .01],[.05 .01],[.06 .08])

axes(ha(1)); 


dates=datenum('1993-01-01'):datenum('2018-06-10');
datess=datevec(dates);


ind=find(dates==datenum(datenum('2017-01-04')))



cpo_merian=cpo(:,8770:8770+27)

apo_merian=apo(:,8770:8770+27)


pcpo=nsum(cpo_merian,2)/28

acpo=nsum(apo_merian,2)/28


nmean(nmean(pcpo))

nmean(nmean(acpo))


percentage_apo=nsum(apo,2)/length(apo)*100;

percentage_cpo=nsum(cpo,2)/length(cpo)*100;

timeseries_cpo=nsum(cpo,1)/282*100;
timeseries_apo=nsum(apo,1)/282*100;

timeseries_both=timeseries_cpo+timeseries_apo;

eddy_index=timeseries_apo-timeseries_cpo;


Time = datetime(datevec(dates));

YourData = table(Time,eddy_index',timeseries_cpo',timeseries_apo',timeseries_both');

TT = table2timetable(YourData);

eddy_index_monthly = retime(TT,'monthly','mean');

eddy_index_monthly_std = retime(TT,'monthly','std');

monthofyearaverage = groupsummary(TT, 'Time', 'monthofyear', 'mean')


figure;
subplot(2,1,1)
plot(X(72:353),percentage_apo_summer,'color',[1    0.3    0.3])
plot(X(72:353),percentage_apo_winter,'color',[0.75    0    0])
plot(X(72:353),percentage_apo,'r','linewidth',1);hold on

plot(X(72:353),percentage_cpo_summer,'color',[0.3    0.3    1])
plot(X(72:353),percentage_cpo_winter,'color',[0   0    0.75 ])
plot(X(72:353),percentage_cpo,'b','linewidth',1);hold on



load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian/stations.mat')
 
ae_in=nansum(ae_in);ce_in=nansum(ce_in);
ae_inn=find(ae_in==1);ce_inn=find(ce_in==1);

cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/satelite_Merian
% tsg=importdata('msm_060_1_no_headers.tsg');
load stations
load merian_contours
cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian



figure('Renderer', 'painters', 'Position', [200 200 600 350]) 
ha = tight_subplot(4,1,[.015 .015],[.05 .015],[.06 .015]); hold on

axes(ha(1)); hold on


text(0.98,0.25,'a','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')
text(0.25,0.95,'All days (n=9292)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','none')

plot(X(68:357),percentage_apo,'r');hold on
hold on
plot(X(68:357),percentage_cpo,'b')
axis tight
%ylabel({'% of time spent with ';'an eddy at 34.5ºS'})
legend('Anticyclonic eddies (19%)','Cyclonic eddies (14%)','orientation','horizontal','location','northeast','AutoUpdate', 'off','box','off')
%xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
xticks([-50:5:15]);xticklabels([]);
yticks([0:10:55]);yticklabels([0:10:55]);
%title('All days (n= 9292)','fontweight','normal')
ylim([0 55]);grid on
 box on

% subplot(4,1,2)
axes(ha(2)); hold on

text(0.98,0.25,'b','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

plot(X(68:357),percentage_apo_winter,'r')
hold on
plot(X(68:357),percentage_cpo_winter,'b')
axis tight
%ylabel({'% of time spent with ';'an eddy at 34.5ºS'})
%legend('Anticyclonic eddies (24%)','Cyclonic eddies (18%)','orientation','horizontal','location','best','AutoUpdate', 'off','box','off')
xticks([-50:5:15]);xticklabels([]);
yticks([0:10:55]);yticklabels([0:10:55]);

%xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
%title('Winter','fontweight','normal')
ylim([0 55]);grid on; box on

text(0.25,0.95,'Winter','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','none')

%subplot(4,1,3)

axes(ha(3)); hold on

text(0.98,0.25,'c','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

plot(X(68:357),percentage_apo_summer,'r')
hold on
plot(X(68:357),percentage_cpo_summer,'b')
axis tight
ylabel({'Time spent within an eddy at 34.375ºS (%)'})
%legend('Anticyclonic eddies (24%)','Cyclonic eddies (18%)','orientation','horizontal','location','best','AutoUpdate', 'off','box','off')
xticks([-50:5:15]);xticklabels([]);
%title('Summer','fontweight','normal')
yticks([0:10:55]);yticklabels([0:10:55]);
text(0.25,0.95,'Summer','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','none')

ylim([0 55]);grid on; box on

%subplot(4,1,4)
axes(ha(4)); hold on

text(0.98,0.25,'d','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

% plot(lon_stations,1,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
% plot(lon_stations(ce_inn),1,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
% plot(lon_stations(ae_inn),1,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on
pos_grid=lon;
fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');
axis ij
xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});

axis tight
grid on; box on

yticks([0:1000:5500]);yticklabels({'0','1','2','3','4','5'});
ylabel('Depth (km)')

xlim([X(68) X(357)])


set(findall(gcf,'-property','FontSize'),'FontSize',11.5)
set(findall(gcf,'-property','Linewidth'),'Linewidth',0.6)

set(gcf,'color','w'); 

% savefig('percentage_eddies')
% 
% %print(gcf, '-dpdf','painters','-r300','fig_percentage_eddies')
% print(gcf, '-dpng','-r600','fig_percentage_eddies')
% 






figure('Renderer', 'painters', 'Position', [200 200 600 350]) 
ha = tight_subplot(5,1,[.015 .015],[.05 .015],[.06 .015]); hold on

axes(ha(1)); hold on


text(0.98,0.25,'a','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')
text(0.25,0.95,'All days (n=9292)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','none')

plot(X(68:357),percentage_apo,'r');hold on
hold on
plot(X(68:357),percentage_cpo,'b')
axis tight
%ylabel({'% of time spent with ';'an eddy at 34.5ºS'})
legend('Anticyclonic eddies (19%)','Cyclonic eddies (14%)','orientation','horizontal','location','northeast','AutoUpdate', 'off','box','off')
%xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
xticks([-50:5:15]);xticklabels([]);
yticks([0:10:55]);yticklabels([0:10:55]);
%title('All days (n= 9292)','fontweight','normal')
ylim([0 55]);grid on
 box on

% subplot(4,1,2)
axes(ha(2)); hold on

text(0.98,0.25,'b','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

plot(X(68:357),percentage_apo_winter,'r')
hold on
plot(X(68:357),percentage_cpo_winter,'b')
axis tight
%ylabel({'% of time spent with ';'an eddy at 34.5ºS'})
%legend('Anticyclonic eddies (24%)','Cyclonic eddies (18%)','orientation','horizontal','location','best','AutoUpdate', 'off','box','off')
xticks([-50:5:15]);xticklabels([]);
yticks([0:10:55]);yticklabels([0:10:55]);

%xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});
%title('Winter','fontweight','normal')
ylim([0 55]);grid on; box on

text(0.25,0.95,'Winter','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','none')

%subplot(4,1,3)

axes(ha(3)); hold on

text(0.98,0.25,'c','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

plot(X(68:357),percentage_apo_summer,'r')
hold on
plot(X(68:357),percentage_cpo_summer,'b')
axis tight
ylabel({'Time spent within an eddy at 34.375ºS (%)'})
%legend('Anticyclonic eddies (24%)','Cyclonic eddies (18%)','orientation','horizontal','location','best','AutoUpdate', 'off','box','off')
xticks([-50:5:15]);xticklabels([]);
%title('Summer','fontweight','normal')
yticks([0:10:55]);yticklabels([0:10:55]);
text(0.25,0.95,'Summer','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','none')

ylim([0 55]);grid on; box on

%subplot(4,1,4)
axes(ha(4)); hold on



m_proj('mercator','long', [X(68) X(357)],'lat', [-36 -33]);hold on
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




axes(ha(5)); hold on

text(0.98,0.25,'d','Units', 'Normalized', 'VerticalAlignment', 'Top', 'Edgecolor','k')

% plot(lon_stations,1,'dk','MarkerFaceColor','k','MarkerEdgeColor','k');hold on
% plot(lon_stations(ce_inn),1,'dc','MarkerFaceColor','c','MarkerEdgeColor','c');hold on
% plot(lon_stations(ae_inn),1,'dm','MarkerFaceColor','m','MarkerEdgeColor','m');hold on
pos_grid=lon;
fill([pos_grid(1) pos_grid(1) pos_grid pos_grid(end)],[5500 topo(1) topo 5500],[.65 .65 .65],'linestyle','none');
axis ij
xticks([-50:5:15]);xticklabels({'50°W','45°W','40°W','35°W','30°W','25°W','20°W','15°W','10°W','5°W','0°E','5°E','10°E','15°E'});

axis tight
grid on; box on

yticks([0:1000:5500]);yticklabels({'0','1','2','3','4','5'});
ylabel('Depth (km)')

xlim([X(68) X(357)])


set(findall(gcf,'-property','FontSize'),'FontSize',11)
set(findall(gcf,'-property','Linewidth'),'Linewidth',0.6)

set(gcf,'color','w'); 

savefig('percentage_eddies2')

% print(gcf, '-dpdf','painters','-r300','fig_percentage_eddies2')
 print(gcf, '-dpng','-r600','fig_percentage_eddies2')





figure('Renderer', 'painters', 'Position', [200 200 600 350]) 
subplot(3,1,1)
plot(monthofyearaverage.mean_Var4,'r','linewidth',1); hold on
plot(monthofyearaverage.mean_Var3,'b','linewidth',1); hold on
axis tight
ylabel({'Time spent within' ; 'an eddy at 34.375ºS (%)'})
legend('Anticyclonic','Cyclonic','Edgecolor','none','orientation','horizontal','location','best')
ylim([13 22])
subplot(3,1,2)
plot(monthofyearaverage.mean_Var3+monthofyearaverage.mean_Var4,'k','linewidth',1); hold on
ylabel({'Total Time spent' ; ' within an eddy  (%)'})
axis tight
subplot(3,1,3)
plot(monthofyearaverage.mean_Var4-monthofyearaverage.mean_Var3,'k','linewidth',1); hold on
axis tight
ylabel({'AEs - CEs (%)'})

xlabel('Month')



aa=monthofyearaverage.mean_Var4-monthofyearaverage.mean_Var3


bb=(abs(monthofyearaverage.mean_Var4)+abs(monthofyearaverage.mean_Var4))

cc=aa./bb


plot(monthofyearaverage.mean_Var4-monthofyearaverage.mean_Var3,'k','linewidth',1); hold on

set(findall(gcf,'-property','FontSize'),'FontSize',11)
set(findall(gcf,'-property','Linewidth'),'Linewidth',0.6)

set(gcf,'color','w'); 

savefig('climatology_eddies')

% print(gcf, '-dpdf','painters','-r300','fig_percentage_eddies2')
 print(gcf, '-dpng','-r600','climatology_eddies')


