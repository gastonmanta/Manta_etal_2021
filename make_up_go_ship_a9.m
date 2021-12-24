cd /Users/gaston/Documents/Phd/Datos_CTD_phd/GO_SHIP/A9_King/cchdo_740H20090307/only_a9_profiles
% projectdir2 = '/Users/gaston/Documents/Phd/Datos_CTD_phd/GO_SHIP/A9_King/cchdo_740H20090307/only_a9_profiles';
% 
% info2 = dir( fullfile(projectdir2, '*.nc') );
% n_files2 = length(info2);
% filenames2 = fullfile( projectdir2, {info2.name} );
% ctd2= struct;
% 
% %first 22 ctd profiles are outside the transect
% for i= 1 : n_files2;
% f = filenames2{i};
% ctd2.lat(i)=ncread(f,'latitude');
% ctd2.lon(i)=ncread(f,'longitude');
% ctd2.temp{i}=ncread(f,'temperature');
% ctd2.pres{i}=ncread(f,'pressure');
% ctd2.sal{i}=ncread(f,'salinity');
% ctd2.oxy{i}=ncread(f,'oxygen');
% end
% 
% cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/go_ship
% 
% pos_grid=[-52:0.05:18];
% 
% ppres=1:max(cellfun(@max, ctd2.pres)); stations=1:length(info2);
% 
% wppres=(cellfun(@max, ctd2.pres));
% 
% min_temp=(cellfun(@min, ctd2.temp));
% 
% % temp_new=NaN(max(cellfun(@max, ctd2.pres)),length(stations));
% 
% TEMP=ctd2.temp';
% PRES=ctd2.pres';
% PSAL=ctd2.sal';
% DOX2=ctd2.oxy';
% 
% 
% ppres=1:5817;
% 
% 
%     for kk=1:length(stations);
%      press=(PRES{kk,1});
%     
%     tempp=(TEMP{kk,1});
%       sall=(PSAL{kk,1});
%         oxx=(DOX2{kk,1});
% %         flu2=(FLU2_QC{kk,1});
% %         turb=(TURB_QC{kk,1});
%      
%   for pp=press(1):length(press) ;
%  ppp=pp-press(1)+1;
%        temp_new(pp,kk)=tempp(ppp);
%         sal_new(pp,kk)=sall(ppp); 
%          ox_new(pp,kk)=oxx(ppp);
% %          flu2_new(pp,kk)=flu2(ppp);
% %          turb_new(pp,kk)=turb(ppp);
% %            pres_new(pp,kk)=press(ppp);
% %            
%     end
%     end
% 
% 
% xx=ctd2.lon;xx=double(xx);
% mask=temp_new./temp_new;
% 
% ind=find(temp_new==0);temp_new(ind)=NaN;
% tt=fillmissing(temp_new,'linear',2,'EndValues','nearest');
% 
% ind=find(sal_new==0);sal_new(ind)=NaN;
% ss=fillmissing(sal_new,'linear',2,'EndValues','nearest');
% 
% 
% ind=find(ox_new==0);ox_new(ind)=NaN;
% oo=fillmissing(ox_new,'linear',2,'EndValues','nearest');
% 
% 
% temp_a9=griddata(xx',1:2:5817,tt,pos_grid',1:5817);
% sal_a9=griddata(xx',1:2:5817,ss,pos_grid',1:5817);
% oxy_a9=griddata(xx',1:2:5817,oo,pos_grid',1:5817);
% 
% mask_a9=griddata(xx',1:2:5817,mask,pos_grid',1:5817);
% 
% bucket=mask_a9;
% 
% % kk=8;
% % mask_a9(:,541-kk:543+kk)=1;
% % mask_a9(1:5500,532)=1;mask_a9(1:5500,552)=1;
% 
% 
% 
% 
% ppres=1:length(sal_a9);
% lon=pos_grid;
% z_a9=ppres;x=lon;
% pres_a9=1:5817;pres_a9=pres_a9';pos_grid_a9=pos_grid;
% buckett=bucket(1:end-1);
% 
% theta_a9=sw_ptmp(sal_a9,temp_a9,ppres',1);
% 
% 
% [SA_a9, in_ocean] = gsw_SA_from_SP(sal_a9,ppres',pos_grid,-30);
% 
% CT_a9 = gsw_CT_from_t(SA_a9,temp_a9,ppres');
% 
% gamma_GP_a9 = gamma_GP_from_SP_pt(sal_a9,theta_a9,ppres',pos_grid,-30);
% 
% 
% 
% % gamma=gamma.*bucket;
% gamma_GPa_9=gamma_GP_a9.*bucket;
% 
% 
% % stations_a9=stations
% 

% clearvars -except  theta_a9 sal_a9 ox_a9 mask_a9 gamma_GP_a9 pres_a9 pos_grid_a9 ;
% 
% save goship_a9

cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/go_ship
load goship_a9



lat=(ones(length(lon),1)*-24)';
lon=pos_grid_a9;

geo_strf_dyn_height_0 = gsw_geo_strf_dyn_height(sal_a9,theta_a9,pres_a9,0);

[geostrophic_velocity_0, mid_lat, mid_long]=gsw_geostrophic_velocity(geo_strf_dyn_height_0,lon,lat,0);


geo_strf_dyn_height_3400 = gsw_geo_strf_dyn_height(sal_a9,theta_a9,pres_a9,3400);

geostrophic_velocity_3400=gsw_geostrophic_velocity(geo_strf_dyn_height_3400,lon,lat,0);


%MEAN AND MEDIAN DEPTH OF GAMMA 28.1
for i=1:length(lon)-1
ind=find(gamma_GP_a9(:,i)>28.099 & gamma_GP_a9(:,i)<28.101);
% pres281(i)=nmedian(pres_a9(ind));
ppres281(i)=nmean(pres_a9(ind));
% v281(i)=nmedian(v_new(ind,i));
% vv281(i)=nmean(v_new(ind,i));
end




%CALCLUATES GEOSTROPHIC VEL AT GAMMA 28.1
for i=1:length(lon)-1
ind=find(gamma_GP_a9 (:,i)>28.099 & gamma_GP_a9(:,i)<28.101);
%ladcp281(i)=nmean(v_new(ind,i));
geostrophic_velocity_281m(i)=nmean(geostrophic_velocity_0(ind,i));
end

geostrophic_velocity_281m(200:end-200)=fillmissing(geostrophic_velocity_281m(200:end-200),'linear');

geostrophic_velocity281r=repmat(geostrophic_velocity_281m,length(geostrophic_velocity_0),1);

geostrophic_velocity_281=geostrophic_velocity_0-geostrophic_velocity281r; %clear geostrophic_velocity281r;

% geostrophic_velocity_281(:,1094)=((geostrophic_velocity_281(:,1093)+geostrophic_velocity_281(:,1095))./2);

length_grid_point=sw_dist([-24.5 -24.5],[-30 -30.05],'km')*1000

sw_dist([-30 -30],[-30 -30.05],'km')

sw_dist([-34.5 -34.5],[-30 -30.05],'km')

area= length_grid_point*mask_a9/1e6;buckett=mask_a9;

%ekman is -0.3

geostrophic_velocity_281_ek=geostrophic_velocity_281;
geostrophic_velocity_281_ek(1:100,:)=geostrophic_velocity_281_ek(1:100,:)-0.0016/2;


geostrophic_velocity_281_ek=geostrophic_velocity_281_ek.*buckett(:,1:end-1);
geostrophic_velocity_281_ek_balanced=geostrophic_velocity_281_ek-nsum(nsum(geostrophic_velocity_281_ek))./nsum(nsum((~isnan(geostrophic_velocity_281_ek))));


geostrophic_velocity_281_ek=geostrophic_velocity_281;
geostrophic_velocity_281_ek(1:100,:)=geostrophic_velocity_281_ek(1:100,:)-0.0016/2;


% geostrophic_velocity_3400_ek=geostrophic_velocity_3400;
% geostrophic_velocity_3400_ek(1:100,:)=geostrophic_velocity_3400_ek(1:100,:)-0.0016/2;
% 
% 
% geostrophic_velocity_281_ek=geostrophic_velocity_281_ek.*buckett(:,1:end-1);
% geostrophic_velocity_281_ek_balanced=geostrophic_velocity_281_ek-nsum(nsum(geostrophic_velocity_281_ek))./nsum(nsum((~isnan(geostrophic_velocity_281_ek))));
% 
% 
% geostrophic_velocity_3400_ek=geostrophic_velocity_3400_ek.*buckett(:,1:end-1);
% geostrophic_velocity_3400_ek_balanced=geostrophic_velocity_3400_ek-nsum(nsum(geostrophic_velocity_3400_ek))./nsum(nsum((~isnan(geostrophic_velocity_3400_ek))));
% 
% 
% 
% transportvgos3400=geostrophic_velocity_3400_ek.*area(:,1:end-1);%transport per area
% 
% transportvgos3400_balanced=geostrophic_velocity_3400_ek_balanced.*area(:,1:end-1);%transport per area
% 
% transportvgos34001_balanced=nansum(transportvgos3400_balanced,2);transportvgos34001_balanced=movmean(transportvgos34001_balanced,10);%averaged in z
% 



transportvgos28=geostrophic_velocity_281_ek.*area(:,1:end-1);%transport per area


transportvgos28_balanced=geostrophic_velocity_281_ek_balanced.*area(:,1:end-1);%transport per area

transportvgos281_balanced=nansum(transportvgos28_balanced,2);transportvgos28_balanced=movmean(transportvgos28_balanced,10);%averaged in z


figure;
subplot(2,2,1)
plot(transportvgos281_balanced,1:5817); hold on
% plot(transportvgos34001_balanced,1:5817)

axis ij

subplot(2,2,2)

plot(cumsum(transportvgos281_balanced),1:5817); hold on
% plot(cumsum(transportvgos34001_balanced),1:5817)


axis ij


subplot(2,2,3:4)
tt=nsum(transportvgos28_balanced,1);
plot(lon(1:end-1),ncumsum(tt))


% %following valla et al 2018

g_n_levels=[19 26.35 27.1 27.6 27.9 28.1 28.27 35];
g_n_levels=[19 26.35 27.1 27.6 27.9 28.1 28.27 38 38];


for i=2:length(g_n_levels)
ind=find(gamma_GP_a9(:,1:end-1)<g_n_levels(i) & gamma_GP_a9(:,1:end-1)>g_n_levels(i-1));

% transport_by_layer_hg_geos281(i)=nansum(transportvgos28_ek(ind));

% transport_by_layer_hgb(i)=nansum(transportva_balanced(ind));
% transport_by_layer_hg_geos3400b(i)=nansum(transportvgos3400a_balanced(ind));
% transport_by_layer_hg_geos3400argob(i)=nansum(transportvgos3400argo_balanced(ind));
transport_by_layer_hg_geos281b(i)=nansum(transportvgos28_balanced(ind));
% transport_by_layer_hg_geos3400b(i)=nansum(transportvgos3400_balanced(ind));

end


figure
tthg=nancumsum(transport_by_layer_hg_geos281b);

tthg=nansum(transport_by_layer_hg_geos281b);


plot(transport_by_layer_hg_geos281b,1:length(transport_by_layer_hg_geos281b));

ylabel('transport by layer');axis ij

ylabel('transport by layer')

axis ij;xline(0);

%yticks('manual')
yticklabels({'','SW','SACW','AAIW','UCDW','NADW','LCDW','AABW',''})
grid on


heat_transport=geostrophic_velocity_281_ek_balanced.*theta_a9(:,1:end-1)*1027.5*3900.*length_grid_point;

nsum(nsum(heat_transport))


nmean(nmean(theta_a9.*bucket))







transport_by_layer_hg_geos281b_a9=transport_by_layer_hg_geos281b;
transportvgos28_balanced_a9=transportvgos28_balanced;


heat_transport=geostrophic_velocity_281_ek_balanced.*CT_a9(:,1:end-1)*1027.5*3900.*area(:,1:end-1)*1e6;
nsum(nsum(heat_transport))


lonn=lon(1:end-1)+0.05/2;

Ctt=griddata(lon,pres_a9,CT_a9,lonn,pres_a9);

heat_transport=geostrophic_velocity_281_ek_balanced.*CT_a9(:,1:end-1)*1027.5*3900.*area(:,1:end-1)*1e6;
nsum(nsum(heat_transport))


heat_transport=geostrophic_velocity_281_ek_balanced.*Ctt*1027.5*3900.*length_grid_point;
nsum(nsum(heat_transport))

% sal_anom=1-(SA_a9-nmean(SA_a9.*bucket,2));

salm=nmean(nmean(SA_a9.*bucket));

sal_anom=1-(SA_a9/salm);sal_anom=sal_anom.*buckett;



freshwater_transport_a9=geostrophic_velocity_281_ek_balanced.*sal_anom(:,1:end-1).*area(:,1:end-1);

nsum(nsum(freshwater_transport_a9))

sal_anomm=griddata(lon,pres_a9,sal_anom,lonn,pres_a9);

freshwater_transport_a9=geostrophic_velocity_281_ek_balanced.*sal_anomm.*length_grid_point;

nsum(nsum(freshwater_transport_a9))/1e6

% nmean(nmean(theta_a9.*buckett))

nmean(nmean(CT_a9.*buckett))

nmean(nmean(SA_a9.*buckett))

nmean(nmean(oxy_a9.*buckett))




 %save transport_a9 transport_by_layer_hg_geos281b_a9 transportvgos28_balanced_a9








figure;
subplot(2,2,1)
plot(transportvgos281_balanced,z_a9')
axis ij

subplot(2,2,2)
plot(cumsum(transportvgos281_balanced),z_a9)

axis ij

subplot(2,2,3:4)
tt=nsum(transportvgos28_balanced,1);
plot(lon(1:end-1),ncumsum(tt))


axis ij

transport_by_layer_hg_geos281b_a9=transport_by_layer_hg_geos281b_a9;

transportvgos28_balanced_a9=transportvgos28_balanced_a9;


save transport_a9 transportvgos28_balanced_a9 transport_by_layer_hg_geos281b_a9


% 
% for i=2:8
% plot(lon(1:end-1),transport_by_layer_hg_geos3400b(i,:)*4,'Color', cmap(i-1, :));
% hold on
% 
% end
% %ylim([-2 2])
% 
% plot(lon(1:end-1),nsum(transportvgos3400a_balanced),'k');
% 
% 
% ylabel('Transport (Sv)')
% 
% 
% plot(transport_by_layer_hgb,1:length(transport_by_layer_hg),'k');hold on;
% 
% 
% % plot(transport_by_layer_hg_geos281bar,1:length(transport_by_layer_hg_geos281bar));hold on;
% % plot(transport_by_layer_hg_geos3400bar,1:length(transport_by_layer_hg_geos3400bar));hold on;
% 
% figure;plot(transport_by_layer_hg_geos281b_a9,1:length(transport_by_layer_hg_geos281b_a9))
% 
% axis ij;xline(0);
% 
% %yticks('manual')
% yticklabels({'','SW','SACW','AAIW','UCDW','NADW','LCDW','AABW',''})
% grid on
