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

%%
cd /Users/gaston/Documents/Phd/Database_South_Atl/Nc_files/2017_ene
adt=ncread('dt_global_allsat_phy_l4_20170101_20170530.nc','adt');
sla=ncread('dt_global_allsat_phy_l4_20170101_20170530.nc','sla');
lx=ncread('dt_global_allsat_phy_l4_20170101_20170530.nc','longitude');
ly=ncread('dt_global_allsat_phy_l4_20170101_20170530.nc','latitude');

mdt=adt-sla;
lx(lx>180)=lx(lx>180)-360;
[lx,id_sort]=sort(lx);
mdt=mdt(id_sort,:);ly=ly(163:300);lx=lx(437:836);mdtt=mdt(437:836,163:300);mdtt=mdtt';


cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian


figure
% p = get(h, 'pos');
%[[x0,y0,width,height]

% p(3) = p(3) + 0.8; p(1) = p(1) - 0.1;p(4) = p(4) + 0.1;p(2) = p(2) + 0.15;set(h, 'pos', p);box on; hold on
m_proj('Mercator','lon',[-68 20],'lat',[-46 -20]); 

 %[CS,CH]=m_etopo2('contourf',[-6000:500:0 250:250:2000],'edgecolor','none');
[CS,CH]=m_etopo2('contourf',[-6500:1000:-500 -200 0],'edgecolor','none');

hold on
%m_contour(lx,ly,mdt',[-0.1:0.1:0.8],'color',[.6 .6 .6]);

%m_contour(lx,ly,mdtt',[0 0.4 0.6 0.8],'color',[.5 .5 .5],'linewidth',1,'showtext','on','labelspacing',3400);

[C,hContour] = m_contour(lx,ly,mdtt,[0 0.4 0.6 0.8],'color',[.1 .1 .1],'linewidth',1);

clabel(C,hContour,'FontSize',12,'Color',[.1 .1 .1],'labelspacing',540);


% clabel(C,hContour,'manual','FontSize',12,'Color',[.5 .5 .5]);

 
% m_contour(lx,ly,mdt',[-0.1:0.2:0.6 0.9],'color','y');

 %[CS,CH]=m_etopo2('contourf',[-6500:500:-500 -200 0],'edgecolor','none');

%   m_grid('linestyle','none','tickdir','out');

m_plot(cell2mat(ctda9.lon),cell2mat(ctda9.lat),'color',[ 1.0000 0.8398 0],'linewidth',1)
m_text(-18,-25,'A09','color',[ 1.0000 0.8398 0],'linewidth',1)


m_plot(cell2mat(ctda10.lon),cell2mat(ctda10.lat),'color',[ 1.0000 0.8398 0],'linewidth',1)
m_text(-18,-31,'A10','color',[ 1.0000 0.8398 0],'linewidth',1)

m_plot(cell2mat(ctda11.lon),cell2mat(ctda11.lat),'color',[ 1.0000 0.8398 0],'linewidth',1)
m_text(-18,-44,'A11','color',[ 1.0000 0.8398 0],'linewidth',1)


m_plot(cell2mat(ctda11.lon),cell2mat(ctda11.lat),'color',[ 1.0000 0.8398 0],'linewidth',1)
m_text(-18,-44,'A11','color',[ 1.0000 0.8398 0],'linewidth',1)

m_plot(cell2mat(ctdmerian.lon),cell2mat(ctdmerian.lat),'color',[0.85 0 0],'linewidth',1)

m_text(-22,-35.5,'MSM60','color',[0.85 0 0],'linewidth',1)

%colormap([ m_colmap('blues',80); m_colmap('gland',48)]);
colormap(m_colmap('blues'));
% colormap('gray');

m_gshhs_h('patch',[0.81818 0.77647 0.70909],'edgecolor','none');

m_grid('xtick',[-60:10:20],'linestyle','none')

% brighten(.5);

%ax=m_contfbar([.002 .19],0.75,CS,CH,'endpiece','no','axfrac',.02);

ax=m_contfbar([.002 .19],0.75,CS,CH,'endpiece','no','axfrac',.02);

ax.XLabel.String='Depth (km)';


m_text(-52,-44,'Zapiola Gyre','color',[ 1.0000 0.8398 0],'linewidth',1)

m_text(-42,-29,'Rio Grande Rise','color',[ 1.0000 0.8398 0],'linewidth',1)

m_text(-0,-26,'Walvis Ridge','color',[ 1.0000 0.8398 0],'linewidth',1)


m_text(4,-37,'Cape Basin','color',[ 1.0000 0.8398 0],'linewidth',1)


m_text(-25,-39.5,'Mid-Atlantic Ridge','color',[ 1.0000 0.8398 0],'linewidth',1)

%ax=m_contfbar([.3 .7],.05,CS,CH);
% set(ax,'fontsize',9)

% title(ax,{'Depth (m)'}); % Move up by inserting a blank line
% 
% title(ax,{'Depth (m)'})
  set(findall(gcf,'-property','FontSize'),'FontSize',12)
% set(findall(gcf,'-property','Linewidth'),'Linewidth',.7)
% 
% 

 set(gcf,'color','w');



cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/figures_september
savefig('figure1_merian_set')

 print(gcf, '-dpng','-r600','figure1_merian_set')

 print(gcf, '-djpeg','-r600','figure1_merian_set')