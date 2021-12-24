clear all
close all

input_dir_nc='/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/msm60ctd_nc';

load /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/msm60_dvs.mat

pos_grid=[-52:0.05:18];

kk=find(depth>6000);
depth(kk)=NaN;clear kk
lon=longitude;

for ip=1:length(pos_grid)
    kk=find(lon>=pos_grid(ip)-0.05 & lon<=pos_grid(ip)+0.05);
    topo(ip)=nmean(depth(kk));
end

topoo=fillmissing(topo,'linear');topoo=sw_pres(topoo',-35);topoo=topoo';


%Get names and variables to extract, it can be seen with:
% aa=ncinfo('msm_060_1_ctd_128.nc')
% vars=aa.Variables.Name

vars_to_load =  {'TIME','TIME_QC','LATITUDE','LATITUDE_QC','LONGITUDE','LONGITUDE_QC','PRES','PRES_QC','TEMP','TEMP_QC','PSAL','PSAL_QC','DOX2','DOX2_QC','TURB','TURB_QC','FLU2','FLU2_QC','PAR','PAR_QC','SPAR','SPAR_QC','BAT','BAT_QC','water_depth'};
%creation_date = ncreadatt('msm_060_1_ctd_126.nc','/','rodb_hdr')

projectdir = '/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/msm60ctd_nc';
info = dir( fullfile(projectdir, '*.nc') );
n_files = length(info);
filenames = fullfile( projectdir, {info.name} );

%create a cell to store each variable, QC is the quality, at flu2 they get
%the other way round
TIME= cell(n_files,1);
TIME_QC= cell(n_files,1);
LATITUDE= cell(n_files,1);
LATITUDE_QC= cell(n_files,1);
LONGITUDE= cell(n_files,1);
LONGITUDE_QC= cell(n_files,1);
PRES= cell(n_files,1);
PRES_QC= cell(n_files,1);
TEMP= cell(n_files,1);
TEMP_QC= cell(n_files,1);
PSAL= cell(n_files,1);
PSAL_QC= cell(n_files,1);
DOX2= cell(n_files,1);
DOX2_QC= cell(n_files,1);
TURB= cell(n_files,1);
TURB_QC= cell(n_files,1);
FLU2= cell(n_files,1);
FLU2_QC= cell(n_files,1);
PAR= cell(n_files,1);
PAR_QC= cell(n_files,1);
SPAR= cell(n_files,1);
SPAR_QC= cell(n_files,1);
BAT= cell(n_files,1);
BAT_QC= cell(n_files,1);
wdepth=cell(n_files,1);

%fill the cell with the values

for k = 1 : n_files
  f = filenames{k};
  TIME{k} = ncread(f, vars_to_load{1});
TIME_QC{k} = ncread(f, vars_to_load{2}); 
  LATITUDE{k} = ncread(f, vars_to_load{3});
LATITUDE_QC{k} = ncread(f, vars_to_load{4});
   LONGITUDE{k}= ncread(f, vars_to_load{5});
LONGITUDE_QC{k}= ncread(f, vars_to_load{6});
  PRES{k}= ncread(f, vars_to_load{7});
PRES_QC{k}= ncread(f, vars_to_load{8});
   TEMP{k}= ncread(f, vars_to_load{9});
TEMP_QC{k}= ncread(f, vars_to_load{10});
   PSAL{k}= ncread(f, vars_to_load{11});
PSAL_QC{k}= ncread(f, vars_to_load{12});
   DOX2{k}= ncread(f, vars_to_load{13});
DOX2_QC{k}= ncread(f, vars_to_load{14});
   TURB{k}= ncread(f, vars_to_load{14});
TURB_QC{k}= ncread(f, vars_to_load{15});
   FLU2{k}= ncread(f, vars_to_load{16});
FLU2_QC{k}= ncread(f, vars_to_load{17});
   PAR{k}= ncread(f, vars_to_load{18});
PAR_QC{k}= ncread(f, vars_to_load{19});
   SPAR{k}= ncread(f, vars_to_load{20});
SPAR_QC{k}= ncread(f, vars_to_load{21});
   BAT{k}= ncread(f, vars_to_load{22});
BAT_QC{k}= ncread(f, vars_to_load{23});
BAT_QC{k}= ncread(f, vars_to_load{23});
wdepth{k}= ncread(f, vars_to_load{24});
end

stations=1:length(PRES); ppres=0:max(cellfun(@max, PRES));

max_pres=cellfun(@max, PRES);

%create a empty matrix
temp_new=NaN(length(ppres),length(stations));
sal_new=NaN(length(ppres),length(stations));
ox_new=NaN(length(ppres),length(stations));
flu2_new=NaN(length(ppres),length(stations));
turb_new=NaN(length(ppres),length(stations));
pres_new=NaN(length(ppres),length(stations));

%put values at thier corresponding depth

    for kk=1:length(stations);
    press=(PRES{kk,1});
    
    tempp=(TEMP{kk,1});
     sall=(PSAL{kk,1});
       oxx=(DOX2{kk,1});
        flu2=(FLU2_QC{kk,1});
        turb=(TURB_QC{kk,1});
     
    for pp=press(1):length(press) ;
ppp=pp-press(1)+1;
      temp_new(pp,kk)=tempp(ppp);
       sal_new(pp,kk)=sall(ppp); 
        ox_new(pp,kk)=oxx(ppp);
         flu2_new(pp,kk)=flu2(ppp);
         turb_new(pp,kk)=turb(ppp);
           pres_new(pp,kk)=press(ppp);
           
    end
    end

% 
ppres=ppres';z=repmat(ppres,1,126);

pres_new=sw_pres(ppres,-35);%zp=sw_pres(z,-35);


real_values=temp_new./temp_new;%real values is a nan 1 matrix with 1 when there's a real value

%temp_new(kk,1:min(temp)-1)=min(temp)

temp=temp_new;sal=sal_new;ox=ox_new;pres=pres_new;
% 
%finds the first and the last measured value
temp_up=cellfun(@(x) x(1),TEMP);
temp_down=cellfun(@(x) x(end),TEMP);
sal_up=cellfun(@(x) x(1),PSAL);
sal_down=cellfun(@(x) x(end),PSAL);
ox_up=cellfun(@(x) x(1),DOX2);
ox_down=cellfun(@(x) x(end),DOX2);


for i=1:11
    for j=1:126
    if isnan(temp(i,j))
temp(i,j)= temp_up(j);
ox(i,j)= ox_up(j);
sal(i,j)= sal_up(j);
    end
    end
end

% u=u1db;v=v1db;
pres=ppres;lon=cell2mat(LONGITUDE);lat=cell2mat(LATITUDE);

max_z=sw_dpth(max_pres,-35);



%    [SA_msm60, in_ocean] = gsw_SA_from_SP(sal_new,pres_new,lon,-34.5);
% 
%    CT_msm60 = gsw_CT_from_t(SA,temp_new,pres_new);
% 
% lon_msm60=lon;


tempp=fillmissing(temp,'linear',2,'EndValues','nearest');
temp=griddata(lon',1:5458,tempp,pos_grid',0:5500);
temp=fillmissing(temp,'nearest');temp=fillmissing(temp,'linear',2,'EndValues','nearest');



sall=fillmissing(sal,'linear',2,'EndValues','nearest');
sal=griddata(lon',1:5458,sall,pos_grid',0:5500);
sal=fillmissing(sal,'nearest');sal=fillmissing(sal,'linear',2,'EndValues','nearest');


oxx=fillmissing(ox,'linear',2,'EndValues','nearest');
ox=griddata(lon',1:5458,oxx,pos_grid',0:5500);
ox=fillmissing(ox,'nearest');ox=fillmissing(ox,'linear',2,'EndValues','nearest');

real_values=griddata(lon',1:5458,real_values,pos_grid',0:5500);



bucket=ones(5501,length(topo));
for i=1:length(topo);
    bucket(topoo(i):end,i)=NaN;
end





% figure;
% pcolor(rr);shading interp;hold on; plot(topoo)
% axis ij
% 
% tempp=temp.*bucket;
% 
% figure;pcolor(tempp);shading interp;axis ij

% %44 AND 45 ARE AT THE SAME POINT, KEEP 45 BECAUSE WENT DEEPER
% sal(:,44)=[];temp(:,44)=[];ox(:,44)=[];real_values(:,44)=[];
% lat(44)=[];lon(44)=[];%max_z(44)=[];Waterdepth(44)=[];
% % u(:,44)=[];v(:,44)=[];
% %llat=-34.5;llat=repmat(llat,1,125);sw_dist(lon,llat','km'); find(ans<5)
% 
% %2 STATIONS MORE IN THE SAME LONGITUDE, KEEP THE DEEPEST ONE
% sal(:,83)=[];temp(:,83)=[];ox(:,83)=[];real_values(:,83)=[];
% lat(83)=[];lon(83)=[];%max_z(83)=[];Waterdepth(83)=[];
% %u(:,83)=[];v(:,83)=[];
stations=cell2mat(LONGITUDE);
lon=pos_grid;
clearvars -except  temp sal ox lon time real_values topo bucket stations max_pres ;

% temp=flip(temp);sal=flip(sal);ox=flip(ox);real_values=flip(real_values);lat=flip(lat);lon=flip(lon);
% topo=flip(topo);bucket=flip(bucket);

pres=0:5500;pres=pres';

lat=ones(1,length(lon)); lat=lat.*-34.5;

theta=sw_ptmp(sal,temp,pres,0); 

ro=sw_pden(sal,temp,pres,0);

gamma = gamma_GP_from_SP_pt(sal,theta,pres,lon,lat);

 save /Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/MSM60_ctd_october