%% try in a different location
longitude=-79.7500;
latitude=42.2500;

[LAG,corr]=XCORRi(longitude, latitude,1);
%% load the data
clear all;
lonssh=cell2mat(struct2cell(load('ssha_tpos_v01u.mat','long')));
latssh=cell2mat(struct2cell(load('ssha_tpos_v01u.mat','lat')));

loncl=cell2mat(struct2cell(load('chlo_swtp_v20b.mat','long')));
latcl=cell2mat(struct2cell(load('chlo_swtp_v20b.mat','lat')));

%%
for ilon=1:1:length(loncl)
    for ilat=1:1:length(latcl)
        [LAG,corr]=XCORRi(loncl(ilon), latcl(ilat),0);
        %[LAG,corr]=XCORR(-49, 33);
        LAGi(ilon,ilat,:)=LAG;
        CORRi(ilon,ilat,:)=corr;
    end
end

[LON,LAT]=meshgrid(loncl,latcl);
save('Q5_uchl_4',"CORRi","LAGi",'LON','LAT');
%%
region0=[-80 -5 20 44];
grayColor = [.7 .7 .7];

figure
p=pcolor(LON,LAT,CORRi'); shading interp; cmocean balance,
hold on; colorbar; caxis([-1 1]);
[c,h]=contour(LON,LAT,LAGi',[-28:0.5:25],'k:');
clabel(c,h);
title('Cross correlations');
hold on
borders('countries','facecolor',grayColor);
axis(region0)
xlabel('Longitude'); ylabel('Latitude');
%% 
clear all; close all; clc;

region0=[-80 -5 20 44];
grayColor = [.7 .7 .7];

load('Q5_uchl_4.mat');
LAG_new = inpaintn(LAGi,10);

figure
%subplot(2,1,1)
p=pcolor(LON,LAT,LAG_new'); shading interp; cmocean balance,
hold on; colorbar; caxis([-15 15]);
hold on
[c,h]=contour(LON,LAT,LAG_new',[-28:0.5:25],'w:','LineWidth',0.35);
clabel(c,h);
borders('countries','facecolor',grayColor);
title('SSH anomaly and  log10(chl/(1mg m^{−3})) Cross-correlation LAG','fontsize',14);
axis(region0)
xlabel('Longitude'); ylabel('Latitude');
hold on
a=plot(-50,34,'o','MarkerSize',15,'MarkerFaceColor','g','Color','g');
legend([a],{'34°N/50°W station'})

% % subplot(2,1,2)
p=pcolor(LON,LAT,LAGi'); shading interp; cmocean balance,
hold on; colorbar; caxis([-12 12]);
hold on
[c,h]=contour(LON,LAT,LAGi',[-28:0.5:25],'k:');
clabel(c,h);
borders('countries','facecolor',grayColor);
title('LAG');
axis(region0)
xlabel('Longitude'); ylabel('Latitude');

