cd D:\Maestria\MER\SOTON\clases\data_analysis\assigment;

fnchlor='AQUA_MODIS.20210101_20211231.L3m.YR.CHL.chlor_a.9km.nc';
fnTopo='ETOPO2v2g_f4.nc';

lontopo=double(ncread(fnTopo,'x'));
latopo=double(ncread(fnTopo,'y'));
bathy=double(ncread(fnTopo,'z'));

%a rectangle with corners at (60°S, 160°E) and (40°S, 180°E)
range0=[160 180 -60 -40];
indxlon=find(lontopo>=range0(1) & lontopo<=range0(2));
indxlat=find(latopo>=range0(3) & latopo<=range0(4));

lontopo2=lontopo(indxlon);
latopo2=latopo(indxlat);
etopo=bathy(indxlon,indxlat);


%% chlor

lonchlor=double(ncread(fnchlor,'lon'));
latchlor=double(ncread(fnchlor,'lat'));
chlor=double(ncread(fnchlor,'chlor_a'));

indxlon2=find(lonchlor>=range0(1) & lonchlor<=range0(2));
indxlat2=find(latchlor>=range0(3) & latchlor<=range0(4));

lonchlor2=lonchlor(indxlon2);
latchlor2=latchlor(indxlat2);
chlor2=chlor(indxlon2,indxlat2);
%% plot
xlabs={'160ºE','162ºE','164ºE','166ºE','168ºE','170ºE','172ºE','174ºE',...
    '176ºE','178ºE','180ºE'}; 
ylabs={'40ºS','42ºS','44ºS','46ºS','48ºS','50ºS','52ºS','54ºS','56ºS','58ºS',...
    '60ºE'};

grayColor = [.7 .7 .7];

subplot(1,2,1)
pcolor(lontopo2,latopo2,etopo'); shading flat; colorbar;
hold on
borders('countries','facecolor',grayColor);
title('Region: 60°S, 160°E and 40°S, 180°E');
axis(range0)
set(gca,'ytick',[-60:2:-40],'yticklabel',flip(ylabs));
set(gca,'xtick',[160:2:180],'xticklabel',xlabs);
caxis([-6000 0]);

subplot(1,2,2)
pcolor(lonchlor2,latchlor2,log10(chlor2')); shading flat; colormap jet;
hc=colorbar;
caxis(log10([0.01 20]));
set(hc,'ticks',log10([0.01 0.1 0.5 1 2 3 4 5 10 20]),...
    'ticklabels',[0.01 0.1 0.5 1 2 3 4 5 10 20],'TickDirection',('out'));
ylabel(hc,'Chlorophyll-a mg/L','FontSize',11,'Rotation',270);
hc.Label.Position(1) = 3;
hold on
borders('countries','facecolor',grayColor);
axis(range0)
title('Region: 60°S, 160°E and 40°S, 180°E');
set(gca,'ytick',[-60:2:-40],'yticklabel',flip(ylabs));
set(gca,'xtick',[160:2:180],'xticklabel',xlabs);

%% save topo
save('etopo.mat','lontopo2','latopo2','etopo');
% chlorophyll
save('chlor.mat','latchlor2','lonchlor2','chlor2');
