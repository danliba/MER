cd D:\Maestria\MER\SOTON\clases\data_analysis\assigment
clear all; close all; clc
%% 
load("chlor.mat");
load("etopo.mat");
%% 
% calculate basic descriptive statistics (mean, σ, skewness, kurtosis) for 
% both datasets and portray the datasets using box plots and histograms. 
% Describe each histogram shape using a single word or phrase; why does the data 
% have that shape?
etopo(etopo>0)=NaN;

meanTopo= nanmean(etopo(:)); %lon lat
meanChlor = nanmean(chlor2(:)); %lon lat 

stdTopo= nanstd(etopo(:)); %lon
stdChlor = nanstd(chlor2(:));

skTopo = skewness(etopo(:));
skChlor = skewness(chlor2(:));

ktTopo= kurtosis(etopo(:));
ktChlor = kurtosis(chlor2(:));

modTopo=mode(etopo(:));
modChlor=mode(chlor2(:));

medTopo=nanmedian(etopo(:));
medChlor=nanmedian(chlor2(:));
%% boxplots
subplot(1,2,1)
boxplot(etopo(:)); 
title('Bathymetry [m]');

subplot(1,2,2)
boxplot(log10(chlor2(:)));
title('Log10 (Chlorophyll-a [mg/m^{3}])')
%% Histograms
subplot(1,2,1)
histogram(etopo(:));
xlabel('Bathymetry [m]'); 
title('Bathymetry')
%bimodal , not normally distributed, 2 modes, median and mean are not equal

subplot(1,2,2)
histogram(log10(chlor2(:))); %not normal
title('Chlorophyll-a');
xlabel('Log10(Chlorophyll-a [mg/m^{3}])')
%% normality test
[h,p] = kstest(etopo(:),'Alpha',0.05); %no normally distributed
[hc,pc] = kstest(log10(chlor2(:)),'Alpha',0.05); %no normally distributed

subplot(1,2,1)
normplot(etopo(:));
subplot(1,2,2)
normplot(log10(chlor2(:)));
%%  Solution 5
clear all; clc; 
load("chlor.mat");
load("etopo.mat");

offshore=[162 166 -58 -54 ];
shelfseas=[174 178 -46 -42 ];

%Chlor
lonOff_chlor=find(lonchlor2>= offshore(1) & lonchlor2<=offshore(2));
lonShelf_chlor=find(lonchlor2>= shelfseas(1) & lonchlor2<=shelfseas(2));

latOff_Chlor=find(latchlor2>= offshore(3) & latchlor2<=offshore(4));
latShelf_chlor=find(latchlor2>= shelfseas(3) & latchlor2<=shelfseas(4));

chlorOffshore=chlor2(lonOff_chlor,latOff_Chlor);
chlorShelf=chlor2(lonShelf_chlor,latShelf_chlor);

%topo
lonOff_topo=find(lontopo2>= offshore(1) & lontopo2<=offshore(2));
lonShelf_topo=find(lontopo2>= shelfseas(1) & lontopo2<=shelfseas(2));

latOff_topo=find(latopo2>= offshore(3) & latopo2<=offshore(4));
latShelf_topo=find(latopo2>= shelfseas(3) & latopo2<=shelfseas(4));

%% plot
% Define the coordinates for the offshore square
offshore_plot = [162, -58, 4, 4];  % [left, bottom, width, height]

% Define the coordinates for the shelfseas square
shelfseas_plot = [174, -46, 4, 4];  % [left, bottom, width, height]

grayColor = [.7 .7 .7];

% Create a figure
figure;
pcolor(lontopo2,latopo2,chlor2'); shading flat; colorbar;
hold on
borders('countries','facecolor',grayColor);
hold on
% Plot the offshore square
rectangle('Position', offshore_plot, 'EdgeColor', 'r', 'LineWidth', 2);
hold on;
% Plot the shelfseas square
rectangle('Position', shelfseas_plot, 'EdgeColor', 'b', 'LineWidth', 2);
% Adjust the axis limits to fit the squares
axis([160 180 -60 -40]);
%axis equal
% Add labels
title('Regions');
xlabel('Longitude');
ylabel('Latitude');

% Add a legend
%legend('Offshore', 'Shelf Seas');
%% means
Mean_chlorOffshore=nanmean(chlorOffshore(:));
Mean_chlorShelf=nanmean(chlorShelf(:));
%%  Plots
text1=['Mean: ',num2str(round(Mean_chlorOffshore,2)),' mg/L'];
text2=['Mean: ',num2str(round(Mean_chlorShelf,2)),' mg/L'];

subplot(1,2,1)
pcolor(lonchlor2(lonOff_chlor),latchlor2(latOff_Chlor),log10(chlorOffshore')); shading flat;
hc=colorbar; colormap jet;
caxis(log10([0.01 20]));
set(hc,'ticks',log10([0.01 0.1 0.5 1 2 3 4 5 10 20]),...
    'ticklabels',[0.01 0.1 0.5 1 2 3 4 5 10 20],'TickDirection',('out'));
ylabel(hc,'Chlorophyll-a mg/L','FontSize',11,'Rotation',270);
hc.Label.Position(1) = 3;
box on;
title(['OffShore ', text1]);

subplot(1,2,2)
pcolor(lonchlor2(lonShelf_chlor),latchlor2(latShelf_chlor),log10(chlorShelf')); shading flat;
hc=colorbar;
caxis(log10([0.01 20]));
set(hc,'ticks',log10([0.01 0.1 0.5 1 2 3 4 5 10 20]),...
    'ticklabels',[0.01 0.1 0.5 1 2 3 4 5 10 20],'TickDirection',('out'));
ylabel(hc,'Chlorophyll-a mg/L','FontSize',11,'Rotation',270);
hc.Label.Position(1) = 3;
box on;
title(['OffShore ', text2]);


%% Solution 6 
%Statistical test
% T-test 1 one way, independent

%Ho = The means are equal , the amount of chlorophyll in the Shelf seas is
%the same as the amount of chlorophyll in Off Shore

%Ha = The concentrations of Chlorophyll in the Shelf Seas are greater than
%average on the Shelf Seas than in Off shore.

%[h1,p1]=ttest2(chlorShelf(:),chlorOffshore(:),'Alpha',0.05);
[p1,h1] = ranksum(chlorShelf(:),chlorOffshore(:),'alpha',0.05);

%At a significance level of 5% null hypothesis is rejected
[p2,h2]=ranksum(chlorShelf(:),chlorOffshore(:),'Alpha',0.01);

%At a significance level of 1% null hypothesis is rejected

[p3,h3]=ranksum(chlorShelf(:),chlorOffshore(:),'Alpha',1.0000e-04);
%At a significance level of 0.01% null hypothesis is rejected


