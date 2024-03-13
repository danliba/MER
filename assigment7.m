cd D:\Maestria\MER\SOTON\clases\data_analysis\assigment
clear all; close all; clc
%% 
load("chlor.mat");
load("etopo.mat");

%% Solution 7

[LAT_chl,LON_chl]=meshgrid(latchlor2,lonchlor2);
[LAT_bathy,LON_bathy]=meshgrid(latopo2,lontopo2);

interp_chlor=interp2(LAT_chl,LON_chl,chlor2,LAT_bathy,LON_bathy,'linear');
interp_chlor2=interp2(LAT_chl,LON_chl,chlor2,LAT_bathy,LON_bathy,"nearest");
% interp_chlor3=interp2(LAT_chl,LON_chl,chlor2,LAT_bathy,LON_bathy,"cubic");
% interp_chlor4=interp2(LAT_chl,LON_chl,chlor2,LAT_bathy,LON_bathy,"spline");
% interp_chlor5=interp2(LAT_chl,LON_chl,chlor2,LAT_bathy,LON_bathy,"makima");
%chlorophyll=reshape(interp_chlor,size(etopo));

%% plots
subplot(3,2,1)
pcolor(lonchlor2,latchlor2,log10(chlor2')); shading flat; colormap jet;
hc=colorbar;
caxis(log10([0.01 20]));
set(hc,'ticks',log10([0.01 0.1 0.5 1 2 3 4 5 10 20]),...
    'ticklabels',[0.01 0.1 0.5 1 2 3 4 5 10 20],'TickDirection',('out'));
ylabel(hc,'Chlorophyll-a mg/L','FontSize',8,'Rotation',270);
hc.Label.Position(1) = 3; title('Chlorophyll Data');

subplot(3,2,2)
pcolor(LON_bathy,LAT_bathy,log10(interp_chlor)); shading flat; colormap jet;
hc=colorbar;
caxis(log10([0.01 20]));
set(hc,'ticks',log10([0.01 0.1 0.5 1 2 3 4 5 10 20]),...
    'ticklabels',[0.01 0.1 0.5 1 2 3 4 5 10 20],'TickDirection',('out'));
ylabel(hc,'Chlorophyll-a mg/L','FontSize',8,'Rotation',270);
hc.Label.Position(1) = 3; title('Linear Interpolation');

subplot(3,2,3)
pcolor(LON_bathy,LAT_bathy,log10(interp_chlor2)); shading flat; colormap jet;
hc=colorbar;
caxis(log10([0.01 20]));
set(hc,'ticks',log10([0.01 0.1 0.5 1 2 3 4 5 10 20]),...
    'ticklabels',[0.01 0.1 0.5 1 2 3 4 5 10 20],'TickDirection',('out'));
ylabel(hc,'Chlorophyll-a mg/L','FontSize',8,'Rotation',270);
hc.Label.Position(1) = 3; title('Nearest Neighbor Interpolation');

subplot(3,2,4)
pcolor(interp_chlor3); shading flat; colormap jet;
colorbar;
title('Cubic Interpolation');

subplot(3,2,5)
pcolor(interp_chlor4); shading flat; colormap jet;
hc=colorbar;
title('Spline Interpolation');

subplot(3,2,6)
pcolor(interp_chlor5); shading flat; colormap jet;
hc=colorbar;
title('Makima Interpolation');

%% 
range0=[160 180 -60 -40];
grayColor = [.7 .7 .7];

subplot(1,3,1)
pcolor(lonchlor2,latchlor2,log10(chlor2')); shading flat; colormap jet;
hc=colorbar;
caxis(log10([0.01 20]));
set(hc,'ticks',log10([0.01 0.1 0.5 1 2 3 4 5 10 20]),...
    'ticklabels',[0.01 0.1 0.5 1 2 3 4 5 10 20],'TickDirection',('out'));
ylabel(hc,'Chlorophyll-a mg/L','FontSize',8,'Rotation',270);
hc.Label.Position(1) = 3; title('Chlorophyll Data');

subplot(1,3,2)
pcolor(LON_bathy,LAT_bathy,log10(interp_chlor)); shading flat; colormap jet;
hc=colorbar;
caxis(log10([0.01 20]));
set(hc,'ticks',log10([0.01 0.1 0.5 1 2 3 4 5 10 20]),...
    'ticklabels',[0.01 0.1 0.5 1 2 3 4 5 10 20],'TickDirection',('out'));
ylabel(hc,'Chlorophyll-a mg/L','FontSize',8,'Rotation',270);
hc.Label.Position(1) = 3; title('Linear Interpolation');

subplot(1,3,3)
pcolor(LON_bathy,LAT_bathy,log10(interp_chlor2)); shading flat; colormap jet;
hc=colorbar;
caxis(log10([0.01 20]));
set(hc,'ticks',log10([0.01 0.1 0.5 1 2 3 4 5 10 20]),...
    'ticklabels',[0.01 0.1 0.5 1 2 3 4 5 10 20],'TickDirection',('out'));
ylabel(hc,'Chlorophyll-a mg/L','FontSize',8,'Rotation',270);
hc.Label.Position(1) = 3; title('Nearest Neighbor Interpolation');

%% 
mean_chlor=nanmean(chlor2(:));
mean_interp_chlor=nanmean(interp_chlor(:));
mean_interp_chlor2=nanmean(interp_chlor2(:));
%Linear interpolation is the best in terms of fitting the right amount of
%chlorophyll
%% 
% 
% plot(log10(interp_chlor(:)),'.');
% hold on
% plot(log10(interp_chlor2(:)),'.');


%% 
vec_topo=etopo(:); vec_chlor=interp_chlor(:);

%plot(etopo(:),interp_chlor(:),'.');
%% SOLUTION 8
idx = ~any(isnan(vec_chlor),2);
topo=vec_topo(idx,:); 
chloro=vec_chlor(idx,:);

topo(topo>=-5)=NaN;
chloro(topo>=-5)=NaN;

idx2 =~any(isnan(topo),2);
topo2=abs(topo(idx2,:)); 
topos2=topo(idx2,:);
chloro2=chloro(idx2,:);

sum(isnan(topo2))

sum(isnan(chloro2))


%% 
semilogy(topo2,chloro2,'.'); grid on; axis square;
ylabel('Chlorophyll [mg/m^{3}]'); xlabel('Bathymetry [m]');

%% Logarithm 

coeff=polyfit(log(topo2),log(chloro2),1);

y_regression = polyval(coeff, log(topo2));

plot(log(topo2),log(chloro2),'.'); grid on; axis square;
ylabel('Ln (Chlorophyll [mg/m^{3}]) '); xlabel('Ln (Bathymetry [m])');
hold on
plot(log(topo2), y_regression, 'r', 'LineWidth', 2);
%b2 = robustfit(log(topo2),log(chloro2));
% [b,bint,r,rint] = regress(log(chloro2),log(topo2));
% rcoplot(r,rint)

rb = robustfit(log(topo2),log(chloro2));
y_rb = rb(1) + rb(2) * log(topo2);
%% Robust fit
plot(log(topo2),log(chloro2),'.'); grid on; axis square;
ylabel('Ln (Chlorophyll [mg/m^{3}]) '); xlabel('Ln (Bathymetry [m])');
hold on
plot(log(topo2), y_regression, 'r', 'LineWidth', 2);
hold on
plot(log(topo2), y_rb, 'g', 'LineWidth', 2);
legend('Data','Linear Correlation','Robust fit');
title('Linearized Correlation')
%eq1='$y = -0.2336x + 0.4161'
%R^2 = 0.3089
% Norm of residuals: 221.5

%% now use polytool


%% 
[fitresult, gof] = createFit(topo2, chloro2);

%It is more meaningful to examine the results with a non-linear function
%such as power function R^2= 0.4182 , C=3.928*B^(-0.3556). 
%% 
% Define the equations with LaTeX formatting

equation1 = '$R^2 = 0.4182$';
equation2 = '$C = 3.928 \cdot B^{-0.3556}$';

% Add text annotations to the figure
annotation('textbox', [0.6, 0.6, 0.4, 0.1], 'Interpreter', 'latex', 'String', equation1, 'FitBoxToText', 'on', 'EdgeColor', 'none');
annotation('textbox', [0.6, 0.5, 0.4, 0.1], 'Interpreter', 'latex', 'String', equation2, 'FitBoxToText', 'on', 'EdgeColor', 'none');

% Customize the text properties (optional)

% Adjust the figure size (optional)
set(gcf, 'Position', [100, 100, 600, 300]);  % Example figure size

ylabel('Chlorphyll [mg/m^3]'); xlabel('Bathymetry [m]'); title('Power Relation');
%% solution 8.2 
%Spearman correlation for non-linear 

[rho,pval] = corr(topos2,chloro2,'Type','spearman','Rows','complete');
%pvalue = 0;
%rho = -0.4898

%There's a negative non-linear correlation of 

%% solution 9
B=8000;
C=3.928*B^(-0.3556);
disp(C) % 0.1608

%I trust this at coefficient of determination of 41.82% 

%% Solution 10: Error propagation

% var_C = abs (dC/dh) * var_h;
rcoplot()

