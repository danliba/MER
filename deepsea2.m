cd D:\Maestria\MER\SOTON\clases\deepSeaEco\assigment2;
%% 
fn='ROVtransect_DATA.xlsx';

[status,sheets] = xlsfinfo(fn);

%datafis=readtable(fn1,char(sheets(1)));
[numData1, textData1, raw1] = xlsread(fn, char(sheets(1))); %

data=numData1(1:12:end,3:end); %lat, lon, distance, seafloor-depth, temp,sal
%% 
plot(data(:,2),data(:,1),'.'); title('Location');
xlabel('Lon'); ylabel('Lat');


plot(data(:,3),data(:,5));
yyaxis right
plot(data(:,3),data(:,6));

plot(data(:,5),data(:,4),'.');

%% 
figure
P=get(gcf,'position');
P(3)=P(3)*4;
P(4)=P(4)*1;
set(gcf,'position',P);
set(gcf,'PaperPositionMode','auto');

ax1= axes;
yyaxis right;
plot(data(:,3),data(:,5),'.-');
ylabel('Temperature [ºC]')
pause(0.1) 
ax1.XTickMode = 'manual'; 
ax1.YTickMode = 'manual'; 
ax1.YLim = [min(ax1.YTick), max(ax1.YTick)];  % see [4]
ax1.XLimMode = 'manual'; 
grid(ax1,'on')
ytick = ax1.YTick;  xlabel('Distance along transect (m)')

yyaxis left;
plot(data(:,3),data(:,6),'x-');
ylabel('Salinity [UPS]')
%grid on
% create 2nd, transparent axes
ax2 = axes('position', ax1.Position);
plot(ax2,data(:,3),-data(:,4), 'k')             % see [3]
ax2.Color = 'none'; 
pause(0.1)  
grid(ax2, 'on')
% Horizontally scale the y axis to alight the grid (again, be careful!)
ax2.XLim = ax1.XLim; 
ax2.XTick = ax1.XTick; 
ax2.YLimMode = 'manual'; 
yl = ax2.YLim; 
ylabel('Seafloor-Depth [m]')
% horzontally offset y tick labels
ax2.YTickLabel = strcat(ax2.YTickLabel, {'              '}); 
ax2.YTick = linspace(yl(1), yl(2), length(ytick));      % see [2]

%set(gca,'xtick',[0:1:121],'xticklabel',round(data(:,3)));

% create 3rd transparent axes
% ax3 = axes('position', ax1.Position,'color','magenta');
% plot(ax3, data(:,3), 'm:','LineWidth',2); % Replace "SomeOtherData" with your actual data
% ax3.Color = 'none';
% pause(0.1);
% grid(ax3, 'on');
% ax3.XLim = ax1.XLim;
% ax3.XTick = ax1.XTick;
% ax3.YLimMode = 'manual';
% yl = ax3.YLim;
% ylabel('Distance')
% ax3.YTickLabel = strcat(ax3.YTickLabel, {'                                 '});
% ax3.YTick = linspace(yl(1), yl(2), length(ytick));
% 
% title('Fisical Variables');
% %set(gca,'xtick',[1:1:12],'xticklabel',meses);
% %xlim([1 12]);
% legend('Sedimentos')
% pause(1)
