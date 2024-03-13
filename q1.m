clear all; close all;
%----------- SET UP DATA -----------
signal_make_grids % sets up time and frequency grids

%% 

load('ssha_tpos_v01u.mat');
indxlon=find(long==-58);
indxlat=find(lat==40);

sshanom=ssha(:,indxlat,indxlon);

tlen=1993:mean(diff(time)):2002;

tnew=linspace(1993,2002,length(time(time>=1993 & time<=2002)));
tnew1 = tnew(1:end-1)';
tnew1(end)+mean(diff(tnew1))-tnew1(1)

sshanom1=interp1(time,sshanom,tnew1);

t1m=tnew1;
s1m=detrend(sshanom1);

N1     = length(t1m);         % sample size
N21    = ceil(N1/2);          % half the number of observations (rounded up)
ts1    = t1m(1);              % first data point - Note: Phases will be
                              % calculated relative to this point, not 0
tlast1 = t1m(end);            % last data point
ti1    = (tlast1-ts1)/(N1-1); % mean sampling interval
te1    = tlast1+ti1;          % end plus one increment (i.e., hypothetical next time point)

% frequency vector
fmin1 = 0;                    % the minimum frequency is always zero
fi1   = 1/(te1-ts1);          % the lowest non-zero frequency in the data
                              % is always 1 divided by the length of the
                              % data set (but note how te was defined). This
                              % is also the frequency increment.
fs1   = 1/ti1;                % the sampling frequency
fmax1 = fs1-fi1;              % before applying the fftshift, the maximum
                              % frequency is always the sampling frequency
                              % minus the frequency increment
fpre1 = (fmin1:fi1:fmax1+fi1/5)';   % frequency grid before applying the fftshift
f1    = [fpre1(N21+1:N1)-fs1; fpre1(1:N21)]; % this moves the second half of the

%%
s = s1m; % Option 18 (loads some different data)

%----------- DO THE FOURIER TRANSFORM -----------
Spre = fft(s); % The spectrum of s (before the fftshift is applied)
S = fftshift(Spre);
%S = Spre; f = fpre; % Uncomment this line to understand what the fftshift did

% Get real and imagniary components
Sreal = real(S);
Simag = imag(S);


% Get |S| and phase information
Sabs = abs(S)./N1;
Senergy = Sabs.*conj(Sabs);
%Sabs=2*Sabs/length(Sabs);
Sphs = angle(S);

% Create alternative phase vector with phase set equal to NaN at
% frequencies where Sabs is very small, useful for visualisation.
phsthr = 1e-1; % This sets the amplitude (as a fraction of the maximum
% amplitude) below which phase is set equal to NaN;
SphsNaN = Sphs;
SphsNaN( Sabs < (phsthr*max(Sabs)) ) = NaN;
%-------------------------------------
flag=1;
%%
if flag==1

%----------- PLOTTING -----------
% Setup figure to plot signal and output
figure(1) % Change figure number to avoid losing earlier figure
clf % Clears figure
set(gcf,'Position',[100 100 1300 1100]) % Sets figure position; adjust to taste


% THE SIGNAL
subplot(2,2,1), box on
plot(time,sshanom);hold on;
plot(t1m,s1m,'ko-','markerfacecolor','r','MarkerSize',3.5);ylim([-0.6 0.65]);
xlabel('Time (year)'); title('Sea surface height anomaly observed at Lat.: 40°N and Long.: 58°W');
ylabel('Sea surface height anomaly (m)');
xlim([1992.5 2002.5]);
legend('SSH Anomaly','SSH detrended anomaly');


%pause % Waits for you to press any key (comment this out if you like)


% THE SPECTRUM
subplot(2,2,2), box on
plot(f1,Sreal,'b-x','linewidth',2), hold on
plot(f1,Simag,'r-.x','linewidth',2)
xlabel('f (s^{-1})')
ylabel('S (m)')
title('Spectrum (real part blue; imaginary part red)','fontsize',14)

% Plot grid lines at x=0, y=0 only
% THE MAGNITUDE OF THE SPECTRUM
subplot(2,2,3), box on
plot(f1,Sabs,'m-x','linewidth',2), hold on
xlabel('f (s^{-1})')
ylabel('|S| (m)')
title('Magnitude of the Spectrum','fontsize',14)

% Plot grid line at y=0 only
ax=axis;
plot([ax(1) ax(2)],[0 0],'k--')
plot([0 0],[ax(3) ax(4)],'k--')


% THE PHASE OF THE SPECTRUM
% subplot(2,2,4), box on
% plot(f1,Sphs,'m-','linewidth',1), hold on
% plot(f1,SphsNaN,'mx','markersize',10,'linewidth',3)
% xlabel('f (s^{-1})')
% ylabel('ang(S)')
% title('Phase of the Spectrum','fontsize',14)
% set(gca,'ytick',pi*([-1:0.5:1]),'yticklabel',{'-pi' '-pi/2' '0' 'pi/2' 'pi'})

% Plot grid line at y=0 only
% ax=axis; ax(3)=-pi; ax(4)=pi; axis(ax);
% plot([ax(1) ax(2)],[0 0],'k--')
% plot([0 0],[ax(3) ax(4)],'k--')

subplot(2,2,4), box on
plot(f1,Senergy,'k-','linewidth',2), hold on
xlabel('1/wavelength, k (m^{-1})')
ylabel('Energy Density (cm^2 s^{-2}/m')
title('Energy Spectral Density','fontsize',14)

%[f Sreal Simag Sabs Sphs];

end

%%

subplot(2,2,1), box on
plot(time,sshanom);hold on;
plot(t1m,s1m,'ko-','markerfacecolor','r','MarkerSize',3.5);ylim([-0.6 0.65]);
xlabel('Time (year)'); title('Sea surface height anomaly observed at Lat.: 40°N and Long.: 58°W','fontsize',12);
ylabel('Sea surface height anomaly (m)');
xlim([1992.5 2002.5]);
legend('SSH Anomaly','SSH detrended anomaly');

subplot(2,2,3), box on
plot(f1,Senergy,'k-','linewidth',2), hold on
xlabel('1/wavelength, k (m^{-1})')
ylabel('Energy Density (cm^2 s^{-2}/m')
title('Energy Spectral Density','fontsize',12)
grid on

subplot(2,2,2), plot(t1m, s1m,'b--.', t1m, s1filt, 'r-')
    xlabel('time (yrs)') 
    ylabel('sea surface height anomaly (m)')
    legend('ssha unfiltered signal','Anual cycle')
    grid on 
    box  on
    hold on
    title('SSH anomaly Anual Cycle','fontsize',12)
    
% freq domain plot    
subplot(2,2,4), plot(f1, S1mag/N1, 'b--.', f1, S1filtmag/N1, 'r-')
    xlabel('frequency (years)'), 
    ylabel('amplitude (m)')
    legend('unfiltered','filtered')
    grid on
    box  on
    hold on
    title('Band-pass Anual cycle filter','fontsize',12)
