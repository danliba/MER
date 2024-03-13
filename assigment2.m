load('ssha_tpos_v01u.mat');
%shows SSHA as a function of time at the location: 40°N, 52°W

indxlon=find(long==-52);
indxlat=find(lat==40);

sshanom=ssha(:,indxlat,indxlon);

plot(time,sshanom,'ko-','markerfacecolor','r');ylim([-0.6 0.65]);
xlabel('Time (year)'); title('Sea surface height anomaly observed at Lat.: 40°N and Long.: 58°W');
ylabel('Sea surface height anomaly (m)');
xlim([1992.5 2002.5]);
%% Problem 1 (25 points for SOES3042, 20 points for SOES6025): Determine the
% amplitude of the annual cycle in the SSHA data at 40°N, 58°W using the Fast Fourier 
% Transform.
clear indxlon indxlat sshanom

indxlon=find(long==-58);
indxlat=find(lat==40);

sshanom=ssha(:,indxlat,indxlon);
plot(diff(time));

%% now we create a new time that has to be evenly spaced finding the delta t
tlen=1993:mean(diff(time)):2002;

tnew=linspace(1993,2002,length(time(time>=1993 & time<=2002)));
tnew1 = tnew(1:end-1)';
tnew1(end)+mean(diff(tnew1))-tnew1(1)

sshanom1=interp1(time,sshanom,tnew1);
%%
% evenly distribute new time axis
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
                              % frequency grid to the start and subtracts the
                              % sampling frequency (The minimum frequency is
                              % essentially 1/length-of-record; the longest
                              % timescale that we can usefully get. However
                              % (N1-1)/N1 is needed, because the periodic
                              % function must return to the start the timestep
                              % AFTER the end of the record.)

% Now we apply the Fourier transform
%
S1=fft(s1m,N1);            % N1 can be changed to a different value. What
                           % happens?
%
% Then we apply fftshift to S1 and S2:
%
S1shift=fftshift(S1);

%
% and compute the magnitude
%
S1mag=fftshift(abs(S1));

%% Filtering
df1   = fi1/2;

f0 = 1; %anual cycle

% find freqs component
flt1 = (fpre1 > f0 - df1       & fpre1 < f0 + df1)'; % select f0
flt2 = (fpre1 > fs1 - f0 - df1 & fpre1 < fs1 - f0 + df1)'; % select fs-f0
flt  = flt1 + flt2;

% apply the filter in freq domain
S1filt = S1.*flt';

S1filtshift = fftshift(S1filt);
S1filtmag   = abs(S1filtshift);

indxFreq=find(f1>=1-df1 & f1<=1+df1);
indxFreq2=find(f1>=-1-df1 & f1<=-1+df1);
Sd1=S1mag(indxFreq);
Sd2=S1mag(indxFreq2);
amp1=(Sd1+Sd2)/N1;
% Apply inverse Fourier transform
s1filt = ifft(S1filt, N1);
amp2=(max(s1filt)-min(s1filt))/2;
disp([num2str(amp1) ' ' num2str(amp2)])

%%
figure(1)
clf

% time domain plot
subplot(2,1,1), plot(t1m, s1m,'b--.', t1m, s1filt, 'r-')
    xlabel('time (yrs)') 
    ylabel('sea surface height anomaly (m)')
    legend('ssha unfiltered signal','Anual cycle')
    grid on 
    box  on
    hold on
    title('SSH anomaly Anual Cycle')
    
% freq domain plot    
subplot(2,1,2), plot(f1, S1mag/N1, 'b--.', f1, S1filtmag/N1, 'r-')
    xlabel('frequency (years^{-1})'), 
    ylabel('amplitude (m)')
    legend('unfiltered','filtered')
    grid on
    box  on
    hold on
    title('Band-pass Anual cycle filter')

%% 
% Problem 1 (25 points for SOES3042, 20 points for SOES6025): Determine the
% amplitude of the annual cycle in the SSHA data at 40°N, 58°W using the Fast Fourier
% Transform.

%Amplitude of the Anual cycle
disp([num2str(amp1) ' ' num2str(amp2)])
%0.1786     %0.1786 not detrended
%0.1869     %0.1869 detrended
%% 
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
subplot(2,2,4), box on
plot(f,Sphs,'m-','linewidth',1), hold on
plot(f,SphsNaN,'mx','markersize',10,'linewidth',3)
xlabel('f (s^{-1})')
ylabel('ang(S)')
title('Phase of the Spectrum','fontsize',14)
set(gca,'ytick',pi*([-1:0.5:1]),'yticklabel',{'-pi' '-pi/2' '0' 'pi/2' 'pi'})

% Plot grid line at y=0 only
ax=axis; ax(3)=-pi; ax(4)=pi; axis(ax);
plot([ax(1) ax(2)],[0 0],'k--')
plot([0 0],[ax(3) ax(4)],'k--')


















