%% Problem 2 (25 points for SOES3042, 20 points for SOES6025): Filter the data at 40°N,
% 58°W to produce two new time series: one showing the annual cycle only, and the
% other showing the signal with the annual cycle and higher frequencies removed.
% Include information on how you constructed your filters and how you applied them etc.
% I will remove frequencies lesser than 1 year.
load('ssha_tpos_v01u.mat');

indxlon=find(long==-58);
indxlat=find(lat==40);

sshanom=ssha(:,indxlat,indxlon);
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
fs1   = 1/ti1;                % the sampling frequency
fmax1 = fs1-fi1;              % before applying the fftshift, the maximum
                              % frequency is always the sampling frequency
                              % minus the frequency increment
fpre1 = (fmin1:fi1:fmax1+fi1/5)';   % frequency grid before applying the fftshift
f1    = [fpre1(N21+1:N1)-fs1; fpre1(1:N21)]; % this moves the second half of the

% Now we apply the Fourier transform
%
S1=fft(s1m,N1);            % N1 can be changed to a different value. What
                           % happens?
%
% Then we apply fftshift to S1 and S2:
S1shift=fftshift(S1);
% and compute the magnitude
S1mag=fftshift(abs(S1));

%% Filtering combined
df1   = fi1/2;

f0 = 1; %anual cycle

% find freqs component
flt1 = (fpre1 > f0 - df1       & fpre1 < f0 + df1)'; % select f0
flt2 = (fpre1 > fs1 - f0 - df1 & fpre1 < fs1 - f0 + df1)'; % select fs-f0

flta= (fpre1 < f0 - df1 | fpre1 > fs1 - f0 + df1)';
% apply the filter in freq domain
fltb=flt1+flt2+flta;
S1filt = S1.*fltb';
    

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
% time domain plot
subplot(2,2,1), plot(t1m, s1m,'b--.', t1m, s1filt, 'r-')
    xlabel('time (yrs)') 
    ylabel('signal (m)')
    legend('unfiltered','filtered')
    grid on 
    box  on
    hold on
    title('SSH anomaly filter')
  
% freq domain plot    
subplot(2,2,3), plot(f1, S1mag/N1, 'b--.', f1, S1filtmag/N1, 'r-')
    xlabel('frequency (years^{-1})'), 
    ylabel('amplitude (m)')
    legend('unfiltered','filtered')
    grid on
    box  on
    hold on
    title('Anual cycle + Higher frequencies filter')

%% filter Anual cycle
clear fltb
df1   = fi1/2;

f0 = 1; %anual cycle

% find freqs component
 flt1 = (fpre1 > f0 - df1       & fpre1 < f0 + df1)'; % select f0
 flt2 = (fpre1 > fs1 - f0 - df1 & fpre1 < fs1 - f0 + df1)'; % select fs-f0

%fltb= (fpre1 < f0 - df1 | fpre1 > fs1 - f0 + df1)';
% apply the filter in freq domain
 fltb=flt1+flt2;
S1filt = S1.*fltb';

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

figure(1)
subplot(2,2,2), plot(t1m, s1m,'b--.', t1m, s1filt, 'r-')
    xlabel('time (yrs)') 
    ylabel('signal (m)')
    legend('unfiltered','filtered')
    grid on 
    box  on
    hold on
    title('SSH anomaly annual cycle')
% freq domain plot    
subplot(2,2,4), plot(f1, S1mag/N1, 'b--.', f1, S1filtmag/N1, 'r-')
    xlabel('frequency (years^{-1})'), 
    ylabel('amplitude (m)')
    legend('unfiltered','filtered')
    grid on
    box  on
    hold on
    title('Band-pass Anual cycle filter')