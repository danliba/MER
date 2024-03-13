% Problem 4 (25 points for SOES3042, 20 points for SOES6025): You will notice that at 
% 34°N, 50°W the annual signal is large in both the SSHA and the log10(chl) data. Plot 
% the two time series from this location in the same plot (yyaxis is useful for this). Then 
% estimate the lag in months between SSHA and log10(chl) using 1-D cross-correlation 
% of the two time series.

%% load the data
clear all;
ssh=cell2mat(struct2cell(load('ssha_tpos_v01u.mat','ssha')));
lonssh=cell2mat(struct2cell(load('ssha_tpos_v01u.mat','long')));
latssh=cell2mat(struct2cell(load('ssha_tpos_v01u.mat','lat')));
timessh=cell2mat(struct2cell(load('ssha_tpos_v01u.mat','time')));

indxlonssh=find(lonssh==-50);
indxlatssh=find(latssh==34);

chlor=cell2mat(struct2cell(load('chlo_swtp_v20b.mat','chlo')));
loncl=cell2mat(struct2cell(load('chlo_swtp_v20b.mat','long')));
latcl=cell2mat(struct2cell(load('chlo_swtp_v20b.mat','lat')));
timecl=cell2mat(struct2cell(load('chlo_swtp_v20b.mat','time')));

indxloncl=find(loncl>=-50.5 & loncl<=-49.5);
indxlatcl=find(latcl>=33.5 & latcl<=34.5);
%% 
chloro=chlor(:,indxlatcl(1),indxloncl(1));
sshanom=ssh(:,indxlatssh,indxlonssh);


%we find a common start and end, 1998 to 2002
%now we crop
indxcl=find(timecl>=1998 & timecl<=2002);
indxssh=find(timessh>=1998 & timessh<=2002);

ntimecl=timecl(indxcl);
ntimessh=timessh(indxssh);

nchlor=chloro(indxcl);
nssh=sshanom(indxssh);
%everything is cropped , now we detrend and use the fourier transform

%% 
%chlor-----1
tlencl=1998:mean(diff(timecl)):2002;

tnew2=linspace(1998,2002,length(timecl(timecl>=1998 & timecl<=2002)));
tnew2 = tnew2(1:end-1)';
tnew2(end)+mean(diff(tnew2))-tnew2(1)

chlor1=interp1(ntimecl,nchlor,tnew2);

%ssh-------2
tlen1=1998:mean(diff(timessh)):2002;

tnew1=linspace(1998,2002,length(timessh(timessh>=1998 & timessh<=2002)));
tnew1 = tnew1(1:end-1)';
tnew1(end)+mean(diff(tnew1))-tnew1(1)

ssh1=interp1(ntimessh,nssh,tnew1);

%% fourier CHLOR
%t1= ssh1; t2=chlor1;
t1m=tnew1(2:end); %remove the first NaN
s1m=detrend(ssh1(2:end),'omitnan');
t2m=tnew2(2:end);
s2m=detrend(chlor1(2:end),'omitnan');

[s1filt,S1mag,N1,f1,S1filtmag,ti1]=fourier_t(t1m,s1m);
[s2filt,S2mag,N2,f2,S2filtmag,ti2]=fourier_t(t2m,s2m);


% time domain plot
figure
subplot(2,1,1)
plot(t1m, s1m,'b--.'); hold on; plot(t1m, s1filt, 'b-','LineWidth',2);
hold on
plot(t2m, s2m,'r--.');hold on; plot(t2m, s2filt, 'r-','LineWidth',2);
title('Anual cycles Band-pass')
xlabel('time (yrs)')
ylabel('signal (m)')
legend('SSHa signal','SSHa Band-pass filtered','Chla signal','Chla Band-pass filtered')
grid on
box  on

subplot(2,1,2)
plot(timecl,chloro); ylabel('log10(chla)');
yyaxis right
plot(timessh,sshanom);ylabel('SSA [m]');
legend('Log10(chla)','SSA')
title('Chla and SSA raw signal');
%% now the cross correlation
% xcorr is most easily applied to datasets that satisfy these conditions.
%
[x12pre,lag12pre]=xcorr(s1m,s2m); % lag12pre stores the lag of s1 behind s2,
                                  % in timesteps, not in hours (xcorr does
                                  % not know anything about the time vector)
%
% Convert lag12pre into values of unit hours.
%
lag12=lag12pre*ti1*12; %mean sampling interval * 
%
% Next we normalise the cross-correlation to fall between -1 and 1,
% for which we use the autocorrelation with zero lag (see
% lecture). This could also be done directly by xcorr (see help xcorr).
%
[a11,lag11]=xcorr(s1m);
[a22,lag22]=xcorr(s2m);
lagind11=find(lag11==0);
lagind22=find(lag22==0);
x12=x12pre/sqrt(a11(lagind11).*a22(lagind22));
%
% Now calculate the cross-spectrum by Fourier transforming the
% cross-correlation. We could also do this using cpsd. Because fft
% assumes the first element to always be at the origin of the
% corresponding axis, we need to make sure the zero lag component is
% the first element of the data we pass over to the fft if we want to
% directly interpret the phase computed from the spectrum as the phase
% difference of the signals. Remember that the signal is assumed to be
% periodic, so we can simply take the data from less than zero and
% append it to the end. That is, we can use an inverse fftshift
% (ifftshift).
%
x12shift=ifftshift(x12);
lag12shift=lag12-min(lag12);
%
% Next we need to create a new frequency vector, corresponding to the
% time-lag axis of the cross-correlation of the signals. lag12 has
% replaced t1 or t2
%
N12=length(lag12);         % number of observations
NT12=lag12(end)-lag12(1);  % length of record (hours)
dT12=NT12/(N12-1);         % mean sample space (hours).
fs12=1/dT12;               % sampling frequency (max frequency)
df12=(1/NT12)*(N12-1)/N12; % minimum non-zero frequency.
%
f12pre=0:df12:fs12-df12;   % the frequency vector
f12=fftshift(f12pre);      % and the vector that corresponds to the
                           % fftshifted signal
f12(f12>(fs12-df12/2)/2)=f12(f12>(fs12-df12/2)/2)-fs12;
%
% Now we calculate the Fourier transform, and fftshift it:
%
X12=fft(x12shift);
X12shift=fftshift(X12);
%
% if we were interested in amplitudes, we would need to normalise the
% magnitude by the length of the record.
%
X12mag=abs(X12shift);
%
% and now we calculate the angle
X12phs=angle(X12shift);
%% 
if 1
   figure(3), clf
   subplot(3,1,1), grid on, box on, hold on
   plot(lag12,x12,'b--.')
   xlabel('lag (yrs)'), ylabel('correlation (m^2)')
   subplot(3,1,2), grid on, box on, hold on
   plot(f12,X12mag,'b--.')
   xlabel('frequency (yrs^{-1})'), ylabel('magnitude of xcorr (m^2)')
   subplot(3,1,3), grid on, box on, hold on
   plot(f12,X12phs,'b--.')
   xlabel('frequency (yrs^{-1})'), ylabel('lag (rads)')
end
%% 
indxLAG=find(X12mag==max(X12mag));
%lagRad=f12(indxLAG(1));
lagRad=X12phs(indxLAG(1));
disp(lagRad)
f_lag=f12(indxLAG(1));
disp(f_lag)
corr=x12(indxLAG(1));
% lagRad=-2.70661;
% f_lag=-1.01031;

LAG=lagRad/(2*pi*f_lag);

disp(LAG)
disp(corr)

%% try in a different location
longitude=-50;
latitude=34;

[LAG,corr]=XCORRi(longitude, latitude,1);


