%this codes helps to find Crosscor using fourier transform, give a
%longitude and latitude for the ssh and chlor data
%Flag is for showing the graphs (1) and dont (0)
function [LAG,corr]=XCORRi(longitude, latitude,flag)
dflonlat=0.3;
ssh=cell2mat(struct2cell(load('ssha_tpos_v01u.mat','ssha')));
lonssh=cell2mat(struct2cell(load('ssha_tpos_v01u.mat','long')));
latssh=cell2mat(struct2cell(load('ssha_tpos_v01u.mat','lat')));
timessh=cell2mat(struct2cell(load('ssha_tpos_v01u.mat','time')));

indxlonssh=find(lonssh>=longitude-dflonlat & lonssh<=longitude+dflonlat);
indxlatssh=find(latssh>=latitude-dflonlat & latssh<=latitude+dflonlat);

chlor=cell2mat(struct2cell(load('chlo_swtp_v20b.mat','chlo')));
loncl=cell2mat(struct2cell(load('chlo_swtp_v20b.mat','long')));
latcl=cell2mat(struct2cell(load('chlo_swtp_v20b.mat','lat')));
timecl=cell2mat(struct2cell(load('chlo_swtp_v20b.mat','time')));

indxloncl=find(loncl>=longitude-dflonlat & loncl<=longitude+dflonlat);
indxlatcl=find(latcl>=latitude-dflonlat & latcl<=latitude+dflonlat);
%%
chloro=chlor(:,indxlatcl(1),indxloncl(1));
sshanom=ssh(:,indxlatssh(1),indxlonssh(1));
%we find a common start and end, 1998 to 2002
%now we crop
indxcl=find(timecl>=1998 & timecl<=2002);
indxssh=find(timessh>=1998 & timessh<=2002);

ntimecl=timecl(indxcl);
ntimessh=timessh(indxssh);

nchlor=chloro(indxcl);
nssh=sshanom(indxssh);
%remove nans
idx2=~any(isnan(nssh),2);
idx=~any(isnan(nchlor),2);
if sum(idx2)<10
    nssh=NaN(length(indxssh),1); ntimessh=NaN(length(indxssh),1);
    disp('SSH is empty')
    LAG=NaN; corr=NaN;
    return
else
    nssh=nssh(idx2); ntimessh=ntimessh(idx2);
end

if sum(idx)<18
    nchlor=NaN(length(indxcl),1);
    ntimecl=ntimessh;
    disp('CHL is empty')
    LAG=NaN; corr=NaN;
    return
else
    nchlor=nchlor(idx); ntimecl=ntimecl(idx);
end
%everything is cropped , now we detrend and use the fourier transform
%% chlor-----1
tlencl=1998:mean(diff(timecl)):2002;
tnew2=linspace(1998,2002,length(timecl(timecl>=1998 & timecl<=2002)));
tnew2 = tnew2(1:end-1)';
tnew2(end)+mean(diff(tnew2))-tnew2(1);
chlor1=interp1(ntimecl,nchlor,tnew2);
%ssh-------2
tlen1=1998:mean(diff(timessh)):2002;
tnew1=linspace(1998,2002,length(timessh(timessh>=1998 & timessh<=2002)));
tnew1 = tnew1(1:end-1)';
tnew1(end)+mean(diff(tnew1))-tnew1(1);
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
if flag==1
figure
plot(t1m, s1m,'b--.', t1m, s1filt, 'b-')
hold on
plot(t2m, s2m,'r--.', t2m, s2filt, 'r-');
title(['Anual cycles at ',num2str(longitude),'ºW & ',num2str(latitude),'ºN']);
xlabel('time (yrs)')
ylabel('signal (m)')
legend('SSh signal','Anual cycle ssh','Chla signal','Anual cycle Chla')
grid on
box  on
hold on
end
%% now the cross correlation
% xcorr is most easily applied to datasets that satisfy these conditions.
[x12pre,lag12pre]=xcorr(s1m,s2m); % lag12pre stores the lag of s1 behind s2,
lag12=lag12pre*ti1*12; %mean sampling interval in months
[a11,lag11]=xcorr(s1m);
[a22,lag22]=xcorr(s2m);
lagind11=find(lag11==0);
lagind22=find(lag22==0);
x12=x12pre/sqrt(a11(lagind11).*a22(lagind22));
x12shift=ifftshift(x12);
lag12shift=lag12-min(lag12);
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
% Now we calculate the Fourier transform, and fftshift it:
X12=fft(x12shift);
X12shift=fftshift(X12);
% if we were interested in amplitudes, we would need to normalise the
% magnitude by the length of the record.
X12mag=abs(X12shift);
% and now we calculate the angle
X12phs=angle(X12shift);
%%
if flag==1
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
if isempty(indxLAG)==1
    disp('EMPTY')
    LAG=NaN;
    corr=NaN;
else
    lagRad=X12phs(indxLAG(1));
    f_lag=f12(indxLAG(1));

    corr=x12(indxLAG(1));

    LAG=lagRad/(2*pi*f_lag);
    disp(['LAG ',num2str(LAG)])
    disp(['CORR ',num2str(corr)])
end
end


