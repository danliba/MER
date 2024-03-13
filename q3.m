%Problem 3 (25 points for SOES3042, 20 points for SOES6025): Create two maps: 
% one showing the amplitude of the annual cycle in the SSHA data at all available 
% locations in the region from 5°W – 80°W in longitude and from 20°N – 44°N in latitude 
% in the North Atlantic, and the other showing the amplitude of the annual cycle in the 
% log10(chl) data within the same region (using the chlorophyll data set in the data file 
% chlo_swtp_v20b.mat). The maps should be associated with appropriate colour 
% scales. Include comments on issues you encountered processing the data and how 
% you dealt with them.

chloro=cell2mat(struct2cell(load('chlo_swtp_v20b.mat','chlo')));
loncl=cell2mat(struct2cell(load('chlo_swtp_v20b.mat','long')));
latcl=cell2mat(struct2cell(load('chlo_swtp_v20b.mat','lat')));
timecl=cell2mat(struct2cell(load('chlo_swtp_v20b.mat','time')));


%% We select the region
%5°W – 80°W | 20°N – 44°N
region0=[-80 -5 20 44];

%chlorophyll
indxlon=find(loncl>=region0(1) & loncl<=region0(2));
indxlat=find(latcl>=region0(3) & latcl<=region0(4));

%new values of chlor
nloncl=loncl(indxlon);
nlatcl=latcl(indxlat);
nchl=chloro(:,indxlat,indxlon);

clear loncl latcl
%icl=interp1(timecl,nchl,tnew1);
%% Find NANs
nchl2=nchl(:);
masknan=isnan(nchl2);
mkdouble=double(masknan);

%reshape to create a nan mask
nanchlor=reshape(mkdouble,[168,48,150]);
nanchlor=permute(nanchlor,[3 2 1]);%lon lat time
nanchlor2=sum(nanchlor,3);
nanchlor3=nanchlor2/168*100;
%NANplot map
subplot(1,2,1)
pcolor(nloncl,nlatcl,nanchlor3'); shading interp;
title('Chlorophyll-a NAN percentage (%)'); colorbar; 
%% new time since 1998
tlen=1998:mean(diff(timecl)):2002;

tnew=linspace(1998,2002,length(timecl(timecl>=1998 & timecl<=2002)));
tnew1 = tnew(1:end-1)';
tnew1(end)+mean(diff(tnew1))-tnew1(1)

%% delete NAN
nchl2=permute(nchl,[3 2 1]);%lon lat time

%figure
for ilon=1:1:length(nloncl)
    for ilat=1:1:length(nlatcl)
        nchltry=permute(nchl2(ilon,ilat,:),[3 2 1]);
        % plot(timecl,nchltry,'.:');
        idx = ~any(isnan(nchltry),2);
        nonanCHL=nchltry(idx);
        nonantime=timecl(idx);
        if sum(idx)<18 %168 -150 threshold
            icl = NaN(length(tnew1),1);
            disp(['Land - sum ' ,num2str(sum(icl))])
        else
        %hold on
        %we interpolate chlor
        icl=interp1(nonantime,nonanCHL,tnew1,'linear');
        if sum(icl)==0
            disp('WHY?')
            break
        end
        %hold on
        % plot(tnew1,icl)
        % pause(0.5)
        % clf
        newchlor(ilon,ilat,:)=icl; %lon lat time
        end
    end
end

%% now we calculate the amplitude with fourier transform
t1m=tnew1;

for ilon=1:1:length(nloncl)
    for ilat=1:1:length(nlatcl)
        chlorf=permute(newchlor(ilon,ilat,:),[3 2 1]); %lon lat time
        [amp1,amp2]=fourier_amp(t1m,chlorf);
        AMPcl(ilon,ilat)=amp1;
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% NOW SSH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ssh=cell2mat(struct2cell(load('ssha_tpos_v01u.mat','ssha')));
lonssh=cell2mat(struct2cell(load('ssha_tpos_v01u.mat','long')));
latssh=cell2mat(struct2cell(load('ssha_tpos_v01u.mat','lat')));
timessh=cell2mat(struct2cell(load('ssha_tpos_v01u.mat','time')));

%% We select the region
%5°W – 80°W | 20°N – 44°N
region0=[-80 -5 20 44];

%ssh
% clear indxlon indxlat
indxlon=find(lonssh>=region0(1) & lonssh<=region0(2));
indxlat=find(latssh>=region0(3) & latssh<=region0(4));

%new values for ssh
nlonssh=lonssh(indxlon);
nlatssh=latssh(indxlat);
nssh=ssh(:,indxlat,indxlon);

%% Find NANs
nssh2=nssh(:);
masknan=isnan(nssh2);
mkdouble=double(masknan);

%reshape to create a nan mask
nanssh=reshape(mkdouble,[350,25,76]);%time lat lon
nanssh=permute(nanssh,[3 2 1]);%lon lat time
nanssh2=sum(nanssh,3);
nanssh3=nanssh2/350*100;% percentaje of NaN
%NANplot map
subplot(1,2,2)
pcolor(nlonssh,nlatssh,nanssh3'); shading interp;
title('Sea Surface Height NAN Percentage (%)'); colorbar; 

%% new time since 1998
tlen=1993:mean(diff(timecl)):2002;

tnew=linspace(1993,2002,length(timecl(timecl>=1993 & timecl<=2002)));
tnew1 = tnew(1:end-1)';
tnew1(end)+mean(diff(tnew1))-tnew1(1)

%%
nssh2=permute(nssh,[3 2 1]);%lon lat time

%figure
for ilon=1:1:length(nlonssh)
    for ilat=1:1:length(nlatssh)
        nsshtry=permute(nssh2(ilon,ilat,:),[3 2 1]);
        % plot(timecl,nchltry,'.:');
        idx = ~any(isnan(nsshtry),2);
        nonanSSH=nsshtry(idx);
        nonantime=timessh(idx);
        if sum(idx)<10 %350 -340 threshold
            issh = NaN(length(tnew1),1);
            disp(['Land - sum ' ,num2str(sum(issh))])
        else
        %hold on
        %we interpolate chlor
        issh=interp1(nonantime,nonanSSH,tnew1,'linear');
        if sum(issh)==0
            disp('WHY?')
            break
        end
        %hold on
        % plot(tnew1,icl)
        % pause(0.5)
        % clf
        newssh(ilon,ilat,:)=issh; %lon lat time
        end
    end
end

%% now we calculate the amplitude with fourier transform
t1m=tnew1;

for ilon=1:1:length(nlonssh)
    for ilat=1:1:length(nlatssh)
        sshf=permute(newssh(ilon,ilat,:),[3 2 1]); %lon lat time
        [amp1,amp2]=fourier_amp(t1m,sshf);
        AMPssh(ilon,ilat)=amp2;
    end
end

%% plot amplitudes
nanchlor2(nanchlor2>150)=NaN; %LAND
nanchlor2(nanchlor2<150)=1;

AMPcl=AMPcl.*nanchlor2;

nanssh2(nanssh2>340)=NaN; %LAND
nanssh2(nanssh2<340)=1;
AMPssh=AMPssh.*nanssh2;
% AMPcl(AMPcl==0)=NaN; % because LAND
% AMPssh(AMPssh==0)=NaN;
%% 
grayColor = [.7 .7 .7];

figure
subplot(1,2,1)
pcolor(nloncl,nlatcl,AMPcl'); shading interp; 
colorbar; colormap jet;
hold on
borders('countries','facecolor',grayColor);
title('Amplitude of Chlorophyll-a from 1998-2002');
axis(region0)
xlabel('Longitude'); ylabel('Latitude');

subplot(1,2,2)
pcolor(nlonssh,nlatssh,AMPssh'); shading interp; 
colorbar; colormap jet;
hold on
borders('countries','facecolor',grayColor);
title('Amplitude of Sea Surface Height Anomaly from 1998-2002');
axis(region0)
xlabel('Longitude'); ylabel('Latitude');
