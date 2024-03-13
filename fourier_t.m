    function [s1filt,S1mag,N1,f1,S1filtmag,ti1]=fourier_t(t1m,s1m)
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
    
    f0 = 1/4; %anual cycle
    
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
end