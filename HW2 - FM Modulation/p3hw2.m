%% HW02 - Communication Systems Matlab Assignment 2
%% Cousrse Instructor : Dr.behrouzi
%% Student Name     Student ID
%% Mehrsa Pourya    95101247
%% Deadline : 16/9/97 
%% Dec2018
%% Problem 3
%% Q1 
clear 
clc
close
[mysignal,fs]=audioread('Audio.wav');              % Open audio                               
dt=1/fs;                                           % each t interval length
t = 0: dt : length(mysignal)/fs-dt ;               % time vector
f = linspace(-fs/2,fs/2,length(mysignal));         % freq. vetcor
fftsig=fft(mysignal);                              % fft
plot(f,abs(fftshift(fftsig)/length(fftsig)))
grid on 
title('Spectrum of audio')
xlabel('frequency(hz)')
ylabel('amplitude')
%%
sound(mysignal,fs)
%% Q2 snrlow 
clear 
clc
close
[mysignal,fs]=audioread('Audio.wav');            % open audio
dt=1/fs;
f = linspace(-fs/2,fs/2,length(mysignal));       % freq. vector
t = 0: dt : length(mysignal)/fs-dt ;             % time vector
nlowvariace=1e-1                                 % noise variance 
noiselow=sqrt(nlowvariace)*randn(length(mysignal),1);
nsiglow=mysignal+noiselow;                       % noised message
nmlft=fftshift(fft(nsiglow)); % calculate fourier transform 
figure
plot(f,abs(nmlft)/length(nmlft))
grid on 
title('Spectrum of noised audio')
xlabel('frequency(hz)')
ylabel('amplitude')
myfilter=zeros(1,length(nsiglow));      % this is gonna be our filter
                                   % k *[ 0 ... 0 ... 1 1 1 ... 0 0 0]
my1337indexlow=max(find(-8000>f)); % calculating our needed 
                                   % index to start 1 indexes considering
                                   % chonsen cut-off freq = 1340
my1337indextop=min(find(8000<f));  % calculating our needed 
                                   % index to stop 1 indexes considering
                                   % chonsen cut-off freq = 1340
myfilter(my1337indexlow:my1337indextop)=1;  % we bulid our ideal filter
filteredml=nmlft.*myfilter';                  % we filter our signal in freq domain 
filmessl=real(ifft(fftshift(filteredml))); % first we use fftshift to cancel 
                                         % effect of previous fftshift used
figure
plot(f,abs(fftshift(fft(filmessl)))/length(nmlft))
grid on 
title('Spectrum of filterd noised audio , SNRlow')
xlabel('frequency(hz)')
ylabel('amplitude')
figure
subplot(3,1,1)
plot(t,mysignal)
grid on 
title('time domain of audio, SNRlow ')
xlabel('time(s)')
ylabel('amplitude')
subplot(3,1,2)
plot(t,nsiglow)
grid on 
title('time domain of noised audio , SNRlow')
xlabel('time(s)')
ylabel('amplitude')
subplot(3,1,3)
plot(t,filmessl)
grid on 
title('time domain of filterd noised audio , SNRlow')
xlabel('time(s)')
ylabel('amplitude')
SNRlow=sum((abs(fft(mysignal)/length(mysignal))).^2)/...
    (sum(abs((fft(noiselow)/length(noiselow))).^2))
%%
sound(filmessl,fs)
%% Q2 snr high
clear 
close
[mysignal,fs]=audioread('Audio.wav');          %open audio
dt=1/fs;
f = linspace(-fs/2,fs/2,length(mysignal));
t = 0: dt : length(mysignal)/fs-dt ;           % time vector
nhighvariance=1e-7                             % noise variance
noisehigh=sqrt(nhighvariance)*randn(length(mysignal),1);
nsighigh=mysignal+noisehigh;
nmhft=fftshift(fft(nsighigh)); % calculate fourier transform 
figure
plot(f,abs(nmhft)/length(nmhft))
grid on 
title('Spectrum of noised audio')
xlabel('frequency(hz)')
ylabel('amplitude')
myfilter=zeros(1,length(nsighigh));      % this is gonna be our filter
                                   % k *[ 0 ... 0 ... 1 1 1 ... 0 0 0]
my1337indexlow=max(find(-8000>f)); % calculating our needed 
                                   % index to start 1 indexes considering
                                   % chonsen cut-off freq = 1340
my1337indextop=min(find(8000<f));  % calculating our needed 
                                   % index to stop 1 indexes considering
                                   % chonsen cut-off freq = 1340
myfilter(my1337indexlow:my1337indextop)=1;  % we bulid our ideal filter
filteredmh=nmhft.*myfilter';                  % we filter our signal in freq domain 
filmessh=real(ifft(fftshift(filteredmh))); % first we use fftshift to cancel 
                                         % effect of previous fftshift used
figure
plot(f,abs(fftshift(fft(filmessh)))/length(nmhft))
grid on 
title('Spectrum of filterd noised audio , SNRhigh')
xlabel('frequency(hz)')
ylabel('amplitude')
figure
subplot(3,1,1)
plot(t,mysignal)
grid on 
title('time domain of audio , SNRhigh')
xlabel('time(s)')
ylabel('amplitude')
subplot(3,1,2)
plot(t,nsighigh)
grid on 
title('time domain of noised audio , SNRhigh')
xlabel('time(s)')
ylabel('amplitude')
subplot(3,1,3)
plot(t,filmessh)
grid on 
title('time domain of filterd noised audio , SNRhigh')
xlabel('time(s)')
ylabel('amplitude')
SNRhigh=sum((abs(fft(mysignal)/length(mysignal))).^2)/...
    (sum(abs((fft(noisehigh)/length(noisehigh))).^2))
%%
sound(filmessh,fs)
%%  Q3 noise of demodulation
clear 
clc
close
[mysignal,fs]=audioread('Audio.wav');     
fc=50e3;
fsmin=400e3;
n=ceil(fsmin/fs);                    % new fs = fs *n
newfs=n*fs;                          % new sampling rate
dt=1/newfs;                          % new time interval length
b=3;                                 % beta
fd=b*8e3/max(mysignal);              % f delta 
upsampledsig=interp(mysignal,n);     % umsampled signal
t=0:dt:length(upsampledsig)*dt-dt;   % time vector
f = linspace(-newfs/2,newfs/2,length(upsampledsig));  % freq. vector
integme=cumtrapz(upsampledsig*dt);   % message integral
xct=cos(2*pi*fc*t'+2*pi*fd*integme); % modulated signal
hxct1=[0 ; xct(1:end-1)];          
hxct2=[xct(2:end); 0 ];
xctprime=(hxct2-hxct1)/(2*dt);       % derivative of xct
y2=abs(xctprime);                    % abs block
Y2=fftshift(fft(y2));                % we filter in freq. domain
g=1/4/fd;                            % gain of filter
myfilter=zeros(1,length(y2));      % this is gonna be our filter
                                   % k *[ 0 ... 0 ... 1 1 1 ... 0 0 0]
my1337indexlow=max(find(-8100>f)); % calculating our needed 
                                   % index to start 1 indexes considering
                                   % chonsen cut-off freq = 1340
my1337indextop=min(find(8100<f));  % calculating our needed 
                                   % index to stop 1 indexes considering
                                   % chonsen cut-off freq = 1340
myfilter(my1337indexlow:my1337indextop)=g;  % we bulid our ideal filter
filteredmh=Y2.*myfilter';                  % we filter our signal in freq domain 
yd=real(ifft(fftshift(filteredmh))); % first we use fftshift to cancel 
                                         % effect of previous fftshift used
output=yd-mean(yd); 
y3=downsample(output,n);                 % downsample
figure
dt=1/fs;
t2 = 0: dt : length(mysignal)/fs-dt ;  
subplot(2,1,1)
plot(t2',mysignal)
title('orginal signal')
subplot(2,1,2)
plot(t2',y3)
title('detected downsampled signal')
snr=sum((abs(fft(mysignal))/length(mysignal)).^2)/...
    sum((abs(fft(mysignal-y3))/length(mysignal)).^2)    % wanted snr
%%
sound(y3,fs)
%% Q 4 noise of demodulation based on beta
clear 
clc
close
[mysignal,fs]=audioread('Audio.wav');     % open audio
fc=50e3;                                  % carrier freq
fsmin=400e3;                              
n=ceil(fsmin/fs);                         % new fs = n * fs
newfs=n*fs;                               % new sampling rate
dt=1/newfs;
c=0;
for b=0.1 : 0.5 : 6;        % values of beta
    c=c+1;
fd=b*8e3/max(mysignal);     % fdelta
% same as previous parts , we modulate , demodulate and calculate snr
upsampledsig=interp(mysignal,n);  
t=0:dt:length(upsampledsig)*dt-dt;
f = linspace(-newfs/2,newfs/2,length(upsampledsig));
integme=cumtrapz(upsampledsig*dt);
xct=cos(2*pi*fc*t'+2*pi*fd*integme);
hxct1=[0 ; xct(1:end-1)];
hxct2=[xct(2:end); 0 ];
xctprime=(hxct2-hxct1)/(2*dt);
y2=abs(xctprime);
Y2=fftshift(fft(y2));
g=1/4/fd;
myfilter=zeros(1,length(y2));      % this is gonna be our filter
                                   % k *[ 0 ... 0 ... 1 1 1 ... 0 0 0]
my1337indexlow=max(find(-8100>f)); % calculating our needed 
                                   % index to start 1 indexes considering
                                   % chonsen cut-off freq = 1340
my1337indextop=min(find(8100<f));  % calculating our needed 
                                   % index to stop 1 indexes considering
                                   % chonsen cut-off freq = 1340
myfilter(my1337indexlow:my1337indextop)=g;  % we bulid our ideal filter
filteredmh=Y2.*myfilter';                  % we filter our signal in freq domain 
yd=real(ifft(fftshift(filteredmh))); % first we use fftshift to cancel 
                                         % effect of previous fftshift used
output=yd-mean(yd);   
y3=downsample(output,n);
snr(c)=sum((abs(fft(mysignal))/length(mysignal)).^2)/...
    sum((abs(fft(mysignal-y3))/length(mysignal)).^2);
end
plot(0.1:0.5:6,snr)
grid on 
ylabel('Output SNR')
xlabel('\beta')
title('\beta effect on Output SNR')
%% q 5 awgn added
clear 
clc
close
% same cms as before just we added noise
[mysignal,fs]=audioread('Audio.wav');     
fc=50e3;  
fsmin=400e3;
n=ceil(fsmin/fs);
newfs=n*fs;
dt=1/newfs;
c=0;
b=4;
    c=c+1;
fd=b*8e3/max(mysignal);
upsampledsig=interp(mysignal,n);  
t=0:dt:length(upsampledsig)*dt-dt;
f = linspace(-newfs/2,newfs/2,length(upsampledsig));
integme=cumtrapz(upsampledsig*dt);
xct=cos(2*pi*fc*t'+2*pi*fd*integme);
noisevariance=0.01      % noise variance
noise=sqrt(noisevariance)*randn(length(t),1);
xctn=xct+noise;         % noised signal
hxct1=[0 ; xctn(1:end-1)];
hxct2=[xctn(2:end); 0 ];
xctprime=(hxct2-hxct1)/(2*dt);
y2=abs(xctprime);%-mean(abs(xctprime));
Y2=fftshift(fft(y2));
g=1/4/fd;
myfilter=zeros(1,length(y2));      % this is gonna be our filter
                                   % k *[ 0 ... 0 ... 1 1 1 ... 0 0 0]
my1337indexlow=max(find(-8100>f)); % calculating our needed 
                                   % index to start 1 indexes considering
                                   % chonsen cut-off freq = 1340
my1337indextop=min(find(8100<f));  % calculating our needed 
                                   % index to stop 1 indexes considering
                                   % chonsen cut-off freq = 1340
myfilter(my1337indexlow:my1337indextop)=g;  % we bulid our ideal filter
filteredmh=Y2.*myfilter';                  % we filter our signal in freq domain 
yd=real(ifft(fftshift(filteredmh))); % first we use fftshift to cancel 
                                         % effect of previous fftshift used
output=yd-mean(yd);   
y3=downsample(output,n);
snroutput=sum((abs(fft(mysignal))/length(mysignal)).^2)/...
    sum((abs(fft(mysignal-y3))/length(mysignal)).^2)
snrinput=sum((abs(fft(xct))/length(xct)).^2)/...
    sum((abs(fft(noise))/length(xct)).^2)
%% Q6 SNRO_SNRI
% same cms as before
clear 
clc
close
[mysignal,fs]=audioread('Audio.wav');     
fc=50e3;
fsmin=400e3;
n=ceil(fsmin/fs);
newfs=n*fs;
dt=1/newfs;
b=4;
nval=logspace(-1,-4,10);
fd=b*8e3/max(mysignal);
upsampledsig=interp(mysignal,n);  
t=0:dt:length(upsampledsig)*dt-dt;
f = linspace(-newfs/2,newfs/2,length(upsampledsig));
integme=cumtrapz(upsampledsig*dt);
xct=cos(2*pi*fc*t'+2*pi*fd*integme);
for c=1 : 10
nv=nval(c);
noise=sqrt(nv)*randn(length(t),1);
xctn=xct+noise;
hxct1=[0 ; xctn(1:end-1)];
hxct2=[xctn(2:end); 0 ];
xctprime=(hxct2-hxct1)/(2*dt);
y2=abs(xctprime);
Y2=fftshift(fft(y2));
g=1/4/fd;
myfilter=zeros(1,length(y2));      % this is gonna be our filter
                                   % k *[ 0 ... 0 ... 1 1 1 ... 0 0 0]
my1337indexlow=max(find(-8100>f)); % calculating our needed 
                                   % index to start 1 indexes considering
                                   % chonsen cut-off freq = 1340
my1337indextop=min(find(8100<f));  % calculating our needed 
                                   % index to stop 1 indexes considering
                                   % chonsen cut-off freq = 1340
myfilter(my1337indexlow:my1337indextop)=g;  % we bulid our ideal filter
filteredmh=Y2.*myfilter';                  % we filter our signal in freq domain 
yd=real(ifft(fftshift(filteredmh))); % first we use fftshift to cancel 
                                         % effect of previous fftshift used
output=yd-mean(yd);   
y3=downsample(output,n);
snroutput(c)=sum((abs(fft(mysignal))/length(mysignal)).^2)/...
    sum((abs(fft(mysignal-y3))/length(mysignal)).^2);
snrinput(c)=sum((abs(fft(xct))/length(xct)).^2)/...
    sum((abs(fft(noise))/length(xct)).^2);
end
% from Q2
snrlow=0.067768532113011;
snrhigh=6.766518433338010e+04;
plot(snrinput,snroutput)
grid on 
xlabel('input SNR')
ylabel('output SNR')
title('output SNR based on input SNR')



