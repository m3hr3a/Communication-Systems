%% COMMUNICATION SYSTEMS - MATLAB ASSIGNMENT 1
%% Instructor    : Dr. H. Behrouzi
%% Deadline      : 18 Aban 1397
%% Student name  : Mehrsa Pourya
%% Student ID    : 95101247
%% Part 2 DSB Modulation
%% Q 1 spectrum of modulated signal
clear
clc
close all
[xt,fs] = audioread('signal.wav');           % message x(t)
t = 0 : 1/fs : (1/fs)*(length(xt)-1) ;       % time vector
f = -fs/2 : fs/(length(xt)-1) : fs/2;        % frequency vector    
Xf = fft(xt);                                % message spectrum 
vt=xt'.*cos(2*pi*75e3*t);                    % modulated signal
Vf=fft(vt);                                  % modulated signal spectrum
subplot(2,1,1)
% double sided signal spectrum
plot(f,fftshift(abs(Xf))/(length(t)-1))
title('|X(f)| Double Sided Spectrum')
xlabel('frequency')
ylabel('amplitude')
grid on
subplot(2,1,2)
% modulated signal spectrum  - we see an alised signal so we need to
% upsample
plot(f,fftshift(abs(Vf))/(length(t)-1))
title('|V(f)| Double Sided Spectrum')
xlabel('frequency')
ylabel('amplitude')
grid on
n=ceil(400000/fs);                            % calculate n for upsampling
fs=n*fs;                                      % new sampling rate
xt=interp(xt,n);                              % upsampled signal
t = 0 : 1/fs : (1/fs)*(length(xt)-1);         % new time vector
f = -fs/2 : fs/(length(xt)-1) : fs/2;         % new frequency vector    
Xf = fft(xt);                                 % upsamppled message spectrum 
vt=xt'.*cos(2*pi*75e3*t);                     % new modulated signal
Vf=fft(vt);                                   % upsampled modulated s spectrum
figure
subplot(2,1,1)
% double sided signal spectrum
plot(f,fftshift(abs(Xf))/(length(t)-1))
title('|X(f)| Double Sided Spectrum upsampled-interp')
xlabel('frequency')
ylabel('amplitude')
grid on
subplot(2,1,2)
% modulated signal spectrum
plot(f,fftshift(abs(Vf))/(length(t)-1))
title('|V(f)| Double Sided Spectrum upsampled-interp')
xlabel('frequency')
ylabel('amplitude')
grid on
%% Q 2 DSB detetecting with semi-synchronized and envelope detector
clear
clc
close all
[xt,fs] = audioread('signal.wav');           % read message x(t)
n=ceil(400000/fs);                           % n used for umsampling
fs=n*fs;                                     % new sampling rate
xt=interp(xt,n);                             % upsampled message
t = 0 : 1/fs : (1/fs)*(length(xt)-1) ;       % time vector
f = -fs/2 : fs/(length(xt)-1) : fs/2;        % frequency vector    
Xf = fft(xt);
mt=xt'.*cos(2*pi*75e3*t);
% semi-synchronized method
pos=find(mt>=0);             % hard limiter +
neg=find(mt<0);              % hard limiter -
HLout(pos)=1;
HLout(neg)=-1;
YH=fftshift((fft(HLout)));
myfilter1=zeros(1,length(mt));      % this is gonna be our filter
                                   % k *[ 0 ... 0 ... 1 1 1 ... 0 0 0]
my1337indexlow=max(find(-76000>f)); % calculating our needed 
                                   % index to start 1 indexes considering
                                   % chonsen cut-off freq = 1340
my1337indextop=min(find(76000<f));  % calculating our needed 
                                   % index to stop 1 indexes considering
                                   % chonsen cut-off freq = 1340
myfilter1(my1337indexlow:my1337indextop)=0.5/(2*sin(pi/2)/pi);  ...
    % we bulid our ideal filter
YHf=YH.*myfilter1;        % we filter our signal in freq domain 

ydt1=real(ifft(fftshift(YHf)));
y1t=mt.*ydt1;    % multiply recieved signal by extracted carrier
Yf=fftshift(fft(y1t));        % calculate fourier transform 
                              % we use fftshift because we want 
                              % a symmetrical spectrum to be multipled
                              % by our ideal rect filter in frequency
                              % domain
myfilter=zeros(1,length(mt));      % this is gonna be our filter
                                   % k *[ 0 ... 0 ... 1 1 1 ... 0 0 0]
my1337indexlow=max(find(-5000>f)); % calculating our needed 
                                   % index to start 1 indexes considering
                                   % chonsen cut-off freq = 1340
my1337indextop=min(find(5000<f));  % calculating our needed 
                                   % index to stop 1 indexes considering
                                   % chonsen cut-off freq = 1340
myfilter(my1337indexlow:my1337indextop)=4;  % we bulid our ideal filter
Ydf=Yf.*myfilter;        % we filter our signal in freq domain 
ydt=real(ifft(fftshift(Ydf)));   % first we use fftshift to cancel 
                                 % effect of previous fftshift used
                                 % then we calculate inverse fourier 
                                 % transform ydt is our detected signal
ll=2000;                % used for ploted signal length
figure
plot(t(length(f)/2:length(f)/2+ll)*1000,xt(length(f)/2:length(f)/2+ll))
hold on
plot(t(length(f)/2:length(f)/2+ll)*1000,ydt(length(f)/2:length(f)/2+ll))
xlim([t(length(f)/2)*1000 t(length(f)/2+ll)*1000])  % setting x axis limits
title('orginal and detected message in a noisless channel-semi synchrozied')
xlabel('t(ms)')
ylabel('message signal')
legend('orginal message','detected message')
grid on
t2=downsample(t,n);
f2=downsample(f,n);
xt2=downsample(xt,n);
ydt2=downsample(ydt,n);
figure
ll=200;
plot(t2((length(f2)+1)/2:(length(f2)+1)/2+ll)*1000,xt2((length(f2)+1)/2:...
    (length(f2)+1)/2+ll))
hold on
plot(t2((length(f2)+1)/2:(length(f2)+1)/2+ll)*1000,ydt2((length(f2)+1)/2:...
    (length(f2)+1)/2+ll))
xlim([t2((length(f2)+1)/2)*1000 t2((length(f2)+1)/2+ll)*1000])  % setting x axis limits
title('orginal and detected message in a noisless channel downsampled -semi synchronized')
xlabel('t(ms)')
ylabel('message signal')
legend('orginal message','detected message ')
grid on
% envelope detector
% cm of all steps is in part 1 code
R=100e3; %1e2 ;  %1e3;
C=1e-6;
vin=[0 mt];
ti=1/fs*(0: 100);
v0=0;
j=1;
for i = 2 :length(mt)-1

    if(vin(i+1)-vin(i))<0
    if(vin(i)-vin(i-1))>=0
          v0=vin(i);
          j=1;
    end
    end
          aux(i)=v0*exp(-ti(j)/(R*C));
          j=j+1;
end
for i = 1 : length(mt)-1
    if(vin(i)-aux(i))>0
        aux(i)=vin(i);
    end
end
ll=2000;
figure
plot(t(length(f)/2:length(f)/2+ll)*1000,xt(length(f)/2:length(f)/2+ll))
hold on
plot(t(length(f)/2:length(f)/2+ll)*1000,aux(length(f)/2:length(f)/2+ll))
xlim([t(length(f)/2)*1000 t(length(f)/2+ll)*1000])
xlabel('t(ms)')
ylabel('signal')
title('detected signal with envelope detector R=100k C=1\mu- envelope detector')
legend('detected ','x(t)')
grid on
% downsample
t2=downsample(t,n);
f2=downsample(f,n);
xt2=downsample(xt,n);
aux2=downsample(aux,n);
figure
ll=200;
plot(t2((length(f2)+1)/2:(length(f2)+1)/2+ll)*1000,xt2((length(f2)+1)/2:...
    (length(f2)+1)/2+ll))
hold on
plot(t2((length(f2)+1)/2:(length(f2)+1)/2+ll)*1000,aux2((length(f2)+1)/2:...
    (length(f2)+1)/2+ll))
xlim([t2((length(f2)+1)/2)*1000 t2((length(f2)+1)/2+ll)*1000])
xlabel('t(ms)')
ylabel('signal')
title('detected signal with envelope detector R=100k C=1\mu downsampled - envelope detector')
legend('detected ','x(t)')
grid on
%% Q 3 synchronized detection
close all
clear
clc
[xt,fs] = audioread('signal.wav');           % read message x(t)
n=ceil(400000/fs);                           % use n for upsampling
fs=n*fs;                                     % new sampling rate
xt=interp(xt,n);
t = 0 : 1/fs : (1/fs)*(length(xt)-1) ;       % time vector
f = -fs/2 : fs/(length(xt)-1) : fs/2;        % frequency vector    
Xf = fft(xt);
mt=xt'.*cos(2*pi*75e3*t);
y1t=mt.*cos(2*pi*75e3*t);    % multiply recieved signal by carrier
Yf=fftshift(fft(y1t));        % calculate fourier transform 
                              % we use fftshift because we want 
                              % a symmetrical spectrum to be multipled
                              % by our ideal rect filter in frequency
                              % domain
myfilter=zeros(1,length(mt));      % this is gonna be our filter
                                   % k *[ 0 ... 0 ... 1 1 1 ... 0 0 0]
myindexlow=max(find(-20000>f)); % calculating our needed 
                                   % index to start 1 indexes considering
                                   % chonsen cut-off freq = 1340
myindextop=min(find(20000<f));  % calculating our needed 
                                   % index to stop 1 indexes considering
                                   % chonsen cut-off freq = 1340
myfilter(myindexlow:myindextop)=2;  % we bulid our ideal filter
Ydf=Yf.*myfilter;        % we filter our signal in freq domain 
ydt=real(ifft(fftshift(Ydf)));   % first we use fftshift to cancel 
                                 % effect of previous fftshift used
                                 % then we calculate inverse fourier 
                                 % transform ydt is our detected signal
ll=2000;                % used for ploted signal length
figure
plot(t(length(f)/2:length(f)/2+ll)*1000,xt(length(f)/2:length(f)/2+ll))
hold on
plot(t(length(f)/2:length(f)/2+ll)*1000,ydt(length(f)/2:length(f)/2+ll))
xlim([t(length(f)/2)*1000 t(length(f)/2+ll)*1000])  % setting x axis limits
title('orginal and detected message in a noisless channel')
xlabel('t(ms)')
ylabel('message signal')
legend('orginal message','detected message')
grid on
t2=downsample(t,n);
f2=downsample(f,n);
xt2=downsample(xt,n);
ydt2=downsample(ydt,n);
figure
ll=200;
plot(t2((length(f2)+1)/2:(length(f2)+1)/2+ll)*1000,xt2((length(f2)+1)/2:...
    (length(f2)+1)/2+ll))
hold on
plot(t2((length(f2)+1)/2:(length(f2)+1)/2+ll)*1000,ydt2((length(f2)+1)/2:...
    (length(f2)+1)/2+ll))
xlim([t2((length(f2)+1)/2)*1000 t2((length(f2)+1)/2+ll)*1000])  % setting x axis limits
title('orginal and detected message in a noisless channel downsampled ')
xlabel('t(ms)')
ylabel('message signal')
legend('orginal message','detected message ')
grid on
%% Q 4 Comparing DSB , AM sent power
clear 
clc
[xt,fs] = audioread('signal.wav');           % read message x(t)
n=ceil(400000/fs);                           % calculate n for upsampling
                                             % to increase sampling rate
                                             % to more than 75k*2=150 khz
fs=n*fs;                                     % new sampling rate
xt=interp(xt,n);                           % upsampled message
t = 0 : 1/fs : (1/fs)*(length(xt)-1) ;       % time vector
f = -fs/2 : fs/(length(xt)-1) : fs/2;        % frequency vector        
mtDSB=xt'.*cos(2*pi*75e3*t);                 % DSB Modulation
mtAM=(1+0.5*xt').*cos(2*pi*75e3*t);          % AM Modulation
MFDSB=fft(mtDSB)/length(t);                  % DSB - Spectrum
MFAM=fft(mtAM)/length(t);                    % AM - spectrum
P_DSB=sum(abs(MFDSB).^2)                     % DSB power
P_AM=sum(abs(MFAM).^2)                       % AM  powr
disp('P AM / P DSB =')
disp(P_AM/P_DSB)                             % AM power / DSB power ratio
%% Q 5 Calculating sent signal amplitude for SNR = 30 DB
clear 
clc
[xt,fs] = audioread('signal.wav');           % read message x(t)
n=ceil(400000/fs);                           % calculate n for upsampling
                                             % to increase sampling rate
                                             % to more than 75k*2=150 khz
fs=n*fs;                                     % new sampling rate
amp=750                                   % controlled signal amplitude
xt=amp*interp(xt,n);                         % upsampled message
t = 0 : 1/fs : (1/fs)*(length(xt)-1) ;       % time vector
f = -fs/2 : fs/(length(xt)-1) : fs/2;        % frequency vector    
mtDSB=xt'.*cos(2*pi*75e3*t);                 % DSB modulated signal
noise=randn(1,length(xt));                   % noise signal
MFDSB=fft(mtDSB)/length(t);                  % DSB modulated signal spectrum
Fn=fft(noise)/length(t);                     % noise signal spectrum
pdsb=sum(abs(MFDSB).^2)                      % signal power
pnoise=sum(abs(Fn).^2)                       % noise power
SNR=pdsb/pnoise                              % snr
SNR_DB=10*log(SNR)/log(10)                   % snr in DB












