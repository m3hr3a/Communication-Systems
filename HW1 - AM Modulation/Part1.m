%% COMMUNICATION SYSTEMS - MATLAB ASSIGNMENT 1
%% Instructor    : Dr. H. Behrouzi
%% Deadline      : 18 Aban 1397
%% Student name  : Mehrsa Pourya
%% Student ID    : 95101247
%% Part 1 AM Modulation
%% Q 1 Message's spectrum
close all
clear
clc
fs  = 3000 ;                                  % sampling frequency
t = -5 : 1/fs : 5 ;                           % time vector
f = -fs/2 : fs/(length(t)-1) : fs/2;          % frequency vector
xt = sin(2* pi*697*t)+ sin(2*pi*1336*t);      % message x(t)
Xf = fft(xt);
subplot(2,1,1)
% double sided spectrum
plot(f,fftshift(abs(Xf))/(length(t)-1))
title('|X(f)| Double Sided Spectrum')
xlabel('frequency')
ylabel('amplitude')
grid on
subplot(2,1,2)
% single sided spectrum
xfshifted=fftshift(abs(Xf));
plot(f((length(t)+1)/2-1:end),2*xfshifted((length(t)+1)/2-1:end)/(length(t)-1))
title('|X(f)| Single Sided Spectrum , L = 30001')
xlabel('frequency')
ylabel('amplitude')
xlim([0 1500])
grid on
% solve amplitude error with signal lentgh variation
figure
clear
clc
fs  = 3000 ;                                  % sampling frequency
t = -5 : 1/fs : 5 -(1/fs);                    % time vector
f = -fs/2 : fs/(length(t)-1) : fs/2;          % frequency vector
xt = sin(2* pi*697*t)+ sin(2*pi*1336*t);      % message x(t)
Xf = fft(xt);
subplot(2,1,1)
% double sided spectrum
plot(f,fftshift(abs(Xf))/(length(t)))
title('|X(f)| Double Sided Spectrum')
xlabel('frequency')
ylabel('amplitude')
grid on
subplot(2,1,2)
% single sided spectrum
xfshifted=fftshift(abs(Xf));
plot(f((length(t))/2:end),2*xfshifted((length(t))/2:end)/(length(t)))
title('|X(f)| Single Sided Spectrum , L=30000')
xlabel('frequency')
ylabel('amplitude')
grid on
xlim([0 1500])
%% Q 2 AWGN 
% no code needed
%% Q 3 Add AWGN to signal
close all
clear
clc
fs  = 120000 ;                                % sampling frequency
t = -0.5 : 1/fs : 0.5 -(1/fs);                % time vector     
xt = sin(2* pi*697*t)+ sin(2*pi*1336*t);      % message x(t)
mio=0.5;
mt=(1+mio*xt).*cos(2*pi*50000*t);             % Modulated message
noisedmt=mt+sqrt(7)*randn(1,length(mt));      % add noise to message
ll=200;                % used for ploted signal length
figure
plot(t(6e4:6e4+ll)*1000,noisedmt(6e4:6e4+ll))
hold on
plot(t(6e4:6e4+ll)*1000,mt(6e4:6e4+ll))
legend('noised modulated message','modulated message')
xlim([t(6e4)*1000 t(6e4+ll)*1000])
title('noise in time domain')
xlabel('t(ms)')
ylabel('amp')
grid on
%% Q 4 noised signal spectrum
close all
clear
clc
fs  = 120000 ;                                % sampling frequency
t = -0.5 : 1/fs : 0.5 -(1/fs);                % time vector     
xt = sin(2* pi*697*t)+ sin(2*pi*1336*t);      % message x(t)
mio=0.5;
mt=(1+mio*xt).*cos(2*pi*50000*t);             % Modulated message
noisedmt=mt+sqrt(7)*randn(1,length(mt));      % add noise to message
f = -fs/2 : fs/(length(t)-1) : fs/2;          % frequency vector
subplot(2,1,1)
MF=fft(mt);
plot(f,fftshift(abs(MF))/(length(t)))
grid on
title('2-sided , Amplitude of M(f)')
xlabel('frequency')
ylabel('|M(f)|')
subplot(2,1,2)
NM=fft(noisedmt);
plot(f,fftshift(abs(NM))/(length(t)))
grid on
title('2-sided , Amplitude of M(f)+AWGN')
xlabel('frequency')
ylabel('|M(f)|')
figure
subplot(2,1,1)
MF=fftshift(fft(mt));
plot(f(length(t)/2:end),2*(abs(MF(length(t)/2:end)))...
     /(length(t)))
grid on
title('One-sided , Amplitude of M(f)')
xlabel('frequency')
ylabel('|M(f)|')
xlim([0 60000])
subplot(2,1,2)
NM=fftshift(fft(noisedmt));
plot(f(length(t)/2:end),2*(abs(NM(length(t)/2:end)))...
     /(length(t)))
grid on
title('One-sided , Amplitude of M(f)+AWGN')
xlabel('frequency')
ylabel('|M(f)|')
xlim([0 60000])
%% Q 5 synchronized detector 
close all
clear
clc
fs  = 2000000;                                % sampling frequency
t = -0.5 : 1/fs : 0.5 -(1/fs);                % time vector     
xt = sin(2* pi*697*t)+ sin(2*pi*1336*t);      % message x(t)
f = -fs/2 : fs/(length(t)-1) : fs/2;  
mio=0.5;    % chosen modulation index
mt=(1+mio*xt).*cos(2*pi*50000*t);   % modulated message
% procedure of detection in a channel without noise
y1t=mt.*cos(2*pi*50000*t);    % multiply recieved signal by carrier
Yf=fftshift(fft(y1t));        % calculate fourier transform 
                              % we use fftshift because we want 
                              % a symmetrical spectrum to be multipled
                              % by our ideal rect filter in frequency
                              % domain
myfilter=zeros(1,length(mt));      % this is gonna be our filter
                                   % k *[ 0 ... 0 ... 1 1 1 ... 0 0 0]
my1337indexlow=max(find(-1339>f)); % calculating our needed 
                                   % index to start 1 indexes considering
                                   % chonsen cut-off freq = 1340
my1337indextop=min(find(1339<f));  % calculating our needed 
                                   % index to stop 1 indexes considering
                                   % chonsen cut-off freq = 1340
myfilter(my1337indexlow:my1337indextop)=4;  % we bulid our ideal filter
figure
% plot(f,abs(Yf)/length(mt))
% hold on
plot(f,myfilter)
xlim([-2000 2000])
title('My Chosen LPF for synchronized detector')
xlabel('frequency')
ylabel('amplitude')
grid on
Ydf=Yf.*myfilter;                % we filter our signal in freq domain 
% figure
% plot(f,abs(Ydf)/length(mt))
xlim([-2000 2000])               % setting x axis limits
ydt=real(ifft(fftshift(Ydf)));   % first we use fftshift to cancel 
                                 % effect of previous fftshift used
                                 % then we calculate inverse fourier 
                                 % transform ydt is our detected signal
ll=20000;                        % used for ploted signal length
figure
plot(t(1e6:1e6+ll)*1000,xt(1e6:1e6+ll))
hold on
plot(t(1e6:1e6+ll)*1000,ydt(1e6:1e6+ll))
xlim([t(1e6) t(1e6+ll)*1000])  % setting x axis limits
title('orginal and detected message in a noisless channel synchronized detector')
xlabel('t(ms)')
ylabel('message signal')
legend('orginal message','detected message' )
grid on
% plot spectrums in diffrent steps 
figure
subplot(4,1,1)
plot(f,fftshift(abs(fft(xt)))/length(f))
xlim([-1.5e5 1.5e5])   
title('orginal signal - synchronized detector')
xlabel('f(Hz)')
ylabel('amp')
grid on
ylim([0 1])
subplot(4,1,2)
plot(f,fftshift(abs(fft(mt)))/length(f))
title('modulated signal  - synchronized detector')
xlabel('f(Hz)')
ylabel('amp')
grid on
xlim([-1.5e5 1.5e5]) 
ylim([0 1])
subplot(4,1,3)
plot(f,fftshift(abs(fft(y1t)))/length(f))
title('output of mixer - synchronized detector  ')
xlabel('f(Hz)')
ylabel('amp')
grid on
xlim([-1.5e5 1.5e5])  
ylim([0 1])
subplot(4,1,4)
plot(f,fftshift(abs(fft(ydt)))/length(f))
title('detected signal (output of LPF)- synchronized detector ')
xlabel('f(Hz)')
ylabel('amp')
grid on
xlim([-1.5e5 1.5e5])
ylim([0 1])
noisedmt=mt+sqrt(70)*randn(1,length(mt));  % noised message
% procedure of detection in a channel with AWGN noise
y1t=noisedmt.*cos(2*pi*50000*t);  % multiply recieved signal by carrier
Yf=fftshift(fft(y1t));            % calculate fourier transform 
Ydf=Yf.*myfilter;                 % we filter our signal in freq domain
% figure 
% plot(f,abs(Ydf)/length(mt))
% xlim([-2000 2000])
ydt=real(ifft(fftshift(Ydf)));    % inversed ft of filtered signal
figure
plot(t(1e6:1e6+ll)*1000,xt(1e6:1e6+ll))
hold on
plot(t(1e6:1e6+ll)*1000,ydt(1e6:1e6+ll))
xlim([t(1e6) t(1e6+ll)*1000])     % setting x axis limits
title('orginal and detected message in channel with AWGN synchronized detector')
xlabel('t(ms)')
ylabel('message signal')
legend('sent message','orginal message')  
grid on
figure
% plot spectrums in diffrent steps 
subplot(4,1,1)
plot(f,fftshift(abs(fft(xt)))/length(f))
xlim([-1.5e5 1.5e5])   
title('orginal signal - synchronized detector')
xlabel('f(Hz)')
ylabel('amp')
grid on
ylim([0 1])
subplot(4,1,2)
plot(f,fftshift(abs(fft(noisedmt)))/length(f))
title('modulated signal + noise - synchronized detector')
xlabel('f(Hz)')
ylabel('amp')
grid on
xlim([-1.5e5 1.5e5]) 
ylim([0 1])
subplot(4,1,3)
plot(f,fftshift(abs(fft(y1t)))/length(f))
title('output of mixer - synchronized detector + AWGN ')
xlabel('f(Hz)')
ylabel('amp')
grid on
xlim([-1.5e5 1.5e5])  
ylim([0 1])
subplot(4,1,4)
plot(f,fftshift(abs(fft(ydt)))/length(f))
title('detected signal (output of LPF)- synchronized detector + AWGN')
xlabel('f(Hz)')
ylabel('amp')
grid on
xlim([-1.5e5 1.5e5])
ylim([0 1])
%% Q 6 effect of shift in phase and frecuency of reciever mixer carrier 
close all
clear
clc
fs  = 2000000;                                % sampling frequency
t = -0.5 : 1/fs : 0.5 -(1/fs);                % time vector     
xt = sin(2* pi*697*t)+ sin(2*pi*1336*t);      % message x(t)
f = -fs/2 : fs/(length(t)-1) : fs/2;  
mio=0.5;    % chosen modulation index
mt=(1+mio*xt).*cos(2*pi*50000*t);   % modulated message
% procedure of detection in a channel without noise
ALO=1;
st=ALO*cos(2*pi*(50000+500)*t+5/180*pi);   % shifted carrier
y1t=mt.*st;    % multiply recieved signal by shifted carrier
Yf=fftshift(fft(y1t));        % calculate fourier transform 
                              % we use fftshift because we want 
                              % a symmetrical spectrum to be multipled
                              % by our ideal rect filter in frequency
                              % domain
myfilter=zeros(1,length(mt));      % this is gonna be our filter
                                   % k *[ 0 ... 0 ... 1 1 1 ... 0 0 0]
my1337indexlow=max(find(-1840>f)); % calculating our needed 
                                   % index to start 1 indexes considering
                                   % chonsen cut-off freq = 1840
my1337indextop=min(find(1840<f));  % calculating our needed 
                                   % index to stop 1 indexes considering
                                   % chonsen cut-off freq = 1840
myfilter(my1337indexlow:my1337indextop)=8;  % we bulid our ideal filter
%figure
% plot(f,abs(Yf)/length(mt))
% hold on
plot(f,myfilter)
xlim([-2000 2000])
title('My Chosen LPF')
xlabel('frequency')
ylabel('amplitude')
Ydf=Yf.*myfilter;        % we filter our signal in freq domain 
figure
plot(f,abs(Ydf)/length(mt))
xlim([-2000 2000])       % setting x axis limits
xlabel('freq')
ylabel('amp')
title('Yd(f) spectrum')
ydt=real(ifft(fftshift(Ydf)));   % first we use fftshift to cancel 
                                 % effect of previous fftshift used
                                 % then we calculate inverse fourier 
                                 % transform ydt is our detected signal
ll=20000;                % used for ploted signal length
figure
plot(t(1e6:1e6+ll)*1000,xt(1e6:1e6+ll))
hold on
plot(t(1e6:1e6+ll)*1000,ydt(1e6:1e6+ll))
xlim([t(1e6) t(1e6+ll)*1000])  % setting x axis limits
title('orginal and detected message in a noiseless channel')
xlabel('t(ms)')
ylabel('message signal')
legend('orginal message','detected message')
noisedmt=mt+sqrt(7)*randn(1,length(mt));  % noised message
% procedure of detection in a channel with AWGN noise
y1t=noisedmt.*st;  % multiply recieved signal by carrier
Yf=fftshift(fft(y1t));            % calculate fourier transform 
Ydf=Yf.*myfilter;                 % we filter our signal in freq domain
% figure 
% plot(f,abs(Ydf)/length(mt))
% xlim([-2000 2000])
ydt=real(ifft(fftshift(Ydf)));   % inversed ft of filtered signal
figure
plot(t(1e6:1e6+ll)*1000,xt(1e6:1e6+ll))
hold on
plot(t(1e6:1e6+ll)*1000,ydt(1e6:1e6+ll))
xlim([t(1e6) t(1e6+ll)*1000])  % setting x axis limits
title('orginal and detected message in channel with AWGN')
xlabel('t(ms)')
ylabel('message signal')
legend('sent message','orginal message')  
%% Q 7 semi-synchronized detector
close all
clear
clc
fs  = 2000000 ;                               % sampling frequency
t = -0.5 : 1/fs : 0.5 -(1/fs);                % time vector     
xt = sin(2* pi*697*t)+ sin(2*pi*1336*t);      % message x(t)
f = -fs/2 : fs/(length(t)-1) : fs/2;  
mio=0.5;
% procedure of detection in a channel without noise
mt=(1+mio*xt).*cos(2*pi*50000*t);             % Modulated message
pos=find(mt>=0);             % hard limiter +
neg=find(mt<0);              % hard limiter -
HLout(pos)=1;
HLout(neg)=-1;
%plot output of hard limiter
plot(t(length(t)/2-250:length(t)/2+250), ...
    HLout(length(t)/2-250:length(t)/2+250))
hold on
plot(t(length(t)/2-250:length(t)/2+250),zeros(501,1))
xlim([t(length(t)/2-250) t(length(t)/2+250)])
title('Output of hard limiter')
xlabel('time (s)')
ylabel('amp')
figure
plot(f,fftshift(abs(fft(HLout)))/length(mt))
YH=fftshift((fft(HLout)));
title('Output of hard limiter spectrum')
xlabel('f')
ylabel('amp')
myfilter1=zeros(1,length(mt));      % this is gonna be our filter
                                   % k *[ 0 ... 0 ... 1 1 1 ... 0 0 0]
my1337indexlow=max(find(-55000>f)); % calculating our needed 
                                   % index to start 1 indexes considering
                                   % chonsen cut-off freq = 1340
my1337indextop=min(find(55000<f));  % calculating our needed 
                                   % index to stop 1 indexes considering
                                   % chonsen cut-off freq = 1340
myfilter1(my1337indexlow:my1337indextop)=0.5/(2*sin(pi/2)/pi);  ...
    % we bulid our ideal filter
YHf=YH.*myfilter1;        % we filter our signal in freq domain 
figure
plot(f,abs(YHf)/length(mt))
title('Out put of first filer spectrum')
xlabel('f')
ylabel('amp')
ydt1=real(ifft(fftshift(YHf)));
y1t=mt.*ydt1;    % multiply recieved signal by extracted carrier
Yf=fftshift(fft(y1t));        % calculate fourier transform 
                              % we use fftshift because we want 
                              % a symmetrical spectrum to be multipled
                              % by our ideal rect filter in frequency
                              % domain
myfilter=zeros(1,length(mt));      % this is gonna be our filter
                                   % k *[ 0 ... 0 ... 1 1 1 ... 0 0 0]
my1337indexlow=max(find(-1339>f)); % calculating our needed 
                                   % index to start 1 indexes considering
                                   % chonsen cut-off freq = 1340
my1337indextop=min(find(1339<f));  % calculating our needed 
                                   % index to stop 1 indexes considering
                                   % chonsen cut-off freq = 1340
myfilter(my1337indexlow:my1337indextop)=4;  % we bulid our ideal filter
%figure
% plot(f,abs(Yf)/length(mt))
% hold on
Ydf=Yf.*myfilter;        % we filter our signal in freq domain 
% figure
% plot(f,abs(Ydf)/length(mt))
%xlim([-2000 2000])       % setting x axis limits
ydt=real(ifft(fftshift(Ydf)));   % first we use fftshift to cancel 
                                 % effect of previous fftshift used
                                 % then we calculate inverse fourier 
                                 % transform ydt is our detected signal
ll=20000;                % used for ploted signal length
figure
plot(t(1e6:1e6+ll)*1000,xt(1e6:1e6+ll))
hold on
plot(t(1e6:1e6+ll)*1000,ydt(1e6:1e6+ll))
xlim([t(1e6) t(1e6+ll)*1000])  % setting x axis limits
title('orginal and detected message in a noisless channel')
xlabel('t(ms)')
ylabel('message signal')
legend('orginal message','detected message')

% procedure of detection in a channel with AWGN noise
noisedmt=mt+sqrt(70)*randn(1,length(mt));  % noised message          
pos=find(noisedmt>=0);             % hard limiter +
neg=find(noisedmt<0);              % hard limiter -
HLout(pos)=1;
HLout(neg)=-1;
YHf=YH.*myfilter1;        % we filter our signal in freq domain 
ydt1=real(ifft(fftshift(YHf)));
y1t=noisedmt.*ydt1;    % multiply recieved signal by extracted carrier
Yf=fftshift(fft(y1t));        % calculate fourier transform 
                              % we use fftshift because we want 
                              % a symmetrical spectrum to be multipled
                              % by our ideal rect filter in frequency
                              % domain
Ydf=Yf.*myfilter;        % we filter our signal in freq domain 
ydt=real(ifft(fftshift(Ydf)));   % first we use fftshift to cancel 
                                 % effect of previous fftshift used
                                 % then we calculate inverse fourier 
                                 % transform ydt is our detected signal
ll=20000;                % used for ploted signal length
figure
plot(t(1e6:1e6+ll)*1000,xt(1e6:1e6+ll))
hold on
plot(t(1e6:1e6+ll)*1000,ydt(1e6:1e6+ll))
xlim([t(1e6) t(1e6+ll)*1000])  % setting x axis limits
title('orginal and detected message in a channel with AWGN noisel')
xlabel('t(ms)')
ylabel('message signal')
legend('orginal message','detected message')
noisedmt=mt+sqrt(7)*randn(1,length(mt));  % noised message
%% Q 8 abs detector
close all
clear
clc
fs  = 2000000;                                % sampling frequency
t = -0.5 : 1/fs : 0.5 -(1/fs);                % time vector     
xt = sin(2* pi*697*t)+ sin(2*pi*1336*t);      % message x(t)
f = -fs/2 : fs/(length(t)-1) : fs/2;  
mio=0.5;    % chosen modulation index
mt=(1+mio*xt).*cos(2*pi*50000*t);   % modulated message
% procedure of detection in a channel without noise
y1t=abs(mt);    % abs xc(t)
Yf=fftshift(fft(y1t));        % calculate fourier transform 
                              % we use fftshift because we want 
                              % a symmetrical spectrum to be multipled
                              % by our ideal rect filter in frequency
                              % domain
myfilter=zeros(1,length(mt));      % this is gonna be our filter
                                   % k *[ 0 ... 0 ... 1 1 1 ... 0 0 0]
my1337indexlow=max(find(-1339>f)); % calculating our needed 
                                   % index to start 1 indexes considering
                                   % chonsen cut-off freq = 1340
my1337indextop=min(find(1339<f));  % calculating our needed 
                                   % index to stop 1 indexes considering
                                   % chonsen cut-off freq = 1340
myfilter(my1337indexlow:my1337indextop)=pi;  % we bulid our ideal filter
figure
% plot(f,abs(Yf)/length(mt))
% hold on
plot(f,myfilter)
xlim([-2000 2000])
title('My Chosen LPF')
xlabel('frequency')
ylabel('amplitude')
Ydf=Yf.*myfilter;        % we filter our signal in freq domain 
% figure
% plot(f,abs(Ydf)/length(mt))
xlim([-2000 2000])       % setting x axis limits
ydt=real(ifft(fftshift(Ydf)));   % first we use fftshift to cancel 
                                 % effect of previous fftshift used
                                 % then we calculate inverse fourier 
                                 % transform ydt is our detected signal
ll=20000;                % used for ploted signal length
figure
plot(t(1e6:1e6+ll)*1000,xt(1e6:1e6+ll))
hold on
plot(t(1e6:1e6+ll)*1000,ydt(1e6:1e6+ll))
xlim([t(1e6) t(1e6+ll)*1000])  % setting x axis limits
title('orginal and detected message in a noisless channel')
xlabel('t(ms)')
ylabel('message signal')
legend('orginal message','detected message')
noisedmt=mt+sqrt(7)*randn(1,length(mt));  % noised message
% procedure of detection in a channel with AWGN noise
y1t=abs(noisedmt);  % abs xc(t)
Yf=fftshift(fft(y1t));            % calculate fourier transform 
Ydf=Yf.*myfilter;                 % we filter our signal in freq domain
% figure 
% plot(f,abs(Ydf)/length(mt))
% xlim([-2000 2000])
ydt=real(ifft(fftshift(Ydf)));   % inversed ft of filtered signal
figure
plot(t(1e6:1e6+ll)*1000,xt(1e6:1e6+ll))
hold on
plot(t(1e6:1e6+ll)*1000,ydt(1e6:1e6+ll))
xlim([t(1e6) t(1e6+ll)*1000])  % setting x axis limits
title('orginal and detected message in channel with AWGN')
xlabel('t(ms)')
ylabel('message signal')
legend('sent message','orginal message') 
%% Q 9 Envelope detector
close all
clear 
clc
R=1e5; %1e2 ;  %1e3;
C=1e-6;
fs  = 2000000;                                % sampling frequency
t = -0.5 : 1/fs : 0.5 -(1/fs);                % time vector     
xt = sin(2* pi*697*t)+ sin(2*pi*1336*t);      % message x(t)
f = -fs/2 : fs/(length(t)-1) : fs/2;  
mio=0.5;                                      % chosen modulation index
mt=(1+mio*xt).*cos(2*pi*50000*t);             % modulated message
vin=[0 mt];                                   % we concat this zero because
                                              % we will need an index
                                              % before our first index
                                              % in our upcoming for loop
ti=1/fs*(0: 100);                             % time interval vector
                                              % we will use it to calculate
                                              % captivator dechargin 
                                              % voltage its length will
                                              % be limited by code steps
                                              % so 100 can be any bigger
                                              % number but it cannot be
                                              % less than 40 ~ f'=0 distanses                                       
for i = 2 :length(mt)-1

    if(vin(i+1)-vin(i))<0
    if(vin(i)-vin(i-1))>=0 
        % this condition implies change of sign of derivation of voltage
        % from + to - in other words when diode goes into its off area
          v0=vin(i);
          j=1;
    end
    end 
         % calculation of decharching values 
          aux(i)=v0*exp(-ti(j)/(R*C));
          j=j+1;
end
for i = 1 : length(mt)-1
    if(vin(i)-aux(i))>0
        % this condion implies that if input voltage is greater than output
        % voltage diode turns on and output = input
        aux(i)=vin(i); % aux is our output voltage
    end
end
ll=20000;
plot(t(1e6:1e6+ll)*1000,aux(1e6:1e6+ll))
xlim([t(1e6) t(1e6+ll)*1000])
xlabel('t(ms)')
ylabel('signal')
% title('detected signal with envelope detector R=100 C=1\mu')
title('detected signal with envelope detector R=100 C=1\mu')
%title('detected signal with envelope detector R=10 C=1\mu')
hold on
plot(t(1e6:1e6+ll)*1000,(xt(1e6:1e6+ll)))
legend('detected 1+\mu x(t)','x(t)')
grid on
figure
noisedmt=mt+sqrt(7)*randn(1,length(mt));  % noised message
R=1e5; %1e ;  %1e3;
C=1e-6;
vin=[0 noisedmt];
ti=1/fs*(0: 100);
% same cm s as above for here
for i = 2 :length(mt)-1

    if(vin(i+1)-vin(i))<0
    if(vin(i)-vin(i-1))>0
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
ll=20000;
plot(t(1e6:1e6+ll)*1000,aux(1e6:1e6+ll))
xlim([t(1e6) t(1e6+ll)*1000])
xlabel('t(ms)')
ylabel('signal')
% title('detected signal with envelope detector R=1k C=1 micro')
title('detected signal with envelope detector R=100k C=1\mu')
%title('detected signal with envelope detector R=10 C=1 micro')
hold on
plot(t(1e6:1e6+ll)*1000,(xt(1e6:1e6+ll)))
legend('detected 1+\mu x(t)','x(t)')
grid on







































