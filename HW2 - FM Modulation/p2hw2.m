%% HW02 - Communication Systems Matlab Assignment 2
%% Cousrse Instructor : Dr.behrouzi
%% Student Name     Student ID
%% Mehrsa Pourya    95101247
%% Deadline : 16/9/97 
%% Dec2018
%% Problem 2
%% Q1 
% no code needed
%% Q2 power noise based on fm
clear 
clc
close all
fs=9e5;                                     % sampling frequency
fc=100e3;                                   % carrier frequency 
c=1;                                        % counter
dt=1/fs;
t = -0.5: dt : 0.5-dt ;                         % time vector
f = -fs/2 : fs/(length(t)-1) : fs/2;        % frequency vector
Am=1;                                       % message amp.
B=10;                                       % beta
for fm =50:50:1000                          % message freq
fd=B*fm/Am;                                 % calculate fdelta
xt=Am*cos(2*pi*fm*t);                       % message
intxt=cumtrapz(xt*dt);                      % integrated message
xct=cos(2*pi*fc*t+2*pi*fd*intxt);           % modulated message
hxct1=[0  xct(1:end-1)];                    % used to calculate derivative
hxct2=[xct(2:end) 0 ];                      % used to calculate derivative
xctprime=(hxct2-hxct1)/(2*dt);              % d/dt block
y2=abs(xctprime);                           % abs block
% filter lpf
g=1/(4*fd);                                 % filter gain
Y2=fftshift(fft(y2));                       % calculate fourier transform 
myfilter=zeros(1,length(y2));      % this is gonna be our filter
                                   % k *[ 0 ... 0 ... 1 1 1 ... 0 0 0]
my1337indexlow=max(find(-fm-100>f)); % calculating our needed 
                                   % index to start 1 indexes considering
                                   % chonsen cut-off freq = 1340
my1337indextop=min(find(fm+100<f));  % calculating our needed 
                                   % index to stop 1 indexes considering
                                   % chonsen cut-off freq = 1340
myfilter(my1337indexlow:my1337indextop)=g;  % we bulid our ideal filter
filteredml=Y2.*myfilter;                  % we filter our signal in freq domain 
yd=real(ifft(fftshift(filteredml))); % first we use fftshift to cancel 
                                         % effect of previous fftshift used
output=yd-mean(yd);
output=output/max(output);               % normilizing
noisepw(c)=sum((abs(fft(xt-output))/length(xt)).^2);  % noise power
c=c+1;
end
fv=50:50:1000  ;
plot(fv,noisepw)
grid on 
title('noise power based on fm')
xlabel('fm(Hz)')
ylabel('noise power')
%% Q3 power noise based on beta
clear 
clc
close all
fs=5e5;                                 % sampling frequency
fc=100e3;                                   % carrier frequency 
c=1;                                        % counter
dt=1/fs;
t = -0.5: dt : 0.5-dt ;                     % time vector
f = linspace(-fs/2,fs/2,length(t));         % frequency vector
Am=1;                                       % message amp.
fm=1e3;                                      % message freq
for B=1:0.5:20                   
fd=B*fm/Am;                                 % calculate fdelta
xt=Am*cos(2*pi*fm*t);                       % message
intxt=cumtrapz(xt*dt);                      % integrated message
xct=cos(2*pi*fc*t+2*pi*fd*intxt);           % modulated message
hxct1=[0  xct(1:end-1)];                    % used to calculate derivative
hxct2=[xct(2:end) 0 ];                      % used to calculate derivative
xctprime=(hxct2-hxct1)/(2*dt);              % d/dt block
y2=abs(xctprime);                           % abs block
% filter lpf
g=1/(4*fd);                                 % filter gain
Y2=fftshift(fft(y2));                       % calculate fourier transform 
myfilter=zeros(1,length(y2));      % this is gonna be our filter
                                   % k *[ 0 ... 0 ... 1 1 1 ... 0 0 0]
my1337indexlow=max(find(-fm-100>f)); % calculating our needed 
                                   % index to start 1 indexes considering
                                   % chonsen cut-off freq = 1340
my1337indextop=min(find(fm+100<f));  % calculating our needed 
                                   % index to stop 1 indexes considering
                                   % chonsen cut-off freq = 1340
myfilter(my1337indexlow:my1337indextop)=g;  % we bulid our ideal filter
filteredml=Y2.*myfilter;                  % we filter our signal in freq domain 
yd=real(ifft(fftshift(filteredml))); % first we use fftshift to cancel 
                                         % effect of previous fftshift used

output=yd-mean(yd);
output=output/max(output);               % normalizing
noisepw(c)=sum((abs(fft(xt-output))/length(xt)).^2);   % noise power
c=c+1;
wantedi(c,:)=xt;
wantedo(c,:)=output;
end
bv=1 : 0.5 :20;
plot(bv,noisepw)
grid on 
title('noise power based on \beta')
xlabel('\beta')
ylabel('noise power')
figure
a=wantedi(4,:);
b=wantedo(4,:);
plot(t,a)
hold on 
plot(t,b)
xlabel('t(s)')
title('FM Modulation Detector for \beta = 2.5')
xlim([0 1e-2])
ylim([-1.25 1.5])
legend('orginal message','detected message','Location','best')
grid on

                                         
                                         
                                         