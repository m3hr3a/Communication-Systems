%% HW02 - Communication Systems Matlab Assignment 2
%% Cousrse Instructor : Dr.behrouzi
%% Student Name     Student ID
%% Mehrsa Pourya    95101247
%% Deadline : 16/9/97 
%% Dec2018
%% Problem 1 
%% Q 1 BesselBW Function
% wanted fuction is attached 
%% Q 2 Bessel BW and Carson approxiamation comparison
clc
clear
close
beta=0:100;                   % beta values , you can change upper bound 
fm=1e3;                       % fm=1kHz
for i = 1 : length(beta)
bBW(i)=BesselBW(beta(i),fm);  % Bessel bandwidth
end
cBW=2*(beta+1)*fm;            % Carson bandwidth
figure
plot(beta,bBW)
hold on
plot(beta,cBW)
legend('BesselBW','CarsonBW','Location','Best')
xlabel('\beta')
ylabel('bandwidth')
title('Bessel BW and Carson approxiamation comparison for a tone')
grid on
%% Q 3
clc
clear
fs=2000000;                                  % sampling frequency
fc=100e3;                                    % carrier frequency
dt=1/fs;                                   
t = -0.1: dt : 0.1-dt ;                     % time vector
f = -fs/2 : fs/(length(t)-1) : fs/2;        % frequency vector
Am=1;
fm=1000;
xt=Am*cos(2*pi*fm*t);                       % primary message
noise=2*rand(1,length(t))-1;                % uniform noise
noisedmessage=xt+noise;                     % FM modulation input
% filtering
nmft=fftshift(fft(noisedmessage)); % calculate fourier transform 
myfilter=zeros(1,length(xt));      % this is gonna be our filter
                                   % k *[ 0 ... 0 ... 1 1 1 ... 0 0 0]
my1337indexlow=max(find(-1000>f)); % calculating our needed 
                                   % index to start 1 indexes considering
                                   % chonsen cut-off freq = 1340
my1337indextop=min(find(1000<f));  % calculating our needed 
                                   % index to stop 1 indexes considering
                                   % chonsen cut-off freq = 1340
myfilter(my1337indexlow:my1337indextop)=1;  % we bulid our ideal filter
filteredm=nmft.*myfilter;                  % we filter our signal in freq domain 
filmess=real(ifft(fftshift(filteredm))); % first we use fftshift to cancel 
                                         % effect of previous fftshift used
filmess=filmess/max(filmess);   % normalizing
DD=1:0.1:10;                    % you can change d upper bound to 20
for c = 1: length(DD)
D=DD(c);
fd=fm*D;
intfilmess=cumtrapz(filmess.*dt);
xct=cos(2*pi*fc*t+2*pi*fd*intfilmess);
l=0.5*0.01;             %our limit max{abs(fft(x(t))}=0.5 , 0.5*0.01=0.005
yy=abs(fftshift(fft(xct)))/length(xct);
k=0;   % k counts next coefs less than 0.01max{abs(fft(x(t))}
p=0;   % detrmining low freq 
q=0;   % helper counter
for i = length(yy)/2 : length(yy)
if (abs(yy(i))<l)
    k=k+1; % if jn(b)<0.01 then k = k +1 
    p=0;
end
if (abs(yy(i))>=l)
    k=0; % when jn(b)>0.01 k will set to 0 again
    p=p+1;
end
if (k ==50000 )         % this means enough next coefs are less than l
    nchosenMAX=i-50000;    % we extract n used to calculate bandwidth
    break 
end
if (p ==1 && q==0)   
    q=q+1;
    % this means enough next coefs are less than l
    nchosenMIN=i;    % we extract n used to calculate bandwidth
end
end
flow=f(nchosenMIN);
if (f(nchosenMIN)<0)
    nchosenMIN=0;     % flow cannot be negative 
end
alf(c)=f(nchosenMAX)-flow;   % simulated bandwidth
carf(c)=2*(D+1)*1e3;         % carson appr.
end                                        
figure
plot(DD,alf)
hold on 
plot(DD,carf)
legend('Simulated BW','Carson BW')
xlabel('D')
ylabel('BW(Hz)')
grid on
title('Simulated BW and Carson approxiamation comparison for multi-tone signal')











