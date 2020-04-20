%% Matlab HW03 - Communication Systems
%% Dr.behrouzi 
%% Fall 2018 
%% Student : Mehrsa Pourya 95101247
%% Part 3
clear
clc
fs = 1e3 ; 
n = fs * 10; % number of samples = 10 s * 1khz
for realization = 1
x(:,realization)= wgn(n,1,0);
b=[1 -0.8 1.2];
a=1;
y(:,realization)= filter(b,a,x(:,realization));
end
subplot(2,1,1)
plot(1:n,x(:,1))
grid on 
xlabel('n')
title('wgn , system input')
subplot(2,1,2)
plot(1:n,y(:,1))
grid on 
xlabel('n')
title('system output')
%% mean 
meanofprocess=mean(y,1)
%% auto-correlation  and psd
figure
n1=20;
b3=xcorr(y,n1)/fs/10;
subplot(2,1,1)
stem(-n1:n1,b3)
xlabel('n')
title('system output auto-correlation')
grid on 
subplot(2,1,2)
plot(-n1:n1,fftshift(abs(fft(b3))))
grid on 
title('system output PSD')
xlabel('f')
%% 10 realization
clear
clc
fs = 1e3 ; 
n = fs * 10; % number of samples = 10 s * 1khz
b=[1 -0.8 1.2];
a=1;
for realization = 1 : 10
x(:,realization)= wgn(n,1,0);
y(:,realization)= filter(b,a,x(:,realization));
end
yavg=mean(y,2);
subplot(3,1,1)
plot(1:n,y(:,1))
grid on 
xlabel('n')
title('avg. of 10 realization system output')
%mean 
meanofprocess=mean(yavg,1)
% auto-correlation and psd
n1=20;
b3=xcorr(yavg,n1)/fs/10*10;
subplot(3,1,2)
stem(-n1:n1,b3)
xlabel('n')
title('avg. of 10 realization system output auto-correlation')
grid on 
subplot(3,1,3)
plot(-n1:n1,fftshift(abs(fft(b3))))
grid on 
title('avg. of 10 realization system output PSD')
xlabel('f')
%% 100 realization
clear
clc
fs = 1e3 ; 
n = fs * 10; % number of samples = 10 s * 1khz
b=[1 -0.8 1.2];
a=1;
for realization = 1 : 100
x(:,realization)= wgn(n,1,0);
y(:,realization)= filter(b,a,x(:,realization));
end
yavg=mean(y,2);
subplot(3,1,1)
plot(1:n,y(:,1))
grid on 
xlabel('n')
title('avg. of 100 realization system output')
%mean 
meanofprocess=mean(yavg,1)
% auto-correlation and psd
n1=20;
b3=xcorr(yavg,n1)/fs/10*100;
subplot(3,1,2)
stem(-n1:n1,b3)
xlabel('n')
title('avg. of 100 realization system output auto-correlation')
grid on 
subplot(3,1,3)
plot(-n1:n1,fftshift(abs(fft(b3))))
grid on 
title('avg. of 100 realization system output PSD')
xlabel('f')
%% 1000 realization
clear
clc
fs = 1e3 ; 
n = fs * 10; % number of samples = 10 s * 1khz
b=[1 -0.8 1.2];
a=1;
for realization = 1 : 1000
x(:,realization)= wgn(n,1,0);
y(:,realization)= filter(b,a,x(:,realization));
end
yavg=mean(y,2);
subplot(3,1,1)
plot(1:n,y(:,1))
grid on 
xlabel('n')
title('avg. of 1000 realization system output')
%mean 
meanofprocess=mean(yavg,1)
% auto-correlation and psd
n1=20;
b3=xcorr(yavg,n1)/fs/10*1000;
subplot(3,1,2)
stem(-n1:n1,b3)
xlabel('n')
title('avg. of 1000 realization system output auto-correlation')
grid on 
subplot(3,1,3)
plot(-n1:n1,fftshift(abs(fft(b3))))
grid on 
title('avg. of 1000 realization system output PSD')
xlabel('f')



















