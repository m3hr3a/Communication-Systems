%% Matlab HW03 - Communication Systems
%% Dr.behrouzi 
%% Fall 2018 
%% Student : Mehrsa Pourya 95101247
%% Part 1 
%% Q1
clc
clear
pv=[0.1 0.3 0.5 0.7 0.9];      % p given values
for k = 1 : length(pv)         % loop for differnet p values
    for j = 1 : 100            % loop for diffrent sample functions
    rng();                     % reproducity of uniform random number
    temp=rand(1,1);            % used to simulate probability
    if ( temp > 0.5 )          % p(x[1] = 1) = 0.5
    x(1)=1;
    end
    if ( temp <= 0.5 )         % p(x[1] = -1) = 0.5
    x(1)=-1;
    end
    p = pv(k);                 % p value used in each kth iteration
    n=99;                      % n+1=100 is length of each sample function
        for i = 1 : n          % loop for each sample function
        rng();                 % reproducity of uniform random number              
        temp(i)=rand(1,1);     % used to simulate probability
        % all conditions below are based on conditional prob. gieven in HW
        if (x(i)==1)           
            if ( temp(i) > p ) % p{x[n]=1 | x[n-1]=1 } = 1 - p 
                x(i+1)=1;
            end
            if ( temp(i) <= p ) 
                x(i+1)=-1;     % p{x[n]=-1 | x[n-1]=1 } = p 
            end
        end
        if (x(i)==-1)
            if ( temp(i) > p ) 
                x(i+1)=-1;    % p{x[n]=-1 | x[n-1]=-1 } = 1 - p 
            end
            if ( temp(i) <= p ) 
                x(i+1)=1;     % p{x[n]=1 | x[n-1]=-1 } = p 
            end
        end
end
mysamples(k,j,:)=x;           % samples(p=pi,iteration,n)
if ( j < 6)                   % plot 5 rsample func.
subplot(5,1,j)
plot(1:n+1,x,'LineWidth',2)
xlim([1  n+1])
ylim([-1.5 1.5])
xlabel('n')
ylabel('X[n]')
title(['sample function ', num2str(j),' for p =',num2str(p)])
end
end
if ( k < 5)
figure
end
end
%% Q 2 
%% NOTE : please run Q1 first
% mean
for k = 1 :5   % for 5 different p values
samplemeans(k,:)=mean(mysamples(k,:,:),2);  % calculating mean of each pro.
subplot(5,1,k)
plot(1:n+1,samplemeans(k,:),'LineWidth',1.5)
title(['mean of process with p =',num2str(pv(k))])
xlim([1  n+1])
ylim([-1.5 1.5])
grid on
end
% auto corrolation
% method 1 , using definition
for k=1 : 5
    for t = 1 : n+1
        for s = 1 :n+1
        R(t,s)=mean(squeeze(mysamples(k,:,t)).*squeeze(mysamples(k,:,s)));
        end
    end
figure
surface(R)
colorbar
xlabel('t')
ylabel('s')
title(['autocorrelation of process with p =',num2str(pv(k))])
%method 2 cal. covariance matrix 
for i2 = 1 : n+1
    avCor(k,i2)=mean(diag(R,-(n+1)/2+i2));         % autocorrelation fucnc;
end
C=cov(squeeze(mysamples(k,:,:)));
figure
surface(C)
colorbar
title(['covariance matrix of process with p =',num2str(pv(k))])
end
figure
for k = 1 : 5
    subplot(5,1,k)
    plot(-(n-1)/2:(n+1)/2,avCor(k,:),'LineWidth',2)
    title(['R(m) with p =',num2str(pv(k))])
    xlabel('m')
end
%% theortical autocorrelation
for k = 1 : 5
    subplot(5,1,k)
    Rt(k,:)=((1-pv(k))*2-1).^abs(-(n-1)/2:(n+1)/2);
    plot(-(n-1)/2:(n+1)/2,Rt(k,:),'LineWidth',2)
    title(['theoritical R(m) with p =',num2str(pv(k))])
    xlabel('m')
end
%% Q3 
% no code neede 
%% Q 4 PSD ,Wiener–Khinchin theorom , please run previous parts first
for k = 1 : 5
for i = 1:j
    Si(:,i)=((abs(fft(squeeze(mysamples(k,i,:))))).^2)/(n+1);
end
Sir=mean(Si,2);        % average PSD s of each sample func to cal PSD of process
subplot(5,1,k)
plot(-(n-1)/2:(n+1)/2,fftshift(Sir)/max(Sir))
hold on
Sf=(abs(fft(squeeze(avCor(k,:)))).^2);  % Wiener–Khinchin PSD
plot(-(n-1)/2:(n+1)/2,fftshift(Sf)/max(Sf))
xlabel('f')
title(['Normilized Power spectrum density of proccess for p=',num2str(pv(k))])
legend('definition','using Wiener–Khinchin')
grid on
end
%%

    





























