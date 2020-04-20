%% Matlab HW03 - Communication Systems
%% Dr.behrouzi 
%% Fall 2018 
%% Student : Mehrsa Pourya 95101247
%% Part 2
%% Q1
% generate teta
clear 
clc
n=100000;
temp=rand(1,n);
indxteta0=find(temp<=1/6);
indxtetapi_4=find((1/6<temp & temp<=1/3));
indxtetapi_2=find((1/3<temp & temp<=2/3));
indxteta2pi_3=find((2/3<temp & temp<=8/9));
indxtetapi=find(8/9<temp);
teta(indxteta0)=0;
teta(indxtetapi_4)=pi/4;
teta(indxtetapi_2)=pi/2;
teta(indxteta2pi_3)=2*pi/3;
teta(indxtetapi)=pi;
% pmf of teta
pmfteta=[length(indxteta0) length(indxtetapi_4) length(indxtetapi_2) ...
    length(indxteta2pi_3) length(indxtetapi)]/n;
x=[0 pi/4 pi/2 2*pi/3 pi];
stem(x,pmfteta)
xlabel('\theta')
title('pmf of \theta')
grid on
% genrate X
figure
X=cos(0.2*pi+teta);
scatter(1:100,X(1:100),'fill')
xlabel('n')
ylabel('X')
title('100 realization of X')
grid on
% pmf of X
figure
pmfX=[length(find(X==cos(0.2*pi)))/n length(find(X==cos(0.2*pi+pi/4)))/n...
    length(find(X==cos(0.2*pi+pi/2)))/n length(find(X==cos(0.2*pi+2*pi/3)))/n...
    length(find(X==cos(0.2*pi+pi)))/n];
Xval=[cos(0.2*pi) cos(0.2*pi+pi/4) cos(0.2*pi+pi/2) cos(0.2*pi+2*pi/3)...
    cos(0.2*pi+pi)];
stem(Xval,pmfX)
grid on 
xlabel('X')
ylabel('Probability')
title('pmf of X')
%% Q 2
clear 
clc
for real =1 : 1000
n=100;
% generate theta
temp=rand(1,n);
indxteta0=find(temp<=1/6);
indxtetapi_4=find((1/6<temp & temp<=1/3));
indxtetapi_2=find((1/3<temp & temp<=2/3));
indxteta2pi_3=find((2/3<temp & temp<=8/9));
indxtetapi=find(8/9<temp);
teta(indxteta0)=0;
teta(indxtetapi_4)=pi/4;
teta(indxtetapi_2)=pi/2;
teta(indxteta2pi_3)=2*pi/3;
teta(indxtetapi)=pi;
nval=1:1:n; 
% process x(n)
X(real,:)=cos(0.2*pi*nval+teta);
    if( real <4 ) 
        subplot(3,1,real)
        plot(X(real,:),'LineWidth',2)
        xlabel('n')
        title('realization of process X(n)=cos(0.2\pi+\theta)')
        grid on
    end
end
% mean of X(n) 
figure
plot(mean(X,1))
xlabel('n')
title('mean of proccess X(n)')
grid on
syms n1 
% theoritical autocorrelation
y=simplify(cos(0.2*pi*n1)*1/6+cos(0.2*pi*n1+pi/4)*...
    1/6+cos(0.2*pi*n1+pi/2)/3+cos(0.2*pi*n1+2*pi/3)*2/9+cos(0.2*pi*n1+pi)/9);
hold on
plot(subs(y,1:100))
legend('simulation','theoritical')
% autocoreelation
for t = 1 : n
    for s = 1 :n
        R(t,s)=mean(X(:,t).*X(:,s));
    end
end
figure
surface(R)
colorbar
title('autocorrelation of process R[X(n)X(m)]')
xlabel('n')
ylabel('m')
for i2 = 1 : n
    avCor(i2)=mean(diag(R,-(n)/2+i2));            % autocorrelation fucnc;
end
figure
plot(-49:50,avCor,'LineWidth',2)
title('autocorrelation of process X(n)')
xlabel('n')
grid on
%% Q3
clear
clc
for real = 1 : 1000;
    rng();
    uvar=2;
    % generate u[n]
    U(real,:)=sqrt(uvar)*randn(1,1000);
    n=-499:1:500;
    X(real,:)=((0.5).^abs(n)).*U(real,:);   % process x[n]
end
for i = 1:1000
    Si(:,i)=((abs(fft(X(i,:)))).^2)/1000;     % cal spectrum for each real.
end
Sir=mean(Si,2);        % average PSD s of each sample func to cal PSD of process
figure
plot(-499:500,fftshift(Sir))
grid on 
title('2-3 PSD of process')
%% Q 4 
clear
clc
for real = 1 :1e3
% generate U(n)
temp=rand(100,1);
U1=find(temp>=0.5);
Um1=find(temp<0.5);
U(U1)=1;
U(Um1)=-1;
% generate X(n)
for n = 1 : 100
X(real,n)=sum(U(1,1:n));
end
    % plot 5 realization
    if ( real < 6)
        subplot(5,1,real)
        plot(1:100,X(real,:),'LineWidth',2)
        xlabel('n')
        title(['part2 - Q4 realization',num2str(real)])
        grid on 
    end
end
% mean 
figure
plot(1:100,mean(X,1))
title('mean of process X[n] part2 Q4')
xlabel('n')
% covariance matrix
figure
surface(cov(X))
R=cov(X);
title('covariance of process X[n] part2 Q4')
xlabel('n')
ylabel('m')
colorbar













