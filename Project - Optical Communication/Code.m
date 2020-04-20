%% Communication Systems Final Project
%% Instructor : Dr.H.Behroozi 
%% Project : OC , Optical Communications
%% Students                  Student ID
%% Mehrsa Pourya             95101247
%% Mahsa Siadati             95101717
%% Sahar Sattari             95101677
%% Winter 2019 
%% Sharif University of Technology - Electrical Engineering Department
%% 5.1 ook coding
clear
clc
close all
c=0;
for Kb = [0 1 5 10]  % given Kb
    c=c+1;  % counter
    figure 
    for Ks =1:1:60;
        % optimized threshold
        if ( Kb~=0)
            kth =Ks/log(1+Ks/Kb);  
        else
            kth=1;
        end
        ksmaller = 0 : 1 : ceil(kth-1);   % k<vth
        kbigger = ceil(kth) : 1 : 100;    % k >= vth
        % error Probability
        perror(c,Ks) = 0.5*sum((Kb.^kbigger)*exp(-Kb)./factorial(kbigger))+...
        0.5*sum(((Kb+Ks).^ksmaller)*exp(-(Kb+Ks))./factorial(ksmaller));
    end 
    
    ks=1:Ks;
    const(ks)=10^-4;  % desied error
    plot(ks,perror(c,:))
    hold on 
    plot(ks,const)
    desiredks=min(find(perror(c,:)<10^-4));  % chosen Ks
    hold on
    scatter(ks(desiredks),perror(c,desiredks),'fill')
    title(['for p(error)=0.0001 chosen Ks=',[num2str(ks(desiredks))],...
        ' is in ook with Kb=',num2str(Kb)])
    grid on 
    legend('error','0.0001','chosen Ks')
    xlabel('Ks')
    ylabel('perror')
end
%% simulation of ook
%% note : please first run 5.1 : ook coding
c=0;
for Kb=[0,1,5,10]  % given Kb
    c= c+1;
    messagelength=1000;                   % simulated message length
    message = rand(messagelength,1)> 0.5; % simulated message P(1)=P(0)=0.5
    for Ks=1: 60;
        % send message
    for i = 1 : messagelength
        if (message(i)==1)
            codedmessage(i)=poissrnd(Kb+Ks);   % on source
        else
            codedmessage(i)=poissrnd(Kb);      % off sourse
        end
    end
    % optimed threshold
    if(Kb==0)
        xth=1;
    else
        xth=Ks/log(1+Ks/Kb);
    end
    % detection of message using optimed threshold
    decodedm=codedmessage >= xth;
    % simulated error
    perrors(c,Ks)= length(find(((message'-decodedm)~=0)))/messagelength;
    end
    figure
    plot(1:Ks,perrors(c,:))   % plot simulated error
    hold on
    plot(1:Ks,perror(c,:))    % plot theoritical error
    grid on
    legend('simulation','theory')
    title(['Perror Both theoritical and simulation''s result based on Ks'...
        ,' ,Kb=',num2str(Kb)])
    xlabel('Ks')
    ylabel('Perror')
end
%% 5.2 Manchester coding
clear 
clc
close all
Ks=1:60;   % range of Ks
c=0;       % counter in for loop
for kb=[0 1 5 10]   % given Kb
    c=c+1;
    for Ks = 1:60
        % calculation of sums reported in report
        sigma2=0;
        sigma3=0;
        for k = 0 : 100
            sigma1=0;
            for koff = k+1 : 100
                sigma1=sigma1+((((kb/2)^koff)*exp(-kb/2))/factorial(koff));
            end
            sigma2=sigma2+sigma1...
                *((((kb+Ks)/2)^k)*exp((-kb-Ks)/2)/factorial(k));
            sigma3=sigma3+((((kb/2)^k)*exp(-kb/2))/factorial(k))...
                *((((kb+Ks)/2)^k)*exp((-kb-Ks)/2)/factorial(k));
        end
        % error Probability
        perror(c,Ks)=sigma2+0.5*sigma3;
    end
        figure
        ks=1:Ks;
        const(ks) = 1e-4;   % desired error
        plot(ks,perror(c,:))
        hold on 
        plot(ks , const);
        desiredks = min(find(perror(c,:)<1e-4));  % find sutibale Ks for
                                                  % given error
        hold on;
        
        scatter(ks(desiredks),perror(c,desiredks),'fill')
        title(['for p(error)=0.0001 chosen Ks='...
            ,num2str(ks(desiredks))...
            ,' is in Manechester coding with kb = ',num2str(kb)])
        grid on
        legend('error' , '0.0001','chosen Ks');
        xlabel('Ks');
        ylabel('perror');
end
%% Manchester Simulation  
%% note : please run part 5.2 : Manchester coding first
messagelength = 1000;    % length of simulated message
message = rand(messagelength,1) <0.5;   % simulated message P(1)=P(0)=0.5
t = 0:0.01:messagelength-0.01;          % t vector using for message plot
signal = zeros(size(t,2),1);
% we upsample message to plot it better
for n=1:messagelength
    for i = 100*(n-1)+1 : 100*n
    signal(i,1) = message(n);
    end
end
plot(t,signal);
xlim([1,20])
ylim([-5,5]);
grid
title('signal');
% code message using Manchester Method
manch = zeros(size(t,2),1);
% code message with logic explained in report
for n=1:messagelength
    for i = 100*(n-1)+1 : 100*n-50
        manch(i,1) = message(n);
    end
    for i = 100*n -49 : 100*n
        if message(n) ==1 
            manch(i,1) = 0;
        end
        if message(n) == 0
            manch(i,1) = 1;
        end
    end
end
manchcodedmessag=downsample(manch,50);  %our final manchester coded message
figure
plot(t,manch);
xlim([1,20])
ylim([-5,5]);
grid
title('manchester encoded signal');
c=0;
% here we send and detect message
for Kb=[0,1,5,10]
    c= c+1;
    for Ks=1: 60;
    for i = 1 : 2*messagelength
        if (manchcodedmessag(i)==1)
            codedmessage(i)=poissrnd(Kb/2+Ks/2);  % on source
        else
            codedmessage(i)=poissrnd(Kb/2);       % off source
        end
    end
    % detection [a b] a>=b -> 1 was sent , a<b -> 0 was sent
    for j = 2 : 2 : 2*messagelength;
        if ( codedmessage(j-1)>=codedmessage(j))
            decodedm(j/2)=1;
        else
            decodedm(j/2)=0;
        end
    end
    % simulated error
    perrors(c,Ks)= length(find(((message'-decodedm)~=0)))...
        /messagelength;
    end
    figure
    plot(1:Ks,perrors(c,:))  % plot simulated error
    hold on
    plot(1:Ks,perror(c,:))   % plot theoritical error
    grid on
    legend('simulation','theory')
    title(['Perror Both theoritical and simulation''s result based on Ks'...
        ,' ,Kb=',num2str(Kb)])
    xlabel('Ks')
    ylabel('Perror')
end
%% 5.3 dispersion effect
clc
clear
close all
% given parameters
Ks=[5,10,50];
deltatb=[0.1 0.2 0.3];
Kb=0.01;
% calculations of Perror
for i =1:3
    for j=1:3
        % optimized threshould
        vth(i*j)=(2*deltatb(j)*Ks(i)-2*Ks(i))/(log(Kb+deltatb(j)*Ks(i))...
            -log(Ks(i)+Kb-deltatb(j)*Ks(i))+log(Kb)-log(Kb+Ks(i)));

   vsmaller=0:1:ceil(vth(i*j)-1);   % smaller indexes in sigma
   vbigger=ceil(vth(i*j)):1:100;    % bigger indexes in sigma
   
      % calculation of sums reported in report
      
       sum1=sum(((Kb+Ks(i)).^vsmaller*exp(-Ks(i)-Kb))...
           ./factorial(vsmaller));
       
       sum2=sum(((Kb+Ks(i)-deltatb(j)*Ks(i)).^vsmaller...
           *exp(-Ks(i)-Kb+deltatb(j)*Ks(i)))./factorial(vsmaller));
  
       sum3=sum(((Kb+deltatb(j)*Ks(i)).^vbigger*...
           exp(-Kb-deltatb(j)*Ks(i)))./factorial(vbigger));
 
       sum4=sum(((Kb).^vbigger*exp(-Kb))./factorial(vbigger));

   
   perror(i*j)=0.25*(sum1+sum2+sum3+sum4);  % error Probability
   
   string=[' Perror for Ks=',num2str(Ks(i)),' , Delta/Tb=',...
       num2str(deltatb(j)),' :'];
       disp(string)
   disp(perror(i*j));
   
end
end



















