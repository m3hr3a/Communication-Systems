function [bandwich]=BesselBW(beta,fm);
% this function calculate bandwidth based on bessel coef
k=0;   % k counts next coefs less than 0.01
for i = 1 : 1e6
jn=besselj(i,beta);   % bessel coef 
if (abs(jn)<0.01)
    k=k+1;            % if jn(b)<0.01 then k = k +1 
end
if (abs(jn)>=0.01)
    k=0;              % when jn(b)>0.01 k will set to 0 again
end
if (k ==1000 )         % this means enough next coefs are less than 0.01
    nchosen=i-1000;    % we extract n used to calculate bandwidth
    break 
end
end
bandwich= 2*nchosen * fm;    % calculate bandwidth
end