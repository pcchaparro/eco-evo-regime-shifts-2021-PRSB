function dx = PopAlleeODE(t,x,bmax,do,d1,K,A,tau,sigma,Eend,eps)

% This function includes the ODE of a population subjected to an Allee
% effect. The ODE has 3 equations:
%       dx(1) --> population density
%       dx(2) --> mean trait
%       dx(3) --> environmental stress

dx=zeros(3,1);

N=x(1,1);
xmean=x(2,1);
E=x(3,1);

bmean = bmax*tau/(sigma^2+tau^2)^(1/2) * exp(-(xmean-E)^2/(2*(sigma^2+tau^2)));
dmean = do+d1*(sigma^2+xmean^2);

meanFit = bmean*(1-N/K)*(N-A)-dmean;
Fitgrad = - 2*d1*xmean - (bmax*tau*exp(-(E - xmean)^2/(2*(sigma^2 + tau^2)))*(A - N)*(K - N)*(2*E - 2*xmean))/(2*K*(sigma^2 + tau^2)^(3/2));

dx(1,1) = N*meanFit;

dx(2,1) = sigma^2*Fitgrad;

if eps>0 && Eend>E
    dx(3,1)= eps;
elseif eps<0 && Eend<E
    dx(3,1)= eps;
else
    dx(3,1)=0;
end
end