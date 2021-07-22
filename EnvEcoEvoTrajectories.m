%
% EnvEcoEvoTrajectories
% Author: Catalina Chaparro
% Last modification: 20/5/2021
%
% This rutine makes figure 3b and 5c of the paper entitled "Fast
% environmental change and eco- evolutionary feedbacks can drive
% regime shifts in ecosystems before tipping points are crossed"

clear

%Parameters----------------------------------------------------

K     = 10;         % Carrying capacity
b     = 1;          % Maximum birth rate
A     = -1;         % Allee threshod
do    = .5;         % Background mortality rate
d13   = 1;          % Mortality rate due to the fecundity-survival tradeoff in fig 3
d15   = .1;         % Mortality rate due to the fecundity-survival tradeoff in fig 5
tau   = 1;          % Degree of especialization
sigma = 0.05;       % Genetic variation
Eend  = 2;          % Maximum environmental stress

eps   =[.02;        % Rate of increase of environmental stress S1 in fig 3b
        .04;        % Rate of increase of environmental stress S2 in fig 3b
        .02;        % Rate of increase of environmental stress S1 in fig 5c
        .04];       % Rate of increase of environmental stress S1 in fig 5c
    

%---------------------------------------------------------------
% Figure 3b

%Calculate the eco-evo-environmental trajectory S1
x0=[10 0 0];
fod=@(t,x) PopAlleeODE(t,x,b,do,d13,K,A,tau,sigma,Eend,eps(1,1));
[tS1,x]=ode23s(fod,[0 500],x0);

NS1     = x(:,1);
xmeanS1 = x(:,2);
ES1     = x(:,3);

%Calculate the eco-evo-environmental trajectory S2
x0=[10 0 0];
fod=@(t,x) PopAlleeODE(t,x,b,do,d13,K,A,tau,sigma,Eend,eps(2,1));
[tS2,x]=ode23s(fod,[0 500],x0);
svecx=size(x);

NS2     = x(:,1);
xmeanS2 = x(:,2);
ES2     = x(:,3);


figure
suptitle('Figure 3b')
subplot(3,1,1)
hold on
plot(tS1,ES1)
plot(tS2,ES2)
ylim([0 2.5])
ylabel('Environmental stress (E)')
legend('S1','S2')
subplot(3,1,2)
hold on
plot(tS1,NS1)
plot(tS2,NS2)
ylabel('Abundance (N)')
ylim([-.8 10])
subplot(3,1,3)
hold on
plot(tS1,xmeanS1)
plot(tS2,xmeanS2)
ylim([-0.05 .5])
ylabel('Mean trait (x)')
xlabel('Time')

%---------------------------------------------------------------
% Figure 5c

%Calculate the eco-evo-environmental trajectory S1
x0=[10 0 0];
fod=@(t,x) PopAlleeODE(t,x,b,do,d15,K,A,tau,sigma,Eend,eps(3,1));
[tS1,x]=ode23s(fod,[0 3000],x0);

NS1     = x(:,1);
xmeanS1 = x(:,2);
ES1     = x(:,3);

%Calculate the eco-evo-environmental trajectory S2
x0=[10 0 0];
fod=@(t,x) PopAlleeODE(t,x,b,do,d15,K,A,tau,sigma,Eend,eps(4,1));
[tS2,x]=ode23s(fod,[0 3000],x0);
svecx=size(x);

NS2     = x(:,1);
xmeanS2 = x(:,2);
ES2     = x(:,3);


figure
suptitle('Figure 5c')
subplot(3,1,1)
hold on
plot(tS1,ES1)
plot(tS2,ES2)
ylim([0 2.5])
ylabel('Environmental stress (E)')
legend('S1','S2')
subplot(3,1,2)
hold on
plot(tS1,NS1)
plot(tS2,NS2)
ylabel('Abundance (N)')
ylim([-.8 10])
subplot(3,1,3)
hold on
plot(tS1,xmeanS1)
plot(tS2,xmeanS2)
ylim([-0.05 1.7])
ylabel('Mean trait (x)')
xlabel('Time')