%
% Fig6
% Author: Catalina Chaparro
% Last modification: 20/5/2021
%
% This rutine makes figure 6 of the paper entitled "Fast
% environmental change and eco- evolutionary feedbacks can drive
% regime shifts in ecosystems before tipping points are crossed"

clear

%Parameters----------------------------------------------------

K     = 10;         % Carrying capacity
b     = 1;          % Maximum birth rate
do    = .5;         % Background mortality rate
d1    = .1;         % Mortality rate due to the fecundity-survival tradeoff
A     = -1;         % Allee threshod
tau   = 1;          % Degree of especialization
sigma = 0;          % Genetic variation

E     = 1:.01:3;    % Environmental stress
x     = 0:.01:3;    % Trait value

lowdens=1E-6;

%Calculate the per capita growth rate (fig 6a)

siE = size(E);
siX = size(x);

percapitagrowth=zeros(siE(1,2),siX(1,2));
 
for i=1:siE(1,2)
    for j=1:siX(1,2)
        bt = b*tau/(sigma^2+tau^2)^(1/2) * exp(-(x(1,j)-E(1,i))^2/(2*(sigma^2+tau^2)));
        dt = do+d1*(sigma^2+x(1,j)^2);

        percapitagrowth(i,j) = bt*(1-lowdens/K)*(lowdens-A)-dt;
    end
end

%Calculate the invasion threshold (TPi)

xmean=0:.01:3;
svecx=size(xmean);

Tipping = zeros(svecx(1,2),1);
for i=1:svecx(1,2)
    if(xmean(1,i)>2)
        Tipping(i,1) = NaN;
    else
        Tipping(i,1) = findTPinv(b,do,d1,K,A,tau,sigma,x(1,i));
    end
end

%Calculate ecological trajectories in fig 6b

traits = 0:.5:2; %Trait values to evaluate
Ef = 2; %Fixed environmental stress
str = size(traits);

for i=1:str(1,2)
    x0=[1E-3 traits(1,i) Ef];
    fod=@(t,x) PopAlleeODE(t,x,b,do,d1,K,A,tau,0,Ef,0);
    [t,x]=ode23s(fod,[0 500],x0);
    ti{i}=t;
    Ni{i}=x(:,1);
end

% Plot --------------------------------------------------

numcolors = 11;
colors=redblue(numcolors);

if rem(numcolors,2)>0
    remove=ceil(numcolors/2);
else
    remove=[numcolors/2 numcolors/2+1];
end

colors(remove,:)=[];

figure
suptitle('Figure 6')
subplot(2,1,1)
contourf(xmean,E,percapitagrowth,'edgecolor','none')
colormap(colors)
caxis([-1 1])
lc = colorbar;
lc.Label.String = 'per capita growth rate at low density';
hold on
plot(xmean,Tipping,'k','LineWidth',3)
xlabel('Trait value (x)')
ylabel('Environmental stress (E)')

subplot(2,1,2)
for i=1:str(1,2)
    plot(ti{i},Ni{i})
    hold on
    leg{i}=sprintf('Trait = %1.1f',traits(1,i));
end
xlabel('Time')
ylabel('Population density (N)')
legend(leg)
