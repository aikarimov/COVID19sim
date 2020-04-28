%% main function
function A_Equation_Russia_mSEIRD
% Smiulate COVID-19 epidemic in Russia
% with a modified SEIRD model
% Actual date: 28-04-2020
% Copyrigth: A.I. Karimov, Youth Research Institute, ETU "LETI"
% You can use this code in accordance with CC BY-ND license
name = 'Russia';
name_rus = 'Россия';
global tmax Npop

Npop = 144.5e6;

t0 = datetime(2020,4,2); %initial date: 2 apr 2020
%total
n = [3548,4149,4731,5389,6343,7497,8672,10131,11917,13584,15770,18328,21102,24490,27938,32008,36793,42853, 47121, 52763, 57999, 62773, 68662, 74588, 80949, 87147, 93558];
%active
n_a = [3283,3834,4355,4989,5890,6945,8029,9357,11028,12433,14349,16710,19238,22306,25402,29145,33423,39201, 43270, 48434, 53066, 57327, 62439, 67657,73435, 79007, 93558 - 867 - 8456];
%death
d = [30,34,43,45,47,58,63,76,94,106,130,148,170,198,232,273,313,361,405,456,513,555,615,681,747, 794, 867];

%daily cases
n_d = n(2:end) - n(1:end-1);
%daily death
d_d = d(2:end) - d(1:end-1);

%time span 
t = 0:length(n)-1;
tmax = max(t);

Ntrials = 100;%trials to find the best fit
ObjLim = 1e8;

%Find fractional-order integration dim2 optimization problem solution
%try to find the best fit
errMod = realmax;
for i = 1:Ntrials
    %initial point
    bet = rand;
    gam = 0.1*rand;
    a1 = 1e-6*rand;
    q = 1e-3*rand;
    mu = 1e-3*rand;
    K0 = 1e-3*rand;
    x0 = [bet gam a1 q mu K0];
    
    %constraints: all values are positive
    A = -eye(6);
    b = [0; 0; 0; 0; 0; 0];
    if abs(ObjFun(x0, t, n, n_a, d)) < ObjLim
        %use fmincon function for fitting the model
        [xopt_t, errMod_t] = fmincon(@(x) ObjFun(x, t, n, n_a, d),x0, A, b);
        if errMod_t < errMod
            errMod = errMod_t;
            xoptMod = xopt_t;
        end
    end
end
%sumulation of the final model
[tspan, y, yact, yd] = SimCoronamSEIRD(xoptMod,n(1),n_a(1), d(1));

%plot fits in lin and log scales
figure(1);
subplot(1,2,1);
plot(t,n,'o',tspan,y,'-',t,n_a,'s',tspan,yact,'-',t,d,'d',tspan,yd,'-');
xlabel('$t$, days','interpreter','latex');
ylabel('$n$, total cases','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
title(name,'interpreter','latex');

subplot(1,2,2);
semilogy(t,n,'o',tspan,y,'-',t,n_a,'s',tspan,yact,'-',t,d,'d',tspan,yd,'-');
xlabel('$t$, days','interpreter','latex');
ylabel('$n$, total cases','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
title(name,'interpreter','latex');
legend('Total cases',...
    ['mSEIRD, $\beta$ = ',num2str(xoptMod(1),'%0.4f'),...
    ' $\gamma$ = ',num2str(xoptMod(2),'%.4f'),...
    ' $\alpha_1$ = ',num2str(xoptMod(3),'%.2e'),...
    ' $q$ = ',num2str(xoptMod(4),'%.2f'),...
    ' $\mu$ = ',num2str(xoptMod(5),'%.2e'),...
    ' $K_0$ = ',num2str(xoptMod(6),'%.2e'),...
    ],'Active cases','mSEIRD active','Deaths','mSEIRD deaths','%.2f','lcn','southeastoutside','interpreter','latex');

disp('Error, mSEIRD:');
disp(errMod);

%now, make a prediction!
tmax = length(n) + 60;%60 days
[tspanp, yp2, zp2, yd2, in] = SimCoronamSEIRD(xoptMod,n(1),n_a(1), d(1));

t = t0 + caldays(t);
tspanp = t0 + caldays(tspanp);

%plot variables n, I, D
figure(3);
subplot(1,2,1);
plot(t,n,'o', tspanp,yp2,'-',t,n_a,'s', tspanp,zp2,'-',t,d,'d',tspanp,yd2,'-');
xlabel('Дни');
ylabel('Число заражений COVID-19');
title([name_rus ,', лин. масштаб']);
xtickformat('dd-MM');
xlim([tspanp(1) tspanp(end)]);
grid;
subplot(1,2,2);
semilogy(t,n,'o', tspanp,yp2,'-',t,n_a,'s',tspanp,zp2,'-',t,d,'d',tspanp,yd2,'-');
xlabel('Дни');
ylabel('Число заражений COVID-19');
title([name_rus ,', лог. масштаб']);
xtickformat('dd-MM');
grid;
legend('Всего, данные','Всего, mod.SEIR','Активных случаев, данные','Активных, mod.SEIRD','Смертей, данные','Смертей, mod.SEIRD');
xlim([tspanp(1) tspanp(end)]);
set(gcf,'position',[250 250 1130 300]);
%plot derivatives dI+, dD
figure(4);
subplot(1,2,1);
t = t(2:end);
tspanp = tspanp(2:end);
plot(t,n_d,'s', tspanp,in,t,d_d,'d',tspanp,yd2(2:end) - yd2(1:end-1),'.-');
xlabel('Дни');
ylabel('Ежедневный прирост');
title([name_rus ,', лин. масштаб']);
xtickformat('dd-MM');
xlim([tspanp(1) tspanp(end)]);
grid;
subplot(1,2,2);
semilogy(t,n_d,'s', tspanp,in,t,d_d,'d',tspanp,yd2(2:end) - yd2(1:end-1),'.-');
xlabel('Дни');
ylabel('Ежедневный прирост');
title([name_rus ,', лог. масштаб']);
xtickformat('dd-MM');
grid;
legend('Новых случаев, данные','Новых случаев, mod.SEIRD','Новых смертей, данные','Новых смертей, mod.SEIRD');
xlim([tspanp(1) tspanp(end)]);
set(gcf,'position',[250 250 1130 300]);
end

%% Objective function for optimization with mean quadratic error estimation
function err = ObjFun(x, t, n, n_a, d)
% x - parameters
% t - time span
% n - data for total cases
% n_a - data for active cases
% d - data for total deaths

    err = 0; %total err
    nf = 0; %number of data points
    [tspan, y, yact, yd, ~] = SimCoronamSEIRD(x,n(1),n_a(1),d(1));
    N = length(tspan);
    M = length(t);
    
    %compare! this code is suitable for any rational simulation stepsize
    for i = 1:N
        for j = 1:M
            if(t(j) == tspan(i))
                err = err + ((y(i) - n(j))^2) + ((yact(i) - n_a(j))^2) + 100*((yd(i) - d(j))^2); %weighted sum of squares
                nf = nf + 1;
            end
        end
    end
    err = sqrt(err)/ nf;
end

%% Simulate COVID-19 epidemic with Euler method 
function [tspan, y, yact, yd, in] = SimCoronamSEIRD(x,n0,n_a0,d)
global tmax Npop
h = 1; %stepsize

%parameters
bet = x(1);
gam = x(2);
a1 = x(3);
q = x(4);
mu = x(5);
K0 = x(6);

tspan = 0:h:tmax;
N = length(tspan);

%data sets for stats
y = zeros(1,N);
yact = zeros(1,N);
yd = zeros(1,N);
in = zeros(1,N-1);

n = n0; %total cumulative I
S = Npop - n0; %susceptible
I = n_a0; %infected
E = q*I; %assumption: Exposed = Infected
R = n0 - n_a0; %actually resistant
D = d;
sigma = 1/3;
theta = 0.6;

%set initial values
y(1) = n;
yact(1) = I;
yd(1) = d;

for i = 2:N %start simulation
    O = exp(-a1*(I + theta*E)^K0); %saturation function
    
    dS = ( - bet*S*(I + theta*E)/Npop*O );
    dE = (   bet*S*(I + theta*E)/Npop*O - sigma*E );
    dI = ( sigma*E - gam*I - mu*I);
    dR = gam*I;
    dD = mu*I;
    
    in(i - 1) = sigma*E; %new cases
    
    S = S + h*dS;
    E = E + h*dE;
    I = I + h*dI;
    R = R + h*dR;
    D = D + h*dD;
    
    n = I + R + D; %total cases
    y(i) = n;
    yact(i) = I; %active cases
    yd(i) = D; %deaths
end
end

