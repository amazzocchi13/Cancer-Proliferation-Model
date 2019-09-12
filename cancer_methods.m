% AR Mazzocchi
% Methods for Determining Cancer Ratios for Prediction
% Updated 23-April-2018

% No Treatment
%tspan = [0 72]; % hours
%np = [670 780 930 1190 1560 1980 2470]; % number of cells
%dnp = [0 110 150 260 370 420 490]; % rate of change
%nq = [56 77 105 138 180 223 272]; % number of cells
%dnq = [0 21 27 33 42 43 49]; % rate of change


% Treatment 1 - G phase inhibitor
%tspan = [0 72]; % hours
%np = [600 680 812 974 1203 1590 2064]; % number of cells
%dnp = [0 80 132 162 229 387 474]; % rate of change
%nq = [54 72 100 132 184 247 312]; % number of cells
%dnq = [0 18 28 32 52 63 65]; % rate of change

% Treatment 2 - M phase inhibitor
tspan = [0 72]; % hours
np = [600 645 580 536 500 476 430]; % number of cells
dnp = [0 45 -65 -44 -36 -24 -46]; % rate of change
nq = [59 77 79 76 70 68 70]; % number of cells
dnq = [0 18 2 -3 -6 -2 2]; % rate of change

%% Function
%cancer_methods(@(z)tspan,np,nq,dnp,dnq)
%function z = cancer_methods(tspan, np, nq, dnp, dnq)

%% Initial Conditions
tint = [tspan(1):12: tspan(2)]; % Works for this but needs adjusted
nq = [nq]/(np(1));
dnp = [dnp]/(np(1))/(tint(2));
dnq = [dnq]/(np(1))/(tint(2));
np = [np]/(np(1));
k1 = 0.125; % cells doubling per hour
k2 = 0.0075*k1;
k3 = 0.005*k1;
yzero = [np(1) nq(1)];

%% Calculating Constants and Determining Regression Error
i = 1;

for i = 1: length(np)
    k4(i) = ((k1 - k3)*np(i)-dnp(i))/(np(i));
    k5(i) = ((k2+2*k3)*np(i)-dnq(i))/(nq(i)); 
end

k4 = mean([k4]);
k5 = mean([k5]);

% Regression Error - k4
% Estimate  
dnphat = (k1-k3-k4)*[np];
% Residuals
r = dnp-dnphat; % Actual minus when average k4 is used  
% SSR
ssr = sum(r.^2);
% Coefficient of Determination
dnpbar = mean(dnp);
R2 = (1-ssr/sum((dnp-dnpbar).^2))*100;

% Plot
figure(1); 
plot(np,dnp,'*')
hold on
plot(np,dnphat,'r')
title(['Linear Model, R2=' num2str(R2)])

% Regression Error - k5
% Estimate 
dnqhat = (k2+2*k3)*[np]-k5*[nq];
% Residuals
r = dnq-dnqhat;
% SSR
ssr = sum(r.^2);
% Coefficient of Determination
dnqbar = mean(dnq);
R2 = (1-ssr/sum((dnq-dnqbar).^2))*100;

% plot

figure(2); 
plot(nq,dnq,'*')
hold on
plot(np,dnqhat,'r')
title(['Linear Model, R2=' num2str(R2)])

%% ODE to Model Behavior

% ODE to Approximate Cancer Cell Behavior
[t,y]=ode45(@cancer_kinetics,tspan,yzero,[],k1,k2,k3,k4,k5);
npint = interp1(t,y(:,1),tint); % interpolate points for np at times matching exact
nqint = interp1(t,y(:,2),tint); % interpolates points for nq at times matching exact

% Error Calculations
error_np = np - npint(:,1)';
error_nq = nq - nqint(:,1)';
maxerr (1) = max(abs(error_np));
maxerr(2) = max(abs(error_nq));
maxerr

% Plot Cancer Cell Behavior from ODE45 
figure(3); 
plot(t,y(:,1),'-',t,y(:,2),'-.')
title('Cancer Cell Populations over Time - ODE45','FontSize',12)
ylabel('Cell Count','FontSize',12); 
legend('Proliferating Cancer Cells','Quiescient Cancer Cells');

%% Exact Data vs ODE Model

% Plot Cancer Cell Behavior from Exact Data
figure(4);
plot([0:12:72],np,'-',tint,nq,'-.')
title('Cancer Cell Populations over Time - ACTUAL','FontSize',12)
ylabel('Cell Count','FontSize',12); 
legend('Proliferating Cancer Cells','Quiescient Cancer Cells');

% Plot of Both ODE45 and Exact Data
figure(5); 
plot(t,y(:,1),'-',t,y(:,2),'-.',[0:12:72],np,'-',[0:12:72],nq,'-.')
title('Cancer Cell Populations over Time - ODE45 and Exact','FontSize',12)
ylabel('Cell Count','FontSize',12); 
legend('ODE - Proliferating Cancer Cells','ODE - Quiescient Cancer Cells','Exact - Proliferating Cancer Cells','Exact - Quiescient Cancer Cells');

%% Outputs of Importance

% Consider ratios at day 1 (24 hrs), day 3 (72 hrs)
% From ODE:

np12 = interp1(t,y(:,1),tint(2)); 
nq12 = interp1(t,y(:,2),tint(2)); 

np72 = interp1(t,y(:,1),tint(7)); 
nq72 = interp1(t,y(:,2),tint(7)); 

np_ratio = np72/np12;
nq_ratio = nq72/nq12;

fprintf ('Np Ratio = %.4f\n', np_ratio);
fprintf ('Nq Ratio = %.4f\n', nq_ratio);









