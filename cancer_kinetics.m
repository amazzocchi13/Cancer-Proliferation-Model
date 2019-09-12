% AR Mazzocchi
% Cancer Kinetics for ODE45
% Numerical Methods Project
% Updated 23-April-2018
% Notes:
% Cancer kinetics function that inputs k values from regression model into
% rate equations related to change in cancer cell populations


function dy = cancer_kinetics(t,y, k1,k2,k3,k4,k5)
% cancer_kinetics.m
% Contains equations for proliferating and quiescent cancer cells
% Variables
np = y(1);
nq = y(2);

% Equations
dy =[(k1-k3-k4)*np
    (k2+2*k3)*np-k5*nq];
  
