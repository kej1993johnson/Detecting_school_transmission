function [f] = seirsc_model_rhs(t, y, params)

SK = y(1);
EK = y(2);
IK = y(3);
RK = y(4);
cumIK = y(5);
SA = y(6);
EA = y(7);
IA = y(8);
RA = y(9);
cumIA = y(10);


f = zeros(10,1); % need to return a column vector
alpha = params(1);
ps = params(2); 
pc = params(3);
gamma = params(4);
muKK = params(5);
muAA = params(6);
muAK = params(7);
muKA = params(8);
N = SA+EA+IA+RA+SK+EK+IK+RK;
NA = SA+EA+IA+RA;
NK  = SK+EK+IK+RK;


f(1) = -SK*((ps*muKK*IK/NK) + (pc*muKA*IA/NA));
f(2) = SK*((ps*muKK*IK/NK) + (pc*muKA*IA/NA)) - alpha*EK;
f(3) = alpha*EK - gamma*IK;
f(4) = gamma*IK;
f(5) = alpha*EK; % cumulative child infections

f(6) = -SA*((pc*muAA*IA/NA) + (pc*muAK*IK/NK));
f(7) = SA*((pc*muAA*IA/NA) + (pc*muAK*IK/NK)) - alpha*EA;
f(8) = alpha*EA - gamma*IA;
f(9) = gamma*IA;
f(10) = alpha*EA; % cumulative community infections 
end

