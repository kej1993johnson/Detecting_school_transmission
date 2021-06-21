function [IK_t, IA_t, cumIK_t, cumIA_t, detIK_t, detIA_t, detcumIK_t, detcumIA_t, tfirst, Ifirst, tclose, Iclose] = sim_school_com2(I0K, I0A, NK,NA,tvec, R0s,R0c,muKK, muAA, muAK, muKA, pSEIR, psymK, psymA, tseek,tdelay,sens,close_thresh)




% epidemiological parameters
alpha = pSEIR(1);
gamma = pSEIR(2);
ps = R0s*gamma/9.7; % set the probability of infection based on R0 and the polymod contact rate
pc = R0c*gamma/9.0;

params = [alpha,ps, pc, gamma, muKK, muAA, muAK, muKA];

y0 = [NK-I0K 0 I0K 0 0 NA-I0A 0 I0A 0 0]; % corresponding to:
%SK0, EK0, IK0, RK0, cumIK,
%SA0, EA0, IA0, RA0, cumIA,

[t,y] = ode23s(@seirsc_model_rhs,tvec,y0,[],params);
    IK_t = y(:,3);
    cumIK_t = y(:,5);
    IA_t = y(:,8);
    cumIA_t = y(:,10);


detIK_t = zeros(length(tvec),1);
detcumIK_t = zeros(length(tvec),1);
dts = diff(tvec);
dt = dts(1);
tseek = tseek./dt;
tdelay = tdelay./dt;

for i = tseek+tdelay+1:length(tvec)
    % detected infections in the school
    detcumIK_t(i) = psymK*cumIK_t(i-tseek-tdelay)*sens; % detected infections are those symptomatic, found X days later,
    % where X is the time to seek a test and then the delay in test results
    detIK_t(i) = psymK*IK_t(i-tseek-tdelay)*sens; 
    % detected infections in the community
    detcumIA_t(i) = psymA*cumIA_t(i-tseek-tdelay)*sens; % detected infections are those symptomatic, found X days later,
    % where X is the time to seek a test and then the delay in test results
    detIA_t(i) = psymA*IA_t(i-tseek-tdelay)*sens; 
end



ifirst = find(detcumIK_t>1,1);
tfirst = tvec(ifirst);
Ifirst = cumIK_t(ifirst);


indclose = find((100*detcumIK_t./NK)>close_thresh,1);
tclose = tvec(indclose);
Iclose = cumIK_t(indclose);


if isempty(Ifirst)
    Ifirst = NaN;
end
if isempty(tfirst)
    tfirst = NaN;
end


if isempty(Iclose)
    Iclose = NaN;
end
if isempty(tclose)
    tclose = NaN;
end















end