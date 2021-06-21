% School reopening + transmission in the community.
% This script runs sim_school_com which takes epidemiological parameters of
% a school and a community and returns the infections in each and the
% detected infections in each (from symptoms).

close all; clear all; clc
%% Start with calling a function with:
% INPUTS: 
% initial condition (I0A, I0K) 
% within group R0 at contact rate estimated by POLYMOD(R0s, R0c)
% contact rates between adults and children (within and between - muAA,
% muKK, muKA, muAK)
% percent symptomatic (psymA, psymK)
% time to seek testing if symptomatic (tseek) 
% time to get test results back (tdelay)
% threshold for closing a school based on detected case (close_thresh)
% OUTPUTS:
% infected over time (IA_t, IC_t)
% detected cases over time (detIA_t, detIA_t)
% time to detect first case (tfirstA, tfirstK,)
% outbreak size at time of detection (IfirstA, IfirstK)
% time to close (tclose)
% outbreak size at time of closing (Iclose)
 
% Set parameters 
NK= 1000;
NA = 4000;
N0 = NA + NK;
I0K = (3/1000)*NK;
I0A = (3/1000)*NA;
R0s = 2.5; 
R0c = 1.1;
psymK = 0.21;
psymA = 0.7;
tseek = 0;
tdelay = 0;
sens = 0.9;
close_thresh = 1;
tvec = 0:1:120;
tvec = tvec';
alpha = 1/3; % 1/latency
gamma = 1/10; % 1/infection duration 
pSEIR = [alpha, gamma];
muAA = 9;
muAK = 2.5;
muKK = 9.7;
muKA = 5.3;
%% Test the model

 [IK_t, IA_t, cumIK_t, cumIA_t, detIK_t, detIA_t, detcumIK_t, detcumIA_t, tfirst, Ifirst, tclose, Iclose] = sim_school_com2(I0K, I0A, NK,NA,tvec, R0s,R0c,muKK, muAA, muAK, muKA, pSEIR, psymK, psymA, tseek,tdelay,sens, close_thresh);
 figure;
 subplot(2,1,1)
 plot(tvec, 1000*detcumIK_t/NK, 'b-', 'LineWidth', 2)
 hold on
 plot(tvec, 1000*detcumIA_t/NA, 'r-', 'LineWidth',2)
 xlabel('Time(days)')
 ylabel('Detected cases per 1000')
 legend('school-aged children', 'adults', 'Location', 'NorthWest')
 legend boxoff
  %ylim([0 1000])
 set(gca,'FontSize',18,'LineWidth',1.5)
 subplot(2,1,2)
 plot(tvec, 1000*cumIK_t/NK, 'b-', 'LineWidth', 2)
 hold on
 plot(tvec, 1000*cumIA_t/NA, 'r--', 'LineWidth',2)
 xlabel('Time(days)')
 ylabel('Cumulative infections per 1000')
 legend('school-aged children', 'adults', 'Location', 'NorthWest')
 legend boxoff
 %ylim([0 1000])
 set(gca,'FontSize',18,'LineWidth',1.5)
 
 figure;
 subplot(2,1,1)
 plot(tvec, 1000*detIK_t/NK, 'b-', 'LineWidth', 2)
 hold on
 plot(tvec, 1000*detIA_t/NA, 'r-', 'LineWidth',2)
 xlabel('Time(days)')
 ylabel('Detected cases per 1000')
 legend('school-aged children', 'adults', 'Location', 'NorthWest')
 legend boxoff
 % ylim([0 300])
 set(gca,'FontSize',18,'LineWidth',1.5)
 subplot(2,1,2)
 plot(tvec, 1000*IK_t/NK, 'b-', 'LineWidth', 2)
 hold on
 plot(tvec, 1000*IA_t/NA, 'r--', 'LineWidth',2)
 xlabel('Time(days)')
 ylabel('True infections per 1000')
 legend('school-aged children', 'adults', 'Location', 'NorthWest')
 legend boxoff
 %ylim([0 300])
 set(gca,'FontSize',18,'LineWidth',1.5)
%%  True and detected infections in schools and communities under 3 scenarios
% Set up the 3 scenarios: 
% 1. School closed (school R0 = community R0, u = 1)
% 2. school open, mixing on weekends (school R0 = 2.5, u = 2/7)
% 3. school open, completely well-mixed (school R0 = 2.5, u = 1)
N0=5000;
R0svec = [R0c, R0s, R0s];
muAKvec = [muAK, muAK, 3*muAK];
muAAvec = [muAA, muAA, muAA]; % shift so that kids in contact with adults twice as much
muKAvec = [muKA, muKA, 3*muKA];
muKKvec = [muKK, muKK, muKK]; % shift so that adults in contact with kids twice as much 




for i= 1:length(R0svec)

    R0si = R0svec(i);
    muAKi = muAKvec(i);
    muAAi = muAAvec(i);
    muKAi = muKAvec(i);
    muKKi = muKKvec(i);

% 
    [IK_t, IA_t, cumIK_t, cumIA_t, detIK_t, detIA_t, detcumIK_t, detcumIA_t, tfirst, Ifirst, tclose, Iclose] = sim_school_com2(I0K, I0A, NK,NA,tvec, R0si,R0c,muKKi, muAAi, muAKi, muKAi, pSEIR, psymK, psymA, tseek,tdelay,sens, close_thresh);
    IK_ti(:,i)=IK_t;
    cumIK_ti(:,i) = cumIK_t;
    detIK_ti(:,i)=detIK_t;
    detcumIK_ti(:,i) = detcumIK_t;
    IA_ti(:,i)=IA_t;
    cumIA_ti(:,i) = cumIA_t;
    detIA_ti(:,i)=detIA_t;
    detcumIA_ti(:,i) = detcumIA_t;
end
tvec = tvec/7
figure;
subplot(1,3,1)
plot(tvec, 1000*IA_ti(:,1)/NA, '-','color', [0.494, 0.184, 0.556], 'LineWidth', 3)
hold on
plot(tvec, 1000*IA_ti(:,2)/NA, '-','color', [0 0.6 0.3], 'LineWidth', 3)
plot(tvec, 1000*IK_ti(:,1)/NK, '--', 'color', [0.494, 0.184, 0.556],'LineWidth', 3)
hold on 
plot(tvec, 1000*IK_ti(:,2)/NK, '--', 'color', [0 0.6 0.3], 'LineWidth', 3)
%plot(tvec, 1000*IK_ti(:,3)/NK,  'color', [0.85, 0.325, 0.098], 'LineWidth', 2)
plot(tvec, 1000*IA_ti(:,1)/NA, '-','color', [0.494, 0.184, 0.556], 'LineWidth', 3)
plot(tvec, 1000*IA_ti(:,2)/NA, '-','color', [0 0.6 0.3], 'LineWidth', 3)
%plot(tvec, 1000*IK_ti(:,3)/NK,  'color', [0.85, 0.325, 0.098], 'LineWidth', 2)
%plot(tvec, 1000*IA_ti(:,3)/NA, '--', 'color', [0.85, 0.325, 0.098], 'LineWidth', 3)
set(gca,'FontSize',18,'LineWidth',1.5)
%legend('school closed', 'school open', 'school open, 2x adult-children contacts', 'Location', 'NorthWest')
%legend boxoff

%legend('closed, school', 'closed, community', 'open partial mixing, school',...
%'open partial mixing, community', 'open well-mixed, school', 'open-well-mixed, community',...
 %'Location', 'NorthWest')
%legend boxoff
ylim([0 250])
xlim([0 tvec(end)])
ylabel('True infections per 1000')
xlabel('Time (days)')
subplot(1,3,2)
plot(tvec, 1000*detIK_ti(:,1)/NK,'--', 'color', [0.494, 0.184, 0.556], 'LineWidth', 3)
hold on 
plot(tvec, 1000*detIK_ti(:,2)/NK, '--', 'color', [0 0.6 0.3], 'LineWidth', 3)
%plot(tvec, 1000*detIK_ti(:,3)/NK, 'color', [0.85, 0.325, 0.098], 'LineWidth', 2)
plot(tvec, 1000*detIA_ti(:,1)/NA, '-','color', [0.494, 0.184, 0.556], 'LineWidth', 3)
plot(tvec, 1000*detIA_ti(:,2)/NA, '-','color', [0 0.6 0.3], 'LineWidth', 3)
%plot(tvec, 1000*detIK_ti(:,3)/NK, 'color', [0.85, 0.325, 0.098], 'LineWidth',2)
%plot(tvec, 1000*detIA_ti(:,3)/NA, '--', 'color', [0.85, 0.325, 0.098], 'LineWidth', 3)
set(gca,'FontSize',18,'LineWidth',1.5)
%legend('school closed', 'school open, partial interaction', 'school open, well-mixed', 'Location', 'NorthWest')
%legend boxoff
ylim([0 80])
xlim([0 tvec(end)])
ylabel('Detected cases per 1000')
xlabel('Time (days)')

ratio_school_closed = (detcumIK_ti(:,1)/NK)./(detcumIA_ti(:,1)/NA);
ratio_school_opened = (detcumIK_ti(:,2)/NK)./(detcumIA_ti(:,2)/NA);

subplot(1,3,3 )
plot(tvec, ratio_school_closed,'-', 'color', [0.494 0.184 0.556], 'LineWidth', 3) % purple
hold on
plot(tvec, ratio_school_opened,'-', 'color', [0 0.6 0.3], 'LineWidth', 3) % green
plot([0 tvec(end)], [(3453/5466) (3453/5466)], 'k-', 'LineWidth', 3)
set(gca,'FontSize',18,'LineWidth',1.5)
%legend('school closed', 'school open', 'Location', 'NorthWest')
%legend boxoff
xlim([0 tvec(end)])
ylabel('Children vs. adult cumulative cases per capita', 'FontSize', 16)
xlabel('Time (days)')

figure;
subplot(1,3,1)
plot(tvec, 1000*cumIA_ti(:,2)/NA, 'k-', 'LineWidth', 3) % adults black
hold on
plot(tvec, 1000*cumIK_ti(:,1)/NK,'k--',  'LineWidth', 3) % kids black

plot(tvec, 1000*cumIA_ti(:,1)/NA, '-','color', [0.494, 0.184, 0.556], 'LineWidth', 3) % adults open
plot(tvec, 1000*cumIA_ti(:,2)/NA, '-','color', [0 0.6 0.3], 'LineWidth', 3) % adults closed
plot(tvec, 1000*cumIK_ti(:,1)/NK,'--', 'color', [0.494, 0.184, 0.556], 'LineWidth', 3)

plot(tvec, 1000*cumIK_ti(:,2)/NK, '--', 'color', [0 0.6 0.3], 'LineWidth', 3)
%plot(tvec, 1000*detIK_ti(:,3)/NK, 'color', [0.85, 0.325, 0.098], 'LineWidth', 2)
plot(tvec, 1000*cumIA_ti(:,1)/NA, '-','color', [0.494, 0.184, 0.556], 'LineWidth', 3)
plot(tvec, 1000*cumIA_ti(:,2)/NA, '-','color', [0 0.6 0.3], 'LineWidth', 3)
%plot(tvec, 1000*detIK_ti(:,3)/NK, 'color', [0.85, 0.325, 0.098], 'LineWidth',2)
%plot(tvec, 1000*detIA_ti(:,3)/NA, '--', 'color', [0.85, 0.325, 0.098], 'LineWidth', 3)
set(gca,'FontSize',18,'LineWidth',1.5)
%legend('adults', 'test', 'low risk', 'high risk','Location', 'NorthWest')
%legend boxoff
xlim([0 tvec(end)])
ylabel('Infections per 1,000')
xlabel('Time (weeks)')
%title('Cumulative infections')

subplot(1,3,2)
plot(tvec, 1000*detcumIA_ti(:,1)/NA, '-','color', [0.494, 0.184, 0.556], 'LineWidth', 3)
hold on
plot(tvec, 1000*detcumIA_ti(:,2)/NA, '-','color', [0 0.6 0.3], 'LineWidth', 3)
plot(tvec, 1000*detcumIK_ti(:,1)/NK,'--', 'color', [0.494, 0.184, 0.556], 'LineWidth', 3)
hold on 
plot(tvec, 1000*detcumIK_ti(:,2)/NK, '--', 'color', [0 0.6 0.3], 'LineWidth', 3)
%plot(tvec, 1000*detIK_ti(:,3)/NK, 'color', [0.85, 0.325, 0.098], 'LineWidth', 2)
plot(tvec, 1000*detcumIA_ti(:,1)/NA, '-','color', [0.494, 0.184, 0.556], 'LineWidth', 3)
plot(tvec, 1000*detcumIA_ti(:,2)/NA, '-','color', [0 0.6 0.3], 'LineWidth', 3)
%plot(tvec, 1000*detIK_ti(:,3)/NK, 'color', [0.85, 0.325, 0.098], 'LineWidth',2)
%plot(tvec, 1000*detIA_ti(:,3)/NA, '--', 'color', [0.85, 0.325, 0.098], 'LineWidth', 3)
set(gca,'FontSize',18,'LineWidth',1.5)
%legend('school closed', 'school open','NorthWest')
%legend boxoff
xlim([0 tvec(end)])
ylabel('Detected infections per 1,000')
xlabel('Time (weeks)')
%title('Cumulative detected infections')

adult_rate = (3393-133)/(73000-4876);
child_rate = 133/4876;
ratio = child_rate/adult_rate;
subplot(1,3,3 )
%plot(91, ratio, 'k*', 'LineWidth', 5)
plot([tvec(1), tvec(end)], [ratio, ratio], 'k-.', 'LineWidth', 3)
hold on
plot(tvec, ratio_school_closed,'-.', 'color', [0.494 0.184 0.556], 'LineWidth', 3) % purple
plot(tvec, ratio_school_opened,'-.', 'color', [0 0.6 0.3], 'LineWidth', 3) % green
%plot(91, ratio, 'k*', 'LineWidth', 5)
set(gca,'FontSize',18,'LineWidth',1.5)
%legend('Wood County', 'Location', 'NorthWest')
%legend boxoff
xlim([0 tvec(end)])
ylim([0 1.5])
ylabel('Child-adult detected infection ratio')
xlabel('Time (weeks)')
%title('Children to adult detected infections ratio')

community_rate = (3393-191)/(73000-4876-654)

%% Find the impact of mixing factor on the total cumulative infections and the disparity between infections
tvec = 0:1:120;
tvec = tvec';


intvec = linspace(0, 5, 20);
colors = varycolor(length(intvec));
figure;
for i= 1:length(intvec)
    int=intvec(i);
    muKAi = muKA*int;
    muAKi = muAK*int;
    [IK_t, IA_t, cumIK_t, cumIA_t, detIK_t, detIA_t, detcumIK_t, detcumIA_t, tfirst, Ifirst, tclose, Iclose] = sim_school_com2(I0K, I0A, NK,NA,tvec, R0s,R0c,muKK, muAA, muAKi, muKAi, pSEIR, psymK, psymA, tseek,tdelay,sens, close_thresh);
    IK_ti(:,i)=IK_t;
    cumIK_ti(:,i) = cumIK_t;
    detIK_ti(:,i)=detIK_t;
    detcumIK_ti(:,i) = detcumIK_t;
    IA_ti(:,i)=IA_t;
    cumIA_ti(:,i) = cumIA_t;
    detIA_ti(:,i)=detIA_t;
    detcumIA_ti(:,i) = detcumIA_t;
    total_inf(i)=100*(cumIK_ti(end,i) +cumIA_ti(end,i))/N0;
    total_inf_comm(i) = cumIA_ti(end,i);
    K_A_ratio(i) = (cumIK_ti(end,i)/NK)/(cumIA_ti(end,i)/NA); % ratio of infection prevalence in school vs in community 
    detK_A_ratio(i)=(detcumIK_ti(end,i)/NK)/(detcumIA_ti(end,i)/NA);
subplot(1,3,1)
plot(tvec, 1000*(IA_ti(:,i) + IK_ti(:,1))/(NK+ NA),'color',  colors(i,:), 'LineWidth', 2)
hold on
xlabel('Time(days)')
ylabel('Total infections per 1000')
set(gca,'FontSize',18,'LineWidth',1.5)
legend('0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1','Location', 'NorthWest')
legend boxoff
xlim([0, tvec(end)])
    
    
subplot(1,3,2)
plot(tvec, 1000*IA_ti(:,i)/NA,'color',  colors(i,:), 'LineWidth', 2)
hold on
xlabel('Time(days)')
ylabel('Adult infections per 1000')
set(gca,'FontSize',18,'LineWidth',1.5)
legend('0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1','Location', 'NorthWest')
legend boxoff
xlim([0, tvec(end)])
subplot(1,3,3)
plot(tvec, 1000*IK_ti(:,i)/NK, 'color',  colors(i,:), 'LineWidth', 2)
hold on
xlabel('Time(days)')
ylabel('School-aged children infections per 1000')
set(gca,'FontSize',18,'LineWidth',1.5)
legend('0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1','Location', 'NorthWest')
legend boxoff
xlim([0, tvec(end)])



end

figure;
subplot(1,2,1)
plot(intvec, total_inf, 'k-', 'LineWidth',2)
xlabel('Factor increase in adult-children contacts')
ylabel('True infection rate (%)')
set(gca,'FontSize',18,'LineWidth',1.5)


subplot(1,2,2)
plot(intvec, K_A_ratio, 'r-', 'LineWidth', 2)
hold on
plot(intvec, linspace(1,1,length(intvec)), 'k--', 'LineWidth', 2)
xlabel('Factor increase in adult-children contacts')
ylabel('Ratio of child to adult infection rate')

set(gca,'FontSize',18,'LineWidth',1.5)
ylim([0.75 2])



figure;
plot(intvec, total_inf_comm, 'r-', 'LineWidth', 3)
xlabel('Factor increase in adult-children contacts')
ylabel('Total infections in the community')
set(gca,'FontSize',18,'LineWidth',1.5)
%% Use NGM to find the group R0 as a function of mu and in-school R0
intvec =linspace(0, 1,40);
Reffsvec = linspace(1.1, 5,40);
phi = Ns/Nc;

[REFFS,INT] = meshgrid(Reffsvec, intvec); % big ol grid of parameters
REFFSflat = reshape(REFFS,1,[]);
INTflat = reshape(INT,1, []);

[R0, K, K1, lambda1, lambda2] = NGM(2.5, 1.1, gamma, 2/7, alpha, Ns, Nc);

for i = 1:length(REFFSflat)
    Reffsi = REFFSflat(i);
    inti = INTflat(i);
    [R0,K, K1, lambda1, lambda2]= NGM(Reffsi, Reffc, gamma, inti, alpha, Ns, Nc);
    lambdas = eig(K);
    groupR0(i) = lambdas(1);
    %groupR0(i) = max([real(lambda1),real(lambda2)]);
end


 GROUPR0 = reshape(groupR0, size(REFFS));
figure;
minmax = @(x)([min(x) max(x)]);
imagesc(minmax(1*(Reffsvec)),minmax(intvec),GROUPR0);
V = linspace(0,7,15);
%V2 = horzcat(linspace(0,15,4), linspace(20,100,9));
[C,h]=contourf(REFFS,INT,GROUPR0, V); clabel(C,h);
%[C2,h2]=contourf(REFF,INIT,TCLOSE, V2); 
hold on
%plot( [Reffc Reffc], [0 1], 'w--', 'LineWidth', 1)

colorbar
colormap(jet); 
%caxis([7 100])
xlabel('School R_0');
ylabel('School-Community Interaction (\mu)');
title('R_0 Students + Community')
set(gca,'FontSize',20,'LineWidth',1.5)











%% Remake figures from original manuscript with linked dynamic model
% First vary the initial prevalence
I0vec = linspace(0.4, 25,6);
I0vec = [0.4, 1, 3, 5,10,25]; % corresponding to low, medium, and high prevalence
colorsets2 = colormap(prism(6));
colorsets2= flip(colorsets2,1);
Ns = 1000;
Nc = 4000;
Reffs = 2.5;
psyms = 0.21;
psymc = 0.8;
for i = 1:length(I0vec)
    I0s=I0vec(i);
    I0c = (I0s/Ns)*Nc;
    [Is_t, Ic_t, cumIs_t, cumIc_t, detIs_t, detIc_t, detcumIs_t, detcumIc_t, tfirst, Ifirst, tclose, Iclose] = sim_school_com(I0s, I0c, Ns,Nc, tvec, Reffs,Reffc,int, pSEIR, psyms, psymc, tseek,tdelay,sens, close_thresh);
    %[I_t,cumI_t, detI_t, detcumI_t, tfirst, Ifirst, tclose, Iclose] = sim_school(I0i,N,tvec, Reff,pSEIR, psym, tseek,tdelay,sens, close_thresh); 
    I_ti(:,i)=Is_t;
    cumI_ti(:,i) = cumIs_t;
    detI_ti(:,i)=detIs_t;
    detcumI_ti(:,i) = detcumIs_t;
    tfirsti(i) = tfirst;
    Ifirsti(i) = Ifirst;
    tclosei(i) = tclose;
    Iclosei(i) = Iclose;
end

% record continuous iterating through I0
I0veccont = linspace(0.4,25,60);
for i = 1:length(I0veccont)
    I0s=I0veccont(i);
    I0c = (I0s/Ns)*Nc;
    [Is_t, Ic_t, cumIs_t, cumIc_t, detIs_t, detIc_t, detcumIs_t, detcumIc_t, tfirst, Ifirst, tclose, Iclose] = sim_school_com(I0s, I0c, Ns,Nc, tvec, Reffs,Reffc,int, pSEIR, psyms, psymc, tseek,tdelay,sens, close_thresh);
    %[I_t,cumI_t, detI_t, detcumI_t, tfirst, Ifirst, tclose, Iclose] = sim_school(I0i,N,tvec, Reff,pSEIR, psym, tseek,tdelay,sens, close_thresh); 
    tfirstci(i) = tfirst;
    Ifirstci(i) = Ifirst;
    tcloseci(i) = tclose;
    Icloseci(i) = Iclose;
end

% calculate upperbound (1 in 10 reporting rate)

for i = 1:length(I0veccont)
    I0s=I0veccont(i);
    I0c = (I0s/Ns)*Nc;
    Reffsi = 3.5;
    [Is_t, Ic_t, cumIs_t, cumIc_t, detIs_t, detIc_t, detcumIs_t, detcumIc_t, tfirst, Ifirst, tclose, Iclose] = sim_school_com(I0s, I0c, Ns,Nc, tvec, Reffsi,Reffc,int, pSEIR, psyms, psymc, tseek,tdelay,sens, close_thresh);
    %[I_t,cumI_t, detI_t, detcumI_t, tfirst, Ifirst, tclose, Iclose] = sim_school(I0i,N,tvec, Reffsi,pSEIR, psym, tseek,tdelay,sens, close_thresh); 
    tfirstupi(i) = tfirst;
    Ifirstupi(i) = Ifirst;
    tcloseupi(i) = tclose;
    Icloseupi(i) = Iclose;
end
% calculate lowerbound (1 in 3 reporting rate)
I0veclo = .6.*linspace(0.4,10,60);
for i = 1:length(I0veccont)
    I0s=I0veccont(i);
    I0c = (I0s/Ns)*Nc;
    Reffsi = 2.2;
    [Is_t, Ic_t, cumIs_t, cumIc_t, detIs_t, detIc_t, detcumIs_t, detcumIc_t, tfirst, Ifirst, tclose, Iclose] = sim_school_com(I0s, I0c, Ns,Nc, tvec, Reffsi,Reffc,int, pSEIR, psyms, psymc, tseek,tdelay,sens, close_thresh);
    %[I_t,cumI_t, detI_t, detcumI_t, tfirst, Ifirst, tclose, Iclose] = sim_school(I0i,N,tvec, Reffi,pSEIR, psym, tseek,tdelay,sens, close_thresh); 
    tfirstloi(i) = tfirst;
    Ifirstloi(i) = Ifirst;
    tcloseloi(i) = tclose;
    Icloseloi(i) = Iclose;
end
%% Vary prevalence and transmission rate
% First take varied prevalences from above run
figure;
subplot(2,2,1)
for i = 1:length(I0vec)
plot(tvec, 100.*cumI_ti(:,i)./Ns, '-', 'Linewidth',2, 'color', colorsets2(i,:))
%plot(tvec, 100.*10.*detI_ti(:,i)./N, '-', 'Linewidth',2, 'color', colorsets2(i,:))
hold on 
end
xlabel('time (days)')
xlim([0 tvec(end)])
i30 = find(tvec==30);
%plot(30, 100*cumI_ti(i30,2)./N, 'r*' ,'LineWidth', 8)
%ylim([0 50])
%ylabel('% of school quarantined')
ylabel('% of school infected')
legend('4 in 10,000','1 in 1000','3 in 1000', '5 in 1000', '10 in 1000','25 in 1000', 'Location', 'NorthWest')
legend boxoff
set(gca,'FontSize',18,'LineWidth',1.5)
%title('Percent quarantined')
title('Effect of prevalence on percent infected')

subplot(2,2,2)
plot(I0veccont, tcloseci, 'k-', 'LineWidth', 2)
hold on
plot(I0veccont, tcloseupi, 'k--', 'LineWidth',2)
plot(I0veccont, tcloseloi, 'k--', 'LineWidth',2)
for i =1:length(I0vec)
plot(I0vec(i), tclosei(i), '*', 'LineWidth', 10, 'color', colorsets2(i,:))
hold on
end
xlabel('Initial infections (out of 1000)')
ylabel('Time to close (days)')
set(gca,'FontSize',18,'LineWidth',1.5, 'Xscale', 'log')
title('Effect of prevalence on time to close')
ylim([5 max(tcloseloi)+1])
xlim([I0veccont(1) I0veccont(end)])
 

% Separate temporarily

% Then vary the transmission rate
I0s0 = 5; % national average prevalence
Reffvec = [1.5, 2.5, 3, 3.5, 4, 5 ]; % corresponding to low, medium, and high prevalence
for i = 1:length(Reffvec)
    Reffsi=Reffvec(i);
    I0c = (I0s0/1000)*Nc;
    I0s = I0s0;
    [Is_t, Ic_t, cumIs_t, cumIc_t, detIs_t, detIc_t, detcumIs_t, detcumIc_t, tfirst, Ifirst, tclose, Iclose] = sim_school_com(I0s, I0c, Ns,Nc, tvec, Reffsi,Reffc,int, pSEIR, psyms, psymc, tseek,tdelay,sens, close_thresh);
    %[I_t,cumI_t, detI_t, detcumI_t, tfirst, Ifirst, tclose, Iclose] = sim_school(I0,N,tvec, Reffi,pSEIR, psym, tseek,tdelay,sens, close_thresh); 
    I_ti(:,i)=Is_t;
    cumI_ti(:,i) = cumIs_t;
    detI_ti(:,i)=detIs_t;
    detcumI_ti(:,i) = detcumIs_t;
    cumcases(:,i) = detcumIs_t;
    tfirsti(i) = tfirst;
    Ifirsti(i) = Ifirst;
    tclosei(i) = tclose;
    Iclosei(i) = Iclose;
end
 %record continuous iterating through Reff
Reffveccont = linspace(1.5,5,10);
tfirstci = [];
Ifirstci = [];
tcloseci = [];
Icloseci = [];
for i = 1:length(Reffveccont)
    Reffsi=Reffveccont(i);
    I0c = (I0s0/1000)*Nc;
    I0s = I0s0;
    [Is_t, Ic_t, cumIs_t, cumIc_t, detIs_t, detIc_t, detcumIs_t, detcumIc_t, tfirst, Ifirst, tclose, Iclose] = sim_school_com(I0s, I0c, Ns,Nc, tvec, Reffsi,Reffc,int, pSEIR, psyms, psymc, tseek,tdelay,sens, close_thresh);
    %[I_t,cumI_t, detI_t, detcumI_t, tfirst, Ifirst, tclose, Iclose] = sim_school(I0,N,tvec, Reffi,pSEIR, psym, tseek,tdelay,sens, close_thresh); 
    tfirstci(i) = tfirst;
    Ifirstci(i) = Ifirst;
    tcloseci(i) = tclose;
    Icloseci(i) = Iclose;
end
 %calculate upperbound (1 in 10 reporting rate)
tfirstupi = [];
Ifirstupi = [];
tcloseupi = [];
Icloseupi = [];
for i = 1:length(Reffveccont)
    Reffsi=Reffveccont(i);
    I0s = 2*I0s0;
    I0c = (I0s0/1000)*Nc;
    [Is_t, Ic_t, cumIs_t, cumIc_t, detIs_t, detIc_t, detcumIs_t, detcumIc_t, tfirst, Ifirst, tclose, Iclose] = sim_school_com(I0s, I0c, Ns,Nc, tvec, Reffsi,Reffc,int, pSEIR, psyms, psymc, tseek,tdelay,sens, close_thresh);
    tfirstupi(i) = tfirst;
    Ifirstupi(i) = Ifirst;
    tcloseupi(i) = tclose;
    Icloseupi(i) = Iclose;
end
% calculate lowerbound (1 in 3 reporting rate)
 tfirstloi = [];
 Ifirstloi = [];
 tcloseloi = [];
 Icloseloi = [];
for i = 1:length(Reffveccont)
    Reffsi=Reffveccont(i);
    I0s = 0.6*I0s0;
    I0c = (I0s0/1000)*Nc;
    [Is_t, Ic_t, cumIs_t, cumIc_t, detIs_t, detIc_t, detcumIs_t, detcumIc_t, tfirst, Ifirst, tclose, Iclose] = sim_school_com(I0s, I0c, Ns,Nc, tvec, Reffsi,Reffc,int, pSEIR, psyms, psymc, tseek,tdelay,sens, close_thresh);
    tfirstloi(i) = tfirst;
    Ifirstloi(i) = Ifirst;
    tcloseloi(i) = tclose;
    Icloseloi(i) = Iclose;
end


% 2 panels below varying Reff

subplot(2,2,3)
for i = 1:length(Reffvec)
%plot(tvec, 100.*10.*detI_ti(:,i)./N, '-', 'Linewidth',2, 'color', colorsets2(i,:))
plot(tvec, 100.*cumI_ti(:,i)./Ns, '-', 'Linewidth',2, 'color', colorsets2(i,:))
hold on 
end
xlabel('time (days)')
xlim([0 tvec(end)])
%ylim([0 50])
ylabel('% of school infected')
%ylabel('% of school quarantined')
legend('R_0=1.5','R_0=2.5', 'R_0=3', 'R_0=3.5', 'R_0=4', 'R_0=5', 'Location', 'NorthWest')
legend boxoff
set(gca,'FontSize',18,'LineWidth',1.5)
title('Effect of R_0 on percent infected')
%title('Percent quarantined')
ylim([0 100])

subplot(2,2,4)
plot(Reffveccont, tcloseci, 'k-', 'LineWidth', 2)
hold on
plot(Reffveccont, tcloseupi, 'k--', 'LineWidth',2)
plot(Reffveccont, tcloseloi, 'k--', 'LineWidth',2)
for i =1:length(I0vec)
plot(Reffvec(i), tclosei(i), '*', 'LineWidth', 10, 'color', colorsets2(i,:))
hold on
end
xlabel('R_0')
ylabel('Time to close (days)')
set(gca,'FontSize',18,'LineWidth',1.5)
title('Effect of R_0 on time to close')
%ylim([5 max(tcloseloi)+1])
xlim([Reffveccont(1) Reffveccont(end)])
%% Remake next figure - discrepancy between detected and true infections

% Then vary the transmission rate
I0s0 = 5; % national average prevalence
Reffvec = [1.5, 2.5, 3, 3.5, 4, 5 ]; % corresponding to low, medium, and high prevalence
for i = 1:length(Reffvec)
    Reffsi=Reffvec(i);
    I0c = (I0s0/1000)*Nc;
    I0s = I0s0;
    [Is_t, Ic_t, cumIs_t, cumIc_t, detIs_t, detIc_t, detcumIs_t, detcumIc_t, tfirst, Ifirst, tclose, Iclose] = sim_school_com(I0s, I0c, Ns,Nc, tvec, Reffsi,Reffc,int, pSEIR, psyms, psymc, tseek,tdelay,sens, close_thresh);
    %[I_t,cumI_t, detI_t, detcumI_t, tfirst, Ifirst, tclose, Iclose] = sim_school(I0,N,tvec, Reffi,pSEIR, psym, tseek,tdelay,sens, close_thresh); 
    I_ti(:,i)=Is_t;
    cumI_ti(:,i) = cumIs_t;
    detI_ti(:,i)=detIs_t;
    detcumI_ti(:,i) = detcumIs_t;
    cumcases(:,i) = detcumIs_t;
    tfirsti(i) = tfirst;
    Ifirsti(i) = Ifirst;
    tclosei(i) = tclose;
    Iclosei(i) = Iclose;
end

i = 2;

figure;
subplot(1,2,1)
plot(tvec, cumI_ti(:,i), '-', 'Linewidth',2, 'color', colorsets2(2,:))
hold on
plot(tvec, detcumI_ti(:,i), '--', 'Linewidth',2, 'color', colorsets2(2,:))
plot([tvec(1) tvec(end)], [1 1], 'k--', 'LineWidth',2)
plot([tvec(1) tvec(end)], [close_thresh*Ns/100, close_thresh*Ns/100], 'k-', 'LineWidth',2)
plot([tfirsti(i) tfirsti(i)], [0 Ifirsti(i)], 'k-', 'LineWidth',2)
plot([tclosei(i) tclosei(i)], [0 Iclosei(i)], 'k-', 'LineWidth',2)
xlabel('time (days)')
xlim([0 tvec(end)])
ylim([0 135])
ylabel('Infections')
legend('True infections', 'Detected cases', 'Location', 'NorthWest')
legend boxoff
title('Relationship between true infections and detected cases')
set(gca,'FontSize',18,'LineWidth',1.5)

ifirst = find(detcumIs_ti(:,i)>1,1);
subplot(1,2,2)
plot(tvec(ifirst:end)-tfirsti(i), cumI_ti(ifirst:end,i)./detcumI_ti(ifirst:end,i), '-', 'Linewidth',2, 'color', colorsets2(2,:))
hold on
xlabel('time since first detected case (days)')
xlim([0 tvec(end)-tfirsti(i)])
%ylim([0 6])
%ylim([0 (close_thresh*N/100)+35])
ylabel('True infections per detected case')
title('True infections per case detected')
set(gca,'FontSize',18,'LineWidth',1.5)
%% Last figure...
% Heatmaps of transmission vs initial prevalence colored by time to close
I0vec =linspace(0.5, 10,40);
Reffvec = linspace(1.1, 5,40);
tvec = 0:0.1:120;
[REFF,INIT] = meshgrid(Reffvec, I0vec); % big ol grid of parameters
REFFflat = reshape(REFF,1,[]);
INITflat = reshape(INIT,1, []);
Nc = 4000; 
Ns = 1000;

for i = 1:length(REFFflat)
    Reffsi = REFFflat(i);
    I0si = INITflat(i);
    I0c = (I0si/1000)*Nc;
    [Is_t, Ic_t, cumIs_t, cumIc_t, detIs_t, detIc_t, detcumIs_t, detcumIc_t, tfirst, Ifirst, tclose, Iclose] = sim_school_com(I0si, I0c, Ns,Nc, tvec, Reffsi,Reffc,int, pSEIR, psyms, psymc, tseek,tdelay,sens, close_thresh);
    tfirstj(i) = tfirst;
    if isnan(tfirst);
        tfirstj(i) = 180;
    end
    Ifirstj(i) = Ifirst;
    if isnan(Ifirst)
        ifirstj(i)=1;
    end
    tclosej(i) = tclose;
    if isnan(tclose)
        tclosej(i)=180;
    end
    Iclosej(i) = Iclose;
    if isnan(Iclose)
        Iclosej(i) = cumIs_t(end);
    end
end


 TFIRST = reshape(tfirstj, size(REFF));
 IFIRST = reshape(Ifirstj, size(REFF));
 TCLOSE = reshape(tclosej, size(REFF));
 ICLOSE = reshape(Iclosej, size(REFF));
%% Plot the heatmaps
figure;
subplot(1,2,2)
minmax = @(x)([min(x) max(x)]);
imagesc(minmax(1*(Reffvec)),minmax(I0vec),TCLOSE);
V = linspace(0,100,11);
V2 = horzcat(linspace(0,15,4), linspace(20,100,9));
[C,h]=contourf(REFF,INIT,TCLOSE); clabel(C,h);
[C2,h2]=contourf(REFF,INIT,TCLOSE, V2); 
clabel(C2,h2);
[x,y,z] = C2xyz(C);
j=find(z==100);
R0crit=cell2mat(x(j));
I0crit =cell2mat(y(j));
%imagesc(minmax(1*(Reffvec)),minmax(I0vec),TCLOSE);
hold on
plot(x{j}, y{j}, 'w-', 'LineWidth', 5)
colorbar
colormap(flipud(jet)); 
caxis([7 100])
xlabel('School R_0');
ylabel('Initial infected per 1000');
title('Time to close (days)')
set(gca,'FontSize',20,'LineWidth',1.5)

subplot(1,2,1)
imagesc(minmax(1*(Reffvec)),minmax(I0vec),TFIRST);
%imagesc(minmax(1*(alphav)),minmax(tremv),PCTINF);
V = linspace(0, 100, 11);
[C,h]=contourf(Reffvec,I0vec,TFIRST, V2); clabel(C,h);
%imagesc(minmax(1*(Reffvec)),minmax(I0vec),TFIRST);
colorbar
colormap(flipud(jet)); 
caxis([7 100])
xlabel('School R_0');
ylabel('Initial infected per 1000');
title('Time to first detected case (days)')
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
subplot(1,2,1)
imagesc(minmax(1*(Reffvec)),minmax(I0vec),IFIRST);
%imagesc(minmax(1*(alphav)),minmax(tremv),PCTINF);
[C,h]=contourf(Reffvec,I0vec,IFIRST); clabel(C,h);
colorbar
colormap(jet); 
xlabel('School R_0');
ylabel('Initial infected per 1000');
title('Number infected in school per 1000 at first detected case')
set(gca,'FontSize',18,'LineWidth',1.5)

% Make heatmap of number of additional community cases


%

Rcvec = linspace(0.8, 1.5,17);
Rsvec= linspace(2, 5, 17);
tvec = 0:1:120;
[RS, RC] = meshgrid(Rsvec,Rcvec); % big ol grid of parameters
RSflat = reshape(RS,1,[]);
RCflat = reshape(RC,1, []);
I0s = 3;
I0c = (I0s/1000)*Nc;

for i = 1:length(RSflat)
     Reffsi = RSflat(i);
     Reffci = RCflat(i);
     % cases with school open
    [Is_t, Ic_t, cumIs_t, cumIc_t, detIs_t, detIc_t, detcumIs_t, detcumIc_t, tfirst, Ifirst, tclose, Iclose] = sim_school_com(I0s, I0c, Ns,Nc, tvec, Reffsi,Reffci,int, pSEIR, psyms, psymc, tseek,tdelay,sens, close_thresh);
    totIc_schooli(i) = 1000*cumIc_t(end)./Nc;
    
      % cases with school closed (R0s = R0c)--
    [Is_t, Ic_t, cumIs_t, cumIc_t, detIs_t, detIc_t, detcumIs_t, detcumIc_t, tfirst, Ifirst, tclose, Iclose] = sim_school_com(I0s, I0c, Ns,Nc, tvec, Reffci,Reffci,int, pSEIR, psyms, psymc, tseek,tdelay,sens, close_thresh);
    totIc_closedi(i) = 1000*cumIc_t(end)./Nc;
    addcasesi(i) = totIc_schooli(i)-totIc_closedi(i); 
    
    % cases with school open for 45 days and then closed...
    % take the sum of the cases with school open for 45 days + school
    % closed for the remainder
    % cases from schools open
    tvec1 = 0:1:45;
    [Is_t, Ic_t, cumIs_t, cumIc_t, detIs_t, detIc_t, detcumIs_t, detcumIc_t, tfirst, Ifirst, tclose, Iclose] = sim_school_com(I0s, I0c, Ns,Nc, tvec1, Reffsi,Reffci,int, pSEIR, psyms, psymc, tseek,tdelay,sens, close_thresh);
    totIc_school45(i) = 1000*cumIc_t(end)./Nc;
    I0s2 = Is_t(end);
    I0c2 = Ic_t(end);
    tvec2 = 0:1:74;
    [Is_t, Ic_t, cumIs_t, cumIc_t, detIs_t, detIc_t, detcumIs_t, detcumIc_t, tfirst, Ifirst, tclose, Iclose] = sim_school_com(I0s2, I0c2, Ns,Nc, tvec2, Reffci,Reffci,int, pSEIR, psyms, psymc, tseek,tdelay,sens, close_thresh);
    totIc_close75(i) = 1000*cumIc_t(end)./Nc;
    totIcopenclose(i) = totIc_school45(i) + totIc_close75(i);
end

ADDCASES = reshape(addcasesi, size(RS));
CASESOPEN = reshape(totIc_schooli, size(RS));
CASESOPENCLOSE = reshape(totIcopenclose, size(RS));
CASESCLOSED = reshape(totIc_closedi, size(RS));

%%







subplot(1,2,2)
minmax = @(x)([min(x) max(x)]);
imagesc(minmax(Rsvec),minmax(Rcvec),ADDCASES);
V = linspace(0,500,11);
[C,h]=contourf(Rsvec,Rcvec, ADDCASES,V); 
clabel(C,h);
colorbar
colormap(jet); 
xlabel('School R_0');
ylabel('Community R_0');
title('Additional community cases per 1000')
set(gca,'FontSize',18,'LineWidth',1.5)

figure;
subplot(1,3,1)
minmax = @(x)([min(x) max(x)]);
imagesc(minmax(Rsvec),minmax(Rcvec),CASESCLOSED);
V = linspace(0,600,13);
[C,h]=contourf(Rsvec,Rcvec, CASESCLOSED,V); 
clabel(C,h);
colorbar
caxis([0 600])
colormap(jet); 
xlabel('School R_0');
ylabel('Community R_0');
title('Schools closed entirely')
set(gca,'FontSize',18,'LineWidth',1.5)

subplot(1,3,2)
minmax = @(x)([min(x) max(x)]);
imagesc(minmax(Rsvec),minmax(Rcvec),CASESOPENCLOSE);
V = linspace(0,600,13);
[C,h]=contourf(Rsvec,Rcvec, CASESOPENCLOSE,V); 
clabel(C,h);
colorbar
caxis([0 600])
colormap(jet); 
xlabel('School R_0');
ylabel('Community R_0');
title('School closed at 45 days')
set(gca,'FontSize',18,'LineWidth',1.5)

subplot(1,3,3)
minmax = @(x)([min(x) max(x)]);
imagesc(minmax(Rsvec),minmax(Rcvec),CASESOPEN);
V = linspace(0,600,13);
[C,h]=contourf(Rsvec,Rcvec, CASESOPEN,V); 
clabel(C,h);
colorbar
caxis([0 600])
colormap(jet); 
xlabel('School R_0');
ylabel('Community R_0');
title('Schools remain open')
set(gca,'FontSize',18,'LineWidth',1.5)