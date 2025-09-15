function [CostIndividual, CostTotal, mps_integral_grams] = plot_mpsaCriteria_231006(figNum, t, x_conc, x_mass, ExpData, ...
    CostParameters, RunDuration, EquilibriumDuration, FilePath, FileName, kValues, interventionDuration)
% Plotting function that creates Figure 2a. The plot simulates plasma 
% leucine, plasma insulin, p70S6K, and Protein Balance with associated 
% experimental data and costs.
%   [CostIndividual, CostTotal, mps_integral_grams] = 
%       plot_mpsaCriteria(figNum, t, x_conc, x_mass, ExpData, 
%       CostParameters, RunDuration, EquilibriumDuration, FilePath, 
%       FileName, kValues, interventionDuration)
%   figNum = figure number for the MATLAB plot
%   t = scalar time outputted from the ODE solver
%   x_conc = column vector of all species across the simulated duration in 
%   units of mol/L as outputted from the ODE solver
%   x_mass = column vector of all species across the simulated duration in 
%   units of moles as outputted from the ODE solver
%   ExpData = extracte experimental data that is plotted with the
%   corresponding simualted time-course
%   CostParameters = Parameters with corresponding experimental data 
%   used to calculate RMSE
%   RunDuration = simulated run duration (minutes)
%   EquilibriumDuration = duration of the equilibrium period to allow model
%   to acheive steady state prior to stimulation (minutes)
%   FilePath = file directory where figure will be saved
%   FileName = name of figure to be saved in the FilePath 
%   kValues = k-value vector use to simulate model
%   interventionDuration = the length of the intervention duration from the
%   extracted experimental data protocol, which is used to select the
%   plotted x-axis (time, minutes)

%% Simulation duration to plot
simDuration = 200;
    
%% Setting up subTightPlot function (allows for better plotting)

gap = [0.16 0.045]; % gap between subplots [vert, horiz]
marg_h = [0.1/1.1 0.08/1.3]; % Vertical exterior borders [bottom, top]
marg_w = [0.1/2.6 0.01/4]; % Horizontal exterior borders [left, right]

 subplotTight = @(m,n,p) subtightplot (m, n, p, gap, marg_h, marg_w);
    % @(m,n,p): m&n = the number of subplots across and down, p = location
    % of subplot in overall plot
    
%% Data cleaning   
% Position of species in data frame
AktS473_pos = 21;
Akt_pos = 19;
AktT308_pos = 20;
AktS473T308_pos = 22;
p70S6K_pos = 27;
p70S6KT389_pos = 28;
LeucinePlasma_pos = 4;
LeucineInt_pos = 6;
Insulin_pos = 3;
Protein_pos = 8;
kicPlasma_pos = 5;
kicInt_pos = 7;
F_ma_pos = 41;
F_m0_pos = 42;
MPS_pos = 44; %rate of MPS
MPB_pos = 45;
FSR_pos = 46; %total synthesized leucine (calculated from FSR)

t_endEq = find(t==EquilibriumDuration);
% p-Akt(S) incorporates both Akt species with phospho serine residues
AktS473_iv = sum(x_conc(t_endEq(1), [AktS473_pos, AktS473T308_pos]));
p70S6KT389_iv = x_conc(t_endEq(1), p70S6KT389_pos);

%% Sorting Experimental Data
ExpDataPlot = nan(length(CostParameters), 13, 4);
ExpDataPlot(:,1,1) = CostParameters';

for n = 1:length(CostParameters)
    if ExpDataPlot(n,1,1) == AktS473_pos
        rows = find(ExpData(:,1) == CostParameters(n));
        ExpDataPlot(n, 1:length(rows),2) = ExpData(rows, 2) + EquilibriumDuration; %Adjusted to Equilibrium period
        ExpDataPlot(n, 1:length(rows),3) = ExpData(rows, 3) *AktS473_iv;
        ExpDataPlot(n, 1:length(rows),4) = sqrt(AktS473_iv^2 * ExpData(rows, 4).^2); % Transforming SE data 
            
    elseif ExpDataPlot(n,1,1) == p70S6KT389_pos
        rows = find(ExpData(:,1) == CostParameters(n));
        ExpDataPlot(n, 1:length(rows),2) = ExpData(rows, 2) + EquilibriumDuration; %Adjusted to Equilibrium period
        ExpDataPlot(n, 1:length(rows),3) = ExpData(rows, 3) *p70S6KT389_iv;
        ExpDataPlot(n, 1:length(rows),4) = sqrt(p70S6KT389_iv^2 * ExpData(rows, 4).^2); % Transforming SE data
    
    else
        rows = find(ExpData(:,1) == CostParameters(n));
        ExpDataPlot(n, 1:length(rows),2) = ExpData(rows, 2) + EquilibriumDuration; %Adjusted to Equilibrium period
        ExpDataPlot(n, 1:length(rows),3) = ExpData(rows, 3);
        ExpDataPlot(n, 1:length(rows),4) = ExpData(rows, 4);
    end
end

%% FSR calculation
% FSR/MPS; integral for FSR comparison
mps = x_mass(:,44);
t_endRun = find(t>=EquilibriumDuration+interventionDuration & t<EquilibriumDuration+interventionDuration+20);
mps_integral = cumtrapz(t(t_endEq(1):t_endRun(1)),...
    mps(t_endEq(1):t_endRun(1)));
mps_integral=mps_integral(end);
LeucineMolarMass = 1/131.17; %mol/g
mps_integral_grams = mps_integral/LeucineMolarMass;

%% x_conc vector updated to include mole/min values for 3-pool and muscle balance parameters
% required for the 'x' parameter inputted in the costCalculation function
x_conc_moleMin = x_conc;
x_conc_moleMin(:,41) = x_mass(:,41); % F_ma in units of mole/min
x_conc_moleMin(:,42) = x_mass(:,42); % F_m0 in units of mole/min
x_conc_moleMin(:,43) = nan; % NB: not correct calculation
x_conc_moleMin(:,44) = x_mass(:,44); % MPS in units of mole/min
x_conc_moleMin(:,45) = x_mass(:,45); % MPB in units of mole/min

%% Root mean squares, normalized
[Title, ~, ~, CostIndividual, CostTotal] = CostCalculation_230920(t, x_conc_moleMin, [AktS473_pos, AktS473T308_pos], p70S6KT389_pos, ...
    EquilibriumDuration, CostParameters, ExpData, mps_integral_grams);

%% Simulation plot

x_lim = [EquilibriumDuration-5, EquilibriumDuration+simDuration];
minuteSpacing=60;
xTicksSeq = 0:minuteSpacing:RunDuration;
xTicks = EquilibriumDuration+xTicksSeq;

opacity = 0.25;  
Angle = 0;%45;
jitter = 2;

length_x = 12;
height_y = 2;

expand_x = 1.4;
expand_y = 1.1;

fontSize_title = 20.5;
fontSize_axes = 16;
fontSize_units = 13.5;
fontSize_inPlotText = 14;

capSize_error = 12;
lineWidth_error = 1.75;

figure(figNum(1));

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 30*expand_x 15*expand_y]); %x_width=10cm y_width=15cm

%% Leucine
subplotTight(height_y, length_x, 1:4)
lines = LeucinePlasma_pos;
    plot (t,x_conc(:,lines), '-', 'LineWidth', 1.5, 'Color',[0 0 1])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;    
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 1.4e-3])
    ylabel('Concentration (mol/L)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    yticks([0 : 0.3e-3 : 1.4e-3])
    title('Leucine', 'fontSize', fontSize_title)
    xtickangle(Angle)

hold on 
    plot(t,x_conc(:,LeucineInt_pos), '-', 'LineWidth',1.25, 'Color',[1 0 0])

% --- Plasma Leucine
% Basal Concentration
plLeu_basal = 1.23e-4; plLeu_basal_range = 0.25;
yneg = plLeu_basal*(plLeu_basal_range);
ypos = plLeu_basal*(plLeu_basal_range);

hold on
errorbar(EquilibriumDuration-jitter, plLeu_basal, ...
    yneg, ypos, ...
    'Color',[0.1 0.1 0.4], 'LineWidth', lineWidth_error, ...
    'CapSize',capSize_error);

% Peak concentration and timing
plLeu_peakConc = 8.9e-4; plLeu_peakConc_range = 0.5;
plLeu_peakTime = 45+EquilibriumDuration; plLeu_peakTime_range = [-15 30];
range = (plLeu_peakConc-plLeu_basal)*(1-plLeu_peakConc_range);

hold on
errorbar(plLeu_peakTime, plLeu_peakConc, ...
    range, range, plLeu_peakTime_range(1), plLeu_peakTime_range(2), ...
    'Color',[0.1 0.1 0.4], 'LineWidth', lineWidth_error, ...
    'CapSize',capSize_error);

% Return to basal
plLeu_endConc = 2.31e-4; plLeu_endConc_range = 2;
yneg = plLeu_endConc*(1/plLeu_endConc_range)-plLeu_endConc;
ypos = plLeu_endConc*plLeu_endConc_range-plLeu_endConc;

hold on
errorbar(EquilibriumDuration+interventionDuration-jitter, ...
    plLeu_endConc, yneg, ypos, ...
    'Color',[0.1 0.1 0.4], 'LineWidth', lineWidth_error, ...
    'CapSize',capSize_error);
hold off

% --- Intracellular leucine
% Basal concentration
intLeu_basal = 1.28e-4; intLeu_basal_range = 0.5;
yneg = intLeu_basal*(intLeu_basal_range);
ypos = intLeu_basal*(intLeu_basal_range);

hold on
errorbar(EquilibriumDuration+jitter, intLeu_basal, ...
    yneg, ypos, ...
    'Color',[0.4 0.1 0.1], 'LineWidth', lineWidth_error, ...
    'CapSize',capSize_error);

% Peak concentration and timing
intLeu_peakConc = 2.8e-4; intLeu_peakConc_range = 0.5;
intLeu_peakTime = 120+EquilibriumDuration; intLeu_peakTime_range = [-60 30];
range = (intLeu_peakConc-intLeu_basal)*(1-intLeu_peakConc_range);

hold on
errorbar(intLeu_peakTime, intLeu_peakConc, ...
    range, range, intLeu_peakTime_range(1), intLeu_peakTime_range(2), ...
    'Color',[0.4 0.1 0.1], 'LineWidth', lineWidth_error, ...
    'CapSize',capSize_error);
hold off

legend({'Plasma', 'Intracellular'}, 'Location','northeast');

%% Insulin
subplotTight(height_y, length_x, 5:8)
    plot(t,(x_conc(:,Insulin_pos)), '-', 'LineWidth',1.5, 'Color',[0 0 1])
    grid on 
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 1.3e-10])
    ylabel('Concentration (mol/L)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('Plasma Insulin', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 

% Basal Concentration
plIns_basal = 3.3e-11; plIns_basal_range = 0.25;
yneg = plIns_basal*(plIns_basal_range);
ypos = plIns_basal*(plIns_basal_range);

hold on
errorbar(EquilibriumDuration-jitter, plIns_basal, ...
    yneg, ypos, ...
    'Color',[0.1 0.1 0.4], 'LineWidth', lineWidth_error, ...
    'CapSize',capSize_error);

% Peak concentration and timing
plIns_peakConc = 9.2e-11; plIns_peakConc_range = 0.5;
plIns_peakTime = 30+EquilibriumDuration; plIns_peakTime_range = [-15 30];
range = (plIns_peakConc-plIns_basal)*(1-plIns_peakConc_range);

hold on
errorbar(plIns_peakTime, plIns_peakConc, ...
    range, range, plIns_peakTime_range(1), plIns_peakTime_range(2), ...
    'Color',[0.1 0.1 0.4], 'LineWidth', lineWidth_error, ...
    'CapSize',capSize_error);

% Return to basal
plIns_endConc = 1.7e-11; %simulated value
plIns_endConc_range = 0.5;
ypos = (plIns_peakConc-plIns_basal)*plIns_endConc_range + plIns_basal;

hold on
errorbar(EquilibriumDuration+interventionDuration-jitter, ...
    plIns_endConc, plIns_endConc, ypos-plIns_endConc, ...
    'Color',[0.1 0.1 0.4], 'LineWidth', lineWidth_error, ...
    'CapSize',capSize_error);
hold off

%% 3-pool parameters (Fm,a & Fm,0)
subplotTight(height_y, length_x, 9:12)
    plot (t,x_mass(:,F_ma_pos), '-', 'LineWidth',1.5, 'Color',[0 0 1]) % F_m,a
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([-1E-5 18.5E-5])
    ylabel('Moles/min', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('3-Pool Model', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on
    plot (t,x_mass(:,F_m0_pos), '-', 'LineWidth',1.5, 'Color',[1 0 0]) % F_m,0

% --- Fma    
% Basal Concentration
F_ma_basal = 2.29e-5; F_ma_basal_range = 0.25;
yneg = F_ma_basal*(F_ma_basal_range);
ypos = F_ma_basal*(F_ma_basal_range);

hold on
errorbar(EquilibriumDuration, F_ma_basal, ...
    yneg, ypos, ...
    'Color',[0.1 0.1 0.4], 'LineWidth', lineWidth_error, ...
    'CapSize',capSize_error);

% Peak concentration and timing
F_ma_peakConc = 1.28e-4; F_ma_peakConc_range = 0.5;
F_ma_peakTime = 45+EquilibriumDuration; F_ma_peakTime_range = [-15 30];
range = (F_ma_peakConc-F_ma_basal)*(1-F_ma_peakConc_range);

hold on
errorbar(F_ma_peakTime, F_ma_peakConc, ...
    range, range, F_ma_peakTime_range(1), F_ma_peakTime_range(2), ...
    'Color',[0.1 0.1 0.4], 'LineWidth', lineWidth_error, ...
    'CapSize',capSize_error);

% Return to basal
F_ma_endConc = 2.91e-5 ; %simulated value
F_ma_endConc_range = 0.5;
ypos = (F_ma_peakConc-F_ma_basal)*F_ma_endConc_range + F_ma_basal;

hold on
errorbar(EquilibriumDuration+interventionDuration-jitter, ...
    F_ma_endConc, F_ma_endConc, ypos-F_ma_endConc, ...
    'Color',[0.1 0.1 0.4], 'LineWidth', lineWidth_error, ...
    'CapSize',capSize_error);
hold off

legend({'F_{m,a}', 'F_{m,0}'},'NumColumns',2);

%% Akt
subplotTight(height_y, length_x, 13:16)
Akt_total = x_conc(:,Akt_pos) + x_conc(:,AktT308_pos) + x_conc(:,AktS473_pos) + x_conc(:,AktS473T308_pos);
    plot(t,Akt_total, '-', 'LineWidth', 1.5, 'Color',[0.5 0.5 0.5])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 Akt_total(1,1)*1.1])
    ylabel('Concentration (mol/L)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    yticks(0 : 0.4e-8 : Akt_total(1,1)*1.1)
    title('Akt', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
    % plotting both Akt species with phosphorylated serine
    % [p-Akt(S) & p-Akt(S,T)]
    Akt_pS = x_conc(:,AktS473_pos) + x_conc(:,AktS473T308_pos);
    plot (t,Akt_pS, '-', 'LineWidth', 1.5, 'Color',[0 0 1]);
hold on 
% Basal Concentration
pAkt_S_basal = 0.49e-9; % Update, sum of both p-S species 
pAkt_S_basal_acceptRange = [0.5, 50];
yneg = max(Akt_total)*pAkt_S_basal_acceptRange(1)/100;
ypos = max(Akt_total)*pAkt_S_basal_acceptRange(2)/100;
middle = (yneg+ypos)/2;

hold on
errorbar(EquilibriumDuration, middle, ...
    middle-yneg, ypos-middle, ...
    'Color',[0.1 0.1 0.4], 'LineWidth', lineWidth_error, ...
    'CapSize',capSize_error);

text(EquilibriumDuration+5, middle*1.25, ...
    {'Basal Phos:', '0.5-50% total Akt'}, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left',...
        'FontSize', 10)

% Peak concentration and timing
pAkt_S_peakFC_acceptRange = [1.1, 2.5];
pAkt_S_peakTime = 60+EquilibriumDuration; pAkt_S_peakTime_range = [-15 60];
yneg = pAkt_S_basal*pAkt_S_peakFC_acceptRange(1);
ypos = pAkt_S_basal*pAkt_S_peakFC_acceptRange(2);
middle = (yneg+ypos)/2;

hold on
errorbar(pAkt_S_peakTime, middle, ...
    middle-yneg, ypos-middle, ...
    pAkt_S_peakTime_range(1), pAkt_S_peakTime_range(2), ...
    'Color',[0.1 0.1 0.4], 'LineWidth', lineWidth_error, ...
    'CapSize',capSize_error);
hold off

text(pAkt_S_peakTime+5, middle*2.2, ...
    {'Peak Phos:', '1.1 - 2.5x basal'}, ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left',...
        'FontSize', 10)

legend({'Total Akt', 'p-Akt^{S473}'});

%% p70S6K
subplotTight(height_y, length_x, 17:20)
p70S6K_total = x_conc(:,p70S6K_pos) + x_conc(:,p70S6KT389_pos);
lines = p70S6KT389_pos;
    plot(t, p70S6K_total, '-', 'LineWidth', 1.5, 'Color',[0.5 0.5 0.5])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 p70S6K_total(1,1)*1.1])
    ylabel('Concentration (mol/L)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('p70S6K', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
    plot (t,x_conc(:,lines), '-', 'LineWidth',1.5, 'Color',[0 0 1]);

% Basal Concentration
p70S6K_T_basal = 1.14e-9; p70S6K_T_basal_acceptRange = [1, 40];
yneg = max(p70S6K_total)*p70S6K_T_basal_acceptRange(1)/100;
ypos = max(p70S6K_total)*p70S6K_T_basal_acceptRange(2)/100;
middle = (yneg+ypos)/2;

hold on
errorbar(EquilibriumDuration-jitter, middle, ...
    middle-yneg, ypos-middle, ...
    'Color',[0.1 0.1 0.4], 'LineWidth', lineWidth_error, ...
    'CapSize',capSize_error);

text(EquilibriumDuration+5, middle*1.4, ...
    {'Basal Phos:', '1-40% total', 'p70S6K'}, ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left',...
        'FontSize', 10)

% Peak concentration and timing
p70S6K_T_peakFC_acceptRange = [1.25, 3];
p70S6K_T_peakTime = 60+EquilibriumDuration; p70S6K_T_peakTime_range = [-15 60];
yneg = p70S6K_T_basal*p70S6K_T_peakFC_acceptRange(1);
ypos = p70S6K_T_basal*p70S6K_T_peakFC_acceptRange(2);
middle = (yneg+ypos)/2;

hold on
errorbar(p70S6K_T_peakTime, middle, ...
    middle-yneg, ypos-middle, ...
    p70S6K_T_peakTime_range(1), p70S6K_T_peakTime_range(2), ...
    'Color',[0.1 0.1 0.4], 'LineWidth', lineWidth_error, ...
    'CapSize',capSize_error);
hold off

text(p70S6K_T_peakTime+5, middle*1.3, ...
    {'Peak Phos:', '1.25 - 3x basal'}, ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left',...
        'FontSize', 10)

legend({'Total p70S6K', 'p-p70S6K^{T389}'}, 'Location','northeast');

%% FSR
% uses x_mass vector as MPS/MPB are in units of *moles*/min
subplotTight(height_y, length_x, 21:24)
    plot (t, x_mass(:,MPS_pos), '-', 'LineWidth',1.5, 'Color',[0 0 1]) 
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ymax = 4e-5;
    ylim([-0.1e-5 ymax])
    ylabel('Moles/min', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('Protein Balance', 'fontSize', fontSize_title)
    xtickangle(Angle)
% plot MPB
hold on
    plot(t, x_mass(:, MPB_pos), '-', 'LineWidth',1.25, 'Color',[1 0 0])                         
% plot NB
hold on
    NetBalance = x_mass(:,44) - x_mass(:,45);
    plot(t, NetBalance,  '-', 'LineWidth',1.5, 'Color',[0 1 0])

% MPS AUC shading
hold on
    area(t(t_endEq(1):t_endRun(1)), x_mass(t_endEq(1):t_endRun(1), MPS_pos), ...
        'FaceColor',[0 0 1], 'FaceAlpha',0.05)

% Peak concentration and timing
FSR_peakTime = 60+EquilibriumDuration; FSR_peakTime_range = [-15 60];

hold on
errorbar(FSR_peakTime, max(x_mass(:, MPS_pos)), ...
    FSR_peakTime_range(1), FSR_peakTime_range(2), 'horizontal',...
    'Color',[0.1 0.1 0.4], 'LineWidth', lineWidth_error, ...
    'CapSize',capSize_error);
hold off 

legend({'MPS', 'MPB', 'NB'}, 'Location','northeast', 'NumColumns',3)

% simulated FSR value (grams)
MPS_simPeriod = mps(t_endEq(1):t_endRun(1));
mps_max = max(MPS_simPeriod);
mps_max_t = find(mps == mps_max);    

rectangle('Position', [t(mps_max_t)-32, 0.5*mps_max-0.45e-5, 64, 0.9e-5], ...
    'FaceColor',[1 1 1 0.9], 'EdgeColor',[1 1 1])

text(t(mps_max_t), 0.5*mps_max, ...
    {'Simulated FSR', '0.18-0.54 g Leu'}, ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center',...
        'FontSize', fontSize_inPlotText)

%% Save plot

prompt = 'Save figure? Y/N: ';
SaveFile = input(prompt,'s');

if SaveFile == 'Y'
    figure1 = figure(figNum(1));
    figfile1 = fullfile(FilePath, FileName);
    saveas(figure1, figfile1, 'epsc' );
else
end

end
