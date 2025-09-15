function plot_mpsaPass_allSpecies(figNum, t, x_conc, x_mass, runDur, eqDur, kValues)

%% Species 
Stomach_pos = 1; Gut_pos = 2;
kicPlasma_pos = 5; kicInt_pos = 7;
Protein_pos = 8;
IR_B_pos = 9; IR_B_p_pos = 10; IR_B_ref_pos = 11;
IRS1_pos = 12; IRS1_pY_pos = 13; IRS1_pS_pos = 14;
PI3K_pos = 15; IRS1_PI3K_pos = 16;
PDK1_pos = 17; PDK1_p_pos = 18;
Akt_pos = 19; AktT308_pos = 20; AktS473_pos = 21; AktS473T308_pos = 22;
TSC_pos = 23; TSC_p_pos = 24;
mTORC1_inactive_pos = 25; mTORC1_active_pos = 26;
p70S6K_pos = 27; p70S6KT389_pos = 28;
PI3K_var_pos = 29; PI3K_var_p_pos = 30;
mTORC2_pos = 31; mTORC2_p_pos = 32;

%% 
figure(figNum)

minuteSpacing=60;
xTicksSeq = 0:minuteSpacing:runDur;
% xTicksSeq = -eqDur:minuteSpacing:runDur;
xTicks = eqDur+xTicksSeq;

total_lineWidth = 0.25;
species_lineWidth = 0.5;

%% Digestive tract
hold on
subplot(4,3,1)
plot(t, x_conc(:,Stomach_pos), 'color', [0 0 1])
    grid on
        ylimMax = 3e-2;
    ylim([-0.05*ylimMax ylimMax])
    xlim([0 eqDur+runDur])
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylabel('Concentration (mol/L)')
    xlabel('Time (min)')
    yticks([0 : 1e-2 : 3e-2])
    title('Digestive Tract')
hold on
    plot(t, x_conc(:,Gut_pos), '--', 'color', [1 0 0])

legend({'Stomach', 'Gut'}, 'location', 'northeast');

%% KIC
hold on
subplot(4,3,2)
plot(t, x_conc(:,kicPlasma_pos), 'LineWidth',species_lineWidth) %'color', [0 0 1])
    grid on
    xlim([0 eqDur+runDur])
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ylimMax = 3.9e-4;
    ylim([-0.05*ylimMax ylimMax])
    ylabel('Concentration (mol/L)')
    xlabel('Time (min)')
    title('KIC')
hold on 
    plot(t,x_conc(:,kicInt_pos), 'LineWidth',species_lineWidth, 'LineStyle', '--') %, 'color', [1 0 0])
legend({'Plasma', 'Intracellular'}, 'location', 'northeast');

%% Insulin receptor
IR_total = x_conc(:,IR_B_pos)+x_conc(:,IR_B_ref_pos)+x_conc(:,IR_B_p_pos);

hold on
subplot(4,3,3)
plot(t, IR_total,'Color',[0.5 0.5 0.5], 'LineWidth',total_lineWidth)
    grid on
    xlim([0 eqDur+runDur])
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ylimMax = 2.4e-9;
    ylim([-0.025*ylimMax ylimMax])
    ylabel('Concentration (mol/L)')
    xlabel('Time (min)')
    title('Insulin Receptor')
hold on
    plot(t, x_conc(:,IR_B_p_pos), 'LineWidth',species_lineWidth) %, 'color', [0 0 1])
hold on 
    plot(t, x_conc(:,IR_B_ref_pos), '--', 'LineWidth',species_lineWidth) %'color', [1 0 0])

legend({'IR-\beta total', 'p-IR-\beta^{Y}', 'IR-\beta_{Ref}'}, 'location', 'northeast');

%% IRS1
IRS1_total = x_conc(:,IRS1_pos) + x_conc(:,IRS1_pS_pos) + x_conc(:,IRS1_pY_pos) + x_conc(:,IRS1_PI3K_pos);

hold on
subplot(4,3,4)
plot(t, IRS1_total, 'Color',[0.5 0.5 0.5], 'LineWidth',total_lineWidth)
    grid on
    xlim([0 eqDur+runDur])
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ylimMax = 6.3e-10;
    ylim([-0.05*ylimMax ylimMax])
    ylabel('Concentration (mol/L)')
    xlabel('Time (min)')
    title('IRS1')
hold on
    plot(t, x_conc(:,IRS1_pS_pos), 'LineWidth',species_lineWidth) %'color', [0 0 1])
hold on
    plot(t, x_conc(:,IRS1_pY_pos), '--', 'LineWidth',species_lineWidth) %'color', [1 0 0])

legend({'IRS1 total', 'p-IRS1^S', 'p-IRS1^{Y}'}, 'location', 'northeast');

%% PI3K
PI3K_total = x_conc(:,PI3K_pos) + x_conc(:,IRS1_PI3K_pos);

hold on
subplot(4,3,5)
plot(t, PI3K_total, 'Color',[0.5 0.5 0.5], 'LineWidth',total_lineWidth)
    grid on
    xlim([0 eqDur+runDur])
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ylimMax = 5.3e-10;
    ylim([-0.05*ylimMax ylimMax])
    ylabel('Concentration (mol/L)')
    xlabel('Time (min)')
    title('PI3K')
hold on
    plot(t,x_conc(:,IRS1_PI3K_pos), 'LineWidth',species_lineWidth) %'color', [0 0 1])

legend({'PI3K total', 'PI3K p-IRS1^{Y}'}, 'location', 'northeast');

%% PDK1
PDK1_total = x_conc(:,PDK1_pos) + x_conc(:,PDK1_p_pos);

subplot(4,3,6)
plot(t, PDK1_total, 'Color',[0.5 0.5 0.5], 'LineWidth',total_lineWidth)
    grid on
    xlim([0 eqDur+runDur])
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ylimMax = 2e-9;
    ylim([-0.025*ylimMax ylimMax])
    ylabel('Concentration (mol/L)')
    xlabel('Time (min)')
    title('PDK1')
hold on
    plot(t,x_conc(:,PDK1_p_pos), 'LineWidth',species_lineWidth) %'Color',[0 0 1])

legend({'PDK1 total', 'p-PDK1'}, 'location', 'northeast');

%% Akt
Akt_total = x_conc(:,Akt_pos) + x_conc(:,AktT308_pos) + x_conc(:,AktS473_pos) + x_conc(:,AktS473T308_pos);

subplot(4,3,7)
plot(t, Akt_total, 'Color',[0.5 0.5 0.5], 'LineWidth',total_lineWidth)
    grid on
    xlim([0 eqDur+runDur])
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ylimMax = 3e-8;
    ylim([-0.05*ylimMax ylimMax])
    ylabel('Concentration (mol/L)')
    xlabel('Time (min)')
    title('Akt')
hold on
    plot(t,x_conc(:,AktT308_pos),'LineWidth',species_lineWidth) %'Color', [204 37 41]./255)
hold on
    plot(t,x_conc(:,AktS473_pos), 'LineStyle', '--','LineWidth',species_lineWidth) %'Color', [107 76 154]./255)
hold on
    plot(t,x_conc(:,AktS473T308_pos),'LineStyle', '-.','LineWidth',species_lineWidth) %'Color', [62 150 81]./255)

legend({'Total Akt', 'p-Akt^{T308}', 'p-Akt^{S473}', 'p-Akt^{S,T}'});

%% TSC
TSC_total = x_conc(:,TSC_pos) + x_conc(:,TSC_p_pos);

subplot(4,3,8)
plot(t, TSC_total, 'Color',[0.5 0.5 0.5], 'LineWidth',total_lineWidth)
    grid on
    xlim([0 eqDur+runDur])
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ylimMax = 1.6e-9;
    ylim([-0.025*ylimMax ylimMax])
    ylabel('Concentration (mol/L)')
    xlabel('Time (min)')
    title('TSC')
hold on
    plot(t,x_conc(:,TSC_p_pos),'LineWidth',species_lineWidth) %'Color', [0 0 1])

legend({'TSC total', 'p-TSC'}, 'location', 'northeast');

%% mTORC1
mTORC1_total = x_conc(:,mTORC1_inactive_pos) + x_conc(:,mTORC1_active_pos);

hold on
subplot(4,3,9)
plot(t, mTORC1_total, 'Color',[0.5 0.5 0.5], 'LineWidth',total_lineWidth)
    grid on
    xlim([0 eqDur+runDur])
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ylimMax = 7.2e-10;
    ylim([-0.05*ylimMax ylimMax])
    ylabel('Concentration (mol/L)')
    xlabel('Time (min)')
    title('mTORC1')
hold on
    plot(t,x_conc(:,mTORC1_active_pos),'LineWidth',species_lineWidth) %'Color', [0 0 1])

legend({'mTORC1 total', 'mTORC1_{active}'}, 'location', 'northeast');

%% PI3K var
PI3K_var_total = x_conc(:,PI3K_var_pos) + x_conc(:,PI3K_var_p_pos);

hold on
subplot(4,3,10)
plot(t, PI3K_var_total, 'Color',[0.5 0.5 0.5], 'LineWidth',total_lineWidth)
    grid on
    xlim([0 eqDur+runDur])
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ylimMax = 5e-10;
    ylim([-0.05*ylimMax ylimMax])
    ylabel('Concentration (mol/L)')
    xlabel('Time (min)')
    title('PI3K Variant')
hold on
    plot(t,x_conc(:,PI3K_var_p_pos),'LineWidth',species_lineWidth) %'Color', [0 0 1])

legend({'PI3K_{var} total', 'p-PI3K_{var}'}, 'location', 'northeast');

%% mTORC2
mTORC2_total = x_conc(:,mTORC2_pos) + x_conc(:,mTORC2_p_pos);

hold on
subplot(4,3,11)
plot(t, mTORC2_total, 'Color',[0.5 0.5 0.5], 'LineWidth',total_lineWidth)
    grid on
    xlim([0 eqDur+runDur])
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ylimMax = 10e-10;
    ylim([-0.05*ylimMax ylimMax])
    ylabel('Concentration (mol/L)')
    xlabel('Time (min)')
    title('mTORC2')
hold on
    plot(t,x_conc(:,mTORC2_p_pos),'LineWidth',species_lineWidth) %'Color', [0 0 1])

legend({'mTORC2 total', 'p-mTORC2^{S2481}'}, 'location', 'northeast');

%% Leucine bound to protein
t_endEq = find(t==eqDur);

% proteinPlot = x_conc(:,Protein_pos);
proteinPlot = x_conc(:,Protein_pos) - x_conc(t_endEq(1),Protein_pos);

hold on
subplot(4,3,12)
plot(t,proteinPlot, 'LineWidth',species_lineWidth)
    grid on
    xlim([0 eqDur+runDur])
    xticks(xTicks)
    xticklabels(xTicksSeq)
%     ylim([0 2*3.5e-10])
    ylabel('Concentration (mol/L)')
    xlabel('Time (min)')
    title('Leucine bound to protein (change from basal)')

end
