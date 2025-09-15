% t = tiledlayout(1, 2, "TileSpacing","tight");
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 40 30]);
xlimit = [295 500];
minuteSpacing = 60;
RunDuration = 180;
xTicksSeq = 0:minuteSpacing:RunDuration;
xTicks = EquilibriumDuration+xTicksSeq;

% time = t;

tileOutput = tiledlayout(3,1); % 'TileSpacing','tight');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0,0,1.37*2,1.5*2]);

nexttile
plot(t, x_conc(:,4), LineWidth=2)
grid on
xlim(xlimit)
xticks(xTicks)
xticklabels([])
% xticklabels(xTicksSeq)
yticks(0:0.5e-3:1e-3)
yticklabels([])
ylim([0 1.2e-3])

nexttile
plot(t, x_conc(:,28), LineWidth=2)
grid on
xlim(xlimit)
xticks(xTicks)
xticklabels([])
% xticklabels(xTicksSeq)
yticks(0:1.5e-9:4.5e-9)
yticklabels([])
ylim([0 5e-9])
% hold on
% plot(t, x_conc(:,27)+x_conc(:,28))

nexttile
plot(t,x(:,44), LineWidth=2)
grid on
xlim(xlimit)
xticks(xTicks)
xticklabels([])
% xticklabels(xTicksSeq)
yticks(0:1e-5:3e-5)
yticklabels([])
ylim([0 3e-5])

%%
FileName = sprintf('Graphical Abstract, simulation - %s', datestr(now, 'yy-mm-dd HH-MM-SS'));
FilePath = FilePathCumulative;

prompt = 'Save figure? Y/N: ';
SaveFile = input(prompt,'s');

if SaveFile == 'Y'
    figfile1 = fullfile(FilePath, FileName);
    saveas(tileOutput, figfile1, 'epsc' );
else
end

%%
% t = tiledlayout(tiledLayout(1), tiledLayout(2), "TileSpacing","tight");
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 40 30]);
% 
% for k = 1:size(mpsaValues,1)
%     
%     nexttile
%     cdfplot(passTest_success(k,:))
%     hold on
%     cdfplot(passTest_fail(k,:))


%%
plot(t, x_conc(:,6))

%%
x = linspace(0,30);
y1 = sin(x/2);
y2 = sin(x/3);
y3 = sin(x/4);

% Plot into first tile three times
tiledlayout('flow')
nexttile
plot(x,y1)