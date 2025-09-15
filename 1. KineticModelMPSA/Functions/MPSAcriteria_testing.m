% x and t are the outputs from the calibrated model
% Both are entered directly into the criteria function to ensure the code
% is correct (i.e., the calibrated model passes the criteria)
% 
% I would run the criteria script piecewise (i.e., comment out the
% non-tested criteria) to ensure that each acceptance function runs and
% that the calibrated model passes each criteria (line 23)
% 
% Following that, alter the parameters at specific times to ensure that
% incorrect x values fail the criteria (lines 17-19)

clc

t_test = t;
x_test = x;
x_conc_test = x_conc;
% x_conc_test(4000:end,3) = x_conc_test(4000:end,3)*10;
% x_conc_test(4000:end,4) = x_conc_test(4000:end,4)*10;
% x_conc_test(3001, 44) = 3e-4;

eqDur = 300;

[passTest] = MPSAcriteria_230705(t_test, x_test, x_conc_test, eqDur)

%%
% Assessing the calibrated model against the acceptance criteria (e.g.,
% basal concentrations, peak timing and concentrations, etc)

clc
simDur = 180;

species = 41;
x_conc_exam = x_conc;
x_conc_exam = x;

% basal conc
basalConc = x_conc_exam(3001,species)
% peak conc
peakConc = max(x_conc_exam(3001:end,species))
% peak conc time
peakConc_t = t(x_conc_exam(:,species)==peakConc)-eqDur
% return to basal (conc at end of simulation)
endConc = x_conc_exam(find(t==eqDur+simDur), species)

%% phospho
% phospho-specific criteria

p_p70S6K = 28;
total_p70S6K = [27:28];
x_conc_test = x_conc;

totalProtein_basal = sum(x_conc_test(3001,total_p70S6K));
phosphoProtein_basal = x_conc_test(3001, p_p70S6K);

p70S6K_phosphoToTotal_percent = phosphoProtein_basal/totalProtein_basal*100

p_Akt_S = 21;
total_Akt = [19:22];

totalProtein_basal = sum(x_conc_test(3001,total_Akt));
phosphoProtein_basal = x_conc_test(3001, p_Akt_S);

Akt_phosphoToTotal_percent = phosphoProtein_basal/totalProtein_basal*100


%% 
% likely ignore the rest

peakConc_range = 0.50;

peakConc_rangeLow = (peakConc - basalConc) * (1-peakConc_range) ...
    + basalConc % ranges calculated using the difference between peak and basal levels
peakConc_rangeHigh = (peakConc - basalConc) * (1+peakConc_range) ...
    + basalConc

%%
(peakConc - basalConc)/2+basalConc
