%% Topological Sampler: Sampling possible 3-node connections
clear variables;
close all;

curdir = 'C:\Users\rmish\OneDrive - Johns Hopkins\Schulman Lab\Topological Sampling Project\5.26.20_Robustness_Analysis';
addpath(curdir,'Model_functions');

N = 5000; %number of samples
%to generate LHS random samples:
% Param_mat = LHS_Sampler_allrandom(N);
ortho_nodes = 3;
Param_mat = table2array(readtable('11508_inputparameters.txt'));

Gen_mat = GTEval_V4(Param_mat, ortho_nodes);

for j = 1:length(Gen_mat)
    hold on;
    plot(Gen_mat(j,:))
end

mat_runs = length(Gen_mat);

pulser = 0;
c=1;
thresh = 0.1;
for i = 1:mat_runs
    [pc_min, ~] = pulsecounter(Gen_mat(i,:), thresh);
    if pc_min
        pulser(c) = i;
        c = c+1;
    end
end
figure;
for j = 1:length(pulser)
    hold on;
    plot([0:0.01:3.23],Gen_mat(pulser(j),:))
end

save('Gen_mat_updated_6.3.20.mat', 'Gen_mat');
genmatpulsesoln = Gen_mat(pulser,:);
param_mat_pulse_soln = Param_mat(pulser,:);
save('Gen_mat_pulser_updated_6.3.20.mat','Gen_mat');

save('Param_mat_pulser_updated_6.3.20.mat', 'param_mat_pulse_soln');