% clearvars -except all_top G Params
% isloaded = false;
% 
% if ~isloaded
%     load('C:\Users\rmish\OneDrive - Johns Hopkins\Schulman Lab\Topological Sampling Project\4.27.20_Topological_Sampling\Updated_Scripts\5-10-2020_fullmat\all_topologies.mat')
%     load('C:\Users\rmish\OneDrive - Johns Hopkins\Schulman Lab\Topological Sampling Project\4.27.20_Topological_Sampling\Updated_Scripts\5-10-2020_fullmat\Parameter_Soln.mat')
%     load('C:\Users\rmish\OneDrive - Johns Hopkins\Schulman Lab\Topological Sampling Project\4.27.20_Topological_Sampling\Updated_Scripts\5-10-2020_fullmat\Topology_Soln.mat')
% end

% load('Gen_mat_updated_6.3.20.mat')
% load('Param_mat_6.3.20.mat')
mat_runs = length(Gen_mat);
%     rtop = Param_mat(i, 13:21);

% pc_min = zeros(mat_runs, 1);
% pc_max = zeros(mat_runs,1);
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
%     [pc_min, pc_max] = pulsecounter(Gen_mat, thresh);

for j = 1:length(pulser)
    hold on;
    plot([0:0.01:3.23],Gen_mat(pulser(j),:))
end

gensoln = Gen_mat(pulser,:);
save('training_parameters.mat', 'gensoln')
% paramsoln = Param_mat(pulser,:);
% save('Gen_pulse_mat_6.3.20.mat', 'gensoln');
% save('Param_pulse_mat_6.3.20.mat', 'paramsoln');
% 
% 
% [x, y] = sort(pc_min, 'descend');
% 
% y10 = y(1:10);
% 
% newdir = fullfile(pwd, date);
% 
% if ~isfolder(newdir)
%     newdir = fullfile(pwd, date);
%     mkdir(newdir);
% end
% 
% for j = 1:10
%     i = y10(j);
%     rtop = all_top(i, :);
%     rtop = reshape(rtop,[3 3]);
%     %create a function that takes the topology and associated G mat and
%     %plots it automatically
%     subplot(2,1,1)
%     diplot = topology_plotter_SM2(rtop);
%     subplot(2,1,2)
%     plot(0:0.01:3,G(:,:,i)');
%     set(gcf, 'Position', get(0, 'Screensize'));
%     ylim([-0.05 1.05]);
%     xlabel('Time (hours)', 'FontSize', 20);
%     ylabel('Normalized Concentration', 'FontSize', 15);
%     set(gca, 'FontSize', 15);
%     saveas(gcf, [newdir, '/','Topol ', num2str(i), '.png']);
%     saveas(gcf, [newdir, '/','Topol ', num2str(i), '.fig']);
%     figure;
% end