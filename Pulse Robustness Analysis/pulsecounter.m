function [pc_min, pc_max] = pulsecounter(vec, thresh)
%% Detect number of pulses

%Pulse counter from zero:
pc_min = 0;
pc_max = 0;
% vec = mat(i,:);
initv = vec(1);

ender = length(vec);
finv2 = vec(ender/2);
finv = vec(ender);
if initv == 0
    if finv < thresh && finv2 < thresh && max(vec) > 0.5
        pc_min = true;
    end
elseif initv == 1
    if 1-finv < thresh && min(vec) < 0.5
        pc_max = true;
    end
end
    

% 
% for i = 1:size(mat,1)
%     vec = mat(i,:);
%     initv = vec(1);
%     finv = vec(end);
%     %initv = 0: start off
%     %initv = 1: start on
%     if initv == 0
%         if finv < thresh && max(vec) > 0.5
%             pc_min = pc_min + 1;
%         end
%     elseif initv == 1
%         if 1-finv < thresh && min(vec) < 0.5
%             pc_max = pc_max + 1;
%         end
%     end
% end