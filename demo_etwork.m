% clearvars
% pv = zeros(15,2);
% for I=1:15
%     disp(I);
%     % ODS
%     ODS_X = load(['etwork-', num2str(I), '-adjacency.txt']); %Load I-th ODS
%     pv(I,1) = RunPoisson(ODS_X);
%     % CRDS
%     pvC  = 0;
%     parfor Run=1:50
%         CRDS_X = Generate_CRDS(ODS_X);
%         pv_t = RunPoisson(CRDS_X);
%         pvC = pvC + pv_t;
%     end
%     pv(I,2) = pvC/50;
%     %% Plot ODS
%     figure
%     ODS_X(ODS_X==1)=0;
%     ODS_X(ODS_X==2)=1;
%     plot(graph(ODS_X));
%     title(['etwork-', num2str(I)],['p-value: ',num2str(pv(I,1))]);
% end



clearvars
pv = zeros(10,2);

%for I=1:6
for I=1:10
    disp(I);
    tic;
    % ODS
    ODS_X = load(['configuration_network_', num2str(I), '.txt']); %Load I-th ODS
    pv(I,1) = RunPoisson(ODS_X);
    elapsedTime = toc;
    fprintf('运行时间:%.2f 秒\n',elapsedTime);
    % CRDS
    pvC  = 0;
%     parfor Run=1:50
%         CRDS_X = Generate_CRDS(ODS_X);
%         pv_t = RunPoisson(CRDS_X);
%         pvC = pvC + pv_t;
%     end
    pv(I,2) = pvC/50;
    
    % Print results
    disp(['graph_Gnp_', num2str(I), ' ODS p-value: ', num2str(pv(I,1)), ' CRDS p-value: ', num2str(pv(I,2))]);

    % Plot ODS
%     figure
%     ODS_X(ODS_X==1)=0;
%     ODS_X(ODS_X==2)=1;
%     G = graph(ODS_X);
%     p = plot(G);
%     p.NodeLabel = {};  % Remove node labels
%     p.NodeColor = 'k'; % Set node color to black (k stands for black in MATLAB)
% 
%     highlight(p, G.Edges.EndNodes(:,1), G.Edges.EndNodes(:,2), 'EdgeColor', [0.5 0.5 0.5]); % Change edge color to gray using RGB triplet
%     title(['Network-', num2str(I)],['p-value: ',num2str(pv(I,1))]);


end





% % karate:sX = 0:34
% clearvars
% pv = zeros(34, 2);  % Modified to accommodate 34 sX values
% 
% for I=21:21
%     disp(['Network: ', num2str(I)]);
%     
%     for sX = 0:34
%         % ODS
%         ODS_X = load(['etwork-', num2str(I), '-adjacency.txt']); %Load I-th ODS
%         pv(sX+1,1) = RunPoisson(ODS_X, sX);
%         
%         % CRDS
%         pvC  = 0;
%         parfor Run=1:50
%             CRDS_X = Generate_CRDS(ODS_X);
%             pv_t = RunPoisson(CRDS_X, sX);
%             pvC = pvC + pv_t;
%         end
%         pv(sX+1,2) = pvC/50;
%         
%         % Print results
%         disp(['sX: ', num2str(sX), ' ODS p-value: ', num2str(pv(sX+1,1)), ' CRDS p-value: ', num2str(pv(sX+1,2))]);
%     end
% end




