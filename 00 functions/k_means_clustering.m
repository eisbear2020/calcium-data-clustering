function cl_id = k_means_clustering(data,type,nClusters,stim,trial)
%performs k-means clustering
% -X: data
% -type: "temp", "temp_avg_trials" or "temp_cluster_localization"
%        "temp_single_trial", "temp_across_trials", "pop_states"
%        "temp_trial_parts"
%
% -nClusters: # of clusters
% returns: cl_id --> cluster id for each row entry of the input data


if type == "temp"
    
    X = data.dFF(:,1:data.i_dur_stim{1,length(data.i_dur_stim)}(end));
    
    % perform k-means on cells
    [idx, C] = kmeans(X, nClusters);

    %visualize individual clusters
    inds = zeros(length(idx), 1);
    idxPic = zeros(length(idx), 1);

    ccmap = jet(nClusters);
    cnt = 0; 
    
    for i = 1:nClusters
        cnt1 = cnt + sum(idx == i);
        inds(cnt+1:cnt1) = find(idx == i);
        idxPic(cnt+1:cnt1) = i;
        cnt = cnt1;
    end

    %visualize entire dataset sorted by cluster
    figure;
    subplot(2, 12, 1); 
    imagesc(idxPic); colormap(jet(nClusters));
    title('CLUSTER ID (COLOR)'); ylabel('CELL ID'); xlabel('');
    subplot(2, 12, 2:12);
    imagesc(X(inds,:)); colormap(jet(256));
    title('dFF');  %set(gca,'xtick',[])
    set(gca, 'YTickLabels', {''}, 'YTick', []);
    set(gcf, 'Position', [127, 1253, 1136, 306]);
    subplot(2, 12,14:24);
    for i = 1: size(data.i_dur_stim,2)
        plot(linspace(data.i_dur_stim{1,i}(1),data.i_dur_stim{1,i}(end),length(stim)),stim,'Color','w')
        hold on
    end
    ylabel({'INTENSITY','(NORM)'})
    xlabel('TIME')
    xlim([0.5, size(X, 2)+0.5])
    set(gcf, 'Color', [0 0 0]);
    set(gcf, 'InvertHardCopy', 'off');
    
    figure
    subplot(2,1,1)
    for i = 1:nClusters
        a=rand(1,3);
         plot(X(idx == i, :)'+5*i, 'Color', a); hold on;
         plot(C(i, :)+5*i, 'LineWidth', 2, 'Color', "k");
         set(gcf, 'Position', [-333, 1514, 1919, 281]);
         xlim([1, size(X, 2)])
         ylim([0,max(C(i, :)+5*i)+5])
         title("Cell activation")
         hold on

        cnt1 = cnt + sum(idx == i);
        inds(cnt+1:cnt1) = find(idx == i);
        idxPic(cnt+1:cnt1) = i;
        cnt = cnt1;
    end
    set(gca,'xtick',[])
    subplot(2,1,2)
    for i = 1: size(data.i_dur_stim,2)
        plot(linspace(data.i_dur_stim{1,i}(1),data.i_dur_stim{1,i}(end),length(stim)),stim,'Color','w')
        hold on
    end
    ylabel({'INTENSITY','(NORM)'})
    xlabel('TIME')
    xlim([0.5, size(X, 2)+0.5])
    set(gcf, 'Color', [0 0 0]);
    set(gcf, 'InvertHardCopy', 'off');
    
    cl_id = idx;
end
    
    
if type == "temp_avg_trials"
    X = data.dFF_avg_over_trials;
    
    % perform k-means on cells
    [idx, C] = kmeans(X, nClusters);

    %visualize individual clusters
    inds = zeros(length(idx), 1);
    idxPic = zeros(length(idx), 1);

    ccmap = jet(nClusters);
    cnt = 0; 
    x_axis_data = 0:(size(X, 2)-1);
    offset = 2.5;
    figure
    subplot(2,1,1)
    for i = 1:nClusters
         plot(x_axis_data,X(idx == i, :)'+offset*i,'Color',ccmap(i,:)); hold on;
         plot(x_axis_data,C(i, :)+offset*i, 'LineWidth', 2.5, 'Color', [1 1 1]);
         set(gcf, 'Position', [-333, 1514, 1919, 281]);
         title("dFF")
         hold on
    end
    xlim([0, size(X, 2)-1])
    ylim([0,max(C(i, :)+offset*i)])
    set(gca,'xtick',[])
    %set(gca,'ytick',[])
    ylabel("CLUSTER ID")
    yticks([offset:offset:offset*nClusters])
    % DONT HARDCODE --> find better solution
    yticklabels({"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"})
    subplot(2,1,2)
    plot(linspace(0,size(X, 2)-1,length(stim)),stim,"Color","w")
    xlim([0,size(X, 2)-1])
    xlabel("TIME BINS")
    ylabel({"INTENSITY","(NORM)"}) 
    
    cl_id = idx;    
end   
    
if type == "temp_cluster_localization"
  
    X = data.dFF(:,1:data.i_dur_stim{1,length(data.i_dur_stim)}(end));
    
    % perform k-means on cells
    [idx, C] = kmeans(X, nClusters);

    %visualize individual clusters
    inds = zeros(length(idx), 1);
    idxPic = zeros(length(idx), 1);
    label_pos = zeros(nClusters, 1);

    ccmap = jet(nClusters);
    cnt = 0; 
    
    for i = 1:nClusters
        cnt1 = cnt + sum(idx == i);
        inds(cnt+1:cnt1) = find(idx == i);
        idxPic(cnt+1:cnt1) = i;
        label_pos(i) = cnt+1 +(cnt1-cnt+1)/2;
        cnt = cnt1;
    end
    labels = cellstr(strcat(repmat("CL",nClusters,1),num2str([1:nClusters]'))); 
    %visualize entire dataset sorted by cluster
    figure;
    subplot(2, 12, 1); 
    imagesc(idxPic); colormap(jet(nClusters));
    title('CLUSTER'); ylabel('CELL ID'); xlabel('');
    set(gca, 'XTickLabels', {''}, 'XTick', []);
    text(repmat(0.9,nClusters,1), label_pos, labels, "Color","k",'fontweight','bold');    
    subplot(2, 12, 2:12);
    imagesc(X(inds,:)); colormap(jet(256));
    title('dFF');  set(gca,'xtick',[])
    set(gca, 'YTickLabels', {''}, 'YTick', []);
    set(gcf, 'Position', [127, 1253, 1136, 306]);
    subplot(2, 12,14:24);
    for i = 1: size(data.i_dur_stim,2)
        plot(linspace(data.i_dur_stim{1,i}(1),data.i_dur_stim{1,i}(end),length(stim)),stim,'Color','w')
        hold on
    end
    ylabel({'INTENSITY','(NORM)'})
    xlabel('TIME')
    xlim([0.5, size(X, 2)+0.5])
    set(gcf, 'Color', [0 0 0]);
    set(gcf, 'InvertHardCopy', 'off');
    ccmap = jet(nClusters);
    
    figure
    for i=1:nClusters
        subplot(2,nClusters/2,i)
        scatter3(data.coordinates(:, 2),data.coordinates(:, 1), ...
        -data.coordinates(:, 3),...
        "MarkerEdgeColor",[0.85 0.85 0.85],"MarkerFaceColor","none"); hold on;     
        scatter3(data.coordinates(idx == i, 2),data.coordinates(idx == i, 1), ...
        -data.coordinates(idx == i, 3),"MarkerFaceColor",ccmap(i,:), "MarkerEdgeColor",ccmap(i,:)); hold on;    
        title(strcat("CLUSTER ",num2str(i)))
        zlabel("SLICES")
        xlabel("X")
        ylabel("Y")
    end    
    set(gcf, 'Color', [0 0 0]);
    set(gcf, 'InvertHardCopy', 'off');
    
    % cluster "density" per slice (how many cells per slice for each cluster)
    cluster_depth = zeros(nClusters, 7);
    coord = data.coordinates;

    figure
    for i=1:nClusters
        for slice = 1:7
            cluster_depth(i,slice) = length(coord(coord(idx == i,3) ==slice));
        end
    end
    cluster_depth = bsxfun(@rdivide, cluster_depth', sum(cluster_depth'))';
    imagesc(cluster_depth.')
    newmap = contrast(cluster_depth.');
    colormap(newmap)
    colorbar
    title("FRACTION OF CELLS")
    xlabel("CLUSTER ID")
    ylabel("SLICE"),
    set(gcf, 'Color', [0 0 0]);
    set(gcf, 'InvertHardCopy', 'off');
    cl_id = idx;    
end


if type == "temp_single_trial"
    trial_nr = trial;
    cnt = 0;
    nClusters = 8;
    % go through all trials

    X = data.dFF(:,data.i_dur_stim{1,trial_nr});

    [idx, C] = kmeans(X, nClusters);
    %visualize individual clusters
    inds = zeros(length(idx), 1);
    idxPic = zeros(length(idx), 1);

    ccmap = jet(nClusters);
    cnt = 0;

    for i = 1:nClusters
        cnt1 = cnt + sum(idx == i);
        inds(cnt+1:cnt1) = find(idx == i);
        idxPic(cnt+1:cnt1) = i;
        cnt = cnt1;
    end


    %visualize entire dataset sorted by cluster
    figure;
    subplot(2, 12, 1); 
    imagesc(idxPic); colormap(jet(nClusters));
    title('CLUSTER ID (COLOR)'); ylabel('CELL ID'); xlabel('');
    subplot(2, 12, 2:12);
    imagesc(X(inds,:)); colormap(jet(256));
    title('CELL ACTIVATIONS');  set(gca,'xtick',[])
    set(gca, 'YTickLabels', {''}, 'YTick', []);
    set(gcf, 'Position', [127, 1253, 1136, 306]);
    %xlim([0,136])
    subplot(2, 12,14:24);
    plot(linspace(1,size(X, 2),length(stim)),stim,"Color","b")
    xlim([0.5,size(X, 2)+0.5])
    xlabel("TIME")
    ylabel("intensity, scaled")   
    
    cl_id = idx;
end


if type == "temp_across_trials"
     
    X = data.dFF_all_trials;

    [idx, C] = kmeans(X, nClusters);
    %visualize individual clusters
    inds = zeros(length(idx), 1);
    idxPic = zeros(length(idx), 1);

    ccmap = jet(nClusters);
    cnt = 0;
    label_pos = zeros(nClusters, 1);
    
    for i = 1:nClusters
        cnt1 = cnt + sum(idx == i);
        inds(cnt+1:cnt1) = find(idx == i);
        idxPic(cnt+1:cnt1) = i;
        label_pos(i) = cnt+1 +(cnt1-cnt+1)/2;
        cnt = cnt1;
    end
    labels = cellstr(strcat(repmat("CL",nClusters,1),num2str([1:nClusters]'))); 

    %visualize entire dataset sorted by cluster
    figure;
    subplot(2, 12, 1); 
    imagesc(idxPic); colormap(jet(nClusters));
    title('CLUSTER'); ylabel('TRACE ID'); xlabel('');
    set(gca, 'XTickLabels', {''}, 'XTick', []);
    text(repmat(0.9,nClusters,1), label_pos, labels, "Color","k",'fontweight','bold');    
    set(gca,'Color','k')
    set(gca,'xtick',[])
    subplot(2, 12, 2:12);
    imagesc(X(inds,:)); colormap(jet(256));
    title('dFF');  set(gca,'xtick',[])
    set(gca, 'YTickLabels', {''}, 'YTick', [],'Color','k');
    set(gcf, 'Position', [127, 1253, 1136, 306]);

    subplot(2, 12, 2:12);
    imagesc(X(inds,:)); colormap(jet(256));
    title('dFF');  set(gca,'xtick',[])
    set(gca, 'YTickLabels', {''}, 'YTick', []);
    set(gcf, 'Position', [127, 1253, 1136, 306]);
    subplot(2, 12,14:24);
    plot(linspace(0,size(X, 2)-1,length(stim)),stim,"Color","w")
    xlim([0,size(X, 2)-1])
    xlabel("TIME BINS")
    ylabel({"INTENSITY","(NORM)"}) 
    set(gcf, 'Color', [0 0 0]);
    set(gcf, 'InvertHardCopy', 'off');
    cl_id = idx;
end

if type == "pop_states"
    %preprocess data
    X = data.dFF_raw;

    %remove "abnormal" cells based on arbitrary definition
    indsToReject = mean(X, 2);
    indsToReject = abs(indsToReject) > 2;
    X(indsToReject, :)  = [];

    %rotate matrix
    X = X';

    %z-score
    X = bsxfun(@minus, X', mean(X'))';
    X = bsxfun(@rdivide, X', sqrt(var(X')))';
    
    [idx, C] = kmeans(X, nClusters);

    %visualize individual clusters
    inds = zeros(length(idx), 1);
    idxPic = zeros(length(idx), 1);

    ccmap = jet(nClusters);
    cnt = 0;
    for i = 1:nClusters
        cnt1 = cnt + sum(idx == i);
        inds(cnt+1:cnt1) = find(idx == i);
        idxPic(cnt+1:cnt1) = i;
        cnt = cnt1;
    end
    %visualize entire dataset sorted by cluster
    figure;
    subplot(1, 12, 1); 
    imagesc(idxPic); colormap(jet(nClusters)); 
    title('CLUSTER ID (COLOR)'); ylabel('Time point'); xlabel('');
    subplot(1, 12, 2:12);
    imagesc(X(inds, :)); colormap(jet(256)); 
    title('CELL ACTIVATIONS');  xlabel('CELL ID');
    set(gca, 'YTickLabels', {''}, 'YTick', []);
    set(gcf, 'Position', [127, 1253, 1136, 306]);
    
    
    figure;
    Xr = data.dFF;
    Xr(indsToReject, :) = [];
    subplot(12, 1, 1); imagesc(idx'); colormap(jet(nClusters));
    set(gca, 'YTickLabels', {''}, 'YTick', [], 'XTickLabels', {''});
    subplot(12, 1, 2:12); imagesc(Xr); colormap(jet(256));
    xlabel('TIME'); ylabel("CELLS");
  
    cl_id = idx;
    
end


    if type == "temp_trial_parts"
        
        X = data.dFF_stacked_trial_parts;

        [idx, C] = kmeans(X, nClusters);
        %visualize individual clusters
        inds = zeros(length(idx), 1);
        idxPic = zeros(length(idx), 1);

        ccmap = jet(nClusters);
        cnt = 0;

        for i = 1:nClusters
            cnt1 = cnt + sum(idx == i);
            inds(cnt+1:cnt1) = find(idx == i);
            idxPic(cnt+1:cnt1) = i;
            cnt = cnt1;
        end


        %visualize entire dataset sorted by cluster
        figure;
        subplot(1, 12, 1); 
        imagesc(idxPic); colormap(jet(nClusters));
        title('CLUSTER ID (COLOR)'); ylabel('CELL ID'); xlabel('');
        set(gca, 'XTickLabels', {''}, 'XTick', []);
        subplot(1, 12, 2:12);
        imagesc(X(inds,:)); colormap(jet(256));
        set(gca, 'YTickLabels', {''}, 'YTick', []);
        title('CELL ACTIVATIONS');
        xlabel("TIME BINS")
        set(gcf, 'Position', [127, 1253, 1136, 306]);
        set(gcf, 'Color', [0 0 0]);
        set(gcf, 'InvertHardCopy', 'off');
        cl_id = idx;        
        
    end
    
if type == "SNR_cluster"

    X = data.SNRTimePointsCells_mat;

    [idx, C] = kmeans(X, nClusters);
    %visualize individual clusters
    inds = zeros(length(idx), 1);
    idxPic = zeros(length(idx), 1);

    ccmap = jet(nClusters);
    cnt = 0;

    for i = 1:nClusters
        cnt1 = cnt + sum(idx == i);
        inds(cnt+1:cnt1) = find(idx == i);
        idxPic(cnt+1:cnt1) = i;
        cnt = cnt1;
    end


    %visualize entire dataset sorted by cluster
    figure;
    subplot(2, 12, 1); 
    imagesc(idxPic); colormap(jet(nClusters));
    title('CLUSTER ID (COLOR)'); ylabel('CELL ID'); xlabel('');
    subplot(2, 12, 2:12);
    imagesc(X(inds,:)); colormap(jet(256));
    title('CELL ACTIVATIONS');  set(gca,'xtick',[])
    set(gca, 'YTickLabels', {''}, 'YTick', []);
    set(gcf, 'Position', [127, 1253, 1136, 306]);
    %xlim([0,136])
    subplot(2, 12,14:24);
    plot(linspace(1,size(X, 2),length(stim)),stim,"Color","b")
    xlim([0.5,size(X, 2)+0.5])
    xlabel("TIME")
    ylabel("intensity, scaled")   
    
    cl_id = idx;
end
  

    
end

