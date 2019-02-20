function FunctionalClusterAnalysis(X,type,idx,per_trial_clus)
    % TODO: use shannon entropy for cells in "char_clusters"
    
    if type == "char_clusters"
    % characterizes clusters using Shannon entropy for trials and diversity
    % of cells within the cluster
        
        % how many cells in the data set
        nr_cells = size(X.dFF,1);
        % how many trials
        nr_trials = length(X.dFF_sep_trials);
        % number of clusters
        nClusters = length(unique(idx));

        cell_clust = zeros(nr_cells, nr_trials);

        trial_dist = zeros(nClusters,nr_trials);
        cells_in_clust = zeros(nClusters,1);
        % go through clusters

        for cl_id = 1:nClusters
        % go through cells
            for i = 1:nr_cells
                trial_id = [];
                trial_id =  find(per_trial_clus(i,:) == cl_id);
                if ~isempty(trial_id)
                   cells_in_clust(cl_id) = cells_in_clust(cl_id) +1; 
                end

                trial_dist(cl_id,trial_id) = trial_dist(cl_id,trial_id) +1;
            end

        end

        cells_in_clust = cells_in_clust./sum(trial_dist,2);

        
        % shannon index for #cells/cluster: NOT IMPLEMENTED
        %-----------------------------------------------------------------
        cell_dist_clust = zeros(nr_cells,nClusters);
        for cl_id = 1:nClusters
            for cell_id = 1:nr_cells
                occ_cell =  length(find(per_trial_clus(cell_id,:) == cl_id));
                if ~isempty(occ_cell)
                   cell_dist_clust(cell_id,cl_id) = occ_cell; 
                end


            end
        end
        
        % shannon index cells
        cl_shan_ind_cell = zeros(nClusters,1);
        cell_dist_clust = cell_dist_clust +  0.00001;
        
        for cl_id = 1:nClusters
            shan_ind = 0;
            % go through all cells
            for i = 1:nr_cells
                p_i = cell_dist_clust(i,cl_id)/sum(cell_dist_clust(:,cl_id));
                shan_ind = shan_ind + p_i*log(p_i);
            end

        cl_shan_ind_cell(cl_id) = - shan_ind;
        end

        
        
        % shannon index for clusters
        %------------------------------------------------------------------
        cl_shan_ind = zeros(nClusters,1);

        % shannon for clusters
        for cl_id = 1:nClusters
            shan_ind = 0;
            % go through all trials
            for i = 1:5
                p_i = trial_dist(cl_id,i)/sum(trial_dist(cl_id,:)) + 0.00001;
                shan_ind = shan_ind + p_i*log(p_i);
            end

        cl_shan_ind(cl_id) = - shan_ind;
        end
        
        figure
        labels = cellstr(strcat(repmat("CL",nClusters,1),num2str([1:nClusters]'))); 
        d = colormap(jet(nClusters));

        dx = 0.001; dy = 0.03; % displacement so the text does not overlay the data points
        for i = 1:size(cells_in_clust,1)
            scatter(cells_in_clust(i),cl_shan_ind(i),120,"MarkerFaceColor",d(i,:),"MarkerEdgeColor",d(i,:))
            ylabel("SHANNON ENTROPY FOR TRIALS")
            xlabel("#DIFFERENT CELLS PER CLUSTER/#TRACES PER CLUSTER")
            hold on
        end
        text(cells_in_clust+dx, cl_shan_ind+dy, labels);    
         set(gcf, 'Color', [0 0 0]);
        set(gcf, 'InvertHardCopy', 'off');
        per_trial_clus = [];
        
    elseif type == "transitions_clusters"
    % calculates and plots transition probabilities between clusters
        
        % number of clusters
        nClusters = length(unique(idx));        
        
        trace_per_clust = zeros(nClusters,1);
        % count how many traces per cluster
        for cl_id = 1:nClusters
            trace_per_clust(cl_id) = sum(idx(:,1) == cl_id);        
        end

        % init transition matrix
        trans_matrix = zeros(nClusters,nClusters);
        total_transitions = size(per_trial_clus,1) * (size(per_trial_clus,2)-1);

        % go through all cells
        for cell_id = 1:size(per_trial_clus,1)
            % go check all transitions
           for trans = 1:(size(per_trial_clus,2) - 1)
               before = per_trial_clus(cell_id,trans);
               after = per_trial_clus(cell_id,trans+1);
               trans_matrix(before,after) = trans_matrix(before,after) +1;
           end
        end

        %turn into prob
        %trans_matrix = trans_matrix / total_transitions;
        
        % scale with size of cluster
        
        
        trans_matrix = trans_matrix ./sum(trans_matrix,2);
        
        %trans_matrix = trans_matrix ./trace_per_clust;
        %sum(trans_matrix,2)
        
        figure
        imagesc(trans_matrix)
        newmap = contrast(trans_matrix);
        colormap(newmap)
        colorbar
        xlabel("TO / CLUSTER ID")
        ylabel("FROM / CLUSTER ID")
        title("#TRANSISTION PROBABILITY")    
        set(gcf, 'Color', [0 0 0]);
        set(gcf, 'InvertHardCopy', 'off');
        
    else
        error("Type not known")
    
    end 
    
        
    
end

