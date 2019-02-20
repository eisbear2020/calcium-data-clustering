function per_trial_clus = FunctionalCellAnalysis(X,type,id_all_trials,id_avg_clus,id_parts,sel_avg_clust,cell_id_array)
    % returns:  - per_trial_clus: rows: cells, col: trials, values: cluster ID        

    if type == "simple_switching"
    % calculates to how many clusters a cell belongs and plots the 
    % histogram
    
        % how many cells in the data set
        nr_cells = size(X.dFF,1);
        % how many trials
        nr_trials = length(X.dFF_sep_trials);

        cell_clust = zeros(nr_cells, nr_trials);
        clust_per_cell = zeros(nr_cells,1);

        for cell_id = 1:nr_cells
            for trial=1:nr_trials
                cell_clust(cell_id,trial) = id_all_trials(cell_id + (trial-1)*nr_cells);
            end
            clust_per_cell(cell_id) = length(unique(cell_clust(cell_id,:)));
        end
        figure
        histogram(clust_per_cell)
        xlabel("#CLUSTERS/CELL")
        ylabel("#CELLS")
        set(gcf, 'Color', [0 0 0]);
        set(gcf, 'InvertHardCopy', 'off');

        figure
        for i=1:5
            subplot(2,3,i)
            scatter3(X.coordinates(:, 2),X.coordinates(:, 1), ...
            -X.coordinates(:, 3),...
            "MarkerEdgeColor",[0.5 0.5 0.5]); hold on;     
            scatter3(X.coordinates(clust_per_cell == i, 2),X.coordinates(clust_per_cell == i, 1), ...
            -X.coordinates(clust_per_cell == i, 3),"MarkerFaceColor",[1 0.2 0.2]); hold on;    
            title(strcat("cell in ",num2str(i)," different cluster"),"Color","white")
            zlabel("slices")
            xlabel("x")
            ylabel("y")
            set(gca,'Color','None')
            set(gca,'XColor',[1 1 1]);
            set(gca,'YColor',[1 1 1]);
            set(gca,'ZColor',[1 1 1]);
            fig = gcf;
            fig.Color = 'None';
            fig.InvertHardcopy = 'off';
            hold on
        end    
    per_trial_clus = cell_clust;
    
    
    elseif type == "switching_of_subset"
    % uses cells from selected average clusters (sel_avg_clust) and plots
    % histogram of cluster assignments for different trials
        
        per_trial_clus = FunctionalCellAnalysis(X,"simple_switching",id_all_trials,[],[]);
        
        % plots histogram for cluster
        nr_trials = length(X.dFF_sep_trials);

        % cells in the selected clusters
        cell_id_array = [];
        in_which_clust = [];

        for i = 1:length(sel_avg_clust)
            cell_id_array = vertcat(cell_id_array,find(id_avg_clus == sel_avg_clust(i)));
        end
        for cell_id = 1:length(cell_id_array)
            in_which_clust = vertcat(in_which_clust,per_trial_clus(cell_id,:));
        end
        figure;
        title("Switching of selected cells")
        for i = 1:nr_trials
            subplot(1,nr_trials,i)
            histogram(in_which_clust(:,i))
            title(strcat("trial",num2str(i)))
            xlim([0 10])
        end
     
    elseif type == "plot_switching_cells"
    % cluster assignment of selected cells for different trials

        per_trial_clus = FunctionalCellAnalysis(X,"simple_switching",id_all_trials,[],[]);
        for cell_id = 1:length(cell_id_array)
            in_which_clust = per_trial_clus(cell_id,:);
            % for single cell
            figure
            scatter([1,2,3,4,5],-(10-in_which_clust),120,"MarkerEdgeColor","w","MarkerFaceColor","w");
            xlabel("TRIAL NR.")
            xlim([1 5])
            ylim([-10 0])
            title("SWITCHING CELL")
            grid on;
            ylabel("CLUSTER ID")
            set(gcf, 'Color', [0 0 0]);
            set(gcf, 'InvertHardCopy', 'off');
        end
    
    elseif type == "mutual_info_ent" 
    % uses mutual information between cluster and sub-stimulus and 
    % entropy of clusters to characterize a single cell
    
        % split trial in parts of equal length
        parts = 3;
        % number of clusters
        nClusters = length(unique(id_parts));
        idx = id_parts;
        % split trials into parts and stack them up
        trial_dat = X.dFF_all_trials;
        part_length = floor(size(trial_dat,2)/parts);

        stacked_parts = [];
        for part = 1:parts
           stacked_parts = vertcat(stacked_parts, trial_dat(:,((part-1)*part_length+1):(part*part_length)));  
        end


        % how many cells in the data set
        nr_cells = size(X.dFF,1);
        % how many trials
        nr_trials = length(X.dFF_sep_trials);

        cell_clust = zeros(nr_cells, nr_trials);
        clust_per_cell = zeros(nr_cells,1);


        for cell_id = 1:nr_cells
            for part = 1:parts
                for trial=1:nr_trials
                    cell_clust(cell_id,(part-1)*nr_trials+trial) = idx(cell_id +(part-1)* (nr_trials *nr_cells)+(trial-1)*nr_cells);
                end
                clust_per_cell(cell_id) = length(unique(cell_clust(cell_id,:)));
            end
        end


        char_cells = zeros(nr_cells,nClusters*parts);

        % go through all cells
        for cell_id = 1:nr_cells
            for part = 1:parts
                for trial=1:nr_trials
                    cl_id = cell_clust(cell_id,(part-1)*nr_trials+trial);
                    char_cells(cell_id,(part-1)*nClusters+cl_id) = char_cells(cell_id,(part-1)*nClusters+cl_id)+1;
                end
            end
        end


        % reduce dimension to plot multiple cells:
        % mutual information and entropy

        sum_matrix= zeros(nr_cells,nClusters);
        mutual_info = zeros(nr_cells,1);

        for cell_id = 1:nr_cells
            sum_vec = zeros(1,nClusters);
            hist_2D = [];
            for part = 1:parts
               sum_vec = sum_vec + char_cells(cell_id,((part-1)*nClusters+1):(part*nClusters)); 
               % generate 2D histogram
               hist_2D = vertcat(hist_2D,char_cells(cell_id,((part-1)*nClusters+1):(part*nClusters))); 
            end
            mutual_info(cell_id) = mutualInfoHist(hist_2D);

            sum_matrix(cell_id,1:nClusters) = sum_vec;
        end

        % shannon

        shannon_ind = zeros(nr_cells,1);
        %non_zero_el = zeros(nr_cells,1);

        % no zero elements for log
        sum_matrix = sum_matrix +  0.00001;

        for cell_id = 1:nr_cells
            %non_zero_el(cell_id) = length(nonzeros(char_cells(cell_id,:)));
            shan_ind = 0;
            for i = 1:nClusters
                p_i = sum_matrix(cell_id,i)/sum(sum_matrix(cell_id,:));
                shan_ind = shan_ind + p_i*log(p_i);
            end
            shannon_ind(cell_id) = -shan_ind;
        end
        
        figure
        scatter(shannon_ind,mutual_info,40, "MarkerFaceColor","w", "MarkerEdgeColor","w")
        xlabel("SHANNON ENTROPY: CLUSTERS")
        ylabel("MUTUAL INFORMATION")
        title("CELL CHARACTERISTIC")

        labels = cellstr(strcat(repmat("C",nr_cells,1),num2str([1:nr_cells]'))); 

        dx = 0.001; dy = 0.03; % displacement so the text does not overlay the data points

        text(shannon_ind+dx, mutual_info+dy, labels);    
        set(gcf, 'Color', [0 0 0]);
        set(gcf, 'InvertHardCopy', 'off');
        
        % visualize matrix

        sel_cell = 260;
        char_matrix = [];

        for part = 1:parts
           char_matrix = vertcat(char_matrix,char_cells(sel_cell,((part-1)*nClusters+1):(part*nClusters))); 
        end

        char_matrix = char_matrix.'/5;
        figure
        imagesc(char_matrix)
        newmap = contrast(char_matrix);
        colormap("gray")
        x = [1.5 1.5];
        y = [0 20];
        line(x,y,"Color","r","LineWidth",2)
        x = [2.5 2.5];
        y = [0 20];
        line(x,y,"Color","r","LineWidth",2)
        xticks([1 2 3])
        %caxis([0 0.6])
        colorbar
        title(strcat("CELL ",int2str(sel_cell),": CONDITIONAL PROB"))
        xlabel("SUB-STIM ID")
        ylabel("CLUSTER ID")
        set(gcf, 'Color', [0 0 0]);
        set(gcf, 'InvertHardCopy', 'off');        
        
        
    else
        error("Type not known")
    
    end 
    
        
    
end

