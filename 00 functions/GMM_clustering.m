function cl_id = GMM_clustering(X,type,nClusters,prop_ts,fitting_rep,clust_nr_array,stim)
% evaluates GMM fitting and generates GMM clustering
% returns cluster IDs for each row of the input X

    if type == "cross_val"
    % cross-validation of #clusters using log-likelihood
 

        % split into training and testing set for cross-validation
        cv = cvpartition(size(X,1),'HoldOut',prop_ts);
        idx = cv.test;

        X_train = X(~idx,:);
        X_test = X(idx,:);

        % options for GMM fitting
        options = statset('Display','final','MaxIter',1000);
        %gm = fitgmdist(X,4,'Options',options,'SharedCovariance',true);

        % loglikelihood array:
        % nr. of clusters, number of results for average, average
        nllh_array = zeros(size(clust_nr_array,2),3);

        % test different nr of clusters
        for i = 1:(size(clust_nr_array,2)) 
            nr_clust = clust_nr_array(i);

            % repeat fitting to minimize influence of initialization
            % use average likelihood over repeats
            lh_array_rep = [];

            for rep = 1:fitting_rep
                try
                    % fit gaussian mixture model with diagonal covariance matrix
                    gm = fitgmdist(X_train,nr_clust,'Options',options,"CovarianceType",'diagonal');
                    %gm = fitgmdist(X_train,nr_clust,'Options',options);
                    [P,nlogL] = posterior(gm,X_test); 
                    lh_array_rep = [lh_array_rep nlogL];
                catch
                    % if an error occurs while fitting the model go to next
                    % iteration
                    continue    
                end

            end
            % calculate average
            nllh_array(i,1:(length(lh_array_rep))) = lh_array_rep;  
        end

        % calculate mean & standard error for each cluster number
        mean_likelihood = zeros(size(nllh_array,1),1);
        std_err_mean = zeros(size(nllh_array,1),1);
        nr_succ_fits = zeros(size(nllh_array,1),1);
        for i = 1:size(nllh_array,1)
            mean_likelihood(i) = mean(nonzeros(nllh_array(i,:)));
            std_err_mean(i) = std(nonzeros(nllh_array(i,:)))/sqrt(size(nonzeros(nllh_array(i,:)),1));
            nr_succ_fits(i) = length(nonzeros(nllh_array(i,:)))/fitting_rep;
        end

        figure
        subplot(1,2,1)
        scatter(clust_nr_array,mean_likelihood,'MarkerEdgeColor',"white",'MarkerFaceColor',"white")
        title('CROSS VALIDATION: #CLUSTERS')
        xlabel('#CLUSTERS')
        ylabel('-(LOG-LIKELIHOOD)')
        hold on
        err = std_err_mean;
        errorbar(clust_nr_array, mean_likelihood, err, 'LineStyle','none');
        subplot(1,2,2)
        scatter(clust_nr_array,nr_succ_fits*100,'MarkerEdgeColor',"white",'MarkerFaceColor',"white")
        title('CROSS VALIDATION: SUCCESSFUL FITS')
        xlabel('#CLUSTERS')
        ylabel('SUCCESSFUL FITS / %')
        ylim([0 100])
        
        cl_id = [];
    end



    if type == "fitting"
    % fit GMM to X
        % options for GMM fitting
        options = statset('Display','final','MaxIter',1000);
        gm = fitgmdist(X,nClusters,'Options',options,"CovarianceType",'diagonal');

        idx = cluster(gm,X);

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

        %visualize entire Xset sorted by cluster
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
    
end

