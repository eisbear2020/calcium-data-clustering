classdef CaData < handle
    %CADATA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        coordinates = []
        dFF_raw = []
        dFF = []
        dFF_sep_trials = []
        dFF_avg_over_trials = []
        dFF_all_trials = []
        dFF_stacked_trial_parts = []
        CaFrameTimes = []
        FileStartTime = []
        i_dur_stim = [] % indices that coincide with the stimulus
        dur = [] % length of sequence coinciding with the stimulus
        SNRTrials_arr = []
        SNRTimePointsCells_mat = []
    end
    
    methods
        function obj = CaData(SI_Aux_file,Ca_file)
        %------------------------------------------------------------------
        % constructs data object using data files as input
        %------------------------------------------------------------------
            % data files as input
            Ca_file_content = load(Ca_file);
            obj.coordinates = Ca_file_content.coordinates;
            obj.dFF = Ca_file_content.dFF;
            obj.dFF_raw = Ca_file_content.dFF;
            SI_Aux_file_content = load(SI_Aux_file);
            % adjust format of CaFrameTimes
            CaFrameTimes_in = SI_Aux_file_content.CaFrameTimes;
            % modify so that CaFramTimes can be divided by 8
            if size(CaFrameTimes_in,1) ~= 8
                CaFrameTimes_in=CaFrameTimes_in(1:end-mod(size(CaFrameTimes_in,2),8)); 
                CaFrameTimes_in=reshape(CaFrameTimes_in,8,[]);
            end
            
            obj.CaFrameTimes = CaFrameTimes_in;
            obj.FileStartTime = SI_Aux_file_content.FileStartTime;
        end
       
        function cutOffFilter(obj,threshold)
%             % filter all cells that have values larger than the threshold
%             obj.coordinates(max(abs(obj.dFF),[],2) > threshold,: ) =[];
%             obj.dFF(max(abs(obj.dFF),[],2) > threshold,: ) =[];
            %remove 'abnormal' cells based on arbitrary definition
            indsToReject = mean(obj.dFF, 2);
            indsToReject = abs(indsToReject) > threshold;
            obj.dFF(indsToReject, :)  = [];
            obj.coordinates(indsToReject, :)  = [];
            
        end
        
        function zScoreRows(obj)
        %------------------------------------------------------------------
        % applies row-wise z-score
        %------------------------------------------------------------------ 
            temp = obj.dFF;
            temp = bsxfun(@minus, temp', mean(temp'))';
            temp = bsxfun(@rdivide, temp', sqrt(var(temp')))';
            obj.dFF = temp;
        end       
        
        
        function SyncStimCaData(obj,stim_start_times, stim_end_times)
        %------------------------------------------------------------------
        % Extract indices of CaData that coincide with the stimulus
        % get times of for each imaging step (use first image of z-stack)
        % save data of single trials in dFF_sep_trials
        %------------------------------------------------------------------
            imaging_seq = obj.FileStartTime+seconds(obj.CaFrameTimes(1,1:end));
 
            for i = 1:size(stim_start_times,2)
                i_dur_stim_in{i} = find((imaging_seq >= stim_start_times(i))&(imaging_seq <= stim_end_times(i)));
            end

            % trim durations to same length
            n=min( cellfun(@(c) size(c,2), i_dur_stim_in));
            for i = 1:size(i_dur_stim_in,2)
                i_dur_stim_in{i}((n+1):end) = [];
            end
            obj.dur = n;
            obj.i_dur_stim = i_dur_stim_in;
            
            % shorten dFF to length of stimulus
            % obj.dFF = obj.dFF(:,obj.i_dur_stim{1,1}(1):obj.i_dur_stim{1,length(obj.i_dur_stim)}(end));
            %obj.dFF = obj.dFF(:,1:obj.i_dur_stim{1,length(obj.i_dur_stim)}(end));
        end
        
        function averageOverTrials(obj) 
        %------------------------------------------------------------------
        % averages the CA signal over all trials
        %------------------------------------------------------------------
            cell_comp_matrix = zeros(size(obj.dFF,1),obj.dur);
            % go through all the cells
            for cell = 1:size(obj.dFF,1)
                comb_data = zeros(size(obj.i_dur_stim,2),obj.dur);

                for i = 1:size(obj.i_dur_stim,2)
                   comb_data(i,:) = obj.dFF(cell,obj.i_dur_stim{i});
                end
                % calculate mean over all trials
                cell_comp_matrix(cell,:)= mean(comb_data);
            end          
            obj.dFF_avg_over_trials = cell_comp_matrix;
        end
        
        function SplitIntoTrials(obj)
         % save data of single trials 
            for i = 1:length(obj.i_dur_stim)
                obj.dFF_sep_trials{i} = obj.dFF(:,obj.i_dur_stim{i});
            end
        end
        
        function SplitTrialsIntoParts(obj, parts)
        % split trials into parts of equal length
            trial_dat = obj.dFF_all_trials;
            part_length = floor(size(trial_dat,2)/parts);
            stacked_parts = [];
            for part = 1:parts
               stacked_parts = vertcat(stacked_parts, trial_dat(:,((part-1)*part_length+1):(part*part_length)));  
            end 
            obj.dFF_stacked_trial_parts = stacked_parts;
        end
        
        
        
        function StackUpTrials(obj)
        % vertcat data of different trials: trials from one cell are 
        % at cell_id + (trial-1)*nr_cells
            all_trial_data = [];
            % go through all trials
            for i = 1: size(obj.i_dur_stim,2)
                all_trial_data = vertcat(all_trial_data,obj.dFF(:,obj.i_dur_stim{1,i}));
            end

            obj.dFF_all_trials = all_trial_data;
        end      
        
        
        function SNRTrials(obj)
        %------------------------------------------------------------------
        % calculates SNR for each cells and each trial
        % output: SNRTrials_arr:
        %           - col: trials
        %           - row: cells
        %------------------------------------------------------------------
             % initialize array
            tmp_arr = zeros(size(obj.dFF,1),size(obj.dFF_sep_trials,2));
            tmp_arr_matlab = zeros(size(obj.dFF,1),size(obj.dFF_sep_trials,2));

            % go through all cells
            for cell_nr = 1:size(obj.dFF,1)
                % go through all trials
                for trial_nr = 1: size(obj.dFF_sep_trials,2)
                    trial = obj.dFF_sep_trials{trial_nr}(cell_nr,:);
                    avg = obj.dFF_avg_over_trials(cell_nr,:);
                    [S,] = sumsqr(trial);
                    [N,] = sumsqr((trial-avg));
                    tmp_arr(cell_nr,trial_nr) = 20*log10(S/N);
                end
            end   
            obj.SNRTrials_arr = tmp_arr;
        end
            
        function plot_SNR_cell_summary(obj,stim_)
            % get mean of SNRs for cells over trials
            mean_SNRs = mean(obj.SNRTrials_arr,2);
            [sorted_SNRS, indices] = sort(mean_SNRs,'descend');
            figure
            histogram(mean_SNRs, 40,"FaceColor","w");
            title('AVERAGE SNR OVER TRIALS PER CELL')
            xlabel('SNR (dB)')
            ylabel('CELLS')
            set(gcf, 'Color', [0 0 0]);
            set(gcf, 'InvertHardCopy', 'off');
            figure
            subplot(2,1,1)
            for i = 1:10
                plot(obj.dFF_avg_over_trials(indices(i),:)+3*i)
                hold on
            end
            title('Avg over trials for 10 cells with best SNR')
            set(gca,'xtick',[])
            subplot(2,1,2)
            plot(linspace(1,136,length(stim_)),stim_)
            xlim([0,140])
            xlabel('TIME')
            ylabel('intensity, scaled')    
            set(gcf, 'Color', [0 0 0]);
            set(gcf, 'InvertHardCopy', 'off');
        end  
            
        
        function SNRTimePointsCells(obj)
        %------------------------------------------------------------------
        % calculates SNR for each time point of a cell using the mean SNR
        % from all trials 
        %
        % output: SNRTimePointsCells_mat:
        %           - col: SNR of time points
        %           - row: cells
        %------------------------------------------------------------------
            % initialize array
            tmp_arr = zeros(size(obj.dFF,1),obj.dur);

            % go through all cells
            for cell_nr = 1:size(obj.dFF,1)
                % go through all trials
                for trial_nr = 1: size(obj.dFF_sep_trials,2)
                    trial = obj.dFF_sep_trials{trial_nr}(cell_nr,:);
                    avg = obj.dFF_avg_over_trials(cell_nr,:);
                    tmp_trial = zeros(size(obj.dFF_sep_trials,2),obj.dur);
                    % go through all time points
                    for t_bin = 1:obj.dur
                        tmp_trial(trial_nr,t_bin) = 20*log10(trial(t_bin)^2/(trial(t_bin)-avg(t_bin))^2);
                    end
                end
                tmp_arr(cell_nr,:) = mean(tmp_trial,1);
            end   
            obj.SNRTimePointsCells_mat = tmp_arr;
        end        
        
        function plot_SNR_per_time_point(obj,stim_)
        %plots SNR for each cell and timepoint
            figure
            subplot(2,1,1)
            X =obj.SNRTimePointsCells_mat;
            imagesc(X); colormap(jet(256));
            c = colorbar;
            c.Label.String = 'SNR (dB)';
            title('SNR FOR EACH TIME BIN PER CELL')
            ylabel('CELLS')
            set(gca,'xtick',[])
            subplot(2,1,2)
            plot(linspace(0,size(X,2),length(stim_)),stim_,"Color","w")
            xlim([0.5,size(X, 2)+0.5])
            xlabel('TIME BINS')
            ylabel('intensity, scaled')
            set(gcf, 'Color', [0 0 0]);
            set(gcf, 'InvertHardCopy', 'off');
        end
       
        function plot_mean_all_cells_and_stim(obj,stim)
            % plots average over all cells and stimulus
            figure
            x_axis = [1:obj.dur];
            mean_all_cells = mean(obj.dFF_avg_over_trials);
            subplot(2,1,1)
            plot(x_axis,rescale(mean_all_cells),'Color','k','LineWidth',1.5)
            xlabel('time / s')
            ylabel('calcium signal, scaled: [0,1]')
            title(sprintf('Average, calcium traces, cells = %d',size(obj.dFF_avg_over_trials,1)))

            x_axis_stim = linspace(1,obj.dur,length(stim));
            subplot(2,1,2)
            plot(x_axis_stim, stim)
            ylabel('intensity, scaled: [0,1]')
            xlabel('time / s')
            title('Stimulus')
        end
       
        function plot_all_data_and_stim(obj, stim)
            X = obj.dFF;
            figure, subplot(2,1,1), imagesc(X); colormap(jet(256));
            title('dFF')
            ylabel('CELLS')
            set(gca,'xtick',[])
            subplot(2,1,2)
            for i = 1: size(obj.i_dur_stim,2)
                plot(linspace(obj.i_dur_stim{1,i}(1),obj.i_dur_stim{1,i}(end),length(stim)),stim,'Color','w')
                hold on
            end
            ylabel({'INTENSITY','(NORM)'})
            xlabel('TIME')
            xlim([0, size(X, 2)])
            set(gcf, 'Color', [0 0 0]);
            set(gcf, 'InvertHardCopy', 'off');
        end
        
    end
end

