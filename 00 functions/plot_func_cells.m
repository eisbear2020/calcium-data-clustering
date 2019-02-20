function plot_func_cells(CaData,stim)
%PLOT_FUNC_CELLS Summary of this function goes here
%   Detailed explanation goes here

    % Light sensing cells
    %----------------------------------------------------------------------
    x_axis = linspace(0,30,size(CaData,2));
    x_axis_stim = linspace(0,30,length(stim));
    cell_comp_matrix = CaData;
    figure

    subplot(4,1,1)
    light_sensing_cells = [];
    for cell = 1:size(cell_comp_matrix,1)
        if (mean(cell_comp_matrix(cell,10:23))/mean(cell_comp_matrix(cell,24:37))) > 4 && mean(cell_comp_matrix(cell,10:23)) > 0 
            light_sensing_cells = [light_sensing_cells, cell];
            plot(x_axis,rescale(cell_comp_matrix(cell,:)),'Color',[0.5 0.5 0.5])
            hold on
        end

    end
    plot(x_axis,rescale(mean(cell_comp_matrix(light_sensing_cells,:))),'Color','w','LineWidth',2.5)
    %plot(x_axis,cell_comp_matrix(light_sensing_cells(7),:),'Color','k','LineWidth',1.5)
    ylabel('dFF (NORM)')
    title(sprintf('Subset 3, #cells = %d',length(light_sensing_cells)))
    set(gca,'xtick',[])
    
    
    % ON cells
    %----------------------------------------------------------------------   
    
    on_cells = [];
    subplot(4,1,2)
    for cell = 1:size(cell_comp_matrix,1)
        if mean(cell_comp_matrix(cell,10:15)/mean(cell_comp_matrix(cell,16:21))) > 2 && mean(cell_comp_matrix(cell,10:15)) > 0 

            on_cells = [on_cells, cell];
            plot(x_axis,rescale(cell_comp_matrix(cell,:)),'Color',[0.5 0.5 0.5])
            hold on
        end

    end

    plot(x_axis,rescale(mean(cell_comp_matrix(on_cells,:))),'Color','w','LineWidth',2.5)
    ylabel('dFF (NORM)')
    title(sprintf('Subset 3, #cells = %d',length(on_cells)))
    set(gca,'xtick',[])


 
    % OFF cells
    %----------------------------------------------------------------------   
    
    off_cells = [];
    subplot(4,1,3)
    for cell = 1:size(cell_comp_matrix,1)
        if mean(cell_comp_matrix(cell,23:32)/mean(cell_comp_matrix(cell,16:23))) > 2 && mean(cell_comp_matrix(cell,10:15)) > 0 

            off_cells = [off_cells, cell];
            plot(x_axis,rescale(cell_comp_matrix(cell,:)),'Color',[0.5 0.5 0.5])
            hold on
        end

    end

    plot(x_axis,rescale(mean(cell_comp_matrix(off_cells,:))),'Color','w','LineWidth',2.5)
    ylabel('dFF (NORM)')
    title(sprintf('Subset 3, #cells = %d',length(off_cells)))
    set(gca,'xtick',[])

    subplot(4,1,4)
    plot(x_axis_stim, stim,"Color","w")
    ylabel({'INTENSITY', '(NORM)'})
    xlabel('TIME / s')
    
    
    
    set(gcf, 'Color', [0 0 0]);
    set(gcf, 'InvertHardCopy', 'off');
    
    
    
    
    
    
end

