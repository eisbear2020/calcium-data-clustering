data_SC_new = CaData(SI_Aux_file_SC, Ca_file_SC);
data_SC_new.SyncStimCaData(new_chirp_SC.stim_start_times, new_chirp_SC.stim_end_times)
data_SC_new.cutOffFilter(7)
data_SC_new.averageOverTrials();
data_SC_new.SplitIntoTrials();
%%
CaData = data_SC_new.dFF_sep_trials{1}(39,:);
stim = stim_raw_SC;

x_axis = linspace(0,30,size(CaData,2));
x_axis_stim = linspace(0,30,length(stim));

figure
subplot(2,1,1)
plot(x_axis,CaData,'Color','w','LineWidth',2.5)
%plot(x_axis,cell_comp_matrix(light_sensing_cells(7),:),'Color','k','LineWidth',1.5)
ylabel('dFF')
set(gca,'xtick',[])

subplot(2,1,2)
plot(x_axis_stim, stim,'Color','w')
ylabel({'INTENSITY','(NORM)'})
xlabel('TIME / s')