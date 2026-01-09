function cellYield = getCellYieldInfo(saveNeuralPath)

%% load session data %%
load([saveNeuralPath '\rawDataBySessionNeural.mat'])

%% determine cell yield %%
pyrCellCtr = 0;
narrowCellCtr = 0;

for clu = 1:length(rawDataBySessionNeural.apData)
    if contains(rawDataBySessionNeural.apData(clu).putativeCellType, 'Pyr')
        pyrCellCtr = pyrCellCtr + 1;
    elseif contains(rawDataBySessionNeural.apData(clu).putativeCellType, 'Narrow')
        narrowCellCtr = narrowCellCtr + 1;                        
    end%if contains
end%clu

cellYield.pyr = pyrCellCtr;
cellYield.narrow = narrowCellCtr;

%% save %%
filename = [saveNeuralPath '\' 'cellYield.mat'];
save(filename, 'cellYield', '-v7.3')


% % %after cell yields per session, the code below can be
% % % used to calculate cell yield across sessions
% % 
% % sessCtr = 0;
% % cellYiledAll = [];
% % 
% % for i =  1:size(allindex,1)%1:size(allindex,1) %loop through every session
% %     if allindex(i,5) > 0%recording sessions only
% %         %%%%% session info %%%%%
% %         subj = [params.iden num2str(allindex(i,1))];
% %         sessDate = num2str(allindex(i,2));
% %         sessNum = num2str(allindex(i,3));
% %         trackInfo = num2str(allindex(i,4));
% %         %virmenDataPath = fullfile(dirs.virmenrawdata, ['\' subj '_' sessDate '_' sessNum '\' 'dataWithLickometer.mat']);
% %         virmenSessDataPath = fullfile(dirs.saveoutputstructs, ['Data\Behavior\sessionData\' subj '\' sessDate '_' sessNum '_' trackInfo]);
% %         neuralRawDataPath = fullfile(dirs.rawdata, [subj '_' sessDate]);
% %         saveNeuralPath = fullfile(dirs.saveoutputstructs, ['Data\Neural\sessionData\' subj '\' sessDate '_' sessNum '_' trackInfo]);
% %         postSpikeSort = 1;
% % 
% %         if str2num(sessNum) == 3
% %             sessCtr = sessCtr + 1;
% % 
% %             load([saveNeuralPath '\' 'cellYield.mat'])
% % 
% %             cellYiledAll(sessCtr,1) = cellYield.pyr;
% %             cellYiledAll(sessCtr,2) = cellYield.narrow;
% %         end
% %     end
% % end

end%function