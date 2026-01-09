function CellExplorer_Pipeline_DC(subj, sessDate, processedDataPath, params)

%% Cell type classification via CellExplorer %%
%adapted from CellExplorerSingerLabTesting.m
%Note: CellExplorer uses Kilosort output, and thus uses data from all sessions
% on a given day to classify putative pyramidal cells and interneurons.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Create and Validate Session Struct %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adapted from NovelTaskTemplate.m
session = [];
%%%%% add general metadata %%%%%
session.general.basePath =  [processedDataPath '\kilosort4']; % Full path
session.general.name = [subj '_' sessDate]; % Session name / basename
session.general.version = 5; % Metadata version
session.general.sessionType = 'Acute'; % Type of recording: Chronic, Acute, Unknown
session.general.sessionName = session.general.name;

%%%%% add animal data %%%%%
session.animal.name = session.general.name; % Animal name is inferred from the data path
session.animal.sex = 'Female'; % Male, Female, Unknown
session.animal.species = 'Mouse'; % Mouse, Rat
session.animal.strain = 'C57';
session.animal.geneticLine = '';

%%%%% add extracellular data %%%%%
session.extracellular.fileName = dir([session.general.basePath '\' 'allrecordings.bin']).name;
%session.extracellular.leastSignificantBit = 0.1950; %for Neuropixels
session.extracellular.probeDepths = 0;
session.extracellular.precision = 'int16';
session.extracellular.sr = params.samprate;
session.extracellular.nChannels = size(params.probeChannels{1,1},2);
session.extracellular.nElectrodeGroups = 1;
session.extracellular.electrodeGroups.channels = {1:session.extracellular.nChannels};
session.extracellular.nSpikeGroups = 1;
session.extracellular.spikeGroups.channels = {1:session.extracellular.nChannels};
channelPositions = readNPY([session.general.basePath '\' 'channel_positions.npy']);
session.extracellular.chanCoords.x = channelPositions(:,1);
session.extracellular.chanCoords.y = channelPositions(:,2);
session.extracellular.chanCoords.source = 'Kilosort';
session.extracellular.chanCoords.verticalSpacing = 28;
session.extracellular.chanCoords.layout = 'poly5';% Probe layout: linear,staggered,poly2,edge,poly3,poly5

%%%%% add spikeSorting data %%%%%
session.spikeSorting{1}.relativePath = '';
session.spikeSorting{1}.format = 'Phy';
session.spikeSorting{1}.method = 'Kilosort';
session.spikeSorting{1}.channels = [];
session.spikeSorting{1}.manuallyCurated = 1;
session.spikeSorting{1}.notes = '';

%add brainRegions data
session.brainRegions.HIP.channels = 1:session.extracellular.nChannels;
session.brainRegions.HIP.electrodeGroups = 1;

%%%%% add epochs data %%%%%
session.epochs{1}.name = session.general.name;
session.epochs{1}.startTime = 0;

%%%%% validate session struct %%%%%
validateSessionStruct(session);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Run the cell metrics pipeline %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exclude_metrics = {'monoSynaptic_connections','spatial_metrics','event_metrics','manipulation_metrics'};
ProcessCellMetrics('session', session, 'excludeMetrics', exclude_metrics,'forceReload', false, 'keepCellClassification', ~false);
%Output notes: This function produces the following stucts: cell_metrics, noiseLevel, session, and spikes. 
% The function includes only clusters manually spike sorted as "good". Clusters, clusters_allrec, and
%clustermetrics take "good" clusters from manual spike sorting "good" and removes clusters that do not meet lab-determined criteria.
%Thus, clusters, clusters_allrec, and clustermetrics include fewer clusters than cell_metrics, noiseLevel, session, and spikes structs, and
%clusters_allrec.ID can be used to index those structs ( e.g., spikes.ts(ismember(spikes.cluID, [clusters_allrec.ID])) ).
%  For spike times, cell_metrics.spikes.times == spikes.times in seconds (so multiplying by 30k will convert to samples),
%and spikes.ts are already in samples.

end%function