function allindex = selectindextable_JLK(spreadsheetfile, varargin)
%Filters:
% Track: 1 = original, 2 = update
% Session Type: 1 = passive, 2 = active, 3 = rest

%% get selection info
filters.animal = nan;           %default is to select none
filters.id = nan;               %default is to include all iden
filters.datesincluded = nan;    %default is to include all
filters.datesexcluded = nan;    %default is to exclude none
filters.sessionnum = nan;       %default is to include all sessions
filters.sessiontype = nan;      %default is to include only active
filters.track = nan;            %default is to include all tracks
filters.recording = nan;        %default is to include regardless of recording

%get input arguments
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'rewritefileinfo'
            rewritefileinfo = varargin{option+1};
        case 'animal'
            filters.animal = varargin{option+1};
        case 'id'
            filters.id = varargin{option+1};
        case 'datesincluded'
            filters.datesincluded = varargin{option+1};
        case 'datesexcluded'
            filters.datesexcluded = varargin{option+1};
        case 'SessionNum'
            filters.sessionnum = varargin{option+1};
        case 'SessionType'
            filters.sessiontype = varargin{option+1};
        case 'Track'
            filters.track = varargin{option+1};
        case 'Recording'
            filters.recording = varargin{option+1};
        otherwise
            disp(['Warning: ' varargin{option} ' is an invalid input argument'])
    end
end

%% read in the spreadsheet
opts = detectImportOptions(spreadsheetfile);
ephysT = readtable(spreadsheetfile, opts);
headerNames = ephysT.Properties.VariableNames;

%% apply selected filters
fnames = fieldnames(filters);
filtersApplied = fnames(cellfun(@(x) any(~isnan(filters.(x))), fnames)); % TODO - extract string based filters

% get included rows for each filter
allFilters = zeros(size(ephysT,1),numel(filtersApplied));
for filtIdx = 1:numel(filtersApplied)
    filtName = filtersApplied{filtIdx};        
    if any(contains(headerNames, filtName, 'IgnoreCase', true)) && contains(filtName, 'track', 'IgnoreCase', true) == 0 && contains(filtName, 'sessiontype', 'IgnoreCase', true) == 0%if the filter name applies to a specific column of the table
        filtCol = headerNames{contains(headerNames, filtName, 'IgnoreCase', true)};
        if iscell(ephysT.(filtCol))
            allFilters(:,filtIdx) = any(ismember(ephysT.(filtCol), filters.(filtName)),2);
        else
            allFilters(:,filtIdx) = any(ephysT.(filtCol) == filters.(filtName),2); %gets any matching criteria if multiple filter items
        end
    elseif contains(filtName, 'date', 'IgnoreCase', true) %special case to deal with dates
        if isempty(filters.datesexcluded)
            filters.datesexcluded = nan;
        end
        if isempty(filters.datesincluded)
            filters.datesincluded = nan;
        end
        if isnan(filters.datesincluded) %default is to include all dates
            allFilters(:,filtIdx) = all((ephysT.Date ~= filters.datesexcluded),2);
        else
            allFilters(:,filtIdx) = (any(ephysT.Date == filters.datesincluded, 2) & (ephysT.Date ~= filters.datesexcluded));
        end
    elseif contains(filtName, 'track', 'IgnoreCase', true) %special case to deal with tracks below
            allFilters(:,filtIdx) = 1;
    elseif contains(filtName, 'sessiontype', 'IgnoreCase', true) %special case to deal with sessiontype below
            allFilters(:,filtIdx) = 1;
    else
        disp(['Warning: no corresponding column found for the ' filtName...
            ' filter. Check the spreadsheet header names'])
    end
end

% apply all filters to the table
allindex = ephysT(all(allFilters,2),:);

%% Change Track from strings to numbers
newtrackcol = strrep(table2cell(allindex(:,strcmp(headerNames, 'Track'))), 'TrackA', "1"); %find TrackA (original)
newtrackcol = strrep(newtrackcol, "1'", "2"); %find TrackA' (update)
newtrackcol = strrep(newtrackcol, 'TrackB', "3"); %find TrackB (novel)
newtrackcol = strrep(newtrackcol, 'TrackC', "4"); %find TrackC (novel2)
newtrackcol = table(double(newtrackcol));%convert to double then table
allindex(:,strcmp(headerNames, 'Track')) = []; %delete track column
allindex = [allindex newtrackcol]; %add new track column
allindex.Properties.VariableNames(end) = {'Track'}; %name new track column "Track"

%now filter for tracks of interest
if ~isnan(filters.track)
    trind = zeros(size(allindex,1),1);
    for tr = filters.track
        trind = trind + (table2array(allindex(:,strcmp(allindex.Properties.VariableNames, 'Track')))==tr);
    end
    allindex = allindex(logical(trind),:);
end

%% Change Session Type from strings to numbers
headerNames = allindex.Properties.VariableNames;%changing track column number changed the order of header names, so re-define here
newsessiontypecol = strrep(table2cell(allindex(:,strcmp(headerNames, 'SessionType'))), 'passive', "1"); %find passive sessions
newsessiontypecol = strrep(newsessiontypecol, 'active', "2"); %find active sessions
newsessiontypecol = strrep(newsessiontypecol, 'rest', "3"); %find rest sessions
newsessiontypecol = table(double(newsessiontypecol));%convert to double then table
allindex(:,strcmp(headerNames, 'SessionType')) = []; %delete sessiontype column
allindex = [allindex newsessiontypecol]; %add new sessiontype column
allindex.Properties.VariableNames(end) = {'SessionType'}; %name new track column "SessionType"

%now filter for sessiontype of interest
if ~isnan(filters.sessiontype)
    sind = zeros(size(allindex,1),1);
    for s = filters.sessiontype
        sind = sind + (table2array(allindex(:,strcmp(allindex.Properties.VariableNames, 'SessionType')))==s);
    end
    allindex = allindex(logical(sind),:);
end

end