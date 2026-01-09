function rat = get_and_save_joshs_64lfp_data_jrm(rat)
%rat = get_and_save_each_rats_64lfp_data(rat, secbefore,secafter, samprate, minexplorationtime, sourcefoldername)
%secbefore and secafter specify the amount of lfp (in seconds) to save before and after each >minexplorationtime exploration event
%minexploration time specifies the smallest exploration epoch for which we want to store LFP data

%sourcefoldername is the location of where the raw data are stored, e.g., 'RAWDATA'

% %get parameters for filtering: 1.5-400 Hz bandpass with .5 Hz shoulders
% [n,fo,mo,w] = remezord([1 1.5 400 400.5], [0 1 0], [1 1 1]*0.01, samprate);
% b = remez(n,fo,mo,w);
% a=1;


% JLK: Commented out 8/10/21 after specifiying BLA + HPC (pyramidal layer)
% channels in specify_joshs_probe_variables_jrm.m
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %Note that I currently use the same 8 channels for all rats, and this temporarily changes the original probe info from specify_joshs_probe_variables_jrm
% % %Delete this section if I end up using all channels, or change specify_joshs_probe_variables_jrm to accutately depict the channels used for each rat
% % for r = 1:length(rat)
% %     rat(r).VtoDchanorder = ...
% %         [9  8   10   7  11  6   12  5   25  24  26  23  27  22  28  21  29  20  30  19  31  18  32  17  1   16   2  15  3   14  4   13, ...
% %         41	40  42	39	43	38	44	37	57	56	58	55	59	54	60	53	61	52	62	51	63	50	64	49	33	48	34	47	35	46	36	45
% %         ];
% %     
% %     rat(r).blachannels = [1, 9, 18, 26];
% %     rat(r).hipchannels = [37, 44, 51, 60];
% % end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%other than sampling rate, I don't think any of these parameters
%actually matter to the actual cleaning process in rmlinesmovingwinc,
%which is a Chronux function and the only reason the params
%structure is here.
samprate=1500;
minexploration = 1;%for now, minimum object exploration is 1s
secbefore = 1;

params.fpass=[1 400]; % band of frequencies to be kept
params.tapers=[6 11]; % taper parameters
params.pad=1; % pad factor for fft
params.trialave=1;
params.Fs = samprate;

for r = 1:length(rat)
    fprintf('Working on rat %d (%s) . . .\n', r, rat(r).name)
    cd(rat(r).name)
    
    for s = 2:length(rat(r).set)
        fprintf('\t...set %d\n', s)
        for d = 1:length(rat(r).set(s).day)
            fprintf('\t...day %d\n', d)
            
            if exist(rat(r).set(s).day(d).foldername,'dir') ~= 0
                cd(rat(r).set(s).day(d).foldername)
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %     LOAD ALL 64 LFPS FOR THIS DAY            %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                tmpVtoDorder = rat(r).VtoDchanorder;%copy into temporary variable so the whole rat struct isn't "broadcast" in parfor loop
                %BLA PROBE
                bla = struct('lfp', [], 'timestamps', []);
                tmpblachans = rat(r).blachannels;%copy into temporary variable so the whole rat struct isn't "broadcast" in parfor loop
                parfor ch = 1:length(rat(r).blachannels)
                    thischnum = tmpVtoDorder(tmpblachans(ch));
                    thisfn = sprintf('ds_113_CH%d.continuous.mat',thischnum);%file name for this channel
                    thisdata = load(thisfn);
                    if( round(thisdata.timestamps(samprate+1)-thisdata.timestamps(1),2) ~= 1 )
                        fprintf('Warning: Sampling rate does not appear to be %d for file %s.\n', samprate, thisfn)
                    end%end if
                    
                    %temp comment out bla(ch).timestamps = thisdata.timestamps;
%                     tmpnumsamps = length(thisdata.timestamps);
%                     bla(ch).timestamps = linspace(0,tmpnumsamps/samprate, tmpnumsamps)';
                    bla(ch).timestamps = thisdata.timestamps;
                    
                    %filter, commented out for now
                    %thisdata.lfp = filtfilt(b,a, thisdata.lfp);
                    
                    %for every 4 sec time bin (overlapping by 2 sec), fits a 60 Hz sine wave to data and subtracts it out
                    %p value set very low (.00000001) becuase I couldn't tell if function would try to
                    %automatically take out other significant sine waves
                    %tau set at 10 based on chronux recommendation
                    thisdata.lfp=rmlinesmovingwinc(thisdata.lfp,[4 2],10,params,.00000001,'n', 60);
                    newlength = length(thisdata.lfp);
                    
                    bla(ch).lfp = thisdata.lfp;
                    bla(ch).timestamps = bla(ch).timestamps(1:newlength);
                    
                    
                end%end ch loop
                clear thisdata
                
                %HIPPOCAMPUS PROBE
                hip = struct('lfp', [], 'timestamps', []);
                tmphipchans = rat(r).hipchannels;%copy into temporary variable so the whole rat struct isn't "broadcast" in parfor loop
                parfor ch = 1:length(rat(r).hipchannels)%will be 33-64 for all rats
                    thischnum = tmpVtoDorder(tmphipchans(ch));%e.g., will be 48 for ch=1 for Layla, will be 41 for ch=1 for all others
                    thisfn = sprintf('ds_113_CH%d.continuous.mat',thischnum);%file name for this channel
                    thisdata = load(thisfn);
                    if( round(thisdata.timestamps(samprate+1)-thisdata.timestamps(1),2) ~= 1 )
                        fprintf('Warning: Sampling rate does not appear to be %d for file %s.\n', samprate, thisfn)
                    end%end if
                    
%                     %temp comment out  hip(ch).timestamps = thisdata.timestamps;
%                     tmpnumsamps = length(thisdata.timestamps);
%                     hip(ch).timestamps = linspace(0,tmpnumsamps/samprate, tmpnumsamps)';
                    hip(ch).timestamps = thisdata.timestamps;
                    
                    %filter, commented out for now
                    %thisdata.lfp = filtfilt(b,a, thisdata.lfp);
                    
                    %for every 4 sec time bin (overlapping by 2 sec), fits a 60 Hz sine wave to data and subtracts it out
                    %p value set very low (.00000001) becuase I couldn't tell if function would try to
                    %automatically take out other significant sine waves
                    %tau set at 10 based on chronux recommendation
                    thisdata.lfp=rmlinesmovingwinc(thisdata.lfp,[4 2],10,params,.00000001,'n', 60);
                    newlength = length(thisdata.lfp);
                    
                    hip(ch).lfp = thisdata.lfp;
                    hip(ch).timestamps = hip(ch).timestamps(1:newlength);
                    
                end%end ch loop
                clear thisdata
                
                %at this point bla(1)-bla(32) should be the bla data in ventral to dorsal order
                %at this point hip(1)-hip(32) should be the hip data in ventral to dorsal order
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %     Extract LFPs for Lap Struct            %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                lap = rat(r).set(s).day(d).lap;
                
                for lp = 1:length(lap)
                    
                    %%%%%  Start with context exploration1  %%%%%
                    probe = struct('name', []);
                    probe(1).name = 'BLA';
                    probe(2).name = 'HIP';
                    
                    %tmptsst = lap(lp).contexploration1.start;%temp starting timestamp (15s before obj presentation)
                    tmptsend = lap(lp).contexploration1.end;%temp ending timestamp (obj presentation)
                    %BLA
                    for ch = 1:length(rat(r).blachannels)
                        %we won't assume that the timestamps are exactly aligned across channels
                        %startindex = find(bla(ch).timestamps >= tmptsst, 1, 'first');
                        stopindex = find(bla(ch).timestamps <= tmptsend, 1, 'last');
                        %tmplfp = bla(ch).lfp(startindex:stopindex);%15s before presenting objects to object presentation
                        tmplfp = bla(ch).lfp(stopindex-22499:stopindex-4500);%12s before presenting objects to object presentation (15s to 3s before)
                        probe(1).lfp(ch,:) = tmplfp;
                    end%end
                    
                    %Hippocampus
                    for ch = 1:length(rat(r).hipchannels)
                        %we won't assume that the timestamps are exactly aligned across channels
                        %startindex = find(bla(ch).timestamps >= tmptsst, 1, 'first');
                        stopindex = find(hip(ch).timestamps <= tmptsend, 1, 'last');
                        %tmplfp = hip(ch).lfp(startindex:stopindex);%15s before presenting objects to object presentation
                        tmplfp = hip(ch).lfp(stopindex-22499:stopindex-4500);%12s before presenting objects to object presentation (15s to 3s before)
                        probe(2).lfp(ch,:) = tmplfp;
                    end%end
                    
                    %add probe to struct
                    lap(lp).contexploration1.probe = probe;
                    
                    %clear tmplfp for next loop
                    clear probe
                    
                    %%%%%  Next, context exploration2  %%%%%
                    if ~isempty(lap(lp).contexploration2)
                        probe = struct('name', []);
                        probe(1).name = 'BLA';
                        probe(2).name = 'HIP';
                        
                        
                        %tmptsst = lap(lp).contexploration2.start;%temp starting timestamp (15s before obj presentation)
                        tmptsend = lap(lp).contexploration2.end;%temp ending timestamp (obj presentation)
                        %BLA
                        for ch = 1:length(rat(r).blachannels)
                            %we won't assume that the timestamps are exactly aligned across channels
                            %startindex = find(bla(ch).timestamps >= tmptsst, 1, 'first');
                            stopindex = find(bla(ch).timestamps <= tmptsend, 1, 'last');
                            %tmplfp = bla(ch).lfp(startindex:stopindex);%15s before presenting objects to object presentation
                            tmplfp = bla(ch).lfp(stopindex-22499:stopindex-4500);%12s before presenting objects to object presentation (15s to 3s before)
                            probe(1).lfp(ch,:) = tmplfp;
                        end%end
                        
                        %Hippocampus
                        for ch = 1:length(rat(r).hipchannels)
                            %we won't assume that the timestamps are exactly aligned across channels
                            %startindex = find(bla(ch).timestamps >= tmptsst, 1, 'first');
                            stopindex = find(hip(ch).timestamps <= tmptsend, 1, 'last');
                            %tmplfp = hip(ch).lfp(startindex:stopindex);%15s before presenting objects to object presentation
                            tmplfp = hip(ch).lfp(stopindex-22499:stopindex-4500);%12s before presenting objects to object presentation (15s to 3s before)
                            probe(2).lfp(ch,:) = tmplfp;
                        end%end
                        
                        %add probe to struct
                        lap(lp).contexploration2.probe = probe;
                        
                        %clear tmplfp for next loop
                        clear probe
                        
                    end%if ~isempty for context exploration2
                    
                    %%%%%  Next, object1 exploration during presentation1  %%%%%
                    if ~isempty(lap(lp).obj1pres1)%determine if object1 was explored this presentation
                        for exp = 1:length(lap(lp).obj1pres1)%if yes, do the following for each time the object was explored
                            if lap(lp).obj1pres1(exp).duration > minexploration%but only do so if this exploration surpasses the minimum exploration
                                probe = struct('name', []);
                                probe(1).name = 'BLA';
                                probe(2).name = 'HIP';
                                
                                tmptsst = lap(lp).obj1pres1(exp).start;%temp ending timestamp (obj presentation)
                                %BLA
                                for ch = 1:length(rat(r).blachannels)
                                    %we won't assume that the timestamps are exactly aligned across channels
                                    tmpindex = find(bla(ch).timestamps >= tmptsst, 1, 'first');
                                    startindex = tmpindex - secbefore*samprate;
                                    stopindex  = tmpindex + secbefore*samprate;
                                    tmplfp = bla(ch).lfp(startindex:stopindex-1);%+/- secbefore/secafter centered and aligned at onset
                                    probe(1).lfp(ch,:) = tmplfp;
                                end%end
                                
                                %Hippocampus
                                for ch = 1:length(rat(r).hipchannels)
                                    %we won't assume that the timestamps are exactly aligned across channels
                                    tmpindex = find(hip(ch).timestamps >= tmptsst, 1, 'first');
                                    startindex = tmpindex - secbefore*samprate;
                                    stopindex  = tmpindex + secbefore*samprate;
                                    tmplfp = hip(ch).lfp(startindex:stopindex-1);%+/- secbefore/secafter centered and aligned at onset
                                    probe(2).lfp(ch,:) = tmplfp;
                                end%end
                                
                                %add probe to struct
                                lap(lp).obj1pres1(exp).probe = probe;
                                
                                %clear tmplfp for next loop
                                clear probe
                            end%if >minexploration
                        end%for exp
                    end%if ~isempty for obj1pres1
                    
                    %%%%%  Next, object2 exploration during presentation1  %%%%%
                    if ~isempty(lap(lp).obj2pres1)%determine if object2 was explored this presentation
                        for exp = 1:length(lap(lp).obj2pres1)%if yes, do the following for each time the object was explored
                            if lap(lp).obj2pres1(exp).duration > minexploration%but only do so if this exploration surpasses the minimum exploration
                                probe = struct('name', []);
                                probe(1).name = 'BLA';
                                probe(2).name = 'HIP';
                                
                                tmptsst = lap(lp).obj2pres1(exp).start;%temp ending timestamp (obj presentation)
                                %BLA
                                for ch = 1:length(rat(r).blachannels)
                                    %we won't assume that the timestamps are exactly aligned across channels
                                    tmpindex = find(bla(ch).timestamps >= tmptsst, 1, 'first');
                                    startindex = tmpindex - secbefore*samprate;
                                    stopindex  = tmpindex + secbefore*samprate;
                                    tmplfp = bla(ch).lfp(startindex:stopindex-1);%+/- secbefore/secafter centered and aligned at onset
                                    probe(1).lfp(ch,:) = tmplfp;
                                end%end
                                
                                %Hippocampus
                                for ch = 1:length(rat(r).hipchannels)
                                    %we won't assume that the timestamps are exactly aligned across channels
                                    tmpindex = find(hip(ch).timestamps >= tmptsst, 1, 'first');
                                    startindex = tmpindex - secbefore*samprate;
                                    stopindex  = tmpindex + secbefore*samprate;
                                    tmplfp = hip(ch).lfp(startindex:stopindex-1);%+/- secbefore/secafter centered and aligned at onset
                                    probe(2).lfp(ch,:) = tmplfp;
                                end%end
                                
                                %add probe to struct
                                lap(lp).obj2pres1(exp).probe = probe;
                                
                                %clear tmplfp for next loop
                                clear probe
                            end%if >minexploration
                        end%for exp
                    end%if ~isempty for obj2pres1
                    
                    %%%%%  Next, object1 exploration during presentation2  %%%%%
                    if ~isempty(lap(lp).obj1pres2)%determine if object1 was explored this presentation
                        for exp = 1:length(lap(lp).obj1pres2)%if yes, do the following for each time the object was explored
                            if lap(lp).obj1pres2(exp).duration > minexploration%but only do so if this exploration surpasses the minimum exploration
                                probe = struct('name', []);
                                probe(1).name = 'BLA';
                                probe(2).name = 'HIP';
                                
                                tmptsst = lap(lp).obj1pres2(exp).start;%temp ending timestamp (obj presentation)
                                %BLA
                                for ch = 1:length(rat(r).blachannels)
                                    %we won't assume that the timestamps are exactly aligned across channels
                                    tmpindex = find(bla(ch).timestamps >= tmptsst, 1, 'first');
                                    startindex = tmpindex - secbefore*samprate;
                                    stopindex  = tmpindex + secbefore*samprate;
                                    tmplfp = bla(ch).lfp(startindex:stopindex-1);%+/- secbefore/secafter centered and aligned at onset
                                    probe(1).lfp(ch,:) = tmplfp;
                                end%end
                                
                                %Hippocampus
                                for ch = 1:length(rat(r).hipchannels)
                                    %we won't assume that the timestamps are exactly aligned across channels
                                    tmpindex = find(hip(ch).timestamps >= tmptsst, 1, 'first');
                                    startindex = tmpindex - secbefore*samprate;
                                    stopindex  = tmpindex + secbefore*samprate;
                                    tmplfp = hip(ch).lfp(startindex:stopindex-1);%+/- secbefore/secafter centered and aligned at onset
                                    probe(2).lfp(ch,:) = tmplfp;
                                end%end
                                
                                %add probe to struct
                                lap(lp).obj1pres2(exp).probe = probe;
                                
                                %clear tmplfp for next loop
                                clear probe
                            end%if >minexploration
                        end%for exp
                    end%if ~isempty for obj1pres2
                    
                    %%%%%  Next, object2 exploration during presentation2  %%%%%
                    if ~isempty(lap(lp).obj2pres2)%determine if object2 was explored this presentation
                        for exp = 1:length(lap(lp).obj2pres2)%if yes, do the following for each time the object was explored
                            if lap(lp).obj2pres2(exp).duration > minexploration%but only do so if this exploration surpasses the minimum exploration
                                probe = struct('name', []);
                                probe(1).name = 'BLA';
                                probe(2).name = 'HIP';
                                
                                tmptsst = lap(lp).obj2pres2(exp).start;%temp ending timestamp (obj presentation)
                                %BLA
                                for ch = 1:length(rat(r).blachannels)
                                    %we won't assume that the timestamps are exactly aligned across channels
                                    tmpindex = find(bla(ch).timestamps >= tmptsst, 1, 'first');
                                    startindex = tmpindex - secbefore*samprate;
                                    stopindex  = tmpindex + secbefore*samprate;
                                    tmplfp = bla(ch).lfp(startindex:stopindex-1);%+/- secbefore/secafter centered and aligned at onset
                                    probe(1).lfp(ch,:) = tmplfp;
                                end%end
                                
                                %Hippocampus
                                for ch = 1:length(rat(r).hipchannels)
                                    %we won't assume that the timestamps are exactly aligned across channels
                                    tmpindex = find(hip(ch).timestamps >= tmptsst, 1, 'first');
                                    startindex = tmpindex - secbefore*samprate;
                                    stopindex  = tmpindex + secbefore*samprate;
                                    tmplfp = hip(ch).lfp(startindex:stopindex-1);%+/- secbefore/secafter centered and aligned at onset
                                    probe(2).lfp(ch,:) = tmplfp;
                                end%end
                                
                                %add probe to struct
                                lap(lp).obj2pres2(exp).probe = probe;
                                
                                %clear tmplfp for next loop
                                clear probe
                            end%if >minexploration
                        end%for exp
                    end%if ~isempty for obj2pres2
                    
                    %%%%%  Next, lastexplorationpres1  %%%%%
                    probe = struct('name', []);
                    probe(1).name = 'BLA';
                    probe(2).name = 'HIP';
                    
                    tmptsst = lap(lp).lastexplorationpres1.start;%temp ending timestamp (obj presentation)
                    %BLA
                    for ch = 1:length(rat(r).blachannels)
                        %we won't assume that the timestamps are exactly aligned across channels
                        tmpindex = find(bla(ch).timestamps >= tmptsst, 1, 'first');
                        startindex = tmpindex - secbefore*samprate;
                        stopindex  = tmpindex + secbefore*samprate;
                        tmplfp = bla(ch).lfp(startindex:stopindex-1);%+/- secbefore/secafter centered and aligned at onset
                        probe(1).lfp(ch,:) = tmplfp;
                    end%end
                    
                    %Hippocampus
                    for ch = 1:length(rat(r).hipchannels)
                        %we won't assume that the timestamps are exactly aligned across channels
                        tmpindex = find(hip(ch).timestamps >= tmptsst, 1, 'first');
                        startindex = tmpindex - secbefore*samprate;
                        stopindex  = tmpindex + secbefore*samprate;
                        tmplfp = hip(ch).lfp(startindex:stopindex-1);%+/- secbefore/secafter centered and aligned at onset
                        probe(2).lfp(ch,:) = tmplfp;
                    end%end
                    
                    %add probe to struct
                    lap(lp).lastexplorationpres1.probe = probe;
                    
                    %clear tmplfp for next loop
                    clear probe
                    
                    %%%%%  Last, lastexplorationpres2  %%%%%
                    if ~isempty(lap(lp).lastexplorationpres2)
                        probe = struct('name', []);
                        probe(1).name = 'BLA';
                        probe(2).name = 'HIP';
                        
                        tmptsst = lap(lp).lastexplorationpres2.start;%temp ending timestamp (obj presentation)
                        %BLA
                        for ch = 1:length(rat(r).blachannels)
                            %we won't assume that the timestamps are exactly aligned across channels
                            tmpindex = find(bla(ch).timestamps >= tmptsst, 1, 'first');
                            startindex = tmpindex - secbefore*samprate;
                            stopindex  = tmpindex + secbefore*samprate;
                            tmplfp = bla(ch).lfp(startindex:stopindex-1);%+/- secbefore/secafter centered and aligned at onset
                            probe(1).lfp(ch,:) = tmplfp;
                        end%end
                        
                        %Hippocampus
                        for ch = 1:length(rat(r).hipchannels)
                            %we won't assume that the timestamps are exactly aligned across channels
                            tmpindex = find(hip(ch).timestamps >= tmptsst, 1, 'first');
                            startindex = tmpindex - secbefore*samprate;
                            stopindex  = tmpindex + secbefore*samprate;
                            tmplfp = hip(ch).lfp(startindex:stopindex-1);%+/- secbefore/secafter centered and aligned at onset
                            probe(2).lfp(ch,:) = tmplfp;
                        end%end
                        
                        %add probe to struct
                        lap(lp).lastexplorationpres2.probe = probe;
                        
                        %clear tmplfp for next loop
                        clear probe
                        
                    end%if ~isempty for lastexplorationpres2
                    
                end%lap
                
                %add to return struct
                rat(r).set(s).day(d).lap = lap;
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %     Extract LFPs for Sleep Struct  %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if d == 1
                    sleep = rat(r).set(s).day(d).sleep;
                    minduration = 5;
                    
                    if ~isempty(sleep)
                        for sl = 1:length(sleep)
                            if sleep(sl).duration > minduration%but only do so if this exploration surpasses the minimum exploration
                                probe = struct('name', []);
                                probe(1).name = 'BLA';
                                probe(2).name = 'HIP';

                                tmptsst = sleep(sl).start;%temp ending timestamp (obj presentation)
                                %BLA
                                for ch = 1:length(rat(r).blachannels)
                                    %we won't assume that the timestamps are exactly aligned across channels
                                    tmpindex = find(bla(ch).timestamps >= tmptsst, 1, 'first');
                                    startindex = tmpindex - (minduration/2)*samprate;
                                    stopindex  = tmpindex + (minduration/2)*samprate;
                                    tmplfp = bla(ch).lfp(startindex:stopindex-1);%+/- secbefore/secafter centered and aligned at onset
                                    probe(1).lfp(ch,:) = tmplfp;
                                end%end

                                %Hippocampus
                                for ch = 1:length(rat(r).hipchannels)
                                    %we won't assume that the timestamps are exactly aligned across channels
                                    tmpindex = find(hip(ch).timestamps >= tmptsst, 1, 'first');
                                    startindex = tmpindex - (minduration/2)*samprate;
                                    stopindex  = tmpindex + (minduration/2)*samprate;
                                    tmplfp = hip(ch).lfp(startindex:stopindex-1);%+/- secbefore/secafter centered and aligned at onset
                                    probe(2).lfp(ch,:) = tmplfp;
                                end%end

                                %add probe to struct
                                sleep(sl).probe = probe;

                                %clear tmplfp for next loop
                                clear probe
                            end%if >minduration
                        end%for sleep
                    end%if ~isempty for sleep
                    
                    %add to return struct
                    rat(r).set(s).day(d).sleep = sleep;
                    
                end %if d == 1
                
                cd ..
                
            else %if doesn't exist, skip folder
                fprintf('\tSkipped set %d day %d for %s.\n', s, d, rat(r).name)
            end %if exist
        end%day
    end%set
    cd ..
end%rat
end %function



%%%%JLK: Will need to do something similar to below for days of incorrect
%%%%plug in

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % now fix the data for Layla Day 3 Trial 8, when the hippocampus probe was plugged in the wrong way %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % r=1;d=3;t=8;p=2;
% %
% % %the hippocampus probe Layla was plugged in the wrong way for this trial
% % %but since this script already reordered the channels into dorsal-to-ventral order, we need to find out the new remapping
% % %tmpch below was determined on 11/11/20 by JRM manually--after we figured out that OE ignored remapping info
% % % AND after we figured out that Layla's Intan chips were plugged in differently than originally assumed
% % %tmpch is the indices of the rempped channels
% % %so, after the chip was plugged in the wrong way, the 6tth lfp is actually the most ventral contact
% % % and the 9th lfp is actually the most dorsal
% %
% %
% % tmpch = [6,5,8,7,2,1,4,3,32,31,28,27,18,17,22,21,14,13,26,25,16,15,30,29,20,19,12,11,24,23,10,9];
% %
% % % first, habituation
% % for l = 1:length(rat(r).day(d).trial(t).habituation.lap)
% %  for e = 1%only first
% %    if(isfield(rat(r).day(d).trial(t).habituation.lap(l).exploration(e), 'probe'))
% %      %copy the original lfp for the hippocampus probe
% %      fprintf('l = %d; e  = %d\n', l, e)
% %      orighiplfp = rat(r).day(d).trial(t).habituation.lap(l).exploration(e).probe(p).lfp;
% %      tmplfp = NaN(size(orighiplfp));
% %       %Hippocampus
% %       for ch = 1:32
% %        tmplfp(ch,:) = orighiplfp(tmpch(ch),:);
% %       end%end
% %       %copy it back
% %       rat(r).day(d).trial(t).habituation.lap(l).exploration(e).probe(p).lfp = tmplfp;
% %
% %     end%end if isfield
% %   end%end of exploration
% % end%lap
% % % then, dishabituation
% % for l = 1:length(rat(r).day(d).trial(t).dishabituation.lap)
% %  for e = 1%only first
% %    if(isfield(rat(r).day(d).trial(t).dishabituation.lap(l).exploration(e), 'probe'))
% %      %copy the original lfp for the hippocampus probe
% %      orighiplfp = rat(r).day(d).trial(t).dishabituation.lap(l).exploration(e).probe(p).lfp;
% %      tmplfp = NaN(size(orighilfp));
% %       %Hippocampus
% %       for ch = 1:32
% %        tmplfp(ch,:) = orighiplfp(tmpch(ch),:);
% %       end%end
% %       %copy it back
% %       rat(r).day(d).trial(t).dishabituation.lap(l).exploration(e).probe(p).lfp = tmplfp;
% %
% %     end%end if isfield
% %   end%end of exploration
% % end%lap
% %
% %
% %
% % cd('..')%go back to the top folder level
% %
% % filedescript = sprintf('erics_data_nis%d', numrats);
% % minstring = sprintf('min%dms',minexplorationtime*1000);%convert to milliseconds since we don't want decimals in file names
% % secstring = sprintf('%dsecs', secbefore);
% % sampstring = sprintf('%dHz', samprate);
% % todaystring = date;%e.g.,  '11-Oct-2020'
% %
% % fn = sprintf('%s_%s_%s_%s_%s.mat', filedescript, minstring, secstring, sampstring, todaystring);
% %
% % %save the data
% % save(fn,'rat', '-v7.3')
% %
% %
% % end%end of function

