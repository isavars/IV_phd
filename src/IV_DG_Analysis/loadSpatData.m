function [cellID, meanRate, burstIndex, env, rMap, SI_spat, SpkTs, waveforms, nSpks] = loadSpatData()
%loadSpatData Loads all the data from a chosen getSpatData .mat file and can be used by all
%my figure making functions. keep working on this - for now it only does
%the first output 

 load ('allRats_spatData_CDB.mat', 'spatData')

% obtaining variables from spatData Table 
    
        meanRate = spatData.meanRate;
        burstIndex = spatData.burstIndex;
        env = spatData.env;
        rMap = spatData.rMap;
        %cellNo = spatData.cellNo;
        SI_spat = spatData.SI_spat;
        cellID = spatData.cellID;
        SpkTs = spatData.SpkTs;
        waveforms = spatData.waveforms;
        nSpks = spatData.nSpks;


        meanRate(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
        burstIndex(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
        SI_spat(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
     
    if nargin == 0
        meanRate;
    
    else
        cellID;
        
    end
end

