function rateOverlap()
%RATEOVERLAP is based on the Leutgeb et al., 2007 analysis: 
%   1. Divide (for each cell) the mean firing rate in the less active
%   enclosure by the mean firing rate in the more active enclosure.
%   2. Average the rates of the active cells for the cell populaiton data.
%   The resulting overlap scores are 1 if the cell populaiton is compared
%   with itself and 0 if a completely difrerent cell population is active.
%   3. should probably do the same stats as in the Leutgeb paper or
%   laurenz's number nSpks in each environment to be considered.  

% load data (i want this to be its own function)

load ('allRats_spatData_CDB.mat', 'spatData')

    meanRate = spatData.meanRate;
    burstIndex = spatData.burstIndex;
    env = spatData.env; 
    SI_spat = spatData.SI_spat;
    nSpks = spatData.nSpks;
    
    meanRate(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    burstIndex(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    SI_spat(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    nSpks(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;

%gather age data from cellInfo 

    cellInfo = getCellInfo();
    age = cellInfo(:,3);
    rat = cellInfo(:,1); 
    
%make indexes for environment to use in comparisons with any file 
    maxFam = 2;
    FamInd = nan(length(meanRate),maxFam); 
    FamIndT = [];     
    NovInd = [];
    famCount = 0;

     for itCell= 1: length(meanRate)
        for itTrial = 1: 5
            if contains(cast(spatData.env(itCell,itTrial),'char'),'fam')
                FamIndT(itCell,itTrial) = itTrial;
                famCount = famCount + 1;
                if famCount <= maxFam
                     FamInd(itCell,famCount)= transpose(nonzeros(FamIndT(itCell,itTrial)));
                end
            elseif strcmp(cast(spatData.env(itCell,itTrial),'char'),'nov')
                NovInd(itCell,itTrial) = itTrial;
                NovInd = nonzeros(NovInd);
            end 
        end
        famCount = 0;
     end
%          FamInd1 = FamInd(:,1) %add these next time to simplify a bit?
%          FamInd2 = FamInd(:,2)

      
%make rate change score
    
    for it_gma = 1: length (meanRate)
        awakeMeanRate(it_gma,:) = mean(meanRate(it_gma,(FamInd(1):FamInd(2))));
    end

    sleepMeanRate = meanRate (:,5);
    
    rateChange = ((awakeMeanRate-sleepMeanRate)./(awakeMeanRate+sleepMeanRate));

% animal binning 
   ratBins = [804 889];%this isnt actually working 
   for itRat = 1:size(ratBins,1)
       indRat = rat>=ratBins(itRat,1) & rat<=ratBins(itRat,2);
       meanRate = meanRate(indRat,:); 
       burstIndex = burstIndex(indRat,:);
       
Age =[];
rateOverlaps =[];
Environment = [];
%   age binning 

    ageBins   =  [18 19;20 21; 22 23; 201206 210226];  % list of age bins each spanning from col1:col2;
    
    for itAge=1:size(ageBins,1)
        
        indAge = age>=ageBins(itAge,1) & age<=ageBins(itAge,2); % index for current age bin
        meanRate = meanRate(indAge,:); 
        burstIndex = burstIndex(indAge,:);

        %gather age data from cellInfo 

        cellInfo = getCellInfo();
        corecorded = cellInfo(:,4);
    
        % make clusters (i want this to be its own function)
        meanMeanRate = [];
        meanBurstIndex = [];
        for it_gm = 1: length (meanRate)
            if nSpks(it_gm,FamInd(it_gm,1))  >= 10 || nSpks(it_gm,FamInd(it_gm,2))  >= 10 || nSpks(it_gm,NovInd(it_gm))  >= 10  % nSpks(it_gm,FamInd(1))  >= 1 || nSpks(it_gm,FamInd(2))  >= 1
                meanMeanRate(it_gm) = mean(meanRate(it_gm,5));% indexing here needs to generalize to all datasets 
                meanBurstIndex(it_gm) = mean(burstIndex(it_gm,5));
            end
        end
        
        cluster1 = [];
        cluster2 = [];
        cluster3 = [];
        cluster4 = [];

        for jj = 1: length(meanMeanRate)
            if  (meanMeanRate(jj) >= 10)   %interneurons (meanMeanRate(jj) >= 1) && (meanMeanRate(jj) <= 10) && (meanBurstIndex(jj) <= 0.05)
                cluster1 = [cluster1; jj];
            elseif  (rateChange(jj) <= 0.3) && (meanMeanRate(jj) >= 0.01) && (meanMeanRate(jj) <= 1.5)  && (meanBurstIndex(jj) >= 0.02) && (meanBurstIndex(jj) < 0.15) %potential granule cells
                cluster2 = [cluster2; jj];
            elseif (meanMeanRate(jj) <= 1.5) && (meanBurstIndex(jj) >= 0.15) && (corecorded(jj) > 3) %potential CA3 cells
                cluster4 = [cluster4; jj];
            elseif (meanMeanRate(jj) >= 0.02) && (meanBurstIndex(jj) >= 0.15) %potential mossy cells
                cluster3 = [cluster3; jj];

%             if  (meanMeanRate(jj) > 1) && (meanMeanRate(jj) < 10) && (meanBurstIndex(jj) < 0.05) %mystery cells
%                 cluster1 = [cluster1; jj];
%             elseif  (rateChange(jj) <= 0.3) && (meanMeanRate(jj) > 0.01) && (meanMeanRate(jj) < 1.5)  && (meanBurstIndex(jj) > 0.02) && (meanBurstIndex(jj) < 0.15) %potential granule cells
%                 cluster2 = [cluster2; jj];
%             elseif (meanMeanRate(jj) > 0.02) && (meanBurstIndex(jj) > 0.15) %potential mossy cells
%                 cluster3 = [cluster3; jj];
%             else 
%                 cluster4 = [cluster4; jj];
            end    
        end 
  
        famOverlap1 = [];
        novOverlap1 = [];
         for itC1 = 1: length(cluster1)              
                 if meanRate(cluster1(itC1),FamInd(cluster1(itC1),1)) > meanRate(cluster1(itC1),FamInd(cluster1(itC1),2))
                    famOverlap1(itC1) = meanRate (cluster1(itC1),FamInd(cluster1(itC1),2)) / meanRate(cluster1(itC1),FamInd(cluster1(itC1),1));
                 elseif meanRate(cluster1(itC1),FamInd(cluster1(itC1),1)) < meanRate(cluster1(itC1),FamInd(cluster1(itC1),2))
                    famOverlap1(itC1) = meanRate (cluster1(itC1),FamInd(cluster1(itC1),1)) / meanRate(cluster1(itC1),FamInd(cluster1(itC1),2));
                 end
                 if meanRate(cluster1(itC1),FamInd(cluster1(itC1),2)) > meanRate(cluster1(itC1),NovInd(cluster1(itC1)))
                    novOverlap1(itC1) = meanRate (cluster1(itC1),NovInd(cluster1(itC1))) / meanRate(cluster1(itC1),FamInd(cluster1(itC1),2));
                 elseif meanRate(cluster1(itC1),FamInd(cluster1(itC1),2)) < meanRate(cluster1(itC1),NovInd(cluster1(itC1)))
                    novOverlap1(itC1) = meanRate (cluster1(itC1),FamInd(cluster1(itC1),2)) / meanRate(cluster1(itC1),NovInd(cluster1(itC1)));
                 end
         end
         famOverlap1 = rmmissing(transpose(famOverlap1));
         novOverlap1 = rmmissing(transpose(novOverlap1));
         FAMOverlapC1 = mean(famOverlap1);
         NOVOverlapC1 = mean(novOverlap1);
         FAMErrC1 = std(famOverlap1)/sqrt(length(cluster1));
         NOVErrC1 = std(novOverlap1)/sqrt(length(cluster1));

        famOverlap2 = [];
        novOverlap2 = [];
         for itC1 = 1: length(cluster2)
%             if (nSpks(cluster2(itC1),1) || nSpks(cluster2(itC1),2) || nSpks(cluster2(itC1),3)) <= 1
%              end
             if meanRate(cluster2(itC1),FamInd(cluster2(itC1),1)) > meanRate(cluster2(itC1),FamInd(cluster2(itC1),2))
                famOverlap2(itC1) = meanRate (cluster2(itC1),FamInd(cluster2(itC1),2)) / meanRate(cluster2(itC1),FamInd(cluster2(itC1),1));
             elseif meanRate(cluster2(itC1),FamInd(cluster2(itC1),1)) < meanRate(cluster2(itC1),FamInd(cluster2(itC1),2))
                famOverlap2(itC1) = meanRate (cluster2(itC1),FamInd(cluster2(itC1),1)) / meanRate(cluster2(itC1),FamInd(cluster2(itC1),2));
             end
             if meanRate(cluster2(itC1),FamInd(cluster2(itC1),2)) > meanRate(cluster2(itC1),NovInd(cluster2(itC1)))
                novOverlap2(itC1) = meanRate (cluster2(itC1),NovInd(cluster2(itC1))) / meanRate(cluster2(itC1),FamInd(cluster2(itC1),2));
             elseif meanRate(cluster2(itC1),FamInd(cluster2(itC1),2)) < meanRate(cluster2(itC1),NovInd(cluster2(itC1)))
                novOverlap2(itC1) = meanRate (cluster2(itC1),FamInd(cluster2(itC1),2)) / meanRate(cluster2(itC1),NovInd(cluster2(itC1)));
             end
         end
         famOverlap2 = rmmissing(transpose(famOverlap2));
         novOverlap2 = rmmissing(transpose(novOverlap2));
         FAMOverlapC2 = mean(famOverlap2);
         NOVOverlapC2 = mean(novOverlap2); 

         famOverlap2 = atanh(famOverlap2);
         novOverlap2 = atanh(novOverlap2);
         
         FAMErrC2 = std(famOverlap2)/sqrt(length(cluster2));
         NOVErrC2 = std(novOverlap2)/sqrt(length(cluster2));
         
        famOverlap3 = [];
        novOverlap3 = [];
         for itC1 = 1: length(cluster3)
%              if (nSpks(cluster3(itC1),1) || nSpks(cluster3(itC1),2) || nSpks(cluster3(itC1),3)) <= 1
%              end
             if meanRate(cluster3(itC1),FamInd(cluster3(itC1),1)) > meanRate(cluster3(itC1),FamInd(cluster3(itC1),2))
                famOverlap3(itC1) = meanRate (cluster3(itC1),FamInd(cluster3(itC1),2)) / meanRate(cluster3(itC1),FamInd(cluster3(itC1),1));
             elseif meanRate(cluster3(itC1),FamInd(cluster3(itC1),1)) < meanRate(cluster3(itC1),FamInd(cluster3(itC1),2))
                famOverlap3(itC1) = meanRate (cluster3(itC1),FamInd(cluster3(itC1),1)) / meanRate(cluster3(itC1),FamInd(cluster3(itC1),2)); %this isn't going to work with all data some have nan in position 1 so i have to use 2vs 4 
             end
             if meanRate(cluster3(itC1),FamInd(cluster3(itC1),2)) > meanRate(cluster3(itC1),NovInd(cluster3(itC1)))
                novOverlap3(itC1) = meanRate (cluster3(itC1),NovInd(cluster3(itC1))) / meanRate(cluster3(itC1),FamInd(cluster3(itC1),2));
             elseif meanRate(cluster3(itC1),FamInd(cluster3(itC1),2)) < meanRate(cluster3(itC1),NovInd(cluster3(itC1)))
                novOverlap3(itC1) = meanRate (cluster3(itC1),FamInd(cluster3(itC1),2)) / meanRate(cluster3(itC1),NovInd(cluster3(itC1)));
             end
         end
         famOverlap3 = rmmissing(transpose(famOverlap3));
         novOverlap3 = rmmissing(transpose(novOverlap3));
         FAMOverlapC3 = mean(famOverlap3);
         NOVOverlapC3 = mean(novOverlap3);
         FAMErrC3 = std(famOverlap3)/sqrt(length(cluster3));
         NOVErrC3 = std(novOverlap3)/sqrt(length(cluster3));
         
         famOverlap3 = atanh(famOverlap3);
         novOverlap3 = atanh(novOverlap3);

        famOverlap4 = [];
        novOverlap4 = [];
         for itC1 = 1: length(cluster4)
             if meanRate(cluster4(itC1),FamInd(cluster4(itC1),1)) > meanRate(cluster4(itC1),FamInd(cluster4(itC1),2))
                famOverlap4(itC1) = meanRate (cluster4(itC1),FamInd(cluster4(itC1),2)) / meanRate(cluster4(itC1),FamInd(cluster4(itC1),1));
             elseif meanRate(cluster4(itC1),FamInd(cluster4(itC1),1)) < meanRate(cluster4(itC1),FamInd(cluster4(itC1),2))
                famOverlap4(itC1) = meanRate (cluster4(itC1),FamInd(cluster4(itC1),1)) / meanRate(cluster4(itC1),FamInd(cluster4(itC1),2));
             end
             if meanRate(cluster4(itC1),FamInd(cluster4(itC1),2)) > meanRate(cluster4(itC1),NovInd(cluster4(itC1)))
                novOverlap4(itC1) = meanRate (cluster4(itC1),NovInd(cluster4(itC1))) / meanRate(cluster4(itC1),FamInd(cluster4(itC1),2));
             elseif meanRate(cluster4(itC1),FamInd(cluster4(itC1),2)) < meanRate(cluster4(itC1),NovInd(cluster4(itC1)))
                novOverlap4(itC1) = meanRate (cluster4(itC1),FamInd(cluster4(itC1),2)) / meanRate(cluster4(itC1),NovInd(cluster4(itC1)));
             end
         end
         
         famOverlap4 = rmmissing(transpose(famOverlap4));
         novOverlap4 = rmmissing(transpose(novOverlap4));        
         FAMOverlapC4 = mean(famOverlap4);
         NOVOverlapC4 = mean(novOverlap4);
         FAMErrC4 = std(famOverlap4)/sqrt(length(cluster4));
         NOVErrC4 = std(novOverlap4)/sqrt(length(cluster4));
         
        CellCount1 = length(cluster1);
        CellCount2 = length(cluster2);
        CellCount3 = length(cluster3);
        CellCount4 = length(cluster4);

        figure()
        y = [FAMOverlapC1 NOVOverlapC1; FAMOverlapC2 NOVOverlapC2; FAMOverlapC3 NOVOverlapC3; FAMOverlapC4 NOVOverlapC4];
        errors = [FAMErrC1  NOVErrC1; FAMErrC2 NOVErrC2; FAMErrC3 NOVErrC3; FAMErrC4 NOVErrC4];
        xticks = ({strcat('C1(',num2str(CellCount1),')'),strcat('C2(',num2str(CellCount2),')'),strcat('C3(',num2str(CellCount3),')'),strcat('C4(',num2str(CellCount4),')')});
        errorBars = gra_groupedbars(y, errors);
        errorBars = set(gca,'xticklabels', xticks);
        hold on;
%         errorBars(1).FaceColor = 'flat';
%         errorBars(1).CData(1,:) = [.5 0 0];
%         errorBars(1).CData(2,:) = [0 0 .5];
%         errorBars(1).CData(3,:) = [0 .5 0];
%         errorBars(1).CData(4,:) = [.5 .5 .5];
%         errorBars(2).FaceColor = 'flat';
%         errorBars(2).CData(1,:) = [1 0 0];
%         errorBars(2).CData(2,:) = [0 0 1];
%         errorBars(2).CData(3,:) = [0 1 0];
%         errorBars(2).CData(4,:) = [0 0 0]; % change these on illustrator
%         if this gets too complicated. 
            title(strcat('P',num2str(ageBins(itAge,1)),' to P',num2str(ageBins(itAge,2)))); %r',num2str(ratBins(itRat,1)),': put this back in when you fix rat binning 
            ylabel('Rate Overlap Score')
            xlabel('Cluster(N)')
                    
%         hFig = gra_multiplot(1, 4, 'figborder', [2 1 1 1]);
%         axArr = getappdata(hFig, 'axesHandles' ); % makes the axes 
% 
%            g1 = categorical({'C1F','C1N'});
%            g2 = categorical({'C2F','C2N'});
%            g3 = categorical({'C3F','C3N'});
%            g4 = categorical({'C4F','C4N'});           
%            boxplot(axArr(1,1),[famOverlap1, novOverlap1], g1)
%            boxplot(axArr(1,2),[famOverlap2, novOverlap2], g2) 
%            boxplot(axArr(1,3),[famOverlap3, novOverlap3], g3)
%            boxplot(axArr(1,4),[famOverlap4, novOverlap4], g4)
%            ylabel(axArr(1,1), 'Rate Overlap Score')
%           title(axArr(1,1), strcat('P',num2str(ageBins(itAge,1)),' to P',num2str(ageBins(itAge,2)))); 

        Age = [Age;repmat(ageBins(itAge,1),length(famOverlap2),1)];% *2 removed for KW atempt 
        rateOverlaps = [rateOverlaps;famOverlap2 ,novOverlap2];%rateOverlaps =[rateOverlaps; novOverlap2 - famOverlap2];%rateOverlaps = [famOverlap2 ;novOverlap2];
%        Environment = [repmat('Fam',length(famOverlap2),1);repmat('Nov', length(novOverlap2),1)];
%        Environment = [Environment;repmat('Fam',length(famOverlap2),1);repmat('Nov', length(novOverlap2),1)];
        %Clusters = {'C2';'C2';'C3';'C3';'C4';'C4'};
        meanRate = spatData.meanRate;
        burstIndex = spatData.burstIndex;
        %kruskalwallis(rateOverlaps, Environment)
%        rateOverlaps
%        [p,h,stats] = ranksum(rateOverlaps, Environment);
    end  
   end
   
   %[p,tbl,stats]= kruskalwallis(rateOverlaps, Age);
   data_for_spss = [rateOverlaps, Age];
   granule_data_v1 = 'granule_data_v1.mat';
   save('granule_data_v1', 'data_for_spss');

%    [p,~,stats]= anovan(rateOverlaps, {Environment,Age});
%    [c,m,h] = multcompare(stats,'Estimate','row');

end       

