function populationData()
%makes summary plots for clustered data in bar charts. 


load ('allRats_spatData_CDB.mat', 'spatData')

%gather data from spatData

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
    corecorded = cellInfo(:,4);
    
%make indexes for environment to use in comparisons with any file 
    maxFam = 2;
    FamInd = nan(length(meanRate),maxFam); 
    FamIndT = [];     
    NovInd = [];
    famCount = 0;
    numTrials = 5;

     for itCell= 1: length(meanRate)
        for itTrial = 1: numTrials
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
     FamInd1 = FamInd(:,1); %add these next time to simply a bit?
     FamInd2 = FamInd(:,2);  
     
%make rate change score
    
    for it_gma = 1: length (meanRate)
        awakeMeanRate(it_gma,:) = mean(meanRate(it_gma,(FamInd1(it_gma):FamInd2(it_gma))));
    end
    sleepMeanRate = meanRate (:,5);
    rateChange = ((awakeMeanRate-sleepMeanRate)./(awakeMeanRate+sleepMeanRate));
          
%make parent figure with axes to be filled in
    
    hFig = gra_multiplot(3, 6, 'figborder', [2 1 1 1]);
    axArr = getappdata(hFig, 'axesHandles' ); % makes the axes 
    axColCount= 1;
    
    GCmeanRates = [];
    MCmeanRates = [];
    Age =[];
    CellType = [];
    firingRates = [];
    SI_scores = [];
%   age binning 

    ageBins   = [ 20 21; 22 23; 201206 210226];  %[18 23]; list of age bins each spanning from col1:col2
    
    for itAge=1:size(ageBins,1)
         
        indAge = age>=ageBins(itAge,1) & age<=ageBins(itAge,2); % index for current age bin
        meanRate = meanRate(indAge,:); 
        burstIndex = burstIndex(indAge,:);
        
%         %make rate change score
% 
%             for it_gma = 1: length (meanRate)
%                 awakeMeanRate(it_gma,:) = mean(meanRate(it_gma,(FamInd1(it_gma):FamInd2(it_gma))));
%             end
%             sleepMeanRate = meanRate (:,5);
%             rateChange = ((awakeMeanRate-sleepMeanRate)./(awakeMeanRate+sleepMeanRate));
        
        %make clusters 
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
            elseif (meanMeanRate(jj) <= 1.5) && (meanBurstIndex(jj) >= 0.15) && (corecorded(jj) > 5) %potential CA3 cells
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
  

        % make spatiality score

            for it_SI = 1: length (meanRate)
                meanSI(it_SI,:) = mean(SI_spat(it_SI,rmmissing(FamInd1(it_SI):FamInd2(it_SI))));
            end
            SI_spat = meanSI;

        %making bar graphs or boxplots for everything   
                           
            
            MeanRatec1 = mean(meanMeanRate(cluster1));
            MeanRatec2 = mean(meanMeanRate(cluster2));
            MeanRatec3 = mean(meanMeanRate(cluster3));
            MeanRatec4 = mean(meanMeanRate(cluster4));
% 
%             x = categorical({'C2','C3'});%'C1',,'C4'});
%             y = [MeanRatec2 MeanRatec3];%MeanRatec1  MeanRatec4

%             bar(axArr(1,axColCount),x(1), y(1), 'r')
%             bar(axArr(1,axColCount),x(1), y(1), 'y')
%             hold (axArr(1,axColCount),'on')
%             bar(axArr(1,axColCount),x(2), y(2), 'g')
% %             bar(axArr(1,axColCount),x(4), y(4), 'k')
%             title(axArr(1,axColCount),strcat('P',num2str(ageBins(itAge,1)),' to P',num2str(ageBins(itAge,2))))
%             hold off;
%                 ylabel(axArr(1,1),'Mean Rate (Hz)')

            group = categorical({'Granule cells', 'Mossy Cells'});
            GC = meanMeanRate(cluster2);% meanMeanRate(cluster3)];
            MC = meanMeanRate(cluster3);

            boxplot(axArr(1,axColCount),GC, group(1))
            boxplot(axArr(1,axColCount +1),MC, group(2))
%             title(axArr(1,axColCount),strcat('P',num2str(ageBins(itAge,1)),' to P',num2str(ageBins(itAge,2))))
                ylabel(axArr(1,1),'Mean Firing Rate (Hz)')
                set(axArr(1,axColCount), 'yscale', 'log')
                set(axArr(1,axColCount + 1), 'yscale', 'log')
                ylim(axArr(1,axColCount),[10^-2 10^1]);
                ylim(axArr(1,axColCount+1),[10^-2 10^1]);

            c1 = mean(meanBurstIndex(cluster1));
            c2 = mean(meanBurstIndex(cluster2));
            c3 = mean(meanBurstIndex(cluster3));
            c4 = meanBurstIndex(cluster4);
            c4 = rmmissing(c4);
            c4= mean (c4);
% 
%             x = categorical({'C2','C3'});%,'C1','C4'});%
%             y = [ c2 c3];%c1 c4];
% 
% %             bar(axArr(2,axColCount),x(1), y(1), 'r')
% 
%             bar(axArr(2,axColCount),x(1), y(1), 'y')
%             hold (axArr(2,axColCount),'on')
%             bar(axArr(2,axColCount),x(2), y(2), 'g')
% %             bar(axArr(2,axColCount),x(4), y(4), 'k')
%             ylabel(axArr(2,1),'Burst Index')
%             hold off;


            group = categorical({'Granule cells', 'Mossy Cells'});
            GC = meanBurstIndex(cluster2);
            MC = meanBurstIndex(cluster3);
            figure()
            boxplot(axArr(2,axColCount),GC, group(1))
            hold (axArr(2,axColCount),'on')
            boxplot(axArr(2,axColCount +1),MC, group(2))
%             title(axArr(2,axColCount),strcat('P',num2str(ageBins(itAge,1)),' to P',num2str(ageBins(itAge,2))))
            hold off;
                ylabel(axArr(2,1),'Burst Index')
                ylim(axArr(2,axColCount),[0 0.5]);
                ylim(axArr(2,axColCount+1),[0 0.5]);
                
%             c1= [];
%             for itRC = 1: length(cluster1)
%                 c1(itRC,:) = rateChange(cluster1(itRC));
%             end
%             c1= mean (c1);
            c2= [];
            for itRC = 1: length(cluster2)
                c2(itRC,:) = rateChange(cluster2(itRC));
            end
            %c2= mean (c2);
            c3= [];
            for itRC = 1: length(cluster3)
                c3(itRC,:) = rateChange(cluster3(itRC));
            end
            %c3= mean (c3);

%             x = categorical({'C2','C3'});%'C1',,'C4'});
%             y = [ c2 c3];%c1 c4];
% 
% %             bar(axArr(3,axColCount),x(1), y(1), 'r')
% 
%             bar(axArr(3,axColCount),x(1), y(1), 'y')
%             hold (axArr(3,axColCount),'on')
%             bar(axArr(3,axColCount),x(2), y(2), 'g')
% %             bar(axArr(3,axColCount),x(4), y(4), 'k')
%             ylabel(axArr(3,1),'State Dependent Rate Change')
%             ylim(axArr(3,axColCount),[-1 1])
%             hold off;
            group = categorical({'Granule cells', 'Mossy Cells'});
            GC = c2;
            MC = c3;

            boxplot(axArr(3,axColCount),GC, group(1))
            hold (axArr(3,axColCount),'on')
            boxplot(axArr(3,axColCount +1),MC, group(2))
%             title(axArr(3,axColCount),strcat('P',num2str(ageBins(itAge,1)),' to P',num2str(ageBins(itAge,2))))
            hold off;
                ylabel(axArr(3,1),'State Dependent Rate Change')
                ylim(axArr(3,axColCount),[-1 1]);
                ylim(axArr(3,axColCount+1),[-1 1]);

            
%             c1 = mean(SI_spat(cluster1));
%             c2 = mean(SI_spat(cluster2));
%             c3 = SI_spat(cluster3);
%             c3 = rmmissing(c3);
%             c3= mean(c3);            
%             c4 = SI_spat(cluster4);
%             c4 = rmmissing(c4);
%             c4= mean(c4);
%             
%             x = categorical({'C1','C2','C3','C4'});
%             y = [c1 c2 c3 c4];
% 
%             bar(axArr(4,axColCount),x(1), y(1), 'r')
%             hold (axArr(4,axColCount),'on')
%             bar(axArr(4,axColCount),x(2), y(2), 'b')
%             bar(axArr(4,axColCount),x(3), y(3), 'g')
%             bar(axArr(4,axColCount),x(4), y(4), 'k')
%             ylabel(axArr(4,1),'Spatiality')
%             hold off;
             GCmeanRates = [GCmeanRates; transpose(meanMeanRate(cluster2))];
             MCmeanRates = [MCmeanRates; transpose(meanMeanRate(cluster3))];
%              Age = [Age;repmat(ageBins(itAge,1),length(cluster3),1)];
             CellType = [CellType; repmat(1,length(cluster2),1);repmat(2,length(cluster3),1)];
             Age = [Age;repmat(ageBins(itAge,1),length(cluster2) +length(cluster3),1)];
%              firingRates = [firingRates;transpose(meanMeanRate(cluster2));transpose(meanMeanRate(cluster3))];
             SI_scores = [SI_scores;SI_spat(cluster2);SI_spat(cluster3)];

            
            meanRate = spatData.meanRate;
            burstIndex = spatData.burstIndex;
            SI_spat = spatData.SI_spat;
            axColCount = axColCount + 2;   


    end 
% % anova on men firing rate over age
% GCages = {'1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3'};
 [p,tbl,stats]=anovan(SI_scores,{Age,CellType},'model','interaction', 'varnames',{'Age','Cell Type'});
% [~,~,stats]=anova1(MCmeanRates,Age)
 [c,~,~,gnames] = multcompare(stats, 'Dimension', [1 2]);

end


 %making bar graphs for everything    

%             MeanRatec1 = mean(awakeMeanRate(cluster1));
%             MeanRatec2 = mean(awakeMeanRate(cluster2));
%             MeanRatec3 = mean(awakeMeanRate(cluster3));
%             MeanRatec4 = mean(awakeMeanRate(cluster4));
% 
%             x = categorical({'Cluster 1','Cluster 2','Cluster 3','Cluster 4'});
%             y = [MeanRatec1 MeanRatec2 MeanRatec3 MeanRatec4];
% 
%             bar(x(1), y(1), 'r')
%             hold all;
%             bar(x(2), y(2), 'g')
%             bar(x(3), y(3), 'b')
%             bar(x(4), y(4), 'k')
%             ylabel('Mean Rate (Hz)')
%             hold off;
% 
%            figure()%nan in cluster 4 for some reason 
%             c1 = mean(meanBurstIndex(cluster1));
%             c2 = mean(meanBurstIndex(cluster2));
%             c3 = mean(meanBurstIndex(cluster3));
%             c4 = mean(meanBurstIndex(cluster4));
% 
%             x = categorical({'Cluster 1','Cluster 2','Cluster 3','Cluster 4'});
%             y = [c1 c2 c3 c4];
% 
%             bar(x(1), y(1), 'r')
%             hold all;
%             bar(x(2), y(2), 'g')
%             bar(x(3), y(3), 'b')
%             bar(x(4), y(4), 'k')
%             ylabel('Burst Index')
%             hold off;
% 
%            figure() %if you average the rate changes they cna only be positive so figure out how to solve this
%             c1 = mean(rateChange(cluster1));
%             c2 = mean(rateChange(cluster2));
%             c3 = mean(rateChange(cluster3));
%             c4 = mean(rateChange(cluster4));
% 
%             x = categorical({'Cluster 1','Cluster 2','Cluster 3','Cluster 4'});
%             y = [c1 c2 c3 c4];
% 
%             bar(x(1), y(1), 'r')
%             hold all;
%             bar(x(2), y(2), 'g')
%             bar(x(3), y(3), 'b')
%             bar(x(4), y(4), 'k')
%             ylabel('State Dependent Rate Change')
%             hold off;
% 
%            figure()
%             c1 = mean(SI_spat(cluster1));
%             c2 = mean(SI_spat(cluster2));
%             c3 = mean(SI_spat(cluster3));
%             c4 = mean(SI_spat(cluster4));
% 
%             x = categorical({'Cluster 1','Cluster 2','Cluster 3','Cluster 4'});
%             y = [c1 c2 c3 c4];
% 
%             bar(x(1), y(1), 'r')
%             hold all;
%             bar(x(2), y(2), 'g')
%             bar(x(3), y(3), 'b')
%             bar(x(4), y(4), 'k')
%             ylabel('Spatiality')
%             hold off;
