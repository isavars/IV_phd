function populationData()
%makes summary plots for clustered data in bar charts. 
%1. add adaptable indexes
% 2. this one is not working right now

load ('allRats_DGCA3_spatData.mat', 'spatData')

%gather data from spatData

    meanRate = spatData.meanRate;
    burstIndex = spatData.burstIndex;
    env = spatData.env; 
    SI_spat = spatData.SI_spat;
    
    meanRate(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    burstIndex(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    SI_spat(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    
%gather age data from cellInfo 

    cellInfo = getCellInfo();
    age = cellInfo(:,3);
    
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
     FamInd1 = FamInd(:,1); %add these next time to simply a bit?
     FamInd2 = FamInd(:,2);    
     
%make parent figure with axes to be filled in
    
    hFig = gra_multiplot(4, 5, 'figborder', [2 1 1 1]);
    axArr = getappdata(hFig, 'axesHandles' ); % makes the axes 
    axColCount= 1;
    
%   age binning 

    ageBins   = [18 19;20 21; 22 23];  % list of age bins each spanning from col1:col2
    
    for itAge=1:size(ageBins,1)
         
        indAge = age>=ageBins(itAge,1) & age<=ageBins(itAge,2); % index for current age bin
        meanRate = meanRate(indAge,:); 
        burstIndex = burstIndex(indAge,:);
        
        %make clusters 
            meanMeanRate = [];
            meanBurstIndex = [];
            for it_gm = 1: length (meanRate)
                meanMeanRate(it_gm) = mean(meanRate(it_gm,5));
                meanBurstIndex(it_gm) = mean(burstIndex(it_gm,5));
            end

            cluster1 = [];
            cluster2 = [];
            cluster3 = [];
            cluster4 = [];

            for jj = 1: length(meanMeanRate)
                if  (meanMeanRate(jj) > 1) && (meanMeanRate(jj) < 10) && (meanBurstIndex(jj) < 0.05) %mystery cells
                    cluster1(jj) = jj;
                elseif  (meanMeanRate(jj) > 0.01) && (meanMeanRate(jj) < 0.99)  && (meanBurstIndex(jj) > 0.02) && (meanBurstIndex(jj) < 0.15) %potential granule cells
                    cluster2(jj) = jj;
                elseif (meanMeanRate(jj) > 0.02) && (meanBurstIndex(jj) > 0.15) %potential mossy cells
                    cluster3(jj) = jj;
                else 
                    cluster4(jj) = jj;
                end    
            end 
            cluster1 = nonzeros(cluster1);
            cluster2 = nonzeros(cluster2);
            cluster3 = nonzeros(cluster3);
            cluster4 = nonzeros(cluster4);

        %make rate change score

            for it_gma = 1: length (meanRate)
                awakeMeanRate(it_gma,:) = mean(meanRate(it_gma,(FamInd(1):FamInd(2))));
            end

            sleepMeanRate = meanRate (:,5);

            rateChange = ((awakeMeanRate-sleepMeanRate)./(awakeMeanRate+sleepMeanRate));

        % make spatiality score

            for it_SI = 1: length (meanRate)
                meanSI(it_SI,:) = mean(SI_spat(it_SI,rmmissing(FamInd1(it_SI):FamInd2(it_SI))));
            end
            SI_spat = meanSI;

        %making bar graphs for everything    

            MeanRatec1 = mean(awakeMeanRate(cluster1));
            MeanRatec2 = mean(awakeMeanRate(cluster2));
            MeanRatec3 = mean(awakeMeanRate(cluster3));
            MeanRatec4 = mean(awakeMeanRate(cluster4));

            x = categorical({'C1','C2','C3','C4'});
            y = [MeanRatec1 MeanRatec2 MeanRatec3 MeanRatec4];

            bar(axArr(1,axColCount),x(1), y(1), 'r')
            hold (axArr(1,axColCount),'on')
            bar(axArr(1,axColCount),x(2), y(2), 'b')
            bar(axArr(1,axColCount),x(3), y(3), 'g')
            bar(axArr(1,axColCount),x(4), y(4), 'k')
            title(axArr(1,axColCount),strcat('P',num2str(ageBins(itAge,1)),' to P',num2str(ageBins(itAge,2))))
            hold off;
                ylabel(axArr(1,1),'Mean Rate (Hz)')

            c1 = mean(meanBurstIndex(cluster1));
            c2 = mean(meanBurstIndex(cluster2));
            c3 = mean(meanBurstIndex(cluster3));
            c4 = meanBurstIndex(cluster4);
            c4 = rmmissing(c4);
            c4= mean (c4);

            x = categorical({'C1','C2','C3','C4'});
            y = [c1 c2 c3 c4];

            bar(axArr(2,axColCount),x(1), y(1), 'r')
            hold (axArr(2,axColCount),'on')
            bar(axArr(2,axColCount),x(2), y(2), 'b')
            bar(axArr(2,axColCount),x(3), y(3), 'g')
            bar(axArr(2,axColCount),x(4), y(4), 'k')
            ylabel(axArr(2,1),'Burst Index')
            hold off;

            c1 = mean(rateChange(cluster1));
            c2 = mean(rateChange(cluster2));
            c3 = mean(rateChange(cluster3));
            c4 = mean(rateChange(cluster4));

            x = categorical({'C1','C2','C3','C4'});
            y = [c1 c2 c3 c4];

            bar(axArr(3,axColCount),x(1), y(1), 'r')
            hold (axArr(3,axColCount),'on')
            bar(axArr(3,axColCount),x(2), y(2), 'b')
            bar(axArr(3,axColCount),x(3), y(3), 'g')
            bar(axArr(3,axColCount),x(4), y(4), 'k')
            ylabel(axArr(3,1),'State Dependent Rate Change')
            ylim(axArr(3,axColCount),[-0.2 1])
            hold off;

            c1 = mean(SI_spat(cluster1));
            c2 = mean(SI_spat(cluster2));
            c3 = SI_spat(cluster3);
            c3 = rmmissing(c3);
            c3= mean(c3);            
            c4 = SI_spat(cluster4);
            c4 = rmmissing(c4);
            c4= mean(c4);
            
            x = categorical({'C1','C2','C3','C4'});
            y = [c1 c2 c3 c4];

            bar(axArr(4,axColCount),x(1), y(1), 'r')
            hold (axArr(4,axColCount),'on')
            bar(axArr(4,axColCount),x(2), y(2), 'b')
            bar(axArr(4,axColCount),x(3), y(3), 'g')
            bar(axArr(4,axColCount),x(4), y(4), 'k')
            ylabel(axArr(4,1),'Spatiality')
            hold off;

            meanRate = spatData.meanRate;
            burstIndex = spatData.burstIndex;
            SI_spat = spatData.SI_spat;
            axColCount = axColCount + 1;            
    end 
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
