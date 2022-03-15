function ageSorting()
%AGESORTING is like firingPropertieScatter but by age
%   In the future I just want to be able to call on it when I want any data
%   sorted by age. "loading data" and "clustering data" should happen in
%   two other functions also. 

%16/06/2020 working on adding sleep state change - fix when I have SWS only


load ('r889_all_wf_means.mat', 'spatData')

%gather data from spatData

    meanRate = spatData.meanRate;
    burstIndex = spatData.burstIndex;
    SI_spat = spatData.SI_spat;
    

%gather age data from cellInfo 

    cellInfo = getCellInfo();
    age = cellInfo(:,3);
    
%make parent figure with axes to be filled in
    
    hFig = gra_multiplot(3, 5, 'figborder', [2 1 1 1]);
    axArr = getappdata(hFig, 'axesHandles' ); % makes the axes 
    
%create indexes to sort data by age make into its own function! 
    
    ageCount = 18;
    counter =0;
    idx1 = 1;
    axColCount = 1;
 
    for itSD = 1: length (age)
        if age(itSD) == ageCount 
            counter = counter +1;
        elseif age(itSD) ~= ageCount 
            idx2 = idx1 + counter;
            meanRate = meanRate(idx1:idx2,:);
            burstIndex = burstIndex(idx1:idx2,:);
            SI_spat = SI_spat(idx1:idx2,:);
                % make spatiality score
                for it_SI = 1: length (meanRate)
                    meanSI(it_SI,:) = mean(SI_spat(it_SI,1:3));
                end
                SI_spat = meanSI;
                %make rate change score
    
                for it_gma = 1: length (meanRate)
                    awakeMeanRate(it_gma,:) = mean(meanRate(it_gma,1:3));
                end

                sleepMeanRate = meanRate (:,5);
                rateChange = ((awakeMeanRate-sleepMeanRate)./(awakeMeanRate+sleepMeanRate));
                
                % make clusters % [cluster1,cluster2,cluster3,cluster4]= makeClusters(meanRate, burstIndex);
                for it_gm = 1: length (meanRate)
                    meanMeanRate(it_gm) = mean(meanRate(it_gm,[1:3 5]));
                    meanBurstIndex(it_gm) = mean(burstIndex(it_gm,[1:3 5]));
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
                
                %any scatter plots should be made here
%                 figure()
                
                scatter(axArr(1,axColCount),meanMeanRate(cluster1), meanBurstIndex(cluster1), 'r') 
                hold (axArr(1,axColCount),'on')
                scatter(axArr(1,axColCount),meanMeanRate(cluster2), meanBurstIndex(cluster2),'b')
                scatter(axArr(1,axColCount),meanMeanRate(cluster3), meanBurstIndex(cluster3), 'g')
                scatter(axArr(1,axColCount),meanMeanRate(cluster4), meanBurstIndex(cluster4), 'k')
                hold off
                    title(axArr(1,axColCount),strcat('P',num2str(ageCount)))
%                     title(axArr(1,axColCount),strcat('Burst Index Vs Mean Rate: P',num2str(ageCount)))
                    xlabel(axArr(1,axColCount),'Mean Rate (seconds)')
                    set(axArr(1,axColCount), 'xscale', 'log')
                    ylabel(axArr(1,axColCount),'Burst Index')
                    axis (axArr(1,axColCount),[0.01 10 0 0.3])
%                     legend(axArr(1,axColCount),'Mystery Cells','Putative Granule','Putative Mossy', 'Unclassified');

                scatter(axArr(2,axColCount),SI_spat(cluster1), meanBurstIndex(cluster1), 'r') 
                hold (axArr(2,axColCount),'on')
                scatter(axArr(2,axColCount),SI_spat(cluster2), meanBurstIndex(cluster2),'b')
                scatter(axArr(2,axColCount),SI_spat(cluster3), meanBurstIndex(cluster3), 'g')
                scatter(axArr(2,axColCount),SI_spat(cluster4), meanBurstIndex(cluster4), 'k')
                hold off
%                     title(axArr(2,axColCount),strcat('Burst Index Vs Spatiality: P',num2str(ageCount)))
                    xlabel(axArr(2,axColCount),'Spatial Score')
                    ylabel(axArr(2,axColCount),'Burst Index')
                    axis (axArr(2,axColCount),[0 0.6 0 0.3])
%                     legend(axArr(1,axColCount),'Mystery Cells','Putative Granule','Putative Mossy', 'Unclassified');

                scatter(axArr(3,axColCount),meanMeanRate(cluster1), SI_spat(cluster1), 'r')  
                hold (axArr(3,axColCount),'on')
                scatter(axArr(3,axColCount),meanMeanRate(cluster2), SI_spat(cluster2),'b')
                scatter(axArr(3,axColCount),meanMeanRate(cluster3), SI_spat(cluster3), 'g')
                scatter(axArr(3,axColCount),meanMeanRate(cluster4), SI_spat(cluster4), 'k')
%                 title(axArr(3,axColCount),strcat('Spatiality Vs Mean Rate: P',num2str(ageCount)))
                    xlabel(axArr(3,axColCount),'Mean Rate (seconds)')
                    set(axArr(3,axColCount), 'xscale', 'log')
                    ylabel(axArr(3,axColCount),'Spatial Score')
                    axis (axArr(3,axColCount),[0.01 10 0 0.6])
%                     legend(axArr(3,axColCount),'Mystery Cells','Putative Granule','Putative Mossy', 'Unclassified');
                hold off
                
%                 scatter(axArr(4,axColCount),meanMeanRate(cluster1), rateChange(cluster1), 'r')  
%                 hold (axArr(4,axColCount),'on')
%                 scatter(axArr(4,axColCount),meanMeanRate(cluster2), rateChange(cluster2),'b')
%                 scatter(axArr(4,axColCount),meanMeanRate(cluster3), rateChange(cluster3), 'g')
%                 scatter(axArr(4,axColCount),meanMeanRate(cluster4), rateChange(cluster4), 'k')
% %                 title(axArr(4,axColCount),strcat('Rate Change Vs Mean Rate: P',num2str(ageCount)))
%                     xlabel(axArr(4,axColCount),'Mean Rate (seconds)')
%                     set(axArr(4,axColCount), 'xscale', 'log')
%                     ylabel(axArr(4,axColCount),'State Dependent Rate Change')
%                     axis (axArr(4,axColCount),[0.01 10 -1 1])
% %                     legend(axArr(4,axColCount),'Mystery Cells','Putative Granule','Putative Mossy', 'Unclassified');
%                 hold off
                
            idx1 = idx2;
            idx2= 0;
            ageCount = ageCount +1;
            counter =0;
            meanRate = spatData.meanRate;
            burstIndex = spatData.burstIndex;
            SI_spat = spatData.SI_spat;
            axColCount = axColCount + 1;
        end  
    end 
    meanRate = meanRate(idx1:end,:);
    burstIndex = burstIndex(idx1:end,:);
    SI_spat = SI_spat(idx1:end,:);
    
    % make spatiality score
        meanSI = [];
        for it_SI = 1: length (meanRate)
            meanSI(it_SI,:) = mean(SI_spat(it_SI,1:3));
        end
        SI_spat = transpose(meanSI);
        
    %make other features
        meanMeanRate = [];
        meanBurstIndex = [];
        for it_gm = 1: length (meanRate)
            meanMeanRate(it_gm) = mean(meanRate(it_gm,[1:3 5]));
            meanBurstIndex(it_gm) = mean(burstIndex(it_gm,[1:3 5]));
        end
        meanMeanRate = transpose(meanMeanRate);
        meanBurstIndex = transpose(meanBurstIndex);
        
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
 
        scatter(axArr(1,axColCount),meanMeanRate(cluster1), meanBurstIndex(cluster1), 'r')  
        hold (axArr(1,axColCount),'on')
        scatter(axArr(1,axColCount),meanMeanRate(cluster2), meanBurstIndex(cluster2),'b')
        scatter(axArr(1,axColCount),meanMeanRate(cluster3), meanBurstIndex(cluster3), 'g')
        scatter(axArr(1,axColCount),meanMeanRate(cluster4), meanBurstIndex(cluster4), 'k')
        title(axArr(1,axColCount),strcat('P',num2str(ageCount)))        
%         title(axArr(1,axColCount),strcat('Burst Index Vs Mean Rate: P',num2str(ageCount)))
            xlabel(axArr(1,axColCount),'Mean Rate (seconds)')
            set(axArr(1,axColCount), 'xscale', 'log')
            ylabel(axArr(1,axColCount),'Burst Index')
            axis (axArr(1,axColCount),[0.01 10 0 0.3])
%             legend(axArr(1,axColCount),'Mystery Cells','Putative Granule','Putative Mossy', 'Unclassified');
        hold off
        
        scatter(axArr(2,axColCount),SI_spat(cluster1), meanBurstIndex(cluster1), 'r')  
        hold (axArr(2,axColCount),'on')
        scatter(axArr(2,axColCount),SI_spat(cluster2), meanBurstIndex(cluster2),'b')
        scatter(axArr(2,axColCount),SI_spat(cluster3), meanBurstIndex(cluster3), 'g')
        scatter(axArr(2,axColCount),SI_spat(cluster4), meanBurstIndex(cluster4), 'k')
%         title(axArr(2,axColCount),strcat('Burst Index Vs Spatiality: P',num2str(ageCount)))
            xlabel(axArr(2,axColCount),'Spatial Score')
            ylabel(axArr(2,axColCount),'Burst Index')
            axis (axArr(2,axColCount),[0 0.6 0 0.3])
%             legend(axArr(2,axColCount),'Mystery Cells','Putative Granule','Putative Mossy', 'Unclassified');
        hold off

        scatter(axArr(3,axColCount),meanMeanRate(cluster1), SI_spat(cluster1), 'r')  
        hold (axArr(3,axColCount),'on')
        scatter(axArr(3,axColCount),meanMeanRate(cluster2), SI_spat(cluster2),'b')
        scatter(axArr(3,axColCount),meanMeanRate(cluster3), SI_spat(cluster3), 'g')
        scatter(axArr(3,axColCount),meanMeanRate(cluster4), SI_spat(cluster4), 'k')
%         title(axArr(3,axColCount),strcat('Spatiality Vs Mean Rate: P',num2str(ageCount)))
            xlabel(axArr(3,axColCount),'Mean Rate (seconds)')
            set(axArr(3,axColCount), 'xscale', 'log')
            ylabel(axArr(3,axColCount),'Spatial Score')
            axis (axArr(3,axColCount),[0.01 10 0 0.6])
            legend(axArr(3,axColCount),'Mystery Cells','Putative Granule','Putative Mossy', 'Unclassified');
        hold off
end

