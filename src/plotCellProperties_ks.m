function plotCellProperties_ks()
    spikeTimes = getCellProperties_ks();

    for x = 1:length(cluLabel)
        [xc, lags]=spk_crosscorr(spikeTimes(x),'AC',0.001,0.5,300,'plot',gca);
    end
end