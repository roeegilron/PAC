function plot_pac_from_results(results)
hfig = figure;
for aa = 1:length(results)
    if length(results) ~=1
        subplot(2,2,aa)
        if aa <= 2
            ttlgrp = 'PAC within';
        else
            ttlgrp = 'PAC between';
        end
    else
        ttlgrp = 'PAC within';
    end
    
    Com_reshaped = results(aa).Comodulogram;
    AmpFreq_BandWidth = results(aa).AmpFreq_BandWidth;
    AmpFreqVector = results(aa).AmpFreqVector;
    PhaseFreq_BandWidth  = results(aa).PhaseFreq_BandWidth;
    PhaseFreqVector  = results(aa).PhaseFreqVector;
    
    contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,Com_reshaped',30,'lines','none')
    set(gca,'fontsize',14)
    ttly = sprintf('Amplitude Frequency %s (Hz)','GPi');
    ylabel(ttly)
    ttlx = sprintf('Phase Frequency %s (Hz)','M1');
    xlabel(ttlx)
    title(ttlgrp);
end
end