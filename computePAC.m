function computePAC(data,sr,varargin)
%% Compute Phase Amplitude coupling
%  Written by Roee Gilron roeegilron@gmail.com
%  Code based on Adriano Tort code described in this paper:
%  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2941206/
%
%  Inputs:
%
%  Required:
%
%  data - 1. matrix of size 1 x time points if computing PAC within one area
%         2. matrix of size 2 x time points if computing PAC between two
%         regions
%  sr   - sampling rate of data
%
%  optional arguments ('name', value) format:
%
%  1. PhaseFreqVector - matrix of phase frequencies.
%     example: computePAC(data,sr,'PhaseFreqVector',2:4:50)
%     default: 2:2:50
%  2. AmpFreqVector - matrix of amp frequencies.
%     example: computePAC(data,sr,'PhaseFreqVector',2:2:50,'AmpFreqVector',2:2:100)
%     default: 5:2:100
%  3. PhaseFreq_BandWidth - bandwith in which to compuate phase freq
%     example: computePAC(data,sr,'PhaseFreqVector',2:2:50,'PhaseFreq_BandWidth',4)
%     default: 4
%  4. AmpFreq_BandWidth - bandwith in which to compuate amp freq
%     example: computePAC(data,sr,'PhaseFreqVector',2:2:50,'AmpFreq_BandWidth',8)
%     default: 10
%  5. useparfor - use parfor to compute PAC (faster, but requires parallel
%     computing toolbox
%     example: computePAC(data,sr,'useparfor',1,'AmpFreq_BandWidth',8)
%     default: 0
%  6. computeSurrogates - compute surrogates for statistical purposes
%     example: computePAC(data,sr,'useparfor',0,'computeSurrogates',1)
%     default: 0
%  7. regionnames - charachter, regions names to be computed (only relevant if computing
%     between two regions)
%     example: computePAC(data,sr,'PhaseFreqVector',2:2:50,'regionnames',{'GPi','Motor Cortex'})
%     default: {'a1','a2'};
%
%  Also accpets structure format, example:
%  params.PhaseFreqVector = 2:2:50;
%  params.AmpFreqVector   = 5:1:100;
%  params.useparfor   = 1;
%  computePAC(data,sr,params);
%
%  Outputs:
%
%
%  plot of PAC in either 1 / 4 subplots:
%  if data is 1 x time points - compute PAC within
%
%  if data is 2 x time points, compute 4 subplots:
%  1 plot PAC within area 1
%  2 plot PAC within area 2
%  3 plot PAC between area 1 (phase) and area 2 (amp)
%  4 plot PAC between area 1 (amp) and area 2 (phase)


%% set params (Define the Amplitude- and Phase- Frequencies)
% Parse input arguments
p = inputParser;
p.CaseSensitive = false;

validationFcn = @(x) validateattributes(x,{'double'},{'nonempty'});

paramName = 'PhaseFreqVector';
defaultVal = 2:2:50;
addParameter(p,paramName,defaultVal,validationFcn)

paramName = 'AmpFreqVector';
defaultVal = 5:2:100;
addParameter(p,paramName,defaultVal,validationFcn)

paramName = 'PhaseFreq_BandWidth';
defaultVal = 4;
addParameter(p,paramName,defaultVal,validationFcn)

paramName = 'AmpFreq_BandWidth';
defaultVal = 10;
addParameter(p,paramName,defaultVal,validationFcn)

paramName = 'useparfor';
defaultVal = 0;
addParameter(p,paramName,defaultVal,validationFcn)

paramName = 'computeSurrogates';
defaultVal = 0;
addParameter(p,paramName,defaultVal,validationFcn)

validationFcn = @(x) validateattributes(x,{'cell'},{'nonempty'});
paramName = 'regionnames';
defaultVal = {'a1','a2'};
addParameter(p,paramName,defaultVal,validationFcn)

p.parse(varargin{:});


%%
%Extract values from the inputParser

PhaseFreqVector      = p.Results.PhaseFreqVector;
AmpFreqVector        = p.Results.AmpFreqVector;
PhaseFreq_BandWidth  = p.Results.PhaseFreq_BandWidth;
AmpFreq_BandWidth    = p.Results.AmpFreq_BandWidth;
useparfor            = p.Results.useparfor; % if true, user parfor, requires parallel computing toolbox
regionames           = p.Results.regionnames;
computesurr          = p.Results.computeSurrogates;

%% Load data
lfp           = data;
data_length   = length(lfp);
data_size     = size(data,1);
srate         = sr;
dt            = 1/srate;
t             = (1:data_length)*dt;
hfig          = figure('Visible','off');
hfig.Position = [1000         666        1132         672];

if data_size == 1
    numplots = 1; 
else
    numplots = 4; 
end
for aa = 1:numplots
    if numplots == 1 
        datAmp = lfp; 
        datPha = lfp; 
        ttlAmp = '';
        ttlPha = '';
    elseif numplots == 4 
        switch aa 
            case 1 
                datAmp = lfp(1,:);
                datPha = lfp(1,:);
                ttlAmp = sprintf('%s',regionames{1});
                ttlPha = sprintf('%s',regionames{1});
            case 2 
                datAmp = lfp(2,:);
                datPha = lfp(2,:);
                ttlAmp = sprintf('%s',regionames{2});
                ttlPha = sprintf('%s',regionames{2});
            case 3 
                datAmp = lfp(1,:);
                datPha = lfp(2,:);
                ttlAmp = sprintf('%s',regionames{1});
                ttlPha = sprintf('%s',regionames{2});
            case 4
                datAmp = lfp(2,:);
                datPha = lfp(1,:);
                ttlAmp = sprintf('%s',regionames{2});
                ttlPha = sprintf('%s',regionames{1});
        end
    end
    
    %% Do filtering and Hilbert transform on CPU
    Comodulogram=single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));
    AmpFreqTransformed = zeros(length(AmpFreqVector), data_length);
    PhaseFreqTransformed = zeros(length(PhaseFreqVector), data_length);
    
    for ii=1:length(AmpFreqVector)
        Af1 = AmpFreqVector(ii);
        Af2 = Af1+AmpFreq_BandWidth;
        AmpFreq=eegfilt_pac(datAmp,srate,Af1,Af2); % just filtering
        AmpFreqTransformed(ii, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
    end
    
    for jj=1:length(PhaseFreqVector)
        Pf1 = PhaseFreqVector(jj);
        Pf2 = Pf1 + PhaseFreq_BandWidth;
        PhaseFreq=eegfilt_pac(datPha,srate,Pf1,Pf2); % this is just filtering
        PhaseFreqTransformed(jj, :) = angle(hilbert(PhaseFreq)); % this is getting the phase time series
    end
    
    %% Do comodulation calculation
    start   = tic;
    % precalcluate vars for comodulation calculation that only need to be calculated once
    nbin     = 18;
    position = zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
    winsize  = 2*pi/nbin;
    for j = 1:nbin
        position(j) = -pi+(j-1)*winsize;
    end
    winsize = 2*pi/nbin;
    lognbin = log(nbin);
    % using this indexing scheme allows for more efficient parfor
    pairuse = [];cnt = 1;
    for jj=1:length(AmpFreqVector)
        for ii=1:length(PhaseFreqVector)
            puse1(cnt) = ii;
            puse2(cnt) = jj;
            cnt = cnt + 1;
        end
    end
    
    % create linearlized Comodulogram
    Comodulogram = zeros(size(pairuse,1),1,'single');
    if useparfor
        parfor p = 1:size(puse1,2)
            Comodulogram(p) = ModIndex_v3(PhaseFreqTransformed(puse1(p), :), AmpFreqTransformed(puse2(p), :)', position,nbin,winsize,lognbin);
        end
    else
        for p = 1:size(puse1,2)
            Comodulogram(p) = ModIndex_v3(PhaseFreqTransformed(puse1(p), :), AmpFreqTransformed(puse2(p), :)', position,nbin,winsize,lognbin);
        end
    end
    Coreshaped = reshape(Comodulogram,length(PhaseFreqVector),length(AmpFreqVector));
    fprintf('comod calc done in %f secs \n',toc(start));
    %% plotting
    if numplots ~=1 
        subplot(2,2,aa)
        if aa <= 2 
            ttlgrp = 'PAC within';
        else
            ttlgrp = 'PAC between';
        end
    else
        ttlgrp = 'PAC within';
    end
    contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,Coreshaped',30,'lines','none')
    set(gca,'fontsize',14)
    ttly = sprintf('Amplitude Frequency %s (Hz)',ttlAmp);
    ylabel(ttly)
    ttlx = sprintf('Phase Frequency %s (Hz)',ttlPha);
    xlabel(ttlx)
    title(ttlgrp);
    colorbar
end
hfig.Visible = 'on';


end