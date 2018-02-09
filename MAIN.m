function MAIN()
%% This is a caller routine to evaluate the user of the parfor on our server 
% the purpose is to test parfor performance with a large number of
% surrogates on our server. 
load ExtractHGHFOOpenField.mat; 
% set paramaters: 
params.PhaseFreqVector      = 2:2:50;
params.AmpFreqVector        = 100:5:200;
params.PhaseFreq_BandWidth  = 4;
params.AmpFreq_BandWidth    = 10;
params.computesurr          = 1;
params.numsurrogate         = 100;
params.alphause             = 0.05;
params.plotdata             = 1;
params.useparfor            = 1; % if true, user parfor, requires parallel computing toolbox

computePAC(lfpHFO,srate,params);
end