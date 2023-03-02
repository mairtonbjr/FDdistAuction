%% Parameters
% lambdaul                        = [4 5 6 25];    % Total number of UL users
% lambdadl                        = [4 5 6 25];    % Total number of DL users
% nchunks                         = [4 5 6 25];    % Total number of freq. channels
% beta                            = 110;    	   % SI cancellation [-dB]
% seedMC                          = 1:100;  	   % Monte Carlo iteration

% Define input parameters
par.nchunks = 25;       % Number of freq. channels [F]
par.lambdaul = 25;      % Number of UL users
par.lambdadl = 25;      % Number of DL users
par.beta = db2lin(-110) % SI cancellation in linear scale

% Define important parameters
par.FreqFad = 1;        % =1 if we consider frequency selective channel
no_cell = 1;		    % Number of cells
par.SINR_tgt_ul = 1;    % SINR target for UL in linear scale = 0 dB
par.SINR_tgt_dl = 1;    % SINR target for DL in linear scale = 0 dB
par.pmaxUL = db2lin(24-30);  % Maximum power for UL
par.pmaxDL = db2lin(24-30);  % Maximum power for DL
par.noise = 10^-19.7*180e3; % Noise per frequency channel (180 kHz)
no_usr = par.lambdaul + par.lambdadl % Total number of users

% Define the alphaUL and alphaDL for path-loss compensation 
% UL
alphaUL = gUL.^-1;
alpha_norm_UL = norm(alphaUL,2);
alphaUL = alphaUL./alpha_norm_UL;
par.alphaUL = alphaUL;
% DL
alphaDL = gDL.^-1;
alpha_norm_DL = norm(alphaDL,2);
alphaDL = alphaDL./alpha_norm_DL;
par.alphaDL = alphaDL;

% Size of the channel matrices
% gUL has dimensions [par.lambdaul par.no_cell par.nchunks]
% gDL has dimensions [par.lambdadl par.no_cell par.nchunks]
% gmm has dimensions [par.lambdaul par.lambdadl par.nchunks]

% Get the corner points necessary for the power allocation
corner_points = find_admiss(no_cell,par,gUL,gDL,gmm);

% Auction solution to maximize the alpha sum-SE
[xUL,xDL,pUL,pDL,SINR_UL,SINR_DL,SINR_Min,SINR,auction_itr] = auction_approx_sol(no_cell,no_usr,par,gUL,gDL,corner_points);
