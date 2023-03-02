%% Parameters
% lambdaul                        = [4 5 6 25];    % Total number of UL users
% lambdadl                        = [4 5 6 25];    % Total number of DL users
% nchunks                         = [4 5 6 25];    % Total number of freq. channels
% beta                            = 110;    	   % SI cancellation [-dB]
% seedMC                          = 1:100;  	   % Monte Carlo iteration

% Set the random number generator seed to 0 and to use Mersenne Twister
% generator
clear; clc;
rng(1,'twister');

% Define input parameters
par.nchunks = 25;       % Number of freq. channels [F]
par.lambdaul = 25;      % Number of UL users
par.lambdadl = 25;      % Number of DL users
par.beta = db2lin(-110); % SI cancellation in linear scale

% Define important parameters
par.FreqFad = 0;        % =1 if we consider frequency selective channel
no_cell = 1;		    % Number of cells
par.SINR_tgt_ul = 1;    % SINR target for UL in linear scale = 0 dB
par.SINR_tgt_dl = 1;    % SINR target for DL in linear scale = 0 dB
par.pmaxUL = db2lin(24-30);  % Maximum power for UL
par.pmaxDL = db2lin(24-30);  % Maximum power for DL
par.noise = 10^-19.7*180e3; % Noise per frequency channel (180 kHz)
no_usr = par.lambdaul + par.lambdadl; % Total number of users
% The size of all matrices will be different if we are in a non-freq. sel.
% fading environment
usedFreq = par.nchunks*par.FreqFad + 1*~par.FreqFad;

%% Size of the channel matrices
% The UL users are from 1:par.lambdaul 
% gUL has dimensions [no_usr, par.no_cell, par.nchunks]
% The DL users are from par.lambdaul +1:par.lambdaul + par.lambdadl
% gDL has dimensions [no_usr, par.no_cell, par.nchunks]
% Since all users need to have a UL and DL channel, these channels
% represent that. Similarly, the UE-to-UE channel is between all users
% regardless of UL and DL.
% gmm has dimensions [no_usr, no_usr, par.nchunks]

%% Define simple channels to create a running example. 
% You should generate and import here your own channel gUL, gDL, and gmm.
% Accordingly, comment the definitions below and import your channels here.
% If you simply want to test the function, use the random channels below  (not 3GPP compliant).
gUL = rand(no_usr, no_cell, 1);
gUL = repmat(gUL, 1, no_cell, par.nchunks);
gDL = rand(no_usr, no_cell, 1);
gDL = repmat(gDL, 1, no_cell, par.nchunks);
gmm = rand(no_usr, no_usr, 1);
gmm = repmat(gmm, 1, 1, par.nchunks);

% Define the alphaUL and alphaDL for path-loss compensation 
% In this case, they are defined specifically for UL and DL users. Hence,
% they have dimensions:
% UL - alphaUL = [par.lambdaul, par.no_cell]
alphaUL = gUL(1:par.lambdaul, :, 1).^-1;
alpha_norm_UL = norm(alphaUL,2);
alphaUL = alphaUL./alpha_norm_UL;
par.alphaUL = alphaUL;
% DL - alphaDL = [par.lambdadl, par.no_cell]
alphaDL = gDL(par.lambdaul+1:no_usr, :, 1).^-1;
alpha_norm_DL = norm(alphaDL,2);
alphaDL = alphaDL./alpha_norm_DL;
par.alphaDL = alphaDL;

% Get the corner points necessary for the power allocation
corner_points = find_admiss(no_cell,par,gUL,gDL,gmm);

% Auction solution to maximize the alpha sum-SE
[xUL,xDL,pUL,pDL,SINR_UL,SINR_DL,SINR_Min,SINR,auction_itr] = auction_approx_sol(no_cell,no_usr,par,gUL,gDL,corner_points);