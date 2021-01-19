function [corner_points] = find_admiss(no_cell,par,gUL,gDL,gmm)
%find_admiss Find the admissible betas
% The function returns the admissible betas in a matrix (IxJxF), where 0
% entries show that the pair (i,j) on resource f does not fulfil the
% minimum SINR requirements
% INPUT
%  no_cell  -- Number of cells
%  par      -- Struct with predefined parameters
%  gUL      -- Matrix gain between users and BS in UL
%  gDL      -- Matrix gain between users and BS in DL
%  gmm      -- Matrix gain between each user
% OUTPUT
% beta_admiss-- Matrix with the beta_s
% chosen_CP-- The area and the chosen corner point
% SINR_min -- Minimum SINR of every pair
% SINR_UL  -- SINR of users in UL for every pair
% SINR_DL  -- SINR of users in DL for every pair

% The size of all matrices will be different if we are in a non-freq. sel.
% fading environment
usedFreq = par.nchunks*par.FreqFad + 1*~par.FreqFad;

%% Preallocation of matrices
% Matrix of the maximum beta_s
beta_matrix=zeros(par.lambdaul,par.lambdadl,usedFreq,no_cell);
% Matrix that indicate the chosen corner points
chosen_CP=zeros(par.lambdaul,par.lambdadl,usedFreq,no_cell);

% Matrix of the Powers in UL - Point N
p_UL_N=zeros(par.lambdaul,par.lambdadl,usedFreq,no_cell);
% Matrix of the Powers in UL - Point O
p_UL_O=zeros(par.lambdaul,par.lambdadl,usedFreq,no_cell);
% Matrix of the Powers in DL - Point L
p_DL_L=zeros(par.lambdaul,par.lambdadl,usedFreq,no_cell);
% Matrix of the Powers in DL - Point Q
p_DL_Q=zeros(par.lambdaul,par.lambdadl,usedFreq,no_cell);

% Matrix of the SINR in UL - Point L
SINR_UL_L=zeros(par.lambdaul,par.lambdadl,usedFreq,no_cell);
% Matrix of the SINR in UL - Point O
SINR_UL_O=zeros(par.lambdaul,par.lambdadl,usedFreq,no_cell);
% Matrix of the SINR in UL - Point M
SINR_UL_M=zeros(par.lambdaul,par.lambdadl,usedFreq,no_cell);

% Matrix of the SINR in DL - Point N
SINR_DL_N=zeros(par.lambdaul,par.lambdadl,usedFreq,no_cell);
% Matrix of the SINR in DL - Point Q
SINR_DL_Q=zeros(par.lambdaul,par.lambdadl,usedFreq,no_cell);
% Matrix of the SINR in DL - Point M
SINR_DL_M=zeros(par.lambdaul,par.lambdadl,usedFreq,no_cell);
%% For each cell
for idxCell=1:no_cell
    % Indices of UL and DL users within this cell
    idxUL = (idxCell - 1)*(par.lambdaul + par.lambdadl) + 1;
    idxDL = (idxCell - 1)*(par.lambdaul + par.lambdadl) + par.lambdaul + 1;
    idxDLfinal = idxCell*(par.lambdaul + par.lambdadl);
    % Vector of UL users
    usersUL = idxUL:idxDL-1;
    % Vector of DL users
    usersDL = idxDL:idxDLfinal;
    % Start filling the matrix beta
    for idxFreq = 1:usedFreq
        for idxUL = 1:par.lambdaul
            %% Evaluation of beta
            % Valeu of beta for UL
            beta_u = (gDL(usersDL,idxCell,idxFreq).*...
                (par.pmaxUL.*gUL(usersUL(idxUL),idxCell,idxFreq) - par.SINR_tgt_ul*par.noise))./...
                (par.SINR_tgt_ul*par.SINR_tgt_dl.*(par.pmaxUL.*gmm(usersUL(idxUL),usersDL,idxFreq)'+par.noise));
            % Value of beta for DL
            beta_d = (par.pmaxDL.*gDL(usersDL,idxCell,idxFreq).*gUL(usersUL(idxUL),idxCell,idxFreq)...
                - par.noise*par.SINR_tgt_dl.*(gmm(usersUL(idxUL),usersDL,idxFreq)'.*par.SINR_tgt_ul...
                + gUL(usersUL(idxUL),idxCell,idxFreq)))./(par.SINR_tgt_ul*par.SINR_tgt_dl*par.pmaxDL...
                .*gmm(usersUL(idxUL),usersDL,idxFreq)');
            % Min value of the maximum beta
            max_upper_beta = min([beta_u beta_d],[],2);
            % Fill the beta_matrix with the minimum beta upper bound
            beta_matrix(idxUL,:,idxFreq) = max_upper_beta';
            %% Evaluation of the powers in the corner points
            % Point N - UL Power P_{iN}^u
            p_UL_N(idxUL,:,idxFreq) = par.SINR_tgt_ul*(par.noise + par.pmaxDL*...
                par.beta)/gUL(usersUL(idxUL),idxCell,idxFreq);
            % Point O - UL Power P_{iO}^u
            p_UL_O(idxUL,:,idxFreq) = (par.pmaxDL.*gDL(usersDL,idxCell,idxFreq) -...
                par.SINR_tgt_dl*par.noise)./(gmm(usersUL(idxUL),usersDL,idxFreq)'.*par.SINR_tgt_dl);

            % Point L - DL Power P_{jL}^d
            p_DL_L(idxUL,:,idxFreq) = par.SINR_tgt_dl*(par.noise + par.pmaxUL.*...
                gmm(usersUL(idxUL),usersDL,idxFreq)')./(gDL(usersDL,idxCell,idxFreq));
            % Point Q - DL Power P_{jQ}^d
            p_DL_Q(idxUL,:,idxFreq) = (par.pmaxUL*gUL(usersUL(idxUL),idxCell,idxFreq) -...
                par.SINR_tgt_ul*par.noise)/(par.beta*par.SINR_tgt_ul);
            %% Evaluation of SINRs in the corner points
            % SINR UL - Point L
            SINR_UL_L(idxUL,:,idxFreq) = par.pmaxUL.*gUL(usersUL(idxUL),idxCell,idxFreq)./...
                (par.noise + p_DL_L(idxUL,:,idxFreq).*par.beta);
            % SINR UL - Point O
            SINR_UL_O(idxUL,:,idxFreq) = p_UL_O(idxUL,:,idxFreq).*gUL(usersUL(idxUL),idxCell,idxFreq)./...
                (par.noise + par.pmaxDL*par.beta);
            % SINR UL - Point M
            SINR_UL_M(idxUL,:,idxFreq) = par.pmaxUL.*gUL(usersUL(idxUL),idxCell,idxFreq)/...
                (par.noise + par.pmaxDL*par.beta);

            % SINR DL - Point N
            SINR_DL_N(idxUL,:,idxFreq) = par.pmaxDL.*gDL(usersDL,idxCell,idxFreq)./...
                (par.noise + p_UL_N(idxUL,:,idxFreq)'.*gmm(usersUL(idxUL),usersDL,idxFreq)');
            % SINR DL - Point Q
            SINR_DL_Q(idxUL,:,idxFreq) = p_DL_Q(idxUL,:,idxFreq)'.*gDL(usersDL,idxCell,idxFreq)./...
                (par.noise + par.pmaxDL.*gmm(usersUL(idxUL),usersDL,idxFreq)');
            % SINR DL - Point M
            SINR_DL_M(idxUL,:,idxFreq) = par.pmaxDL.*gDL(usersDL,idxCell,idxFreq)./...
                (par.noise + par.pmaxUL.*gmm(usersUL(idxUL),usersDL,idxFreq)');
        end
    end
end

% Set 1 if the pair is admissible and 0 if it is not - @lt stands for less
% than
beta_admiss = bsxfun(@lt,par.beta,beta_matrix);
% Separating in areas the corner points - Area 1
chosen_CP(SINR_UL_M > par.SINR_tgt_ul & SINR_DL_M >= par.SINR_tgt_dl) = 10;
% Separating in areas the corner points - Area 2
chosen_CP(SINR_UL_M > par.SINR_tgt_ul & SINR_DL_M < par.SINR_tgt_dl) = 20;
% chosen_CP_A2 = bsxfun(@eq,20,chosen_CP);
% chosen_CP_A2 = bsxfun(@times,beta_admiss,chosen_CP_A2);
% Separating in areas the corner points - Area 3
chosen_CP(SINR_UL_M <= par.SINR_tgt_ul & SINR_DL_M >= par.SINR_tgt_dl) = 30;
% chosen_CP_A3 = bsxfun(@eq,30,chosen_CP);
% chosen_CP_A3 = bsxfun(@times,beta_admiss,chosen_CP_A3);
% Selecting only the acceptable pairs - point to point
% product between the matrices beta_admiss and chosen_CP
chosen_CP = bsxfun(@times,beta_admiss,chosen_CP);

% Fill struct
corner_points.chosen_CP = chosen_CP;

corner_points.p_UL_N = p_UL_N;
corner_points.p_UL_O = p_UL_O;
corner_points.p_DL_L = p_DL_L;
corner_points.p_DL_Q = p_DL_Q;

corner_points.SINR_UL_L = SINR_UL_L;
corner_points.SINR_UL_O = SINR_UL_O;
corner_points.SINR_UL_M = SINR_UL_M;

corner_points.SINR_DL_N = SINR_DL_N;
corner_points.SINR_DL_Q = SINR_DL_Q;
corner_points.SINR_DL_M = SINR_DL_M;
% % Nuling the pairs that do not respect the upper bound on beta
% p_UL_N = bsxfun(@times,~~(chosen_CP_A1+chosen_CP_A2),p_UL_N);
% p_UL_O = bsxfun(@times,chosen_CP_A2,p_UL_O);
% 
% p_DL_L = bsxfun(@times,~~(chosen_CP_A1+chosen_CP_A3),p_DL_L);
% p_DL_Q = bsxfun(@times,chosen_CP_A3,p_DL_Q);
end

