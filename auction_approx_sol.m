function [xUL,xDL,pUL,pDL,SINR_UL,SINR_DL,SINR_Min,SINR,auction_itr] = auction_approx_sol(no_cell,no_usr,par,gUL,gDL,corner_points)
%auction_approx_sol Uses auction theory to select the best pairs
% The function returns the powers assigned to users and their SINRs
% INPUT
%  no_cell  -- Number of cells
%  par      -- Struct with predefined parameters
%  gUL      -- Matrix gain between users and BS in UL
%  gDL      -- Matrix gain between users and BS in DL
%  corner_points -- Struct with matrices related to the admissible areas
% OUTPUT
% xUL       -- Assignment matrix for UL
% xDL       -- Assignment matrix for DL
% SINR_UL  -- SINR of users in UL for every pair
% SINR_DL  -- SINR of users in DL for every pair
% SINR_Min -- Minimum SINR of the pair in a resource
% SINR     -- Total SINR perceived

% General vectors
% Assignment
xUL = zeros(no_usr,no_cell,par.nchunks);
xDL = zeros(no_usr,no_cell,par.nchunks);
% Power
pUL = zeros(no_usr,no_cell,par.nchunks);
pDL = zeros(no_usr,no_cell,par.nchunks);
% SINR
SINR_UL = zeros(no_usr,no_cell,par.nchunks);
SINR_DL = zeros(no_usr,no_cell,par.nchunks);
SINR_Min = zeros(par.nchunks,no_cell);

% For each cell
for idxCell = 1:no_cell
    
    % Vector of UL users - in the cell
    usersUL = 1:par.lambdaul;
    % Vector of DL users - in the cell
    usersDL = 1:par.lambdadl;
    % Vector of RBs
    freqAv = 1:par.nchunks;
    
    % Indices of UL and DL users within this cell
    idxUL = (idxCell - 1)*(par.lambdaul + par.lambdadl) + 1;
    idxDL = (idxCell - 1)*(par.lambdaul + par.lambdadl) + par.lambdaul + 1;
    idxDLfinal = idxCell*(par.lambdaul + par.lambdadl);
    % Vector of UL - in the system
    usersULFull = idxUL:idxDL-1;
    % Vector of DL - in the system
    usersDLFull = idxDL:idxDLfinal;
        
    %% Separate CP in the areas
    chosen_CP_A1 = bsxfun(@eq,10,corner_points.chosen_CP);
    chosen_CP_A2 = bsxfun(@eq,20,corner_points.chosen_CP);
    chosen_CP_A3 = bsxfun(@eq,30,corner_points.chosen_CP);
    % Chosen corner points
    chosen_CP = zeros(size(chosen_CP_A1));
    
    % Nuling the pairs that do not respect the upper bound on beta - Powers
    p_UL_N = bsxfun(@times,~~(chosen_CP_A1+chosen_CP_A2),corner_points.p_UL_N);
    p_UL_O = bsxfun(@times,chosen_CP_A2,corner_points.p_UL_O);
    
    p_DL_L = bsxfun(@times,~~(chosen_CP_A1+chosen_CP_A3),corner_points.p_DL_L);
    p_DL_Q = bsxfun(@times,chosen_CP_A3,corner_points.p_DL_Q);
    
    % Nuling the pairs that do not respect the upper bound on beta - SINR
    SINR_UL_L = bsxfun(@times,~~(chosen_CP_A1+chosen_CP_A3),corner_points.SINR_UL_L);
    SINR_UL_O = bsxfun(@times,chosen_CP_A2,corner_points.SINR_UL_O);
    SINR_UL_M = bsxfun(@times,chosen_CP_A1,corner_points.SINR_UL_M);
    
    SINR_DL_N = bsxfun(@times,~~(chosen_CP_A1+chosen_CP_A2),corner_points.SINR_DL_N);
    SINR_DL_Q = bsxfun(@times,chosen_CP_A3,corner_points.SINR_DL_Q);
    SINR_DL_M = bsxfun(@times,chosen_CP_A1,corner_points.SINR_DL_M);
    
    SINR_UL_N = ~~SINR_DL_N.*par.SINR_tgt_ul;
    SINR_UL_Q = ~~SINR_DL_Q.*par.SINR_tgt_ul;
    SINR_DL_L = ~~SINR_UL_L.*par.SINR_tgt_dl;
    SINR_DL_O = ~~SINR_UL_O.*par.SINR_tgt_dl;
    
    %% Generate matrix of alpha's
    alphaUL_mtx = repmat(par.alphaUL,1,size(SINR_UL_N,2));
    alphaDL_mtx = repmat(par.alphaDL,1,size(SINR_UL_N,1))';
    
    %% Select CP for admissible area 1 - Corner Points L, M and N
    % Evaluation of the alpha sum-SE of the pairs for admissible area 1
    SE_L_sum_A1 = alphaUL_mtx.*log2(1+SINR_UL_L.*chosen_CP_A1) + alphaDL_mtx.*log2(1+SINR_DL_L.*chosen_CP_A1);
    SE_M_sum = alphaUL_mtx.*log2(1+SINR_UL_M) + alphaDL_mtx.*log2(1+SINR_DL_M);
    SE_N_sum_A1 = alphaUL_mtx.*log2(1+SINR_UL_N.*chosen_CP_A1) + alphaDL_mtx.*log2(1+SINR_DL_N.*chosen_CP_A1);
     % Select the corner point that maximizes the alpha sum-SE
    CP_L_A1 = bsxfun(@and,bsxfun(@gt,SE_L_sum_A1,SE_M_sum),bsxfun(@gt,SE_L_sum_A1,SE_N_sum_A1));
    chosen_CP(CP_L_A1) = 11;
    CP_M = bsxfun(@and,bsxfun(@gt,SE_M_sum,SE_L_sum_A1),bsxfun(@gt,SE_M_sum,SE_N_sum_A1));
    chosen_CP(CP_M) = 12;
    CP_N_A1 = bsxfun(@and,bsxfun(@gt,SE_N_sum_A1,SE_L_sum_A1),bsxfun(@gt,SE_N_sum_A1,SE_M_sum));
    chosen_CP(CP_N_A1) = 13;
    % Evaluation of the alpha sum-SE of the pair for admissible area 1 
    SE_A1_sum = (chosen_CP == 11).*SE_L_sum_A1 + (chosen_CP == 13).*SE_N_sum_A1 + (chosen_CP == 12).*SE_M_sum;
    
    %% Select CP for admissible area 2 - Corner Points O and N
    % Evaluation of the alpha sum-SE of the pairs for admissible area 2
    SE_O_sum = alphaUL_mtx.*log2(1+SINR_UL_O) + alphaDL_mtx.*log2(1+SINR_DL_O);
    SE_N_sum_A2 = alphaUL_mtx.*log2(1+SINR_UL_N.*chosen_CP_A2) + alphaDL_mtx.*log2(1+SINR_DL_N.*chosen_CP_A2);
    % Select the corner point that maximizes the alpha sum-SE
    chosen_CP(bsxfun(@gt,SE_O_sum,SE_N_sum_A2)) = 21;
    chosen_CP(bsxfun(@gt,SE_N_sum_A2,SE_O_sum)) = 22;
    SE_A2_sum = (chosen_CP == 21).*SE_O_sum + (chosen_CP == 22).*SE_N_sum_A2;
    
    %% Select CP for admissible area 3 - Corner Points L and Q
    % Evaluation of the alpha sum-SE of the pairs for admissible area 3
    SE_L_sum_A3 = alphaUL_mtx.*log2(1+SINR_UL_L.*chosen_CP_A3) + alphaDL_mtx.*log2(1+SINR_DL_L.*chosen_CP_A3);
    SE_Q_sum = alphaUL_mtx.*log2(1+SINR_UL_Q) + alphaDL_mtx.*log2(1+SINR_DL_Q);
    SE_A3_sum = (chosen_CP == 31).*SE_L_sum_A3 + (chosen_CP == 32).*SE_Q_sum;
    % Select the corner point that maximizes the alpha sum-SE
    chosen_CP(bsxfun(@gt,SE_L_sum_A3,SE_Q_sum)) = 31;
    chosen_CP(bsxfun(@gt,SE_Q_sum,SE_L_sum_A3)) = 32;
    
    % Total SE
    % Sum
    SE_sum_total = SE_A1_sum + SE_A2_sum + SE_A3_sum;
    
    
    %% The auction starts here
    % Value of epsilon
    epsilon = 0.1;
    % Number of iterations
    auction_itr = 0;
    % Prices stored at the UL user and BS
    Prices_UE = zeros(size(SE_A1_sum));
    prices_BS = zeros(1,par.lambdadl);
    % Association of UL to DL users
    assign_ul_to_dl = zeros(1,par.lambdaul);
    % Flag to determine whether the UL users are assigned or not
    flag_assign = zeros(1,par.lambdaul);
    % Bid matrix
    bid_mtx = zeros(size(SE_A1_sum));
    
    % Continue until there are users not assigned
    while (sum(flag_assign)~=par.lambdaul)
        % For all UL users
        for idxUL = 1:par.lambdaul
            if flag_assign(idxUL) == 0
                % Evaluate v_i from Eq. 13 and j_i from Eq. 14
                [max_util_vi,idx_max_util_ji] = max(SE_sum_total(idxUL,:)-Prices_UE(idxUL,:)); %#ok<ASGLU>
                % Evaluate w_i from Eq. 15
                sec_max_util_wi = sort(SE_sum_total(idxUL,:)-Prices_UE(idxUL,:),'descend');
                sec_max_util_wi = sec_max_util_wi(2);
                % Evaluate the bids for this user
                bid_mtx(idxUL,idx_max_util_ji) = SE_sum_total(idxUL,idx_max_util_ji)- sec_max_util_wi + epsilon;
                % Selected DL user
                assign_ul_to_dl(idxUL) = idx_max_util_ji;
                % Update the number of iterations
                auction_itr = auction_itr + 1;
            end
        end
        
        % For the selected DL users
        selec_DL = unique(assign_ul_to_dl);
        for idxSelec = 1:length(selec_DL)
            % Get the users that bid on resource selec_DL(idxSelec)
            idxUL_bidders = find(assign_ul_to_dl == selec_DL(idxSelec));
            % If it is unique, then update the price of this DL user
            if length(idxUL_bidders) == 1
                % Update price if necessary
                if prices_BS(selec_DL(idxSelec)) ~= bid_mtx(idxUL_bidders,selec_DL(idxSelec))
                    prices_BS(selec_DL(idxSelec)) = bid_mtx(idxUL_bidders,selec_DL(idxSelec));
                    % Flag necessary to count the number of iterations
                    flag_update = 1;
                else
                    % If I do not update the prices, then there is no need
                    % to count this as an iteration
                    flag_update = 0;
                end
                % Assign this DL user to this UL user
                flag_assign(idxUL_bidders) = 1;
            else % There is competition for this DL user
                % Sort the bids in the descending order
                [sort_prices_BS,idx_max_bid] = sort(bid_mtx(idxUL_bidders,selec_DL(idxSelec)),'descend');
                % Select the UL user with the highest bid
                prices_BS(selec_DL(idxSelec)) = sort_prices_BS(1);
                % Assign this DL user to this UL user
                flag_assign(idxUL_bidders(idx_max_bid(1))) = 1;
                % Remove the assignment for the other UL users
                assign_ul_to_dl(idxUL_bidders(idx_max_bid(2:end))) = 0;
                % Remove from the assignment the other UL users
                flag_assign(idxUL_bidders(idx_max_bid(2:end))) = 0;
                % Flag necessary to count the number of iterations
                flag_update = 1;
            end
            % Update the number of iterations
            auction_itr = auction_itr + flag_update*1;
            % Update the prices for UL users
            Prices_UE(:,selec_DL(idxSelec)) = prices_BS(selec_DL(idxSelec));
        end
    end
    % End of the Auction
    
    %% Start to assign the UL to DL users in any resource (remember that it
    % is non-frequency selective fading)
    % Start with DL, since lambdadl>= lambdaul always
    for idxDL = 1:par.lambdadl
        % UL user to which this DL user is assigned to
        idx_UL_selec = find(assign_ul_to_dl == idxDL, 1);
        if isempty(idx_UL_selec) % This DL user is not sharing the resource
            % Assign resource DL
            xDL(usersDLFull(idxDL),idxCell,freqAv(idxDL)) = 1;
            % Power
            pDL(usersDLFull(idxDL),idxCell,freqAv(idxDL)) = par.pmaxDL;
            % SINR
            SINR_DL(usersDLFull(idxDL),idxCell,freqAv(idxDL)) = par.pmaxDL*...
                gDL(usersDLFull(idxDL),idxCell,freqAv(idxDL))/par.noise;
            % Minimum SINR of this resource
            SINR_Min(freqAv(idxDL),idxCell) = SINR_DL(usersDLFull(idxDL),idxCell,freqAv(idxDL));
        else % This DL user is sharing the resource with a UL user
            % Assign resource to UL and DL
            xUL(usersULFull(idx_UL_selec),idxCell,freqAv(idxDL)) = chosen_CP(usersUL(idx_UL_selec),...
                usersDL(idxDL));
            xDL(usersDLFull(idxDL),idxCell,freqAv(idxDL)) = chosen_CP(usersUL(idx_UL_selec),...
                usersDL(idxDL));
            % Assign the power to users in UL and DL
            switch chosen_CP(usersUL(idx_UL_selec),usersDL(idxDL))
                case 12 % Point M
                    % Powers
                    pUL(usersULFull(idx_UL_selec),idxCell,freqAv(idxDL)) = par.pmaxUL;
                    pDL(usersDLFull(idxDL),idxCell,freqAv(idxDL)) = par.pmaxDL;
                    % SINRs
                    SINR_UL(usersULFull(idx_UL_selec),idxCell,freqAv(idxDL)) = SINR_UL_M(usersUL(idx_UL_selec),...
                        usersDL(idxDL));
                    SINR_DL(usersDLFull(idxDL),idxCell,freqAv(idxDL)) = SINR_DL_M(usersUL(idx_UL_selec),...
                        usersDL(idxDL));
                case 21 % Point O
                    % Powers
                    pUL(usersULFull(idx_UL_selec),idxCell,freqAv(idxDL)) = p_UL_O(usersUL(idx_UL_selec),...
                        usersDL(idxDL));
                    pDL(usersDLFull(idxDL),idxCell,freqAv(idxDL)) = par.pmaxDL;
                    % SINRs
                    SINR_UL(usersULFull(idx_UL_selec),idxCell,freqAv(idxDL)) = SINR_UL_O(usersUL(idx_UL_selec),...
                        usersDL(idxDL));
                    SINR_DL(usersDLFull(idxDL),idxCell,freqAv(idxDL)) = par.SINR_tgt_dl;
                case {13, 22} % Point N - Area 1 or 2
                    % Powers
                    pUL(usersULFull(idx_UL_selec),idxCell,freqAv(idxDL)) = p_UL_N(usersUL(idx_UL_selec),...
                        usersDL(idxDL));
                    pDL(usersDLFull(idxDL),idxCell,freqAv(idxDL)) = par.pmaxDL;
                    % SINRs
                    SINR_UL(usersULFull(idx_UL_selec),idxCell,freqAv(idxDL)) = par.SINR_tgt_ul;
                    SINR_DL(usersDLFull(idxDL),idxCell,freqAv(idxDL)) = SINR_DL_N(usersUL(idx_UL_selec),...
                        usersDL(idxDL));
                case {11, 31} % Point L - Area 1 or 3
                    % Powers
                    pUL(usersULFull(idx_UL_selec),idxCell,freqAv(idxDL)) = par.pmaxUL;
                    pDL(usersDLFull(idxDL),idxCell,freqAv(idxDL)) = p_DL_L(usersUL(idx_UL_selec),...
                        usersDL(idxDL));
                    % SINRs
                    SINR_UL(usersULFull(idx_UL_selec),idxCell,freqAv(idxDL)) = SINR_UL_L(usersUL(idx_UL_selec),...
                        usersDL(idxDL));
                    SINR_DL(usersDLFull(idxDL),idxCell,freqAv(idxDL)) = par.SINR_tgt_dl;
                case 32 % Point Q
                    % Powers
                    pUL(usersULFull(idx_UL_selec),idxCell,freqAv(idxDL)) = par.pmaxUL;
                    pDL(usersDLFull(idxDL),idxCell,freqAv(idxDL)) = p_DL_Q(usersUL(idx_UL_selec),...
                        usersDL(idxDL));
                    % SINRs
                    SINR_UL(usersULFull(idx_UL_selec),idxCell,freqAv(idxDL)) = par.SINR_tgt_ul;
                    SINR_DL(usersDLFull(idxDL),idxCell,freqAv(idxDL)) = SINR_DL_Q(usersUL(idx_UL_selec),...
                        usersDL(idxDL));
                case 0 % These 2 users have not been assigned to any pair and cannot form one
                    % Check which one achieves higher SINR with maximum power
                    % UL
                    sinr_users(1) = par.pmaxUL*gUL(usersULFull(idx_UL_selec),idxCell,freqAv(idxDL))/par.noise;
                    % DL
                    sinr_users(2) = par.pmaxDL*gDL(usersDLFull(idxUL(1)),idxCell,freqAv(idxDL))/par.noise;
                    if sinr_users(1) >= sinr_users(2) % UL is chosen
                        % Assignment
                        xUL(usersULFull(idx_UL_selec),idxCell,freqAv(idxDL)) = 1;
                        xDL(usersDLFull(idxDL),idxCell,freqAv(idxDL)) = inf; % Necessary for I different than J
                        % Power
                        pUL(usersULFull(idx_UL_selec),idxCell,freqAv(idxDL)) = par.pmaxUL;
                        % SINR
                        SINR_UL(usersULFull(idx_UL_selec),idxCell,freqAv(idxDL)) = sinr_users(1);
                    else % DL is chosen
                        % Assignment
                        xDL(usersDLFull(idxDL),idxCell,freqAv(idxDL)) = 1;
                        xUL(usersULFull(idx_UL_selec),idxCell,freqAv(idxDL)) = inf;% Necessary for I different than J
                        % Power
                        pDL(usersDLFull(idxDL),idxCell,freqAv(idxDL)) = par.pmaxDL;
                        % SINR
                        SINR_DL(usersDLFull(idxDL),idxCell,freqAv(idxDL)) = sinr_users(2);
                    end
            end
            % Minimum SINR of this resource
            SINR_Min(freqAv(idxDL),idxCell) = min([SINR_UL(usersULFull(idx_UL_selec),idxCell,freqAv(idxDL))...
            SINR_DL(usersDLFull(idxDL),idxCell,freqAv(idxDL))]);
        end
    end
    
end

% Full SINR
SINR = SINR_UL + SINR_DL;
end

