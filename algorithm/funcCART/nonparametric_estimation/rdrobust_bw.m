function output = rdrobust_bw(Y, X, T, Z, C, W, c, o, nu, o_B, h_V, h_B, scale, vce, nnmatch, kernel, dups, dupsid, covs_drop_coll, ginv_tol)
% MATLAB version of rdrobust_bw from R
% Accepts name-value pairs

    
    dT = 0;
    dZ = 0;
    sC = 0;
    % Step 1: V
    kernFnc = str2func( [ kernel , '2' ] );
    w = feval( kernFnc , ( ( X - c ) ./ h_V ) ) ./ h_V;
    if W ~= 0
        w = W .* w;
    end

    ind_V = (w > 0);
    eY = Y(ind_V);
    eX = X(ind_V);
    eW = w(ind_V);
    n_V = sum(ind_V);
    D_V = eY;
    R_V = zeros(n_V, o + 1);
    for j = 1:(o + 1)
        R_V(:, j) = (eX - c).^(j - 1);
    end
    aux_m = R_V .* repmat( sqrt(eW), [1,o+1]);
    U = chol(aux_m' * aux_m);
    invG_V = inv(U) * inv(U)'; %qrXXinv(R_V .* sqrt(eW));

    e_v = zeros(o + 1, 1);
    e_v(nu + 1) = 1;

    s = 1;
    dT = 0;
    dZ = 0;
    dC = 0;

    if ~isempty(T)
        dT = 1;
        eT = T(ind_V, :);
        D_V = [D_V, eT];
    else
        eT = [];
    end

    if ~isempty(Z)
        dZ = size(Z, 2);
        eZ = Z(ind_V, :);
        D_V = [D_V, eZ];

        U = (R_V' .* eW') * D_V;
        ZWD = (eZ' .* eW') * D_V;
        colsZ = (2 + dT):(1 + dT + dZ);
        UiGU = (U(:, colsZ)') * (invG_V * U);

        ZWZ = ZWD(:, colsZ) - UiGU(:, colsZ);
        ZWY = ZWD(:, 1:(1 + dT)) - UiGU(:, 1:(1 + dT));

        if covs_drop_coll == 1
            gamma = pinv(ZWZ, ginv_tol) * ZWY;
        else
            gamma = ZWZ \ ZWY;
        end

        s = [1; -gamma(:, 1)];
    else
        eZ = [];
    end

    if ~isempty(C)
        dC = 1;
        eC = C(ind_V);
    else
        eC = [];
    end

    beta_V = invG_V * ((R_V' .* eW') * D_V);

    if isempty(Z) && ~isempty(T)
        tau_Y = factorial(nu) * beta_V(nu + 1, 1);
        tau_T = factorial(nu) * beta_V(nu + 1, 2);
        s = [1 / tau_T; -(tau_Y / tau_T^2)];
    end

    if ~isempty(Z) && ~isempty(T)
        s_T = [1; -gamma(:, 2)];
        tau_Y = factorial(nu) * [beta_V(nu + 1, 1); beta_V(nu + 1, colsZ)]' * s;
        tau_T = factorial(nu) * [beta_V(nu + 1, 2); beta_V(nu + 1, colsZ)]' * s_T;
        s = [1 / tau_T;
            -(tau_Y / tau_T^2);
            -(1 / tau_T) * gamma(:, 1) + (tau_Y / tau_T^2) * gamma(:, 2)];
    end

    dups_V = 0; dupsid_V = 0; predicts_V = 0; hii = [];

    if strcmp(vce, 'nn')
        dups_V = dups(ind_V);
        dupsid_V = dupsid(ind_V);
    elseif any(strcmp(vce, {'simple','hce-0', 'hce-1', 'hce-2', 'hce-3'}))
        predicts_V = R_V * beta_V;
        hii = [];
        if any(strcmp(vce, {'hce-2', 'hce-3'}))
            hii = sum((R_V * invG_V) .* (R_V .* eW), 2);
        end
    end

    res_V = rdrobust_res(eX, eY, eT, eZ, predicts_V, hii, vce, ...
                         nnmatch, dups_V, dupsid_V, o + 1);
    aux = rdrobust_vce(dT + dZ, s, R_V .* eW, res_V, eC, 1);

    V_V = (invG_V * aux * invG_V);
    V_V = V_V(nu + 1, nu + 1);

    v = (R_V' .* eW') * ((eX - c) / h_V).^(o + 1);
    Hp = h_V.^((0:o)');
    BConst = (Hp .* (invG_V * v));
    BConst = BConst(nu + 1);

    % Step 2: B
    w = feval( kernFnc , ( ( X - c ) ./ h_B ) ) ./ h_B;
    if W ~= 0
        w = W .* w;
    end

    ind_B = (w > 0);
    eY = Y(ind_B);
    eX = X(ind_B);
    eW = w(ind_B);
    D_B = eY;
    n_B = sum(ind_B);

    R_B = zeros(n_B, o_B + 1);
    for j = 1:(o_B + 1)
        R_B(:, j) = (eX - c).^(j - 1);
    end
    aux_b = R_B .* repmat( sqrt(eW), [1,o_B + 1]);
    Ub = chol(aux_b' * aux_b);
    invG_B = inv(Ub) * inv(Ub)'; %qrXXinv(R_V .* sqrt(eW));

    if ~isempty(T)
        eT = T(ind_B, :);
        D_B = [D_B, eT];
    else
        eT = [];
    end

    if ~isempty(Z)
        eZ = Z(ind_B, :);
        D_B = [D_B, eZ];
    else
        eZ = [];
    end

    if ~isempty(C)
        eC = C(ind_B);
    end

    beta_B = invG_B * ((R_B' .* eW') * D_B);

    BWreg = 0;
    if scale > 0
        e_B = zeros(o_B + 1, 1);
        e_B(o + 2) = 1;
        dups_B = 0; 
        dupsid_B = 0; 
        hii = []; 
        predicts_B = 0;

        if strcmp(vce, 'nn')
            dups_B = dups(ind_B);
            dupsid_B = dupsid(ind_B);
        elseif any(strcmp(vce, {'simple','hce-0', 'hce-1', 'hce-2', 'hce-3'}))
            predicts_B = R_B * beta_B;
            if any(strcmp(vce, {'hce-2', 'hce-3'}))
                hii = sum((R_B * invG_B) .* (R_B .* eW), 2);
            end
        end

        res_B = rdrobust_res(eX, eY, eT, eZ, predicts_B, hii, vce, ...
                             nnmatch, dups_B, dupsid_B, o_B + 1);
        V_B = (invG_B * rdrobust_vce(dT + dZ, s, R_B .* eW, res_B, eC, 1) * invG_B ) ;
        V_B = V_B(o + 2, o + 2);
        BWreg = 3 * BConst^2 * V_B;
    end

    B = sqrt(2 * (o + 1 - nu)) * BConst * (s' * beta_B(o + 2, :)');
    V = (2 * nu + 1) * h_V^(2 * nu + 1) * V_V;
    R = scale * (2 * (o + 1 - nu)) * BWreg;
    rate = 1 / (2 * o + 3);

    output = struct('V', V, 'B', B, 'R', R, 'rate', rate);
end
