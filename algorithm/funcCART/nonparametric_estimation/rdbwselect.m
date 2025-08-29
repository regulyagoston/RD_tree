function out = rdbwselect(y, x, varargin)
% MATLAB translation of rdbwselect from R
% Required inputs: y (outcome), x (running var)
% Optional name-value pairs

% Default parameters
opts = struct(...
    'c', 0, ...
    'fuzzy', [], ...
    'deriv', 0, ...
    'p', 1, ...
    'ginv_tol', 1e-20, ...
    'kernel', 'Triang', ...
    'weights', [], ...
    'bwselect', 'mserd', ...
    'vce', 'nn', ...
    'cluster', [], ...
    'nnmatch', 3, ...
    'scaleregul', 1, ...
    'sharpbw', false, ...
    'bwcheck', [], ...
    'bwrestrict', true ...
    );

covs_drop_coll = [];

% Parse name-value pair arguments
opts = parseInputs(opts, varargin{:});

% Set polynomial order
p = opts.p;
vce = opts.vce;
cluster = opts.cluster;
c = opts.c;
fuzzy = opts.fuzzy;
ginv_tol = opts.ginv_tol;
weights = opts.weights;
bwselect = opts.bwselect;
nnmatch = opts.nnmatch;
scaleregul = opts.scaleregul;
sharpbw = opts.sharpbw;
bwcheck = opts.bwcheck;
bwrestrict = opts.bwrestrict;
deriv = opts.deriv;
kernel = opts.kernel;
q = p + 1;

% Rescaling and bandwidth pilot
x_iq = quantile(x, 0.75) - quantile(x, 0.25);
BWp = min([std(x), x_iq / 1.349]);
N = length(x);

% Split the data based on cutoff c
ind_l = x < c;
ind_r = x >= c;

X_l = x(ind_l);
X_r = x(ind_r);

x_l_min = min(X_l);
x_l_max = max(X_l);
x_r_min = min(X_r);
x_r_max = max(X_r);

range_l = abs(c - x_l_min);
range_r = abs(c - x_r_max);

Y_l = y(ind_l);
Y_r = y(ind_r);

N_l = length(X_l);
N_r = length(X_r);

x_min = min(x);
x_max = max(x);

N = N_l + N_r;

M_l = N_l;
M_r = N_r;

% Check if p, q, deriv, and nnmatch are integers
p_round = round(p) / p;
q_round = round(q) / q;
d_round = round(deriv + 1) / (deriv + 1);
m_round = round(nnmatch) / nnmatch;

if ((p_round ~= 1 && p > 0) || (q_round ~= 1 && q > 0) || d_round ~= 1 || m_round ~= 1)
    warning('p, q, deriv, and matches should be integer numbers');
    return;  % or: error('Invalid input parameters');
end

% Kernel type and constant C_c
if strcmp(kernel, 'epanechnikov') || strcmp(kernel, 'epa')
    kernel_type = 'Epanechnikov';
    C_c = 2.34;
elseif strcmp(kernel, 'uniform') || strcmp(kernel, 'uni')
    kernel_type = 'Uniform';
    C_c = 1.843;
else
    kernel_type = 'Triangular';
    C_c = 2.576;
end

% ***********************************************************************
% Initialize optional inputs
Z_l = [];
Z_r = [];
T_l = [];
T_r = [];
C_l = [];
C_r = [];
g_l = 0;
g_r = 0;

% Handle nearest-neighbor VCE
if strcmp(vce, 'nn')
    nn_l = ones(N_l, 1);
    nn_r = ones(N_r, 1);

    % MATLAB equivalent of ave(..., FUN=sum)
    [~, ~, gidx_l] = unique(X_l);
    dups_l = accumarray(gidx_l, 1);            % Sum
else
    dups_l = [];
    dups_r = [];
    dupsid_l = [];
    dupsid_r = [];
end

% If cluster is provided
if ~isempty(cluster)
    C_l = cluster(ind_l, :);
    C_r = cluster(ind_r, :);
    g_l = length(unique(C_l));
    g_r = length(unique(C_r));
end

% Fixed weights
fw_l = 0;
fw_r = 0;

% ----------------------------------------------------------------------
% Initial bandwidth calculation
c_bw = C_c * BWp * N^(-1/5);

% Apply bandwidth restriction if enabled
if bwrestrict
    bw_max_l = abs(c - x_min);
    bw_max_r = abs(c - x_max);
    bw_max = max(bw_max_l, bw_max_r);
    c_bw = min(c_bw, bw_max);
end

% Bandwidth adjustment if bwcheck is provided
bw_adj = 0;
if ~isempty(bwcheck)
    bwcheck_l = min(bwcheck, M_l);
    bwcheck_r = min(bwcheck, M_r);

    % NOTE: X_uniq_l and X_uniq_r must be defined beforehand as the sorted unique values
    bw_min_l = abs(X_uniq_l(bwcheck_l) - c) + 1e-8;
    bw_min_r = abs(X_uniq_r(bwcheck_r) - c) + 1e-8;

    c_bw = max([c_bw, bw_min_l, bw_min_r]);
    bw_adj = 1;
end

%% First estimation
C_d_l = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c, q+1, q+1, q+2, c_bw, range_l, 0, vce, nnmatch, kernel, dups_l, dupsid_l, covs_drop_coll, ginv_tol);
C_d_r = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c, q+1, q+1, q+2, c_bw, range_r, 0, vce, nnmatch, kernel, dups_r, dupsid_r, covs_drop_coll, ginv_tol);
%% Second round
if ismember(bwselect, ["msetwo", "certwo", "msecomb2", "cercomb2"])
    d_bw_l = (C_d_l.V / C_d_l.B^2)^C_d_l.rate;
    d_bw_r = (C_d_r.V / C_d_r.B^2)^C_d_l.rate;

    if bwrestrict
        d_bw_l = min(d_bw_l, bw_max_l);
        d_bw_r = min(d_bw_r, bw_max_r);
    end

    if ~isempty(bwcheck)
        d_bw_l = max(d_bw_l, bw_min_l);
        d_bw_r = max(d_bw_r, bw_min_r);
    end

    C_b_l = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c, q, p+1, q+1, c_bw, d_bw_l, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l, covs_drop_coll, ginv_tol);
    b_bw_l = (C_b_l.V / (C_b_l.B^2 + scaleregul * C_b_l.R))^C_b_l.rate;

    C_b_r = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c, q, p+1, q+1, c_bw, d_bw_r, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r, covs_drop_coll, ginv_tol);
    b_bw_r = (C_b_r.V / (C_b_r.B^2 + scaleregul * C_b_r.R))^C_b_l.rate;

    if bwrestrict
        b_bw_l = min(b_bw_l, bw_max_l);
        b_bw_r = min(b_bw_r, bw_max_r);
    end

    C_h_l = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c, p, deriv, q, c_bw, b_bw_l, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l, covs_drop_coll, ginv_tol);
    h_bw_l = (C_h_l.V / (C_h_l.B^2 + scaleregul * C_h_l.R))^C_h_l.rate;

    C_h_r = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c, p, deriv, q, c_bw, b_bw_r, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r, covs_drop_coll, ginv_tol);
    h_bw_r = (C_h_r.V / (C_h_r.B^2 + scaleregul * C_h_r.R))^C_h_l.rate;

    if bwrestrict
        h_bw_l = min(h_bw_l, bw_max_l);
        h_bw_r = min(h_bw_r, bw_max_r);
    end
end

if ismember(bwselect, ["msesum", "cersum", "msecomb1", "msecomb2", "cercomb1", "cercomb2"])
    d_bw_s = ((C_d_l.V + C_d_r.V) / (C_d_r.B + C_d_l.B)^2)^C_d_l.rate;

    if bwrestrict
        d_bw_s = min(d_bw_s, bw_max);
    end

    if ~isempty(bwcheck)
        d_bw_s = max([d_bw_s, bw_min_l, bw_min_r]);
    end

    C_b_l = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c, q, p+1, q+1, c_bw, d_bw_s, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l, covs_drop_coll, ginv_tol);
    C_b_r = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c, q, p+1, q+1, c_bw, d_bw_s, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r, covs_drop_coll, ginv_tol);

    b_bw_s = ((C_b_l.V + C_b_r.V) / ((C_b_r.B + C_b_l.B)^2 + scaleregul * (C_b_r.R + C_b_l.R)))^C_b_l.rate;

    if bwrestrict
        b_bw_s = min(b_bw_s, bw_max);
    end

    C_h_l = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c, p, deriv, q, c_bw, b_bw_s, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l, covs_drop_coll, ginv_tol);
    C_h_r = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c, p, deriv, q, c_bw, b_bw_s, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r, covs_drop_coll, ginv_tol);

    h_bw_s = ((C_h_l.V + C_h_r.V) / ((C_h_r.B + C_h_l.B)^2 + scaleregul * (C_h_r.R + C_h_l.R)))^C_h_l.rate;

    if bwrestrict
        h_bw_s = min(h_bw_s, bw_max);
    end
end


if ismember(bwselect, ["mserd", "cerrd", "msecomb1", "msecomb2", "cercomb1", "cercomb2", ""])
    d_bw_d = ((C_d_l.V + C_d_r.V) / (C_d_r.B - C_d_l.B)^2)^C_d_l.rate;

    if bwrestrict
        d_bw_d = min(d_bw_d, bw_max);
    end

    if ~isempty(bwcheck)
        d_bw_d = max([d_bw_d, bw_min_l, bw_min_r]);
    end

    C_b_l = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c, q, p+1, q+1, c_bw, d_bw_d, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l, covs_drop_coll, ginv_tol);
    C_b_r = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c, q, p+1, q+1, c_bw, d_bw_d, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r, covs_drop_coll, ginv_tol);

    b_bw_d = ((C_b_l.V + C_b_r.V) / ((C_b_r.B - C_b_l.B)^2 + scaleregul * (C_b_r.R + C_b_l.R)))^C_b_l.rate;

    if bwrestrict
        b_bw_d = min(b_bw_d, bw_max);
    end

    C_h_l = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c, p, deriv, q, c_bw, b_bw_d, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l, covs_drop_coll, ginv_tol);
    C_h_r = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c, p, deriv, q, c_bw, b_bw_d, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r, covs_drop_coll, ginv_tol);

    h_bw_d = ((C_h_l.V + C_h_r.V) / ((C_h_r.B - C_h_l.B)^2 + scaleregul * (C_h_r.R + C_h_l.R)))^C_h_l.rate;

    if bwrestrict
        h_bw_d = min(h_bw_d, bw_max);
    end
end

x_sd = 1;
if ismember(bwselect, ["mserd", "cerrd", "msecomb1", "msecomb2", "cercomb1", "cercomb2", ""])
    h_mserd = x_sd * h_bw_d;
    b_mserd = x_sd * b_bw_d;
end

if ismember(bwselect, ["msesum", "cersum", "msecomb1", "msecomb2", "cercomb1", "cercomb2"])
    h_msesum = x_sd * h_bw_s;
    b_msesum = x_sd * b_bw_s;
end

if ismember(bwselect, ["msetwo", "certwo", "msecomb2", "cercomb2"])
    h_msetwo_l = x_sd * h_bw_l;
    h_msetwo_r = x_sd * h_bw_r;
    b_msetwo_l = x_sd * b_bw_l;
    b_msetwo_r = x_sd * b_bw_r;
end

if ismember(bwselect, ["msecomb1", "cercomb1"])
    h_msecomb1 = min([h_mserd, h_msesum]);
    b_msecomb1 = min([b_mserd, b_msesum]);
end

if ismember(bwselect, ["msecomb2", "cercomb2"])
    h_msecomb2_l = median([h_mserd, h_msesum, h_msetwo_l]);
    h_msecomb2_r = median([h_mserd, h_msesum, h_msetwo_r]);
    b_msecomb2_l = median([b_mserd, b_msesum, b_msetwo_l]);
    b_msecomb2_r = median([b_mserd, b_msesum, b_msetwo_r]);
end

cer_h = N^(-(p/((3+p)*(3+2*p))));
if (~isempty(cluster))
  cer_h = (g_l+g_r)^(-(p/((3+p)*(3+2*p))));
end
cer_b = 1;
if strcmp(bwselect, "cerrd")
    h_cerrd = h_mserd * cer_h;
    b_cerrd = b_mserd * cer_b;
end
if strcmp(bwselect, "cersum")
    h_cersum = h_msesum * cer_h;
    b_cersum = b_msesum * cer_b;
end
if strcmp(bwselect, "certwo")
    h_certwo_l = h_msetwo_l * cer_h;
    h_certwo_r = h_msetwo_r * cer_h;
    b_certwo_l = b_msetwo_l * cer_b;
    b_certwo_r = b_msetwo_r * cer_b;
end
if strcmp(bwselect, "cercomb1")
    h_cercomb1 = h_msecomb1 * cer_h;
    b_cercomb1 = b_msecomb1 * cer_b;
end
if strcmp(bwselect, "cercomb2")
    h_cercomb2_l = h_msecomb2_l * cer_h;
    h_cercomb2_r = h_msecomb2_r * cer_h;
    b_cercomb2_l = b_msecomb2_l * cer_b;
    b_cercomb2_r = b_msecomb2_r * cer_b;
end

% Initialize bws as NaN array 1x4
bws = nan(1,4);

% Fill bws according to bwselect
if strcmp(bwselect, "mserd") || strcmp(bwselect, "")
    bws = [h_mserd, h_mserd, b_mserd, b_mserd];
elseif strcmp(bwselect, "msetwo")
    bws = [h_msetwo_l, h_msetwo_r, b_msetwo_l, b_msetwo_r];
elseif strcmp(bwselect, "msesum")
    bws = [h_msesum, h_msesum, b_msesum, b_msesum];
elseif strcmp(bwselect, "msecomb1")
    bws = [h_msecomb1, h_msecomb1, b_msecomb1, b_msecomb1];
elseif strcmp(bwselect, "msecomb2")
    bws = [h_msecomb2_l, h_msecomb2_r, b_msecomb2_l, b_msecomb2_r];
elseif strcmp(bwselect, "cerrd")
    bws = [h_cerrd, h_cerrd, b_cerrd, b_cerrd];
elseif strcmp(bwselect, "certwo")
    bws = [h_certwo_l, h_certwo_r, b_certwo_l, b_certwo_r];
elseif strcmp(bwselect, "cersum")
    bws = [h_cersum, h_cersum, b_cersum, b_cersum];
elseif strcmp(bwselect, "cercomb1")
    bws = [h_cercomb1, h_cercomb1, b_cercomb1, b_cercomb1];
elseif strcmp(bwselect, "cercomb2")
    bws = [h_cercomb2_l, h_cercomb2_r, b_cercomb2_l, b_cercomb2_r];
end

% Prepare output struct
out.h_l = bws(1);
out.h_r = bws(2);
out.b_l = bws(3);
out.b_r = bws(4);

end

function opts = parseInputs(defaults, varargin)
% Parses name-value pairs
opts = defaults;
for i = 1:2:length(varargin)
    name = varargin{i};
    val = varargin{i+1};
    if isfield(opts, name)
        opts.(name) = val;
    else
        error('Unknown parameter: %s', name);
    end
end
end
