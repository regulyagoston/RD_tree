%% rdrobust replication
%
% y, x, fuzzy, cluster are Nx1 vectors
% c is a scalar for cutoff value
% paramLP is a local polynomial parameter object

function out = rdrobust( y , x , c , bw_struct , fuzzy , cluster , ord_poly , vce , optTree )

% Order of polynomials
p = ord_poly;
q = p + 1;

% Order X and Y
[x, ord_x] = sort(x);
y = y(ord_x);
if ~isempty(fuzzy)
    fuzzy = fuzzy(ord_x);
end

%% Two subsamples
id_l = x < c;
id_r = x >= c;
x_l = x( id_l , : );
x_r = x( id_r , : );
y_l = y( id_l , : );
y_r = y( id_r , : );


% Get the bandwidths: 
% h is the main bandwidth
h_l = bw_struct.h_l;
h_r = bw_struct.h_r;
b_l = bw_struct.b_l;
b_r = bw_struct.b_r;

%% Kernel weights
kernFnc = str2func( [ optTree.kernel , '2' ] );
w_h_l = feval( kernFnc , ( ( x_l - c ) ./ h_l ) ) ./ h_l;
w_h_r = feval( kernFnc , ( ( x_r - c ) ./ h_r ) ) ./ h_r;
w_b_l = feval( kernFnc , ( ( x_l - c ) ./ b_l ) ) ./ b_l;
w_b_r = feval( kernFnc , ( ( x_r - c ) ./ b_r ) ) ./ b_r;

% New indexing
ind_h_l = w_h_l > 0;
ind_h_r = w_h_r > 0;
ind_b_l = w_b_l > 0;
ind_b_r = w_b_r > 0;
if h_l > b_l
    ind_l = ind_h_l;
else
    ind_l = ind_b_l;
end
if h_r > b_r
    ind_r = ind_h_r;
else
    ind_r = ind_b_r;
end

eN_l = sum( ind_l );
eN_r = sum( ind_r );

% if ( eN_l < 20 || eN_r < 20 )
%     out = struct;
%     out.tau_bc = NaN;
%     out.tau_cl = NaN;
%     out.tau_V_cl  = Inf;
%     out.tau_V_rb  = Inf;
%     out.tau_V_cl_est  = Inf;
%     out.tau_V_rb_est  = Inf;
%     out.tau_se_cl = Inf;
%     out.tau_se_rb = Inf;
%     out.bias_l    = Inf;
%     out.bias_r    = Inf;
%     out.N         = 0;
%     return;
% end

eY_l = y_l( ind_l );
eY_r = y_r( ind_r );
eX_l = x_l( ind_l );
eX_r = x_r( ind_r );
eW_h_l = w_h_l( ind_l );
eW_h_r = w_h_r( ind_r );
eW_b_l = w_b_l( ind_l );
eW_b_r = w_b_r( ind_r );

edups_l   = zeros( eN_l , 1 );
edups_r   = zeros( eN_r , 1 );
edupsid_l = zeros( eN_l , 1 );
edupsid_r = zeros( eN_r , 1 );

%% Matrices for estimations
eX_lc = eX_l - c;
eX_rc = eX_r - c;
u_l = ( eX_lc ) ./ h_l;
u_r = ( eX_rc ) ./ h_r;
R_q_l = ones( eN_l , q + 1 );
R_q_r = ones( eN_r , q + 1 );
for i = 1 : q
    R_q_l( : , i + 1 ) = eX_lc .^ i;
    R_q_r( : , i + 1 ) = eX_rc .^ i;
end
% Matrices for p
R_p_l = R_q_l( : , 1 : end - 1 );
R_p_r = R_q_r( : , 1 : end - 1 );
R_p_W_h_l = R_p_l .* eW_h_l;
R_p_W_h_r = R_p_r .* eW_h_r;
R_p_sW_h_l = R_p_l .* sqrt( eW_h_l );
R_p_sW_h_r = R_p_r .* sqrt( eW_h_r );
L_l = R_p_W_h_l' * u_l.^( p + 1 );
L_r = R_p_W_h_r' * u_r.^( p + 1 );
iGamma_p_l = pinv( R_p_sW_h_l' * R_p_sW_h_l );
iGamma_p_r = pinv( R_p_sW_h_r' * R_p_sW_h_r );
% aux_1 = chol(  R_p_sW_h_l' * R_p_sW_h_l  );
% iGamma_p_l = inv( aux_1 ) * inv( aux_1' );
% aux_2 = chol(  R_p_sW_h_r' * R_p_sW_h_r  );
% iGamma_p_r = inv( aux_2 ) * inv( aux_2' );

% Matrices for q
R_q_sW_b_l = R_q_l .* sqrt( eW_b_l );
R_q_sW_b_r = R_q_r .* sqrt( eW_b_r );
iGamma_q_l = pinv( R_q_sW_b_l' * R_q_sW_b_l );
iGamma_q_r = pinv( R_q_sW_b_r' * R_q_sW_b_r );

e_p1 = zeros( q + 1 , 1 );
e_p1( p + 2 ) = 1;

Q_q_l = ( R_p_W_h_l' - h_l.^( p + 1 ) .* ( L_l * e_p1' ) * ( ( iGamma_q_l * R_q_l' )' .* eW_b_l )' )';
Q_q_r = ( R_p_W_h_r' - h_r.^( p + 1 ) .* ( L_r * e_p1' ) * ( ( iGamma_q_r * R_q_r' )' .* eW_b_r )' )';

D_l = eY_l;
D_r = eY_r;
eC_l = [];
eC_r = [];
eT_l = [];
eT_r = [];
eZ_l = [];
eZ_r = [];
dZ = 0;
dT = 0;

if ~isempty( fuzzy )
    dT = 1;
    T_l = fuzzy( id_l );
    T_r = fuzzy( id_r );
    eT_l = T_l( ind_l );
    eT_r = T_r( ind_r );
    D_l = [ D_l , eT_l ];
    D_r = [ D_r , eT_r ]; 
end

if ~isempty( cluster )
    C_l = cluster( id_l , : );
    C_r = cluster( id_r , : );
    eC_l = C_l( ind_l , : );
    eC_r = C_r( ind_r , : );
end

%% Estimates
% p'th order
beta_p_l = iGamma_p_l * ( R_p_W_h_l' * D_l );
beta_p_r = iGamma_p_r * ( R_p_W_h_r' * D_r );
% q'th order
beta_q_l = iGamma_q_l * ( ( R_q_l .* eW_b_l)' * D_l );
beta_q_r = iGamma_q_r * ( (R_q_r .* eW_b_r)' * D_r );
% Bias corrected
beta_bc_l = iGamma_p_l * ( Q_q_l' * D_l );
beta_bc_r = iGamma_p_r * ( Q_q_r' * D_r );
% Raw taus
beta_p = beta_p_r - beta_p_l;
%beta_q = beta_q_r - beta_q_l;
beta_bc = beta_bc_r - beta_bc_l;

% Auxilary values
fDeriv = factorial( optTree.deriv );

sFd = optTree.scalepar .* fDeriv;
tau_cl = sFd .* beta_p( optTree.deriv + 1 , 1 ); 
tau_Y_cl = tau_cl;
tau_bc = sFd .* beta_bc( optTree.deriv + 1 , 1 ); 
tau_Y_bc = tau_bc;
s_Y = 1;

tau_Y_cl_l = sFd .* beta_p_l( optTree.deriv + 1 , 1 );
tau_Y_cl_r = sFd .* beta_p_r( optTree.deriv + 1 , 1 );
tau_Y_bc_l = sFd .* beta_bc_l( optTree.deriv + 1 , 1 );
tau_Y_bc_r = sFd .* beta_bc_r( optTree.deriv + 1 , 1 );
bias_l = tau_Y_cl_l - tau_Y_bc_l;
bias_r = tau_Y_cl_r - tau_Y_bc_r;

if ~isempty( fuzzy )
    tau_T_cl = fDeriv .* beta_p( optTree.deriv + 1 , 2 );
    tau_T_bc = fDeriv .* beta_bc( optTree.deriv + 1 , 2 );
    tau_cl = tau_Y_cl ./ tau_T_cl;
    s_Y = [ 1./ tau_T_cl , -( tau_Y_cl./tau_T_cl.^2 ) ];
    B_F = [ tau_Y_cl - tau_Y_bc , tau_T_cl - tau_T_bc ];
    tau_bc = tau_cl - s_Y*B_F';
    sV_T = [ 0 , 1 ];

    tau_T_cl_l = fDeriv .* beta_p_l( optTree.deriv + 1 , 2 );
    tau_T_cl_r = fDeriv .* beta_p_r( optTree.deriv + 1 , 2 );
    tau_T_bc_l = fDeriv .* beta_bc_l( optTree.deriv + 1 , 2 );
    tau_T_bc_r = fDeriv .* beta_bc_r( optTree.deriv + 1 , 2 );

    B_F_l = [ tau_Y_cl_l - tau_Y_bc_l , tau_T_cl_l - tau_T_bc_l ];
    B_F_r = [ tau_Y_cl_r - tau_Y_bc_r , tau_T_cl_r - tau_T_bc_r ];
    bias_l = s_Y * B_F_l';
    bias_r = s_Y * B_F_r'; 
end



%% Estimating the variance
hii_l = 0;
hii_r = 0;
predicts_p_l = 0;
predicts_p_r = 0;
predicts_q_l = 0;
predicts_q_r = 0;
if any( strcmpi( vce , {'simple','hce-1','hce-2','hce-3'} ) )
    predicts_p_l = R_p_l * beta_p_l;
    predicts_p_r = R_p_r * beta_p_r;
    predicts_q_l = R_q_l * beta_q_l;
    predicts_q_r = R_q_r * beta_q_r;
    if any( strcmpi( vce , { 'hce-2' , 'hce-3' } ) )
        hii_l = NaN( eN_l , 1 );
        for i = 1 : eN_l
            hii_l( i ) = R_p_l( i , : ) * iGamma_p_l * R_p_W_h_l( i , : )';
        end
        hii_r = NaN( eN_r , 1 );
        for i = 1 : eN_r
            hii_r( i ) = R_p_r( i , : ) * iGamma_p_r * R_p_W_h_r( i , : )';
        end
    end
end

% Main bandwidth
res_h_l = rdrobust_res( eX_l , eY_l , eT_l , eZ_l , predicts_p_l , hii_l ,...
                        vce , [] , edups_l , edupsid_l , p + 1 );
res_h_r = rdrobust_res( eX_r , eY_r , eT_r , eZ_r , predicts_p_r , hii_r ,...
                        vce , [] , edups_r , edupsid_r , p + 1 );

% Secondary                    
res_b_l = rdrobust_res( eX_l , eY_l , eT_l , eZ_l , predicts_q_l , hii_l ,...
                        vce , [] , edups_l , edupsid_l , q + 1 );
res_b_r = rdrobust_res( eX_r , eY_r , eT_r , eZ_r , predicts_q_r , hii_r ,...
                        vce , [] , edups_r , edupsid_r , q + 1 );


Sigma_Y_p_l = rdrobust_vce( dT + dZ , s_Y , R_p_W_h_l , res_h_l , eC_l , 1 );
V_Y_cl_l = iGamma_p_l * Sigma_Y_p_l * iGamma_p_l;
Sigma_Y_p_r = rdrobust_vce( dT + dZ , s_Y , R_p_W_h_r , res_h_r , eC_r , 1 );
V_Y_cl_r = iGamma_p_r * Sigma_Y_p_r * iGamma_p_r;
Sigma_Y_bc_l = rdrobust_vce( dT + dZ , s_Y , Q_q_l , res_b_l , eC_l , 1 );
V_Y_rb_l = iGamma_p_l * Sigma_Y_bc_l * iGamma_p_l;
Sigma_Y_bc_r = rdrobust_vce( dT + dZ , s_Y , Q_q_r , res_b_r , eC_r , 1 );
V_Y_rb_r = iGamma_p_r * Sigma_Y_bc_r * iGamma_p_r;
V_tau_cl = fDeriv.^2.*( V_Y_cl_l + V_Y_cl_r );
V_tau_cl = V_tau_cl( optTree.deriv + 1 , optTree.deriv + 1 );
V_tau_rb = fDeriv.^2.*( V_Y_rb_l + V_Y_rb_r );
V_tau_rb = V_tau_rb( optTree.deriv + 1 , optTree.deriv + 1 );

if ~isempty( fuzzy )
    Sigma_T_p_l = rdrobust_vce( dT + dZ , sV_T , R_p_W_h_l , res_h_l , eC_l , 1 );
    V_T_cl_l = iGamma_p_l * Sigma_T_p_l * iGamma_p_l;
    Sigma_T_p_r = rdrobust_vce( dT + dZ , sV_T , R_p_W_h_r , res_h_r , eC_r , 1 );
    V_T_cl_r = iGamma_p_r * Sigma_T_p_r * iGamma_p_r;
    Sigma_T_rb_l = rdrobust_vce( dT + dZ , sV_T , Q_q_l , res_b_l , eC_l , 1 );
    V_T_rb_l = iGamma_p_l * Sigma_T_rb_l * iGamma_p_l;
    Sigma_T_rb_r = rdrobust_vce( dT + dZ , sV_T , Q_q_r , res_b_r , eC_r , 1 );
    V_T_rb_r = iGamma_p_r * Sigma_T_rb_r * iGamma_p_r;
    V_T_cl = fDeriv.^2 .* ( V_T_cl_l + V_T_cl_r );
    V_T_cl = V_T_cl( optTree.deriv + 1 , optTree.deriv + 1 );
    V_T_rb = fDeriv.^2 .* ( V_T_rb_l + V_T_rb_r );
    V_T_rb = V_T_rb( optTree.deriv + 1 , optTree.deriv + 1 );
end

% Variance
V_tau_cl_h = V_tau_cl;
V_tau_rb_h = V_tau_rb;
if ~isempty( fuzzy )
    V_T_cl_h   = V_T_cl;
    V_T_rb_h   = V_T_rb;
end

%% Output
out = struct;
out.tau_bc = tau_bc;
out.tau_cl = tau_cl;
out.tau_V_cl  = V_tau_cl;
out.tau_V_rb  = V_tau_rb;
out.tau_V_cl_est  = V_tau_cl_h;
out.tau_V_rb_est  = V_tau_rb_h;
out.tau_se_cl = optTree.scalepar .* sqrt( V_tau_cl );
out.tau_se_rb = optTree.scalepar .* sqrt( V_tau_rb );
out.bias_l    = bias_l;
out.bias_r    = bias_r;
out.N         = eN_r + eN_l;
if ~isempty( fuzzy )
    out.tau_Y_bc = tau_Y_bc;
    out.tau_Y_cl = tau_Y_cl;
    out.tau_T_bc = tau_T_bc;
    out.tau_T_cl = tau_T_cl;
    out.tau_V_T_cl  = V_T_cl;
    out.tau_V_T_rb  = V_T_rb;
    out.tau_V_T_cl_est  = V_T_cl_h;
    out.tau_V_T_rb_est  = V_T_rb_h;
end


end