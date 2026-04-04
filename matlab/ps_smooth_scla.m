function []=ps_smooth_scla(use_small_baselines)
%PS_SMOOTH_SCLA Smooth spatially-correlated look angle error (HPC Optimized)
%   Removes outliers in look angle error and master atmosphere/orbit error
%   by comparing against neighbor extremes defined by Delaunay triangulation.
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          February 2026
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization:
%   1. Topology Generation: Replaced legacy 'triangle' C-call and disk I/O with 
%      native 'delaunayTriangulation', drastically improving speed and reliability.
%   2. Vectorized Graph Traversal: Replaced millions of loop iterations with 
%      'accumarray', dropping neighborhood extreme retrieval time to seconds.
%   3. Memory Projection & Expansion: Replaced giant matrix transpose division 
%      and 'repmat' expansions with Streamed Block-IO Projection (H_op) and native 
%      implicit expansion, effectively eliminating >80GB RAM spikes.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, March 2007
%   ======================================================================

logit;
logit('Smoothing spatially-correlated look angle error (HPC Mode)...')

if nargin<1
    use_small_baselines=0;
end

small_baseline_flag=getparm('small_baseline_flag');

load psver
psname=['ps',num2str(psver)];
bpname=['bp',num2str(psver)];
if use_small_baselines==0
    sclaname=['scla',num2str(psver)];
    sclasmoothname=['scla_smooth',num2str(psver)];
else
    sclaname=['scla_sb',num2str(psver)];
    sclasmoothname=['scla_smooth_sb',num2str(psver)];
end

ps=load(psname);
scla=load(sclaname);
K_ps_uw=scla.K_ps_uw;
C_ps_uw=scla.C_ps_uw;
ph_ramp=scla.ph_ramp;
clear scla

n_ps=ps.n_ps;
logit(sprintf('Number of points per ifg: %d', n_ps))

% --- HPC Opt 1: Native Delaunay Triangulation (No I/O) ---
logit('   -> Building Delaunay Triangulation...', 2);
DT=delaunayTriangulation(double(ps.xy(:,2:3)));
edgs=edges(DT);

n_edge=size(edgs,1);
logit(sprintf('Number of arcs per ifg: %d', n_edge))

% --- HPC Opt 2: Vectorized neighbor extreme search ---
logit('   -> Vectorizing neighbor extreme search...', 2);

% Build bidirectional edges for undirected graph
src = [edgs(:,1); edgs(:,2)];
dst = [edgs(:,2); edgs(:,1)];

% Extract property values for destination nodes
K_vals = K_ps_uw(dst);
C_vals = C_ps_uw(dst);

% Dynamic type casting for accumarray fill values
fill_K_inf = cast(inf, class(K_vals));
fill_C_inf = cast(inf, class(C_vals));

% Fast aggregation using accumarray to find neighbor extremes
Kneigh_max = accumarray(src, K_vals, [n_ps, 1], @max, -fill_K_inf);
Kneigh_min = accumarray(src, K_vals, [n_ps, 1], @min,  fill_K_inf);
Cneigh_max = accumarray(src, C_vals, [n_ps, 1], @max, -fill_C_inf);
Cneigh_min = accumarray(src, C_vals, [n_ps, 1], @min,  fill_C_inf);

% --- Smooth outliers (Original logic) ---
ix1 = K_ps_uw > Kneigh_max;
ix2 = K_ps_uw < Kneigh_min;
K_ps_uw(ix1) = Kneigh_max(ix1); % reduce positive outliers
K_ps_uw(ix2) = Kneigh_min(ix2); % increase negative outliers
   
ix1 = C_ps_uw > Cneigh_max;
ix2 = C_ps_uw < Cneigh_min;
C_ps_uw(ix1) = Cneigh_max(ix1); % reduce positive outliers
C_ps_uw(ix2) = Cneigh_min(ix2); % increase negative outliers

% --- HPC Opt 3: Memory-safe projection operator ---
logit('   -> Calculating bperp matrix mapping...', 2);
m_bp = matfile(bpname);
sz_bp = size(m_bp, 'bperp_mat');
is_bp_vector = (sz_bp(1) == 1); 

if use_small_baselines==0
    if strcmpi(small_baseline_flag,'y')
        G = zeros(ps.n_ifg, ps.n_image);
        for i = 1:ps.n_ifg
             G(i, ps.ifgday_ix(i,1)) = -1;
             G(i, ps.ifgday_ix(i,2)) =  1;
        end
        G = G(:, [1:ps.master_ix-1, ps.master_ix+1:end]);
        
        % Pre-calculate pseudo-inverse operator H_op
        H_op = (G' * G) \ G'; 
        
        % Block-flow pseudo-inverse projection
        bperp_mat_raw = zeros(ps.n_ps, size(H_op, 1), 'single');
        block_size = 1000000;
        for b = 1:ceil(ps.n_ps / block_size)
            idx = (b-1)*block_size + 1 : min(b*block_size, ps.n_ps);
            if is_bp_vector
                bp_chunk = repmat(m_bp.bperp_mat, length(idx), 1);
            else
                bp_chunk = m_bp.bperp_mat(idx, :);
            end
            bperp_mat_raw(idx, :) = single(double(bp_chunk) * H_op');
        end
        % Insert zero column for master image
        bperp_mat = [bperp_mat_raw(:, 1:ps.master_ix-1), zeros(ps.n_ps, 1, 'single'), bperp_mat_raw(:, ps.master_ix:end)];
        clear bperp_mat_raw;
    else
        if is_bp_vector
            bp_full = repmat(m_bp.bperp_mat, ps.n_ps, 1);
        else
            bp_full = m_bp.bperp_mat;
        end
        bperp_mat = [bp_full(:, 1:ps.master_ix-1), zeros(ps.n_ps, 1, 'single'), bp_full(:, ps.master_ix:end)];
        clear bp_full;
    end
else
    if is_bp_vector
        bperp_mat = repmat(m_bp.bperp_mat, ps.n_ps, 1);
    else
        bperp_mat = m_bp.bperp_mat;
    end
end
% =============================================================================================

% Calculate smoothed SCLA phase
ph_scla = K_ps_uw .* bperp_mat;

stamps_save(sclasmoothname, K_ps_uw, C_ps_uw, ph_scla, ph_ramp)    
logit(1);