function [dna_constrained, stats] = analyse_dna_constraints(encoded_image, motifs, max_run, max_iters, seq_len)
%ANALYSE_DNA_CONSTRAINTS
%  对 encoded_image 做：
%    1) 约束前：GC、同聚物、禁忌 motif 统计
%    2) 约束编码：限制同聚物长度 + 打碎 forbidden motifs + 轻量 GC 平衡
%    3) 约束后：再做同样统计
%
%  [dna_constrained, stats] = analyse_dna_constraints(encoded_image, motifs, max_run, max_iters, seq_len)
%
%  encoded_image : H×W 或 H×W×C 的 uint8 密文图像（DNA 层后）
%  motifs        : 禁忌 motif 列表（cell数组），可留空用默认
%  max_run       : 允许的最大同聚物长度，默认 3（比原来更严格一点）
%  max_iters     : 约束迭代轮数，默认 5
%  seq_len       : 统计时的序列长度，默认 512
%
%  返回：
%    dna_constrained : 约束后 DNA 序列（与 encoded_image 对应的 char 矩阵，'A','C','G','T'）
%    stats           : 结构体，含 before/after 各种统计

    if nargin < 2 || isempty(motifs)
        motifs = {'AAAAAA','CCCCCC','GGGGGG','TTTTTT', ...
                  'GCGC','CGCG','GAATTC','AAGCTT'};
    end
    if nargin < 3 || isempty(max_run)
        max_run = 3;      % 比原来4更严格一点
    end
    if nargin < 4 || isempty(max_iters)
        max_iters = 5;
    end
    if nargin < 5 || isempty(seq_len)
        seq_len = 512;
    end

    fprintf('================ DNA constraint analysis ================\n');

    %--------------------------------------------------------
    % 1) 把 encoded_image 映射到 DNA 字符矩阵 'A','C','G','T'
    %    优先使用严格 bit-level 映射
    %--------------------------------------------------------
    [H0, W0, C0] = size(encoded_image);
    try
        dna_raw = to_dna_char_matrix_strict(encoded_image, [H0, W0*4, C0]);
    catch
        dna_raw = local_to_dna_chars(encoded_image);
    end

    dna_raw = upper(char(dna_raw));
    [H,W,C] = size(dna_raw);
    total_bases = numel(dna_raw);

    fprintf('[info] DNA matrix size = %d x %d x %d, total_bases = %d\n', H, W, C, total_bases);

    %--------------------------------------------------------
    % 2) 约束前统计
    %--------------------------------------------------------
    fprintf('\n>>> BEFORE CONSTRAINT:\n');
    stats_before.gc    = compute_gc_stats(dna_raw, seq_len, '  ');
    stats_before.homo  = compute_homopolymer_stats(dna_raw, seq_len, '  ');
    stats_before.motif = compute_forbidden_stats(dna_raw, motifs, seq_len, '  ');

    %--------------------------------------------------------
    % 3) 约束编码：限制同聚物 + 打碎 forbidden motif + 轻量 GC 平衡
    %--------------------------------------------------------
    fprintf('\n>>> RUN CONSTRAINT ENCODER (max_run=%d, max_iters=%d)...\n', max_run, max_iters);
    [dna_constrained, cons_stats] = constraint_encode_dna(dna_raw, motifs, max_run, max_iters);
    fprintf('  [constraint_encode_dna] total changes: %d, iters: %d\n', ...
            cons_stats.total_changes, cons_stats.num_iters);

    %--------------------------------------------------------
    % 4) 约束后统计
    %--------------------------------------------------------
    fprintf('\n>>> AFTER CONSTRAINT:\n');
    stats_after.gc    = compute_gc_stats(dna_constrained, seq_len, '  ');
    stats_after.homo  = compute_homopolymer_stats(dna_constrained, seq_len, '  ');
    stats_after.motif = compute_forbidden_stats(dna_constrained, motifs, seq_len, '  ');

    %--------------------------------------------------------
    % 5) 汇总到 stats 结构体
    %--------------------------------------------------------
    stats = struct();
    stats.motifs      = motifs;
    stats.max_run     = max_run;
    stats.seq_len     = seq_len;
    stats.constraint  = cons_stats;
    stats.before      = stats_before;
    stats.after       = stats_after;

    fprintf('================== DONE DNA CONSTRAINT ==================\n\n');
end

% ============================================================
% 辅助函数 1：数值/字符 -> 'A','C','G','T'
% ============================================================
function dna = local_to_dna_chars(X)
    if isstring(X) || ischar(X)
        dna = upper(char(X));
        if any(~ismember(dna(:), 'ACGT'))
            error('local_to_dna_chars: 非 A/C/G/T 字符存在。');
        end
        return;
    end
    X = double(X);
    v = mod(X,4);
    map = 'ACGT';
    dna = map(v+1);
end

% ============================================================
% 辅助函数 2：GC 统计
% ============================================================
function gc_stats = compute_gc_stats(dna, seq_len, prefix)
    if nargin < 3, prefix = ''; end
    dna_vec = dna(:)';
    N = numel(dna_vec);
    num_seq = floor(N / seq_len);
    if num_seq < 1
        error('compute_gc_stats: total_bases(%d) < seq_len(%d)。', N, seq_len);
    end

    G = (dna_vec=='G');
    C = (dna_vec=='C');
    gc_global = (sum(G)+sum(C))/N * 100;

    gc_per_seq = zeros(num_seq,1);
    for s = 1:num_seq
        seg = dna_vec((s-1)*seq_len+1 : s*seq_len);
        gc_per_seq(s) = (sum(seg=='G') + sum(seg=='C')) / seq_len * 100;
    end

    gc_mean = mean(gc_per_seq);
    gc_std  = std(gc_per_seq);

    fprintf('%s[GC stats]\n', prefix);
    fprintf('%s  Total bases          : %d\n', prefix, N);
    fprintf('%s  Num sequences        : %d\n', prefix, num_seq);
    fprintf('%s  Seq length           : %d\n', prefix, seq_len);
    fprintf('%s  GC global            : %.2f %%\n', prefix, gc_global);
    fprintf('%s  GC per-seq (mean±std): %.2f ± %.2f %%\n', prefix, gc_mean, gc_std);

    gc_stats = struct();
    gc_stats.total_bases = N;
    gc_stats.num_sequences = num_seq;
    gc_stats.seq_length = seq_len;
    gc_stats.gc_global  = gc_global;
    gc_stats.gc_mean    = gc_mean;
    gc_stats.gc_std     = gc_std;
end

% ============================================================
% 辅助函数 3：同聚物统计
% ============================================================
function homo_stats = compute_homopolymer_stats(dna, seq_len, prefix)
    if nargin < 3, prefix=''; end
    dna_vec = dna(:)';
    N = numel(dna_vec);
    num_seq = floor(N / seq_len);

    all_runs = [];
    max_run_per_seq = zeros(num_seq,1);

    for s = 1:num_seq
        seg = dna_vec((s-1)*seq_len+1 : s*seq_len);
        runs = local_run_lengths(seg);
        all_runs = [all_runs; runs]; %#ok<AGROW>
        max_run_per_seq(s) = max(runs);
    end

    overall_max_run    = max(max_run_per_seq);
    mean_run_length    = mean(all_runs);
    median_run_length  = median(all_runs);
    mean_max_run_per_seq = mean(max_run_per_seq);
    max_run_per_seq_max  = overall_max_run;

    fprintf('%s[Homopolymer stats]\n', prefix);
    fprintf('%s  Num sequences        : %d\n', prefix, num_seq);
    fprintf('%s  Max run (overall)    : %d\n', prefix, overall_max_run);
    fprintf('%s  Run length (mean/median)  : %.2f / %.2f\n', prefix, ...
            mean_run_length, median_run_length);
    fprintf('%s  Max run per-seq (mean/max): %.2f / %d\n', prefix, ...
            mean_max_run_per_seq, max_run_per_seq_max);

    homo_stats = struct();
    homo_stats.num_sequences        = num_seq;
    homo_stats.overall_max_run      = overall_max_run;
    homo_stats.mean_run_length      = mean_run_length;
    homo_stats.median_run_length    = median_run_length;
    homo_stats.mean_max_run_per_seq = mean_max_run_per_seq;
    homo_stats.max_run_per_seq_max  = max_run_per_seq_max;
end

function runs = local_run_lengths(seq_row)
    if isempty(seq_row)
        runs = [];
        return;
    end
    seq_row = char(seq_row);
    current_base = seq_row(1);
    current_len  = 1;
    runs = [];
    for k = 2:numel(seq_row)
        if seq_row(k) == current_base
            current_len = current_len + 1;
        else
            runs = [runs; current_len]; %#ok<AGROW>
            current_base = seq_row(k);
            current_len  = 1;
        end
    end
    runs = [runs; current_len];
end

% ============================================================
% 辅助函数 4：禁忌 motif 统计
% ============================================================
function motif_stats = compute_forbidden_stats(dna, motifs, seq_len, prefix)
    if nargin < 4, prefix=''; end
    dna_vec = dna(:)';
    N = numel(dna_vec);
    num_seq = floor(N / seq_len);

    fprintf('%s[Forbidden motif stats]\n', prefix);
    fprintf('%s  Num sequences        : %d\n', prefix, num_seq);
    fprintf('%s  Seq length           : %d\n', prefix, seq_len);

    K = numel(motifs);
    total_occ          = zeros(1,K);
    seq_with_motif     = zeros(1,K);
    max_occ_per_seq    = zeros(1,K);

    for k = 1:K
        motif = upper(char(motifs{k}));
        occ_per_seq = zeros(num_seq,1);

        for s = 1:num_seq
            seg = dna_vec((s-1)*seq_len+1 : s*seq_len);
            idxs = strfind(seg, motif);
            occ_per_seq(s) = numel(idxs);
        end

        total_occ(k)       = sum(occ_per_seq);
        seq_with_motif(k)  = sum(occ_per_seq > 0);
        max_occ_per_seq(k) = max(occ_per_seq);

        frac = seq_with_motif(k)/num_seq*100;
        fprintf('%s  Motif "%s": total_occ=%d, seq_with_motif=%d (%.2f%%), max_occ_per_seq=%d\n', ...
            prefix, motifs{k}, total_occ(k), seq_with_motif(k), frac, max_occ_per_seq(k));
    end

    motif_stats = struct();
    motif_stats.motifs              = motifs;
    motif_stats.total_occ           = total_occ;
    motif_stats.seq_with_motif      = seq_with_motif;
    motif_stats.seq_with_motif_frac = seq_with_motif/num_seq;
    motif_stats.max_occ_per_seq     = max_occ_per_seq;
end

% ============================================================
% 辅助函数 5：约束编码器（同聚物 + forbidden + 轻量GC平衡）
% ============================================================
function [dna_fixed, cons_stats] = constraint_encode_dna(dna_raw, motifs, max_run, max_iters)
    dna_raw = upper(dna_raw);
    if any(~ismember(dna_raw(:),'ACGT'))
        error('constraint_encode_dna: dna_raw 中存在非 A/C/G/T 字符');
    end

    s = dna_raw(:)';
    N = numel(s);
    total_changes = 0;

    for it = 1:max_iters
        changes_this_iter = 0;

        % Pass 1: 限制同聚物
        prev = s(1);
        run_len = 1;
        for i = 2:N
            if s(i) == prev
                run_len = run_len + 1;
                if run_len > max_run
                    left_char  = s(max(1, i-1));
                    right_char = s(min(N, i+1));
                    s(i) = pick_alt_base_smart(s(i), left_char, right_char, s, i, motifs, 0.5);
                    changes_this_iter = changes_this_iter + 1;
                    total_changes     = total_changes + 1;
                    prev = s(i);
                    run_len = 1;
                end
            else
                prev = s(i);
                run_len = 1;
            end
        end

        % Pass 2: 打碎 forbidden motifs
        for k = 1:numel(motifs)
            motif = upper(char(motifs{k}));
            Lm = length(motif);
            while true
                idxs = strfind(s, motif);
                if isempty(idxs)
                    break;
                end
                p = idxs(1);
                q = p + floor((Lm-1)/2);
                left_char  = s(max(1, q-1));
                right_char = s(min(N, q+1));
                s(q) = pick_alt_base_smart(s(q), left_char, right_char, s, q, motifs, 0.5);
                changes_this_iter = changes_this_iter + 1;
                total_changes     = total_changes + 1;
            end
        end

        % Pass 3: 轻量 GC 平衡（局部窗口）
        [s, gc_changes] = local_gc_balance(s, max_run, motifs);
        changes_this_iter = changes_this_iter + gc_changes;
        total_changes     = total_changes + gc_changes;

        if changes_this_iter == 0
            break;
        end
    end

    dna_fixed = reshape(s, size(dna_raw));
    cons_stats = struct();
    cons_stats.total_changes = total_changes;
    cons_stats.num_iters     = it;
end

% ============================================================
% 智能替换：尽量避免新run / motif，同时兼顾GC平衡
% ============================================================
function b_new = pick_alt_base_smart(b_old, left_char, right_char, s, pos, motifs, target_gc)
    candidates = 'ACGT';
    candidates(candidates == b_old) = [];

    best_score = inf;
    b_new = 'A';

    current_gc = mean(s == 'G' | s == 'C');

    for k = 1:numel(candidates)
        c = candidates(k);

        score = 0;

        % 1) 避免和左右邻居相同
        if c == left_char
            score = score + 5;
        end
        if c == right_char
            score = score + 5;
        end

        % 2) 尽量靠近目标 GC
        gc_after = current_gc;
        if ismember(b_old, 'GC')
            gc_after = gc_after - 1/numel(s);
        end
        if ismember(c, 'GC')
            gc_after = gc_after + 1/numel(s);
        end
        score = score + abs(gc_after - target_gc) * 10;

        % 3) 避免在局部制造 forbidden motif
        s_tmp = s;
        s_tmp(pos) = c;
        local_seg = get_local_segment(s_tmp, pos, 12);
        for t = 1:numel(motifs)
            if ~isempty(strfind(local_seg, upper(char(motifs{t})))) %#ok<STREMP>
                score = score + 20;
            end
        end

        % 4) 避免制造长 run
        local_max_run = get_local_max_run(local_seg);
        score = score + local_max_run;

        if score < best_score
            best_score = score;
            b_new = c;
        end
    end
end

% ============================================================
% 局部 GC 平衡
% ============================================================
function [s, changes] = local_gc_balance(s, max_run, motifs)
    changes = 0;
    N = numel(s);

    win = 80;
    step = 40;

    for st = 1:step:(N-win+1)
        ed = st + win - 1;
        seg = s(st:ed);
        gc = mean(seg == 'G' | seg == 'C');

        if gc > 0.60
            idx_local = find(seg == 'G' | seg == 'C');
            for kk = 1:numel(idx_local)
                p = st + idx_local(kk) - 1;
                left_char  = s(max(1, p-1));
                right_char = s(min(N, p+1));
                cand = pick_alt_from_set('AT', left_char, right_char, s, p, motifs, max_run);
                if cand ~= s(p)
                    s(p) = cand;
                    changes = changes + 1;
                    break;
                end
            end
        elseif gc < 0.40
            idx_local = find(seg == 'A' | seg == 'T');
            for kk = 1:numel(idx_local)
                p = st + idx_local(kk) - 1;
                left_char  = s(max(1, p-1));
                right_char = s(min(N, p+1));
                cand = pick_alt_from_set('GC', left_char, right_char, s, p, motifs, max_run);
                if cand ~= s(p)
                    s(p) = cand;
                    changes = changes + 1;
                    break;
                end
            end
        end
    end
end

function b_new = pick_alt_from_set(base_set, left_char, right_char, s, pos, motifs, max_run)
    candidates = base_set;
    best_score = inf;
    b_new = s(pos);

    for k = 1:numel(candidates)
        c = candidates(k);
        if c == s(pos)
            continue;
        end

        score = 0;
        if c == left_char
            score = score + 5;
        end
        if c == right_char
            score = score + 5;
        end

        s_tmp = s;
        s_tmp(pos) = c;
        local_seg = get_local_segment(s_tmp, pos, 12);

        local_max_run = get_local_max_run(local_seg);
        if local_max_run > max_run
            score = score + 50;
        end

        for t = 1:numel(motifs)
            if ~isempty(strfind(local_seg, upper(char(motifs{t})))) %#ok<STREMP>
                score = score + 20;
            end
        end

        if score < best_score
            best_score = score;
            b_new = c;
        end
    end
end

function seg = get_local_segment(s, pos, radius)
    L = max(1, pos-radius);
    R = min(numel(s), pos+radius);
    seg = s(L:R);
end

function mr = get_local_max_run(seg)
    runs = local_run_lengths(seg);
    mr = max(runs);
end