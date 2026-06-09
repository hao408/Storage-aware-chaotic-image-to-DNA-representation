function stats = analyse_dna_gc_length(dna_in, seq_len)
% analyse_dna_gc_length
% 统计 DNA 序列/矩阵的 GC 含量。
%
% 当前版本兼容两类输入：
%   1) char DNA矩阵，例如 dna_constrained，元素为 A/C/G/T；
%   2) uint8/double encoded_image，此时仍按 strict 规则转为 DNA。
%
% 注意：
%   对 char DNA 矩阵，采用 MATLAB column-major 顺序展开，
%   与 build_oligos.m 当前导出 payload 的顺序保持一致。

    if nargin < 2
        seq_len = 512;
    end

    % ===== 1. 输入转为 DNA 字符序列 =====
    if ischar(dna_in) || isstring(dna_in)
        dna_char = upper(char(dna_in));
    else
        [H, W, C] = size(dna_in);
        targetSize = [H, W*4, C];
        dna_char = to_dna_char_matrix_strict(dna_in, targetSize);
        dna_char = upper(char(dna_char));
    end

    % column-major 展平，必须与 build_oligos.m 一致
    dna_seq = reshape(dna_char, 1, []);
    total_bases = numel(dna_seq);

    % ===== 2. 按 seq_len 分块统计 =====
    num_sequences = ceil(total_bases / seq_len);
    pad_len = num_sequences * seq_len - total_bases;

    if pad_len > 0
        dna_seq_pad = [dna_seq, repmat('A', 1, pad_len)];
    else
        dna_seq_pad = dna_seq;
    end

    seq_mat = reshape(dna_seq_pad, seq_len, num_sequences).';

    gc_per_seq = zeros(num_sequences, 1);
    len_per_seq = zeros(num_sequences, 1);

    for i = 1:num_sequences
        s = seq_mat(i, :);

        if i == num_sequences && pad_len > 0
            s = s(1:end-pad_len);
        end

        len_per_seq(i) = numel(s);
        gc_per_seq(i) = 100 * sum(s == 'G' | s == 'C') / max(1, numel(s));
    end

    % ===== 3. 全局 GC =====
    nA = sum(dna_seq == 'A');
    nC = sum(dna_seq == 'C');
    nG = sum(dna_seq == 'G');
    nT = sum(dna_seq == 'T');

    gc_global = 100 * (nG + nC) / max(1, total_bases);

    gc_mean_seq = mean(gc_per_seq);
    gc_std_seq  = std(gc_per_seq);

    mean_len = mean(len_per_seq);
    min_len  = min(len_per_seq);
    max_len  = max(len_per_seq);

    % ===== 4. 画 GC 直方图 =====
    figure;
    histogram(gc_per_seq, 20);
    xlabel('GC% per sequence');
    ylabel('Count');
    title(sprintf('GC%% per sequence (mean=%.2f%%, std=%.2f%%)', ...
        gc_mean_seq, gc_std_seq));
    grid on;

    % ===== 5. 输出 stats =====
    stats = struct();
    stats.total_bases         = total_bases;
    stats.num_sequences       = num_sequences;
    stats.mean_seq_length     = mean_len;
    stats.min_seq_length      = min_len;
    stats.max_seq_length      = max_len;
    stats.gc_percent_global   = gc_global;
    stats.gc_percent_mean_seq = gc_mean_seq;
    stats.gc_percent_std_seq  = gc_std_seq;
    stats.count_A             = nA;
    stats.count_C             = nC;
    stats.count_G             = nG;
    stats.count_T             = nT;

    fprintf('\n[analyse_dna_gc_length - column-major]\n');
    fprintf('  Total bases          : %d\n', total_bases);
    fprintf('  Num sequences        : %d\n', num_sequences);
    fprintf('  Seq length (mean/min/max): %.1f / %d / %d\n', ...
        mean_len, min_len, max_len);
    fprintf('  GC global            : %.2f %%\n', gc_global);
    fprintf('  GC per-seq (mean±std): %.2f ± %.2f %%\n', gc_mean_seq, gc_std_seq);
end