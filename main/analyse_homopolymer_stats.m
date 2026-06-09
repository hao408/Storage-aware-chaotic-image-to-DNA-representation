function stats = analyse_homopolymer_stats(dna_in, seq_len)
% analyse_homopolymer_stats
% 统计 DNA 序列/矩阵的 homopolymer run length。
%
% 当前版本兼容：
%   1) char DNA矩阵，例如 dna_constrained；
%   2) uint8/double encoded_image。
%
% 对 char DNA 矩阵，采用 column-major 顺序展开，
% 与 build_oligos.m 当前导出 payload 的顺序一致。

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

    % ===== 2. 按 seq_len 分块 =====
    num_sequences = ceil(total_bases / seq_len);
    pad_len = num_sequences * seq_len - total_bases;

    if pad_len > 0
        dna_seq_pad = [dna_seq, repmat('A', 1, pad_len)];
    else
        dna_seq_pad = dna_seq;
    end

    seq_mat = reshape(dna_seq_pad, seq_len, num_sequences).';

    all_run_lengths = [];
    max_run_per_seq = zeros(num_sequences, 1);

    for i = 1:num_sequences
        s = seq_mat(i, :);

        if i == num_sequences && pad_len > 0
            s = s(1:end-pad_len);
        end

        [max_run_i, runs_i] = local_run_stats(s);
        max_run_per_seq(i) = max_run_i;
        all_run_lengths = [all_run_lengths, runs_i]; %#ok<AGROW>
    end

    if isempty(all_run_lengths)
        mean_run_len   = NaN;
        median_run_len = NaN;
        max_overall    = 0;
    else
        mean_run_len   = mean(all_run_lengths);
        median_run_len = median(all_run_lengths);
        max_overall    = max(all_run_lengths);
    end

    max_run_per_seq_mean = mean(max_run_per_seq);
    max_run_per_seq_max  = max(max_run_per_seq);

    % ===== 3. 画 run length 直方图 =====
    figure;
    if max_overall <= 10
        edges = 0.5:1:(max_overall + 0.5);
    else
        edges = 0.5:1:10.5;
    end

    histogram(all_run_lengths, edges);
    xlabel('Run length (same base)');
    ylabel('Count');

    if max_overall <= 10
        title(sprintf('Homopolymer run lengths (max=%d, mean=%.2f)', ...
            max_overall, mean_run_len));
    else
        title(sprintf('Homopolymer run lengths (max=%d, mean=%.2f; displayed 1-10)', ...
            max_overall, mean_run_len));
    end

    grid on;

    % ===== 4. 输出 stats =====
    stats = struct();
    stats.total_bases            = total_bases;
    stats.num_sequences          = num_sequences;
    stats.max_run_overall        = max_overall;
    stats.mean_run_length        = mean_run_len;
    stats.median_run_length      = median_run_len;
    stats.max_run_per_seq_mean   = max_run_per_seq_mean;
    stats.max_run_per_seq_max    = max_run_per_seq_max;

    fprintf('\n[analyse_homopolymer_stats - column-major]\n');
    fprintf('  Total bases               : %d\n', total_bases);
    fprintf('  Num sequences             : %d\n', num_sequences);
    fprintf('  Max run (overall)         : %d\n', max_overall);
    fprintf('  Run length (mean/median)  : %.2f / %.2f\n', ...
        mean_run_len, median_run_len);
    fprintf('  Max run per-seq (mean/max): %.2f / %d\n', ...
        max_run_per_seq_mean, max_run_per_seq_max);
end


function [max_run, runs] = local_run_stats(seq)
    seq = upper(char(seq));

    if isempty(seq)
        max_run = 0;
        runs = [];
        return;
    end

    cur = 1;
    max_run = 1;
    runs = [];

    for k = 2:numel(seq)
        if seq(k) == seq(k-1)
            cur = cur + 1;
        else
            runs(end+1) = cur; %#ok<AGROW>
            max_run = max(max_run, cur);
            cur = 1;
        end
    end

    runs(end+1) = cur;
    max_run = max(max_run, cur);
end