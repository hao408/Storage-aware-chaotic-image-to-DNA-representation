function [noisy_dna, noise_stats] = analyse_dna_channel_noise(dna_in, error_prob)
% analyse_dna_channel_noise
% 模拟 DNA 通道替换噪声。
%
% 当前推荐输入：
%   dna_in = dna_constrained
%
% 说明：
%   如果输入是 A/C/G/T 字符矩阵，则直接在 DNA 碱基上加 substitution noise。
%   统计时调用 analyse_dna_gc_length 和 analyse_homopolymer_stats，
%   二者均使用 column-major 顺序，与 build_oligos.m 一致。

    if nargin < 2
        error_prob = 0.05;
    end

    % ===== 1. 判断输入类型 =====
    if ischar(dna_in) || isstring(dna_in)
        dna_char = upper(char(dna_in));
        noisy_dna = dna_char;

        bases = 'ACGT';

        noise_mask = rand(size(noisy_dna)) < error_prob;
        noise_pos = find(noise_mask);

        for k = 1:numel(noise_pos)
            pos = noise_pos(k);
            old_base = noisy_dna(pos);

            candidates = bases(bases ~= old_base);
            noisy_dna(pos) = candidates(randi(numel(candidates)));
        end

    else
        % 如果仍然传入 encoded_image，则保留原来的数值型扰动逻辑
        noisy_dna = dna_in;
        [H, W, C] = size(dna_in);

        for i = 1:H
            for j = 1:W
                for c = 1:C
                    if rand() < error_prob
                        original_base = noisy_dna(i, j, c);
                        new_base = randi([0, 3], 1, 1);

                        while new_base == original_base
                            new_base = randi([0, 3], 1, 1);
                        end

                        noisy_dna(i, j, c) = new_base;
                    end
                end
            end
        end
    end

    % ===== 2. 统计噪声前后 GC 和 homopolymer =====
    noise_stats = struct();
    noise_stats.noise_level = error_prob;

    noise_stats.gc_before = analyse_dna_gc_length(dna_in);
    noise_stats.gc_after  = analyse_dna_gc_length(noisy_dna);

    noise_stats.homo_stats_before = analyse_homopolymer_stats(dna_in);
    noise_stats.homo_stats_after  = analyse_homopolymer_stats(noisy_dna);

    noise_stats.gc_percent_before = noise_stats.gc_before.gc_percent_global;
    noise_stats.gc_percent_after  = noise_stats.gc_after.gc_percent_global;

    % ===== 3. 打印摘要 =====
    fprintf('\n[analyse_dna_channel_noise]\n');
    fprintf('Noise level: %.2f\n', error_prob);
    fprintf('GC before noise: %.2f%%\n', noise_stats.gc_percent_before);
    fprintf('GC after noise : %.2f%%\n', noise_stats.gc_percent_after);
    fprintf('Homopolymer before noise: Max run = %d, Mean run = %.2f\n', ...
        noise_stats.homo_stats_before.max_run_overall, ...
        noise_stats.homo_stats_before.mean_run_length);
    fprintf('Homopolymer after noise : Max run = %d, Mean run = %.2f\n', ...
        noise_stats.homo_stats_after.max_run_overall, ...
        noise_stats.homo_stats_after.mean_run_length);
end