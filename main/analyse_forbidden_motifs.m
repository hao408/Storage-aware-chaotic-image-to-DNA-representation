function stats = analyse_forbidden_motifs(encoded_image, motifs)
%ANALYSE_FORBIDDEN_MOTIFS  检测 DNA 序列中的禁忌 motif.
%
%   stats = ANALYSE_FORBIDDEN_MOTIFS(encoded_image, motifs)
%
%   encoded_image : H×W 或 H×W×C 的密文图像（uint8 或 'A''C''G''T'）
%   motifs        : cell 数组，包含若干 DNA 短序列字符串，例如：
%                   {'AAAAAA','CCCCCC','GCGC','CGCG','GAATTC'}
%                   若省略，使用一个默认集合。
%
%   输出 stats 结构体：
%       .motifs            1×K cell，每个元素是 motif 字符串
%       .total_occ         1×K double，总出现次数
%       .seq_with_motif    1×K double，有至少一个该 motif 的序列条数
%       .seq_with_motif_frac 1×K double，上述比例
%       .max_occ_per_seq   1×K double，单条序列中最多出现次数
%
%   并在命令行打印报告.

    fprintf('[analyse_forbidden_motifs]\n');

    if nargin < 2 || isempty(motifs)
        % 一个示例性的默认禁忌集合，你可按需要修改/扩展
        motifs = {'AAAAAA', 'CCCCCC', 'GGGGGG', 'TTTTTT', ...
                  'GCGC', 'CGCG', ...
                  'GAATTC', 'AAGCTT'};  % 常见限制性酶位点 EcoRI/HindIII
    end

    % 1) 转成 DNA 字符矩阵
    dna = local_to_dna_chars(encoded_image);
    [H, W, C] = size(dna);

    num_sequences = H * C;
    seq_length    = W;

    fprintf('  Num sequences(H*C) : %d\n', num_sequences);
    fprintf('  Seq length         : %d\n', seq_length);

    % 2) 把每一行/通道当作一条序列
    sequences = cell(num_sequences,1);
    idx = 0;
    for ch = 1:C
        for i = 1:H
            idx = idx + 1;
            sequences{idx} = dna(i,:,ch);   % 1×W char row
        end
    end

    K = numel(motifs);
    total_occ          = zeros(1,K);
    seq_with_motif     = zeros(1,K);
    max_occ_per_seq    = zeros(1,K);

    % 3) 针对每个 motif 做统计
    for k = 1:K
        motif = upper(char(motifs{k}));
        occ_per_seq = zeros(num_sequences,1);

        for s = 1:num_sequences
            seq = sequences{s};
            % 在这一条序列中查找 motif 的所有出现位置
            idxs = strfind(seq, motif);
            occ_per_seq(s) = numel(idxs);
        end

        total_occ(k)       = sum(occ_per_seq);
        seq_with_motif(k)  = sum(occ_per_seq > 0);
        max_occ_per_seq(k) = max(occ_per_seq);
    end

    seq_with_motif_frac = seq_with_motif / num_sequences;

    % 4) 打印摘要
    for k = 1:K
        fprintf('  Motif "%s": total_occ=%d, seq_with_motif=%d (%.2f%%), max_occ_per_seq=%d\n', ...
            motifs{k}, total_occ(k), seq_with_motif(k), ...
            seq_with_motif_frac(k)*100, max_occ_per_seq(k));
    end

    % 填 stats
    stats = struct();
    stats.motifs              = motifs;
    stats.total_occ           = total_occ;
    stats.seq_with_motif      = seq_with_motif;
    stats.seq_with_motif_frac = seq_with_motif_frac;
    stats.max_occ_per_seq     = max_occ_per_seq;
end

% ======= 与前一个文件类似的辅助函数 =======
function dna = local_to_dna_chars(X)
    if isstring(X) || ischar(X)
        dna = upper(char(X));
        if any(~ismember(dna(:), 'ACGT'))
            error('DNA 矩阵中存在非 A/C/G/T 字符。');
        end
        return;
    end
    X = double(X);
    v = mod(X, 4);
    map = 'ACGT';
    dna = map(v + 1);
end
