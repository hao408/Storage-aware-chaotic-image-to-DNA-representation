function [dna_seq, dna_matrix] = oligos_to_dna(oligos, meta)
%OLIGOS_TO_DNA  结构化 oligo 的逆变换（无噪声版）
%
%   [dna_seq, dna_matrix] = oligos_to_dna(oligos, meta)
%
%   oligos : cell(num_oligos,1) 或 string 数组，每条为 char 行向量
%   meta   : build_oligos 返回的结构体（含 total_bases / params / original_size）
%
%   输出：
%       dna_seq    : 1×total_bases char, 'A','C','G','T'
%       dna_matrix : 按原尺寸恢复后的 DNA char 矩阵

    if iscell(oligos)
        oligos = string(oligos);
    else
        oligos = string(oligos);
    end

    params       = meta.params;
    num_oligos   = numel(oligos);
    payload_len  = params.payload_len;
    index_len    = params.index_len;
    ecc_len      = params.ecc_len;
    primer_left  = upper(char(params.primer_left));
    primer_right = upper(char(params.primer_right));
    LpL          = length(primer_left);
    LpR          = length(primer_right);

    indices  = zeros(num_oligos,1);
    payloads = cell(num_oligos,1);

    for k = 1:num_oligos
        s = upper(char(oligos(k)));

        core_len_expect = index_len + payload_len + ecc_len;
        if length(s) < (LpL + core_len_expect + LpR)
            error('oligo %d 长度不足以拆出 primer/index/payload/ecc。', k);
        end

        % 去掉两端 primer
        core = s(LpL+1 : end-LpR);

        % 切 index / payload / ecc
        idx_nt     = core(1:index_len);
        payload_nt = core(index_len+1 : index_len+payload_len);
        % ecc_nt   = core(index_len+payload_len+1 : index_len+payload_len+ecc_len);

        idx_int = decode_index_nt(idx_nt);
        indices(k) = idx_int;
        payloads{k} = payload_nt;
    end

    % 按 index 排序，恢复原始顺序
    [~, order] = sort(indices);
    payloads_sorted = payloads(order);

    % 拼接所有 payload，截断回 total_bases
    concat = [payloads_sorted{:}];
    dna_seq = concat(1:meta.total_bases);

    % 恢复成原始 DNA 矩阵
    if isfield(meta, 'original_size')
        if isfield(meta, 'flatten_order') && strcmp(meta.flatten_order, 'column-major')
            dna_matrix = reshape(dna_seq, meta.original_size);
        else
            dna_matrix = unflatten_dna_row_major(dna_seq, meta.original_size);
        end
    else
        dna_matrix = dna_seq;
    end
end

%================== 辅助函数: index 解码 ==================%
function idx = decode_index_nt(index_nt)
%DECODE_INDEX_NT  将 A/C/G/T 序列还原为整数 idx（base-4 解码）

    index_nt = upper(char(index_nt));
    n = length(index_nt);
    digits = zeros(1,n);

    for i = 1:n
        switch index_nt(i)
            case 'A'
                digits(i) = 0;
            case 'C'
                digits(i) = 1;
            case 'G'
                digits(i) = 2;
            case 'T'
                digits(i) = 3;
            otherwise
                error('decode_index_nt: 非法碱基 %c（不是 A/C/G/T）。', index_nt(i));
        end
    end

    idx = 0;
    for i = 1:n
        idx = idx * 4 + digits(i);
    end
end

%================== 辅助函数: row-major 逆展平 ==================%
function dna_matrix = unflatten_dna_row_major(dna_seq, original_size)

    if numel(original_size) == 2
        H = original_size(1);
        W = original_size(2);
        dna_matrix = reshape(dna_seq, [W, H]).';
    elseif numel(original_size) == 3
        H = original_size(1);
        W = original_size(2);
        C = original_size(3);

        dna_matrix = repmat('A', H, W, C);
        pos = 1;
        for c = 1:C
            block = dna_seq(pos : pos + H*W - 1);
            dna_matrix(:,:,c) = reshape(block, [W, H]).';
            pos = pos + H*W;
        end
    else
        error('original_size 维度不支持。');
    end
end