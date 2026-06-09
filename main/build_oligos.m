function [oligos, meta] = build_oligos(dna_constrained, params)
%BUILD_OLIGOS  结构化 oligo: primerL | index | payload | ecc | primerR
%
%   [oligos, meta] = build_oligos(dna_constrained, params)
%
%   dna_constrained : H×W×C 的 char 矩阵，元素为 'A','C','G','T'
%
%   params 结构体字段（需要预先设置）：
%       .payload_len   每条 oligo 中 payload 的长度 (nt)
%       .index_len     index 长度 (nt)，需满足 4^index_len > num_oligos
%       .ecc_len       ECC 字段长度 (nt)，当前先占位，用 'A' 填充
%       .primer_left   左引物序列 (char 行向量)
%       .primer_right  右引物序列 (char 行向量)
%
%   输出：
%       oligos : cell(num_oligos,1)，每个元素为 char 行向量（一条 oligo）
%       meta   : 结构体，记录 total_bases / num_oligos / params 等信息

    %---------------- 参数 & 输入检查 ----------------%
    if ~ischar(dna_constrained) && ~isstring(dna_constrained)
        error('dna_constrained 必须是 char 或 string 类型，且元素为 A/C/G/T。');
    end

    dna_constrained = upper(char(dna_constrained));

    if any(~ismember(dna_constrained(:), 'ACGT'))
        error('dna_constrained 中存在非 A/C/G/T 字符。');
    end

    required_fields = {'payload_len','index_len','ecc_len','primer_left','primer_right'};
    for k = 1:numel(required_fields)
        f = required_fields{k};
        if ~isfield(params, f)
            error('params.%s 未设置。', f);
        end
    end

    payload_len  = params.payload_len;
    index_len    = params.index_len;
    ecc_len      = params.ecc_len;
    primer_left  = upper(char(params.primer_left));
    primer_right = upper(char(params.primer_right));

    if any(~ismember(primer_left, 'ACGT')) || any(~ismember(primer_right, 'ACGT'))
        error('primer_left / primer_right 必须只包含 A/C/G/T。');
    end

    %---------------- 统一按 row-major 展平成一维 DNA 串 ----------------%
    dna_seq = reshape(dna_constrained, 1, []);
    total_bases = numel(dna_seq);

    % 需要多少条 oligo
    num_oligos = ceil(total_bases / payload_len);

    % 检查 index_len 是否足够编码 num_oligos
    if 4^index_len < num_oligos
        error('index_len=%d 太短，不能编码 %d 条 oligos（需要 4^index_len >= num_oligos）。', ...
              index_len, num_oligos);
    end

    %---------------- 构造每条 oligo ----------------%
    oligos = cell(num_oligos, 1);
    pos = 1;
    pad_pattern = 'ACGT'; ;  % 最后一块不满 payload_len 时的填充值

    for k = 1:num_oligos
        tail = min(pos + payload_len - 1, total_bases);
        payload = dna_seq(pos:tail);

        if numel(payload) < payload_len
            pad_len = payload_len - numel(payload);
            pad_seq = repmat(pad_pattern, 1, ceil(pad_len / numel(pad_pattern)));
            pad_seq = pad_seq(1:pad_len);
            payload = [payload, pad_seq];
        end

        % index: 用 base-4 编码 (k-1)
        index_nt = encode_index_nt(k-1, index_len);

        % ECC: 当前先占位
        ecc_pattern = 'ACGT';
        ecc_nt = repmat(ecc_pattern, 1, ceil(ecc_len / numel(ecc_pattern)));
        ecc_nt = ecc_nt(1:ecc_len);

        % 拼接 oligo: primerL | index | payload | ecc | primerR
        oligos{k} = [primer_left, index_nt, payload, ecc_nt, primer_right];

        pos = pos + payload_len;
    end

    %---------------- 记录 meta 信息 ----------------%
    meta = struct();
    meta.total_bases   = total_bases;
    meta.num_oligos    = num_oligos;
    meta.params        = params;
    meta.pad_pattern   = pad_pattern;
    meta.original_size = size(dna_constrained);
    meta.flatten_order = 'column-major';
end

%================== 辅助函数: index 编码 ==================%
function index_nt = encode_index_nt(idx, index_len)
%ENCODE_INDEX_NT  将非负整数 idx 编码为长度为 index_len 的 DNA 序列 (A/C/G/T).
%   使用 base-4: 0->A, 1->C, 2->G, 3->T.

    if idx < 0
        error('encode_index_nt: idx 必须为非负整数。');
    end

    digits = zeros(1, index_len);  % base-4 digits
    v = idx;
    for i = index_len:-1:1
        digits(i) = mod(v, 4);
        v = floor(v / 4);
    end

    if v > 0
        warning('encode_index_nt: idx 超出了 index_len 所能表示的范围，将被截断。');
    end

    map = 'ACGT';   % 0->A,1->C,2->G,3->T
    index_nt = map(digits + 1);
end

%================== 辅助函数: row-major 展平 ==================%
function dna_seq = flatten_dna_row_major(dna_matrix)
    dna = upper(char(dna_matrix));

    if ndims(dna) == 3
        [H, W, C] = size(dna);
        buf = cell(1, C);
        for c = 1:C
            buf{c} = reshape(dna(:,:,c).', 1, H * W);
        end
        dna_seq = [buf{:}];
    else
        [H, W] = size(dna);
        dna_seq = reshape(dna.', 1, H * W);
    end
end