close all,  clc

% Import the image.
%original_image = imread('photos\baboon.jpg');
original_image = imread('photos\house.tiff');
subplot(2,4,1); imshow(original_image);title('original image')

% The complementary principle of DNA.
% 正确的DNA互补规则定义
dna_complementary_principle = [
    % 标准DNA互补：A-T, C-G
    84, 71, 67, 65;  % 0(A)->T(84), 1(C)->G(71), 2(G)->C(67), 3(T)->A(65)
    
    % 其他可能的互补规则
    84, 67, 71, 65;  % A->T, C->C, G->G, T->A
    71, 84, 65, 67;  % A->G, C->T, G->A, T->C
    67, 84, 65, 71;  % A->C, C->T, G->A, T->G
    84, 65, 71, 67;  % A->T, C->A, G->G, T->C
    71, 65, 84, 67;  % A->G, C->A, G->T, T->C
    67, 65, 84, 71;  % A->C, C->A, G->T, T->G
    65, 84, 67, 71;  % A->A, C->T, G->C, T->G
];
dna_complementary_principle = uint8(dna_complementary_principle);

%验证 DNA 编码/解码规则是否正确
verify_dna_rules(dna_complementary_principle);

%{
    ===============================================
                Encryption procedure
    ===============================================
%}
% Generate two pseudo-random sequences key_stream_diffusion and
% key_stream_scrambling, which are used to scramble the original image and 
% diffuse the scrambled image respectively.
[key_stream_diffusion, key_stream_scrambling, key_stream_dna] = generate_chaotic_quence_ten(original_image);

% 在加密过程前添加数据类型检查
fprintf('图像尺寸: %d x %d x %d, 数据类型: %s\n', size(original_image), class(original_image));
fprintf('扩散密钥流尺寸: %d x %d x %d, 数据类型: %s\n', size(key_stream_diffusion), class(key_stream_diffusion));
fprintf('置乱密钥流尺寸: %d x %d x %d, 数据类型: %s\n', size(key_stream_scrambling), class(key_stream_scrambling));
fprintf('DNA密钥流尺寸: %d x %d x %d, 数据类型: %s\n', size(key_stream_dna), class(key_stream_dna));

% 确保图像是uint8类型
if ~isa(original_image, 'uint8')
    original_image = uint8(original_image);
    fprintf('已将图像转换为uint8类型\n');
end

% Verify the reversibility of DNA encoding and decoding (small area test)
verify_dna_roundtrip(original_image, key_stream_dna, dna_complementary_principle);

% Scramble the original image
scrambled_image = scramble_img(original_image, key_stream_scrambling);
subplot(2,4,2); imshow(scrambled_image);title('scrambled image')

% Diffuse the scrambled image
diffused_image = diffuse_img(scrambled_image, key_stream_diffusion);
% diffused_image = cut_and_paste_attack(diffused_image);
% diffused_image = noisy_attack(diffused_image);
subplot(2,4,3); imshow(diffused_image);title('diffused image')

% Convert the diffused image with DNA computing
encoded_image = encode_img(diffused_image, dna_complementary_principle, key_stream_dna);
subplot(2,4,4); imshow(encoded_image);title('encoded image')

%{
    ===========================================================
      DNA constraints + oligo generation + FASTA/CSV export
    ===========================================================
%}
% 1) 先得到受约束 DNA 表示
[dna_constrained, dna_stats] = analyse_dna_constraints(encoded_image, [], 3, 5, 512);

% 如果 analyse_dna_constraints 没返回 char DNA 矩阵，就退回严格转换
if isempty(dna_constrained) || ~ischar(dna_constrained)
    dna_constrained = to_dna_char_matrix_strict(encoded_image);
end
dna_constrained = upper(char(dna_constrained));


% 2) 建立结果输出目录
out_dir = 'oligo_results';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% 3) 设置结构化 oligo 参数
oligo_params = struct();
oligo_params.payload_len  = 120;                        % payload长度
oligo_params.index_len    = 8;                          % index长度，4^8=65536，一般够用
oligo_params.ecc_len      = 8;                          % 先占位
oligo_params.primer_left  = 'ACGTACGTACGTACGTACGT';    % 20nt
oligo_params.primer_right = 'TGCATGCATGCATGCATGCA';    % 20nt

% 4) 生成结构化 oligos
[oligos, oligo_meta] = build_oligos(dna_constrained, oligo_params);

fprintf('\n===== Oligo Parameter Check =====\n');
fprintf('total_bases   = %d\n', oligo_meta.total_bases);
fprintf('num_oligos    = %d\n', oligo_meta.num_oligos);
fprintf('payload_len   = %d\n', oligo_params.payload_len);
fprintf('index_len     = %d\n', oligo_params.index_len);
fprintf('ecc_len       = %d\n', oligo_params.ecc_len);
fprintf('primer_left   = %d nt\n', length(oligo_params.primer_left));
fprintf('primer_right  = %d nt\n', length(oligo_params.primer_right));

oligo_lens = cellfun(@length, oligos);
fprintf('mean oligo length = %.2f nt\n', mean(oligo_lens));
fprintf('min oligo length  = %d nt\n', min(oligo_lens));
fprintf('max oligo length  = %d nt\n', max(oligo_lens));

fprintf('expected oligo length = %d nt\n', ...
    length(oligo_params.primer_left) + oligo_params.index_len + ...
    oligo_params.payload_len + oligo_params.ecc_len + ...
    length(oligo_params.primer_right));

fprintf('=================================\n');
% 5) 无噪声 round-trip 验证（证明结构化 oligo 可逆）
[dna_seq_recovered, dna_matrix_recovered] = oligos_to_dna(oligos, oligo_meta);

if isequal(dna_matrix_recovered, dna_constrained)
    fprintf('\n[OK] structured oligo round-trip 成功：无噪声下可完全恢复。\n');
else
    error('structured oligo round-trip 失败：恢复后的 DNA 矩阵与原始不一致。');
end

% 6) 导出前3条 oligo 示例，论文里可直接截图/展示
example_txt = fullfile(out_dir, 'oligo_examples.txt');
fid = fopen(example_txt, 'w');
for i = 1:min(3, numel(oligos))
    fprintf(fid, 'Oligo %d\n', i-1);
    fprintf(fid, '%s\n\n', oligos{i});
end
fclose(fid);

% 7) 导出 FASTA / CSV（直接导出 build_oligos 生成的结构化 176-nt oligo）
fasta_path = fullfile(out_dir, 'dna_cipher_pool.fasta');
csv_path   = fullfile(out_dir, 'dna_cipher_pool.csv');

oligo_len_list = cellfun(@length, oligos);
if numel(unique(oligo_len_list)) ~= 1
    error('导出的结构化 oligo 长度不一致，请检查 build_oligos。');
end

fid = fopen(fasta_path, 'w');
for i = 1:numel(oligos)
    fprintf(fid, '>OLI_%05d|LEN=%d\n', i-1, oligo_len_list(i));
    fprintf(fid, '%s\n', oligos{i});
end
fclose(fid);

T_export = table( ...
    (0:numel(oligos)-1).', ...
    string(oligos(:)), ...
    oligo_len_list(:), ...
    'VariableNames', {'index','sequence','length'} ...
);
writetable(T_export, csv_path);

export_summary = struct();
export_summary.num_oligos = numel(oligos);
export_summary.oligo_length = unique(oligo_len_list);
export_summary.payload_len = oligo_meta.params.payload_len;
export_summary.fasta = fasta_path;
export_summary.csv = csv_path;

primer_left_len  = length(oligo_params.primer_left);
index_len        = oligo_params.index_len;
payload_len      = oligo_params.payload_len;
ecc_len          = oligo_params.ecc_len;

payload_start = primer_left_len + index_len + 1;
payload_end   = primer_left_len + index_len + payload_len;

payload_only = cellfun(@(s) s(payload_start:payload_end), oligos, 'UniformOutput', false);

export_summary.gc_mean_full = mean(cellfun(@calc_gc_percent_local, oligos));
export_summary.gc_mean_payload = mean(cellfun(@calc_gc_percent_local, payload_only));

export_summary.homopolymer_max_full = max(cellfun(@max_homopolymer_local, oligos));
export_summary.homopolymer_max_payload = max(cellfun(@max_homopolymer_local, payload_only));

% 8) 打印结果摘要
fprintf('\n========== Oligo Summary ==========\n');
fprintf('DNA total bases           : %d\n', oligo_meta.total_bases);
fprintf('Structured oligo number   : %d\n', oligo_meta.num_oligos);
fprintf('Payload length per oligo  : %d nt\n', oligo_meta.params.payload_len);
fprintf('Index length              : %d nt\n', oligo_meta.params.index_len);
fprintf('ECC length                : %d nt\n', oligo_meta.params.ecc_len);
fprintf('Primer left length        : %d nt\n', length(oligo_meta.params.primer_left));
fprintf('Primer right length       : %d nt\n', length(oligo_meta.params.primer_right));

fprintf('\n========== Export Summary ==========\n');
fprintf('Exported oligo number     : %d\n', export_summary.num_oligos);
fprintf('Exported oligo length     : %d nt\n', export_summary.oligo_length);
fprintf('Payload length per oligo  : %d nt\n', export_summary.payload_len);
fprintf('Mean GC content full oligo : %.2f %%\n', export_summary.gc_mean_full);
fprintf('Mean GC content payload    : %.2f %%\n', export_summary.gc_mean_payload);
fprintf('Max homopolymer full oligo : %d\n', export_summary.homopolymer_max_full);
fprintf('Max homopolymer payload    : %d\n', export_summary.homopolymer_max_payload);
fprintf('FASTA file                : %s\n', export_summary.fasta);
fprintf('CSV file                  : %s\n', export_summary.csv);
fprintf('Example oligo txt         : %s\n', example_txt);

% 9) 保存mat，后面论文画图、做统计会更方便
save(fullfile(out_dir, 'oligo_workspace.mat'), ...
    'dna_constrained', 'dna_stats', ...
    'oligos', 'oligo_meta', ...
    'dna_seq_recovered', 'dna_matrix_recovered', ...
    'export_summary');

%{
    ===============================================
                Decryption procedure
    ===============================================
%}
% Decode the encoded image
decoded_image = decode_img(encoded_image, dna_complementary_principle, key_stream_dna);
subplot(2,4,5); imshow(decoded_image);title('decoded image')

% Verify intermediate step
diffused_uint8 = uint8(diffused_image);
decoded_uint8 = uint8(decoded_image);
diff_check = sum(abs(double(diffused_uint8(:)) - double(decoded_uint8(:))));
fprintf('Diffused vs Decoded difference: %d\n', diff_check);

% Dediffuse the diffused image
dediffused_image = dediffuse_img(decoded_image, key_stream_diffusion);
subplot(2,4,6); imshow(dediffused_image);title('dediffused image')

% Verify intermediate step  
scrambled_uint8 = uint8(scrambled_image);
dediffused_uint8 = uint8(dediffused_image);
diff_check2 = sum(abs(double(scrambled_uint8(:)) - double(dediffused_uint8(:))));
fprintf('Scrambled vs Dediffused difference: %d\n', diff_check2);

% Descramble the dediffused image.
descrambled_image = descramble_img(dediffused_image, key_stream_scrambling);
subplot(2,4,7); imshow(descrambled_image);title('descrambled image')

% Check if decryption is successful
if isequal(original_image, descrambled_image)
    disp('✅ Decryption successful: Original and decrypted images are identical');
else
    disp('❌ Decryption failed: Images are different');
    
    % Calculate difference
    diff = abs(double(original_image) - double(descrambled_image));
    error_rate = sum(diff(:) > 0) / numel(original_image);
    fprintf('Error rate: %.6f%%\n', error_rate * 100);
    
    % Find where differences occur
    [row, col, ch] = ind2sub(size(original_image), find(diff > 0));
    if ~isempty(row)
        fprintf('First difference at: row=%d, col=%d, channel=%d\n', row(1), col(1), ch(1));
        fprintf('Original value: %d, Decrypted value: %d\n', ...
                original_image(row(1), col(1), ch(1)), ...
                descrambled_image(row(1), col(1), ch(1)));
    end
end
%% 还原成标准的 A/C/G/T 矩阵再可视化
function Mout = to_dna_char_matrix(X, targetSize)
%TO_DNA_CHAR_MATRIX  Convert various representations into 'A','C','G','T' char matrix
% targetSize = [H W C]   (use size(original_image))

    H = targetSize(1); W = targetSize(2);
    C = 1; if numel(targetSize)>=3, C = targetSize(3); end
    total_bases = H*W*C;

    % Case A: already char A/C/G/T
    if ischar(X) || isstring(X)
        M = upper(char(X));
        if numel(M) ~= total_bases
            error('DNA矩阵元素数(%d)与目标尺寸所需碱基数(%d)不一致。', numel(M), total_bases);
        end
        if any(~ismember(M(:),'ACGT'))
            error('发现非 A/C/G/T 字符；这不是解码后的DNA矩阵。');
        end
        Mout = reshape(M, [W,H,C]); Mout = permute(Mout,[2 1 3]); % 列主到行主统一
        return;
    end

    % 强制为数值 uint8
    if ~isnumeric(X), X = uint8(X); end

    % Case B: 0..3 每元素一碱基
    if all(ismember(X(:), 0:3))
        map = 'ACGT';
        if numel(X) ~= total_bases
            error('0..3形式的元素数(%d)与目标碱基数(%d)不一致。', numel(X), total_bases);
        end
        Mout = reshape(map(double(X(:))+1), [W,H,C]); 
        Mout = permute(Mout,[2 1 3]);
        return;
    end

    % Case C: 打包字节（每个字节含4个2bit碱基，顺序：高位到低位：b7b6,b5b4,b3b2,b1b0）
    % 例如:  uint8 -> [d3 d2 d1 d0]  (每个 d∈{0,1,2,3})
    % 需要的字节数：
    need_bytes = ceil(total_bases/4);
    if numel(X) < need_bytes
        error('打包字节数不足：提供了 %d 字节，但至少需要 %d 字节。', numel(X), need_bytes);
    end
    bytes = uint8(X(:));
    bytes = bytes(1:need_bytes);

    % 解包：得到 4*need_bytes 个digits，取前 total_bases 个
    digs = zeros(need_bytes*4,1,'uint8');
    % 提取顺序： (b7b6)->d3, (b5b4)->d2, (b3b2)->d1, (b1b0)->d0
    d3 = bitshift(bitand(bytes, uint8(192)),-6); % 1100 0000
    d2 = bitshift(bitand(bytes, uint8(48)) , -4); % 0011 0000
    d1 = bitshift(bitand(bytes, uint8(12)) , -2); % 0000 1100
    d0 = bitand(bytes, uint8(3));                 % 0000 0011
    digs(1:4:end) = d3;
    digs(2:4:end) = d2;
    digs(3:4:end) = d1;
    digs(4:4:end) = d0;

    digs = digs(1:total_bases);
    map = 'ACGT';
    M = map(double(digs)+1);           % 行向量
    Mout = reshape(M, [W,H,C]); 
    Mout = permute(Mout,[2 1 3]);      % 到 HxWxC
end


%{
    ===============================================
                Analysis procedure
    ===============================================
%}

% Histogram analysis
figure;
analyse_histogram(encoded_image);

% Entropy analysis
figure;
entropy_arr = analyse_entropy(encoded_image);

% Correlation analysis
figure;
red_correlation = analyse_correlation(encoded_image(:,:,1), 5000);
green_correlation = analyse_correlation(encoded_image(:,:,2), 5000);
blue_correlation = analyse_correlation(encoded_image(:,:,3), 5000);

%correlation_3dir
figure;
[rH, rV, rD] = analyse_correlation_3dir(encoded_image);

% NPRC analysis
NPRC = analyse_NPRC(original_image, encoded_image);
fprintf('NPRC: %.6f\n', NPRC);

% UACI analysis
UACI = analyse_UACI(original_image, encoded_image);
fprintf('UACI: %.6f\n', UACI);

% GC content length homopolymer Forbidden motif analysis
% 前面已经完成 DNA constraint analysis，这里不要重新覆盖 dna_constrained
% [dna_constrained, dna_stats] = analyse_dna_constraints(encoded_image);

% DNA storage density
storage_stats = analyse_dna_storage_density(original_image, encoded_image);

% Table 1
% results_table1 = run_experiments_table1(original_image, dna_complementary_principle);

% npcr_uaci
plain2 = original_image;
plain2(1,1,:) = uint8( mod(double(plain2(1,1,:)) + 1, 256) );

scrambled2 = scramble_img(plain2, key_stream_scrambling);
diffused2  = diffuse_img(scrambled2, key_stream_diffusion);
encoded2   = encode_img(diffused2, dna_complementary_principle, key_stream_dna);

results_diff = analyse_uaci_nprc(encoded_image, encoded2, false, 255);
fprintf('【差分攻击意义】两张密文的 NPCR: %.6f\n', results_diff.npcr_score);
fprintf('【差分攻击意义】两张密文的 UACI: %.6f\n', results_diff.uaci_score);

% uaci_npcr
stats = analyse_npcr_uaci(encoded_image, encoded2);

% noisy
error_prob = 0.05;
[noisy_encoded_image, noise_stats] = analyse_dna_channel_noise(dna_constrained, error_prob);


% 选取一个“像素域密文/伪灰度密文”用于原图 vs 密文对比
cipher_like = [];

if exist('cipher_image','var')
    cipher_like = cipher_image;
elseif exist('diffused_image','var')
    cipher_like = diffused_image;
elseif exist('encoded_image','var')
    % 数字 0..3 → 伪灰度 0/85/170/255
    E = encoded_image;
    if isnumeric(E)
        cipher_like = uint8(mod(double(E),4) * 85);
    else
    error('Unsupported encoded_image type for visualization.');
    end
end

if ~isempty(cipher_like)
    stats_pc = analyse_plain_cipher_stats(original_image, cipher_like);
end


%{
    ===============================================
                Oligo 结构参数设置
    ===============================================
%}
% params_oligo = struct();
% params_oligo.payload_len  = 120;   % payload 长度 (nt)，你可以改成你想要的
% params_oligo.index_len    = 8;     % 8nt -> 4^8=65536 > 常见 oligo 数，足够
% params_oligo.ecc_len      = 8;     % 占位 ECC 字段
% params_oligo.primer_left  = 'ACGTACGTACGTACGTACGT';    % 20nt 占位
% params_oligo.primer_right = 'TGCATGCATGCATGCATGCA';    % 20nt 占位
% 
% % 1) dna_constrained -> oligos 
% [oligos, oligo_meta] = build_oligos(dna_constrained, params_oligo);
% 
% % 2) oligos -> dna_rec（无噪声）
% dna_rec_seq = oligos_to_dna(oligos, oligo_meta);
% dna_rec     = reshape(dna_rec_seq, size(dna_constrained));
% 
% % 检查是否完全一致
% diff_count = sum(dna_constrained(:) ~= dna_rec(:));
% fprintf('[Oligo无噪声测试] dna_constrained vs dna_rec 不同碱基数 = %d\n', diff_count);


% 关闭空白图窗（无 axes 或 axes 下无子对象）
allF = findall(0,'Type','figure');
for f = allF.'
    ax = findall(f,'Type','axes');
    if isempty(ax) || all(arrayfun(@(a) isempty(get(a,'Children')), ax))
        close(f);
    end
end

function gc = calc_gc_percent_local(seq)
    seq = upper(char(seq));
    gc = 100 * sum(seq == 'G' | seq == 'C') / numel(seq);
end

function maxhp = max_homopolymer_local(seq)
    seq = upper(char(seq));
    if isempty(seq)
        maxhp = 0;
        return;
    end

    maxhp = 1;
    current_run = 1;

    for i = 2:numel(seq)
        if seq(i) == seq(i-1)
            current_run = current_run + 1;
        else
            maxhp = max(maxhp, current_run);
            current_run = 1;
        end
    end

    maxhp = max(maxhp, current_run);
end

function dna_seq = local_flatten_row_major_check(dna_matrix)
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

function m = local_max_run_check(seq)
    s = char(seq);
    if isempty(s)
        m = 0;
        return;
    end

    cur = 1;
    m = 1;

    for i = 2:length(s)
        if s(i) == s(i-1)
            cur = cur + 1;
        else
            m = max(m, cur);
            cur = 1;
        end
    end

    m = max(m, cur);
end