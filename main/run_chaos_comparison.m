clc; clear; close all;

%% =========================================================
%  Multi-chaotic-system comparison experiment
%  对比：
%  1) CIS-T / ten
%  2) Logistic
%  3) Henon
%  4) 6D
%
%  输出：
%  1) comparison_results/chaos_comparison_results.csv
%  2) comparison_results/chaos_comparison_results.txt
%  3) comparison_results/chaos_comparison_workspace.mat
%
%  注意：
%  本脚本已经改为使用 build_oligos，与主程序 oligo 结构保持一致：
%  primerL(20) | index(8) | payload(120) | ecc(8) | primerR(20)
%  Total oligo length = 176 nt
%% =========================================================

%% =========================
%  输出目录
%% =========================
out_dir = 'comparison_results';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% =========================
%  读取测试图像
%  如果你想换图，就改这里
%% =========================
original_image = imread('photos\house.tiff');
if ~isa(original_image, 'uint8')
    original_image = uint8(original_image);
end

%% =========================
%  DNA互补规则（与你主程序一致）
%% =========================
dna_complementary_principle = [
    84, 71, 67, 65;
    84, 67, 71, 65;
    71, 84, 65, 67;
    67, 84, 65, 71;
    84, 65, 71, 67;
    71, 65, 84, 67;
    67, 65, 84, 71;
    65, 84, 67, 71;
];
dna_complementary_principle = uint8(dna_complementary_principle);

%% =========================
%  Oligo 参数：必须与主程序 build_oligos 保持一致
%% =========================
oligo_params = struct();
oligo_params.payload_len  = 120;
oligo_params.index_len    = 8;
oligo_params.ecc_len      = 8;

% 下面两条引物序列必须是 20 nt。
% 如果你主程序里用的不是这两条，请替换成主程序里的真实 primer_left / primer_right。
oligo_params.primer_left  = 'ACGTACGTACGTACGTACGT';   % 20 nt
oligo_params.primer_right = 'TGCATGCATGCATGCATGCA';   % 20 nt

%% =========================
%  方法列表
%% =========================
method_names = {'CIS-T', 'Logistic', 'Henon', '6D'};
num_methods = numel(method_names);

%% =========================
%  结果初始化
%% =========================
Entropy_list   = zeros(num_methods,1);
CorrH_list     = zeros(num_methods,1);
CorrV_list     = zeros(num_methods,1);
CorrD_list     = zeros(num_methods,1);
NPCR_list      = zeros(num_methods,1);
UACI_list      = zeros(num_methods,1);
GCmean_list    = zeros(num_methods,1);
GCrange_list   = strings(num_methods,1);
HPmean_list    = zeros(num_methods,1);
HPmax_list     = zeros(num_methods,1);
OligoNum_list  = zeros(num_methods,1);
Runtime_list   = zeros(num_methods,1);

%% =========================
%  主循环
%% =========================
for m = 1:num_methods
    method_name = method_names{m};

    fprintf('\n====================================================\n');
    fprintf('Running method: %s\n', method_name);
    fprintf('====================================================\n');

    tStart = tic;

    %% 1. 生成密钥流
    [key_stream_diffusion, key_stream_scrambling, key_stream_dna] = ...
        get_chaotic_streams(method_name, original_image);

    %% 2. 加密主流程
    scrambled_image = scramble_img(original_image, key_stream_scrambling);
    diffused_image  = diffuse_img(scrambled_image, key_stream_diffusion);
    encoded_image   = encode_img(diffused_image, dna_complementary_principle, key_stream_dna);

    %% 3. 安全性指标
    % 3.1 熵
    Entropy_list(m) = calc_entropy_local(encoded_image);

    % 3.2 三方向相关性
    [corrH, corrV, corrD] = calc_correlation_3dir_local(encoded_image, 5000);
    CorrH_list(m) = corrH;
    CorrV_list(m) = corrV;
    CorrD_list(m) = corrD;

    % 3.3 NPCR / UACI
    [npcr_val, uaci_val] = calc_npcr_uaci_local(original_image, encoded_image);
    NPCR_list(m) = npcr_val;
    UACI_list(m) = uaci_val;

    %% 4. DNA约束分析
    % 先尽量调用你原有 analyse_dna_constraints
    try
        [dna_constrained, ~] = analyse_dna_constraints(encoded_image);
        if isempty(dna_constrained) || ~ischar(dna_constrained)
            dna_constrained = to_dna_char_matrix_strict(encoded_image);
        end
    catch
        dna_constrained = to_dna_char_matrix_strict(encoded_image);
    end

    dna_constrained = upper(char(dna_constrained));

    % 4.1 统计整体GC
    dna_seq = flatten_dna_row_major(dna_constrained);
    gc_all = calc_gc_percent_local(dna_seq);

    %% 4.2 使用 build_oligos 生成 oligo，并统计 oligo 级指标
    % 这一段替代原来的 export_oligos，保证和主程序一致。
    export_subdir = fullfile(out_dir, method_name);
    if ~exist(export_subdir, 'dir')
        mkdir(export_subdir);
    end

    [oligos, oligo_meta] = build_oligos(dna_constrained, oligo_params);

    % 检查当前 oligo 参数，确保和主程序一致
    fprintf('\n===== Oligo Parameter Check: %s =====\n', method_name);
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
    fprintf('=====================================\n');

    % 逐条 oligo 统计 GC 和最大 homopolymer
    num_oligos = numel(oligos);
    gc_percent = zeros(num_oligos, 1);
    max_hp = zeros(num_oligos, 1);
    length_nt = zeros(num_oligos, 1);

    for oi = 1:num_oligos
        seq = upper(char(oligos{oi}));
        length_nt(oi) = length(seq);
        gc_percent(oi) = calc_gc_percent_local(seq);
        max_hp(oi) = calc_max_homopolymer_local(seq);
    end

    GCmean_list(m) = mean(gc_percent);
    GCrange_list(m) = sprintf('%.2f-%.2f', min(gc_percent), max(gc_percent));
    HPmean_list(m) = mean(max_hp);
    HPmax_list(m) = max(max_hp);
    OligoNum_list(m) = num_oligos;

    % 保存每种方法对应的 oligo 明细，方便后续画图和核对
    oligo_id = (0:num_oligos-1).';
    sequence = string(oligos);
    T = table(oligo_id, sequence, length_nt, gc_percent, max_hp, ...
        'VariableNames', {'oligo_id','sequence','length_nt','gc_percent','max_homopolymer'});

    csv_file = fullfile(export_subdir, 'dna_cipher_pool.csv');
    writetable(T, csv_file);

    %% 5. 运行时间
    Runtime_list(m) = toc(tStart);

    %% 6. 打印
    fprintf('Entropy              = %.6f\n', Entropy_list(m));
    fprintf('Corr-H / V / D       = %.6f / %.6f / %.6f\n', CorrH_list(m), CorrV_list(m), CorrD_list(m));
    fprintf('NPCR / UACI          = %.6f / %.6f\n', NPCR_list(m), UACI_list(m));
    fprintf('Overall DNA GC       = %.6f\n', gc_all);
    fprintf('Mean GC              = %.6f\n', GCmean_list(m));
    fprintf('GC range             = %s\n', GCrange_list(m));
    fprintf('Mean max homopolymer = %.6f\n', HPmean_list(m));
    fprintf('Max homopolymer      = %d\n', HPmax_list(m));
    fprintf('Number of oligos     = %d\n', OligoNum_list(m));
    fprintf('Runtime (s)          = %.4f\n', Runtime_list(m));
end

%% =========================
%  汇总表
%% =========================
ResultTable = table( ...
    string(method_names(:)), ...
    Entropy_list, ...
    CorrH_list, CorrV_list, CorrD_list, ...
    NPCR_list, UACI_list, ...
    GCmean_list, GCrange_list, ...
    HPmean_list, HPmax_list, ...
    OligoNum_list, Runtime_list, ...
    'VariableNames', ...
    {'Method','Entropy','CorrH','CorrV','CorrD','NPCR','UACI', ...
     'MeanGC','GCRange','MeanMaxHomopolymer','MaxHomopolymer', ...
     'NumOligos','RuntimeSec'});

disp(' ');
disp('================ Chaos Comparison Results ================');
disp(ResultTable);
disp('=========================================================');

%% =========================
%  保存 CSV
%% =========================
csv_out = fullfile(out_dir, 'chaos_comparison_results.csv');
writetable(ResultTable, csv_out);

%% =========================
%  保存 TXT
%% =========================
txt_out = fullfile(out_dir, 'chaos_comparison_results.txt');
fid = fopen(txt_out, 'w');

fprintf(fid, '================ Chaos Comparison Results ================\n');
fprintf(fid, 'Method\tEntropy\tCorrH\tCorrV\tCorrD\tNPCR\tUACI\tMeanGC\tGCRange\tMeanMaxHP\tMaxHP\tNumOligos\tRuntimeSec\n');

for i = 1:height(ResultTable)
    fprintf(fid, '%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%s\t%.6f\t%d\t%d\t%.4f\n', ...
        ResultTable.Method(i), ...
        ResultTable.Entropy(i), ...
        ResultTable.CorrH(i), ...
        ResultTable.CorrV(i), ...
        ResultTable.CorrD(i), ...
        ResultTable.NPCR(i), ...
        ResultTable.UACI(i), ...
        ResultTable.MeanGC(i), ...
        ResultTable.GCRange(i), ...
        ResultTable.MeanMaxHomopolymer(i), ...
        ResultTable.MaxHomopolymer(i), ...
        ResultTable.NumOligos(i), ...
        ResultTable.RuntimeSec(i));
end

fprintf(fid, '=========================================================\n');
fclose(fid);

%% =========================
%  保存 MAT
%% =========================
save(fullfile(out_dir, 'chaos_comparison_workspace.mat'), 'ResultTable');

fprintf('\n已保存：\n');
fprintf('1. CSV 结果表 : %s\n', csv_out);
fprintf('2. TXT 结果表 : %s\n', txt_out);
fprintf('3. MAT 工作区 : %s\n', fullfile(out_dir, 'chaos_comparison_workspace.mat'));

%% =========================================================
%  辅助函数：按方法名生成密钥流
%% =========================================================
function [key_stream_diffusion, key_stream_scrambling, key_stream_dna] = ...
    get_chaotic_streams(method_name, original_image)

    switch lower(method_name)
        case {'cis-t', 'cis_t', 'ten'}
            [key_stream_diffusion, key_stream_scrambling, key_stream_dna] = ...
                generate_chaotic_quence_ten(original_image);

        case {'logistic'}
            [key_stream_diffusion, key_stream_scrambling, key_stream_dna] = ...
                generate_chaotic_quence_logistic(original_image);

        case {'henon'}
            [key_stream_diffusion, key_stream_scrambling, key_stream_dna] = ...
                generate_chaotic_quence_henon(original_image);

        case {'6d'}
            [key_stream_diffusion, key_stream_scrambling, key_stream_dna] = ...
                generate_chaotic_quence_6d(original_image);

        otherwise
            error('未知方法名: %s', method_name);
    end
end

%% =========================================================
%  辅助函数：局部计算熵
%% =========================================================
function ent = calc_entropy_local(img)
    img = uint8(img);
    counts = imhist(img(:));
    p = counts / sum(counts);
    p(p==0) = [];
    ent = -sum(p .* log2(p));
end

%% =========================================================
%  辅助函数：三方向相关性
%% =========================================================
function [corrH, corrV, corrD] = calc_correlation_3dir_local(img, sample_num)

    if nargin < 2
        sample_num = 5000;
    end

    img = double(img);
    if ndims(img) == 3
        img = img(:,:,1);  % 先只取一个通道，避免太复杂
    end

    [H, W] = size(img);

    % Horizontal
    xh = zeros(sample_num,1);
    yh = zeros(sample_num,1);
    for k = 1:sample_num
        i = randi([1,H]);
        j = randi([1,W-1]);
        xh(k) = img(i,j);
        yh(k) = img(i,j+1);
    end
    corrMat = corrcoef(xh, yh);
    corrH = corrMat(1,2);

    % Vertical
    xv = zeros(sample_num,1);
    yv = zeros(sample_num,1);
    for k = 1:sample_num
        i = randi([1,H-1]);
        j = randi([1,W]);
        xv(k) = img(i,j);
        yv(k) = img(i+1,j);
    end
    corrMat = corrcoef(xv, yv);
    corrV = corrMat(1,2);

    % Diagonal
    xd = zeros(sample_num,1);
    yd = zeros(sample_num,1);
    for k = 1:sample_num
        i = randi([1,H-1]);
        j = randi([1,W-1]);
        xd(k) = img(i,j);
        yd(k) = img(i+1,j+1);
    end
    corrMat = corrcoef(xd, yd);
    corrD = corrMat(1,2);
end

%% =========================================================
%  辅助函数：NPCR/UACI
%% =========================================================
function [npcr_val, uaci_val] = calc_npcr_uaci_local(img1, img2)

    img1 = double(img1);
    img2 = double(img2);

    % 如果尺寸一致，直接算；不一致报错
    if ~isequal(size(img1), size(img2))
        error('calc_npcr_uaci_local: 两幅图尺寸不一致');
    end

    D = img1 ~= img2;
    npcr_val = sum(D(:)) / numel(img1) * 100;

    uaci_val = mean(abs(img1(:) - img2(:)) / 255) * 100;
end

%% =========================================================
%  辅助函数：row-major flatten
%% =========================================================
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

%% =========================================================
%  辅助函数：GC 计算
%% =========================================================
function gc = calc_gc_percent_local(seq)
    seq = upper(char(seq));
    gc = 100 * (sum(seq == 'G') + sum(seq == 'C')) / numel(seq);
end

%% =========================================================
%  辅助函数：最大 homopolymer 长度
%% =========================================================
function max_run = calc_max_homopolymer_local(seq)
    seq = upper(char(seq));

    if isempty(seq)
        max_run = 0;
        return;
    end

    max_run = 1;
    current_run = 1;

    for i = 2:numel(seq)
        if seq(i) == seq(i-1)
            current_run = current_run + 1;
            if current_run > max_run
                max_run = current_run;
            end
        else
            current_run = 1;
        end
    end
end