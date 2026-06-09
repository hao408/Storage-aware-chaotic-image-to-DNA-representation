clc; clear; close all;

%% =========================
%  路径设置
%% =========================
out_dir   = 'oligo_results';
mat_file  = fullfile(out_dir, 'oligo_workspace.mat');
csv_file  = fullfile(out_dir, 'dna_cipher_pool.csv');

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

if ~exist(mat_file, 'file')
    error('找不到文件: %s', mat_file);
end

if ~exist(csv_file, 'file')
    error('找不到文件: %s', csv_file);
end

%% =========================
%  读取约束后的 DNA 数据
%% =========================
S = load(mat_file);

if ~isfield(S, 'dna_constrained')
    error('mat文件中未找到变量 dna_constrained');
end

dna_constrained = upper(char(S.dna_constrained));

if any(~ismember(dna_constrained(:), 'ACGT'))
    error('dna_constrained 中存在非 A/C/G/T 字符');
end

%% =========================
%  统计“约束前/导出前”的 DNA 表示
%  这里把 dna_constrained 本身看作导出前的 DNA 表示
%% =========================
before_seq = flatten_dna_row_major(dna_constrained);

before_gc = calc_gc_percent(before_seq);
before_hp = calc_max_homopolymer(before_seq);

% 如果你还想统计 forbidden motifs，可在这里追加
before_len = numel(before_seq);

%% =========================
%  读取“导出后”的 oligo CSV
%% =========================
T = readtable(csv_file, 'TextType', 'string');

if ~ismember("gc_percent", string(T.Properties.VariableNames))
    error('CSV中未找到 gc_percent 列');
end
if ~ismember("max_homopolymer", string(T.Properties.VariableNames))
    error('CSV中未找到 max_homopolymer 列');
end
if ~ismember("sequence", string(T.Properties.VariableNames))
    error('CSV中未找到 sequence 列');
end

after_gc_mean = mean(T.gc_percent);
after_gc_std  = std(T.gc_percent);
after_gc_min  = min(T.gc_percent);
after_gc_max  = max(T.gc_percent);

after_hp_mean = mean(T.max_homopolymer);
after_hp_std  = std(T.max_homopolymer);
after_hp_min  = min(T.max_homopolymer);
after_hp_max  = max(T.max_homopolymer);

after_num_oligos = height(T);
after_len_mean = mean(strlength(T.sequence));
after_len_min  = min(strlength(T.sequence));
after_len_max  = max(strlength(T.sequence));

%% =========================
%  打印结果
%% =========================
fprintf('\n================ Before vs After Constraint/Export ================\n');

fprintf('\n[Before export / sequence-level organization]\n');
fprintf('Total DNA bases            : %d\n', before_len);
fprintf('GC content (%%)             : %.4f\n', before_gc);
fprintf('Maximum homopolymer        : %d\n', before_hp);

fprintf('\n[After export / oligo-level statistics]\n');
fprintf('Number of oligos           : %d\n', after_num_oligos);
fprintf('Mean GC content (%%)        : %.4f\n', after_gc_mean);
fprintf('Std GC content (%%)         : %.4f\n', after_gc_std);
fprintf('GC range (%%)               : [%.2f, %.2f]\n', after_gc_min, after_gc_max);
fprintf('Mean max homopolymer       : %.4f\n', after_hp_mean);
fprintf('Std max homopolymer        : %.4f\n', after_hp_std);
fprintf('Homopolymer range          : [%d, %d]\n', after_hp_min, after_hp_max);
fprintf('Mean oligo length (nt)     : %.2f\n', after_len_mean);
fprintf('Oligo length range (nt)    : [%d, %d]\n', after_len_min, after_len_max);

fprintf('===================================================================\n');

%% =========================
%  生成对比表
%% =========================
Setting = [
    "Before export (DNA representation)";
    "After export (Oligo pool mean)"
    ];

GC_mean = [
    before_gc;
    after_gc_mean
    ];

GC_range = [
    string(sprintf('%.2f', before_gc));
    string(sprintf('%.2f-%.2f', after_gc_min, after_gc_max))
    ];

HP_mean = [
    before_hp;
    after_hp_mean
    ];

HP_range = [
    string(sprintf('%d', before_hp));
    string(sprintf('%d-%d', after_hp_min, after_hp_max))
    ];

Length_info = [
    string(sprintf('%d bases', before_len));
    string(sprintf('%.2f nt (range %d-%d)', after_len_mean, after_len_min, after_len_max))
    ];

CompareTable = table(Setting, GC_mean, GC_range, HP_mean, HP_range, Length_info);

disp(' ');
disp('================ Table for Paper ================');
disp(CompareTable);
disp('=================================================');

%% =========================
%  保存为CSV
%% =========================
csv_out = fullfile(out_dir, 'before_after_constraint_summary.csv');
writetable(CompareTable, csv_out);

%% =========================
%  保存为TXT
%% =========================
txt_out = fullfile(out_dir, 'before_after_constraint_summary.txt');
fid = fopen(txt_out, 'w');

fprintf(fid, '================ Before vs After Constraint/Export ================\n');

fprintf(fid, '\n[Before export / sequence-level organization]\n');
fprintf(fid, 'Total DNA bases            : %d\n', before_len);
fprintf(fid, 'GC content (%%)             : %.4f\n', before_gc);
fprintf(fid, 'Maximum homopolymer        : %d\n', before_hp);

fprintf(fid, '\n[After export / oligo-level statistics]\n');
fprintf(fid, 'Number of oligos           : %d\n', after_num_oligos);
fprintf(fid, 'Mean GC content (%%)        : %.4f\n', after_gc_mean);
fprintf(fid, 'Std GC content (%%)         : %.4f\n', after_gc_std);
fprintf(fid, 'GC range (%%)               : [%.2f, %.2f]\n', after_gc_min, after_gc_max);
fprintf(fid, 'Mean max homopolymer       : %.4f\n', after_hp_mean);
fprintf(fid, 'Std max homopolymer        : %.4f\n', after_hp_std);
fprintf(fid, 'Homopolymer range          : [%d, %d]\n', after_hp_min, after_hp_max);
fprintf(fid, 'Mean oligo length (nt)     : %.2f\n', after_len_mean);
fprintf(fid, 'Oligo length range (nt)    : [%d, %d]\n', after_len_min, after_len_max);

fprintf(fid, '\n================ Table for Paper ================\n');
fprintf(fid, 'Setting\tGC_mean\tGC_range\tHP_mean\tHP_range\tLength_info\n');
for i = 1:height(CompareTable)
    fprintf(fid, '%s\t%.4f\t%s\t%.4f\t%s\t%s\n', ...
        CompareTable.Setting(i), ...
        CompareTable.GC_mean(i), ...
        CompareTable.GC_range(i), ...
        CompareTable.HP_mean(i), ...
        CompareTable.HP_range(i), ...
        CompareTable.Length_info(i));
end
fprintf(fid, '=================================================\n');

fclose(fid);

fprintf('\n已保存：\n');
fprintf('1. 对比表 CSV : %s\n', csv_out);
fprintf('2. 对比表 TXT : %s\n', txt_out);

%% =========================
%  辅助函数
%% =========================
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

function gc = calc_gc_percent(seq)
    seq = upper(char(seq));
    gc = 100 * (sum(seq == 'G') + sum(seq == 'C')) / numel(seq);
end

function maxhp = calc_max_homopolymer(seq)
    seq = upper(char(seq));
    if isempty(seq)
        maxhp = 0;
        return;
    end

    maxhp = 1;
    runlen = 1;
    for i = 2:numel(seq)
        if seq(i) == seq(i-1)
            runlen = runlen + 1;
            if runlen > maxhp
                maxhp = runlen;
            end
        else
            runlen = 1;
        end
    end
end