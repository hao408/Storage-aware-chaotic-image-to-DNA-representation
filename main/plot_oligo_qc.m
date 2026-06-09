clc; clear; close all;

% =========================
%  设置CSV文件路径
% =========================
csv_file = 'oligo_results/dna_cipher_pool.csv';   % 按你当前路径写
out_dir  = 'oligo_results';

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

if ~exist(csv_file, 'file')
    error('找不到文件: %s', csv_file);
end

% =========================
%  读取CSV
% =========================
T = readtable(csv_file, 'TextType', 'string');

% 检查列名
disp('CSV列名如下：');
disp(T.Properties.VariableNames);

% 兼容你的列名
gc_col  = "gc_percent";
hp_col  = "max_homopolymer";
seq_col = "sequence";
len_col = "length";

vars = string(T.Properties.VariableNames);

% 必须至少有 sequence 列，否则无法计算 GC 和同聚物
if ~ismember(seq_col, vars)
    error('CSV中未找到 sequence 列，无法计算 GC content 和 homopolymer。');
end

seq_data_full = upper(string(T.(seq_col)));

% =========================
%  选择统计对象
% =========================
use_payload_only = true;

primer_left_len  = 20;
index_len        = 8;
payload_len      = 120;
ecc_len          = 8;
primer_right_len = 20;

if use_payload_only
    payload_start = primer_left_len + index_len + 1;
    payload_end   = primer_left_len + index_len + payload_len;

    seq_data = extractBetween(seq_data_full, payload_start, payload_end);
    seq_data = string(seq_data);

    fprintf('[plot_oligo_qc] Using payload-only region: %d-%d\n', payload_start, payload_end);
else
    seq_data = seq_data_full;
    fprintf('[plot_oligo_qc] Using full exported oligo sequence.\n');
end

% =========================
%  计算 GC content
% =========================
seq_len = strlength(seq_data);
gc_data = 100 * (count(seq_data, "G") + count(seq_data, "C")) ./ seq_len;

% =========================
%  计算最大同聚物长度
% =========================
hp_data = zeros(height(T), 1);

for i = 1:height(T)
    s = char(seq_data(i));

    if isempty(s)
        hp_data(i) = 0;
    else
        max_run = 1;
        cur_run = 1;

        for j = 2:length(s)
            if s(j) == s(j-1)
                cur_run = cur_run + 1;
            else
                max_run = max(max_run, cur_run);
                cur_run = 1;
            end
        end

        max_run = max(max_run, cur_run);
        hp_data(i) = max_run;
    end
end

% =========================
%  统计长度
% =========================
if use_payload_only
    oligo_lengths = strlength(seq_data);     % 此时其实是 payload length，应为120
else
    if ismember(len_col, vars)
        oligo_lengths = T.(len_col);
    else
        oligo_lengths = strlength(seq_data);
    end
end

% =========================
%  统计量计算
% =========================
num_oligos = height(T);

gc_mean = mean(gc_data);
gc_std  = std(gc_data);
gc_min  = min(gc_data);
gc_max  = max(gc_data);

hp_mean = mean(hp_data);
hp_std  = std(hp_data);
hp_min  = min(hp_data);
hp_max  = max(hp_data);

if ~isempty(oligo_lengths)
    oligo_len_mean = mean(oligo_lengths);
    oligo_len_min  = min(oligo_lengths);
    oligo_len_max  = max(oligo_lengths);
else
    oligo_len_mean = NaN;
    oligo_len_min  = NaN;
    oligo_len_max  = NaN;
end

% =========================
%  打印统计结果
% =========================
fprintf('\n================ Payload Statistics ================\n');
fprintf('Number of payloads        : %d\n', num_oligos);
fprintf('Mean GC content (%%)       : %.4f\n', gc_mean);
fprintf('Std GC content (%%)        : %.4f\n', gc_std);
fprintf('GC range (%%)              : [%.2f, %.2f]\n', gc_min, gc_max);
fprintf('Mean max homopolymer      : %.4f\n', hp_mean);
fprintf('Std max homopolymer       : %.4f\n', hp_std);
fprintf('Homopolymer range         : [%d, %d]\n', hp_min, hp_max);
fprintf('Mean payload length (nt)  : %.2f\n', oligo_len_mean);
fprintf('Payload length range (nt) : [%d, %d]\n', oligo_len_min, oligo_len_max);
fprintf('====================================================\n');

% =========================
%  保存统计结果到txt
% =========================
txt_file = fullfile(out_dir, 'oligo_statistics_summary.txt');
fid = fopen(txt_file, 'w');

fprintf(fid, '================ Payload Statistics ================\n');
fprintf(fid, 'Number of payloads        : %d\n', num_oligos);
fprintf(fid, 'Mean GC content (%%)       : %.4f\n', gc_mean);
fprintf(fid, 'Std GC content (%%)        : %.4f\n', gc_std);
fprintf(fid, 'GC range (%%)              : [%.2f, %.2f]\n', gc_min, gc_max);
fprintf(fid, 'Mean max homopolymer      : %.4f\n', hp_mean);
fprintf(fid, 'Std max homopolymer       : %.4f\n', hp_std);
fprintf(fid, 'Homopolymer range         : [%d, %d]\n', hp_min, hp_max);
fprintf(fid, 'Mean payload length (nt)  : %.2f\n', oligo_len_mean);
fprintf(fid, 'Payload length range (nt) : [%d, %d]\n', oligo_len_min, oligo_len_max);
fprintf(fid, '====================================================\n');

fclose(fid);

% =========================
%  画图：GC直方图
% =========================
fig1 = figure('Position', [100, 100, 800, 600]);
histogram(gc_data, 15);
xlabel('GC content (%)');
ylabel('Number of oligos');
title('Distribution of GC Content in Exported Oligos');
grid on;

saveas(fig1, fullfile(out_dir, 'gc_content_histogram.png'));

% =========================
%  画图：最大同聚物直方图
% =========================
fig2 = figure('Position', [150, 150, 800, 600]);
histogram(hp_data, 'BinMethod', 'integers');
xlabel('Maximum homopolymer length');
ylabel('Number of oligos');
title('Distribution of Maximum Homopolymer Length');
grid on;

saveas(fig2, fullfile(out_dir, 'homopolymer_histogram.png'));

% =========================
%  画图：汇总双图
% =========================
fig3 = figure('Position', [200, 100, 1200, 500]);

subplot(1,2,1);
histogram(gc_data, 15);
xlabel('GC content (%)');
ylabel('Number of oligos');
title('GC Content Distribution');
grid on;

subplot(1,2,2);
histogram(hp_data, 'BinMethod', 'integers');
xlabel('Maximum homopolymer length');
ylabel('Number of oligos');
title('Homopolymer Distribution');
grid on;

saveas(fig3, fullfile(out_dir, 'oligo_qc_summary.png'));

% =========================
%  生成一个统计表CSV，方便论文直接引用
% =========================
summary_names = {
    'num_oligos';
    'gc_mean';
    'gc_std';
    'gc_min';
    'gc_max';
    'hp_mean';
    'hp_std';
    'hp_min';
    'hp_max';
    'oligo_len_mean';
    'oligo_len_min';
    'oligo_len_max'
    };

summary_values = [
    num_oligos;
    gc_mean;
    gc_std;
    gc_min;
    gc_max;
    hp_mean;
    hp_std;
    hp_min;
    hp_max;
    oligo_len_mean;
    oligo_len_min;
    oligo_len_max
    ];

SummaryTable = table(summary_names, summary_values, ...
    'VariableNames', {'Metric', 'Value'});

writetable(SummaryTable, fullfile(out_dir, 'oligo_statistics_summary.csv'));

fprintf('\n已完成：\n');
fprintf('1. GC直方图: %s\n', fullfile(out_dir, 'gc_content_histogram.png'));
fprintf('2. 同聚物直方图: %s\n', fullfile(out_dir, 'homopolymer_histogram.png'));
fprintf('3. 双图汇总: %s\n', fullfile(out_dir, 'oligo_qc_summary.png'));
fprintf('4. 统计文本: %s\n', txt_file);
fprintf('5. 统计CSV : %s\n', fullfile(out_dir, 'oligo_statistics_summary.csv'));