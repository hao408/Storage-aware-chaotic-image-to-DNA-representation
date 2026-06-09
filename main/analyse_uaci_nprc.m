function results = analyse_uaci_nprc(img_a, img_b, need_display, largest_allowed_val)
% analyse_uaci_nprc  计算两张密文图像之间的 NPCR 和 UACI
%
% 用法：
%   results = analyse_uaci_nprc(img_a, img_b);
%   results = analyse_uaci_nprc(img_a, img_b, need_display, largest_allowed_val);
%
% 输入：
%   img_a, img_b        - 两张大小和类型相同的图像（通常是两张密文）
%   need_display        - 是否在命令行显示结果(true/false)，缺省为 true
%   largest_allowed_val - 像素最大值，8bit 图像一般为 255；缺省自动用图像最大值
%
% 输出：
%   results.npcr_score  - NPCR 数值（越接近 1 越好）
%   results.uaci_score  - UACI 数值（8bit 理论期望 ≈ 0.33）
%
% 说明：
%   这是"差分攻击意义下"的 NPCR/UACI：两张密文的差异，
%   而不是"原图 vs 密文"的差异。

    if nargin < 3 || isempty(need_display)
        need_display = true;
    end

    % 转 double 方便计算
    img_a = double(img_a);
    img_b = double(img_b);

    % 尺寸检查
    if any(size(img_a) ~= size(img_b))
        error('analyse_uaci_nprc: 两幅图像尺寸不一致！');
    end

    % 像素总数
    [M, N, C] = size(img_a);
    num_pix = M * N * C;

    % 若未给最大值，则自动取两图中最大像素
    if nargin < 4 || isempty(largest_allowed_val)
        largest_allowed_val = max([img_a(:); img_b(:)]);
        if largest_allowed_val == 0
            largest_allowed_val = 255;  % 兜底
        end
    end

    % ---------- NPCR ----------
    % 不相等的像素数量 / 总像素数
    diff_mask = (img_a ~= img_b);
    results.npcr_score = sum(diff_mask(:)) / num_pix;

    % ---------- UACI ----------
    % 所有像素差的绝对值和 / (像素个数 * 最大像素值)
    diff_val = abs(img_a - img_b);
    results.uaci_score = sum(diff_val(:)) / (num_pix * largest_allowed_val);

    % ---------- 可选显示 ----------
    if need_display
        fprintf('NPCR = %.6f\n', results.npcr_score);
        fprintf('UACI = %.6f\n', results.uaci_score);
    end
end
