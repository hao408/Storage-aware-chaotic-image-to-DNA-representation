function stats = analyse_npcr_uaci(img1, img2)
% analyse_npcr_uaci  计算两幅同尺寸图像的 NPCR 和 UACI
%
% 用法：
%   stats = analyse_npcr_uaci(cipher1, cipher2);
%   fprintf('NPCR = %.6f, UACI = %.6f\n', stats.npcr, stats.uaci);
%
% 输入：
%   img1, img2 : uint8 图像，大小相同，可以是灰度或彩色
%
% 输出：
%   stats.npcr : NPCR 数值（0~1）
%   stats.uaci : UACI 数值（0~1）

    if ~isa(img1, 'uint8') || ~isa(img2, 'uint8')
        error('img1, img2 必须是 uint8 图像');
    end
    if any(size(img1) ~= size(img2))
        error('img1, img2 尺寸必须相同');
    end

    img1 = double(img1);
    img2 = double(img2);

    MAX = 255;   % 8bit 图像
    N = numel(img1);

    % NPCR: 不同像素个数 / 总像素
    D = img1 ~= img2;
    npcr = sum(D(:)) / N;

    % UACI: |C1 - C2| / MAX 的平均
    uaci = sum(abs(img1(:) - img2(:))) / (N * MAX);

    stats = struct('npcr', npcr, 'uaci', uaci);

    fprintf('analyse_npcr_uaci: NPCR = %.6f, UACI = %.6f\n', npcr, uaci);
end
