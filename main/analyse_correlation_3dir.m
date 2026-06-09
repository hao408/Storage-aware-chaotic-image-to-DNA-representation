function [rH, rV, rD] = analyse_correlation_3dir(img)
% analyse_correlation_3dir  计算并可视化密文的水平/垂直/对角相关系数
%
% 用法：
%   figure;
%   [rH, rV, rD] = analyse_correlation_3dir(encoded_image);
%
% 输入：
%   img : 加密后的图像，uint8，HxW 或 HxWx3
%
% 输出：
%   rH, rV, rD : 水平、垂直、对角三个方向相邻像素相关系数
%
% 说明：
%   - 若为彩色图，默认只用第 1 通道（你也可以自己改为 3 通道平均）
%   - 本函数会在当前 figure 里画 3 个子图散点图，并在标题中标出相关系数

    if ndims(img) == 3
        ch = img(:,:,1);   % 只用第1通道
    else
        ch = img;
    end
    ch = double(ch);

    % 水平方向相邻像素
    Xh = ch(:, 1:end-1);
    Yh = ch(:, 2:end);
    Xh = Xh(:);
    Yh = Yh(:);
    C = corrcoef(Xh, Yh);
    rH = C(1,2);

    % 垂直方向相邻像素
    Xv = ch(1:end-1, :);
    Yv = ch(2:end,   :);
    Xv = Xv(:);
    Yv = Yv(:);
    C = corrcoef(Xv, Yv);
    rV = C(1,2);

    % 对角方向相邻像素
    Xd = ch(1:end-1, 1:end-1);
    Yd = ch(2:end,   2:end);
    Xd = Xd(:);
    Yd = Yd(:);
    C = corrcoef(Xd, Yd);
    rD = C(1,2);

    % -------- 画散点图 --------
    % 你可以在主脚本里先 figure; 再调用本函数
    clf;
    subplot(1,3,1);
    scatter(Xh, Yh, 1, '.');
    xlabel('x_i'); ylabel('x_{i+1}');
    title(sprintf('Horizon: r = %.4f', rH));

    subplot(1,3,2);
    scatter(Xv, Yv, 1, '.');
    xlabel('x_i'); ylabel('x_{i+1}');
    title(sprintf('Vertical: r = %.4f', rV));

    subplot(1,3,3);
    scatter(Xd, Yd, 1, '.');
    xlabel('x_i'); ylabel('x_{i+1}');
    title(sprintf('Diagonal: r = %.4f', rD));

end
