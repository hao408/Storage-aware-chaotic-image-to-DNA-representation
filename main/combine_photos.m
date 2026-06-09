clc; clear; close all;

% =========================================================
% make_lightweight_supp_figures.m
% Purpose:
%   Combine the three split Histogram figures into one lightweight JPG.
%   Combine the three split Correlation figures into one lightweight JPG.
%
% Input folder:
%   supplementary_figures
%
% Required input images:
%   Figure_S1a_House_Histograms.png
%   Figure_S1b_Peppers_Histograms.png
%   Figure_S1c_Baboon_Histograms.png
%   Figure_S2a_House_Correlation.png
%   Figure_S2b_Peppers_Correlation.png
%   Figure_S2c_Baboon_Correlation.png
%
% Output images:
%   SuppFig/Figure_S1_Histograms_combined.jpg
%   SuppFig/Figure_S2_Correlation_combined.jpg
%
% Important:
%   This script does NOT delete or overwrite your original split PNG files.
% =========================================================

input_dir = 'supplementary_figures';
output_dir = 'SuppFig';

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

hist_files = {
    fullfile(input_dir, 'Figure_S1a_House_Histograms.png')
    fullfile(input_dir, 'Figure_S1b_Peppers_Histograms.png')
    fullfile(input_dir, 'Figure_S1c_Baboon_Histograms.png')
};

corr_files = {
    fullfile(input_dir, 'Figure_S2a_House_Correlation.png')
    fullfile(input_dir, 'Figure_S2b_Peppers_Correlation.png')
    fullfile(input_dir, 'Figure_S2c_Baboon_Correlation.png')
};

make_combined_jpg(hist_files, fullfile(output_dir, 'Figure_S1_Histograms_combined.jpg'), 1800);
make_combined_jpg(corr_files, fullfile(output_dir, 'Figure_S2_Correlation_combined.jpg'), 1800);

fprintf('\n已生成轻量版合并图：\n');
fprintf('%s\n', fullfile(output_dir, 'Figure_S1_Histograms_combined.jpg'));
fprintf('%s\n', fullfile(output_dir, 'Figure_S2_Correlation_combined.jpg'));
fprintf('\n原始 PNG 文件不会被删除。\n');

function make_combined_jpg(file_list, out_file, target_width)
    imgs = cell(numel(file_list), 1);

    for i = 1:numel(file_list)
        if exist(file_list{i}, 'file') ~= 2
            error('找不到文件：%s', file_list{i});
        end

        img = imread(file_list{i});

        % 如果是透明 PNG，转成白底
        if size(img, 3) == 4
            alpha = double(img(:,:,4)) / 255;
            rgb = double(img(:,:,1:3));
            white = 255 * ones(size(rgb));
            img = uint8(rgb .* alpha + white .* (1 - alpha));
        end

        if size(img, 3) == 1
            img = repmat(img, [1 1 3]);
        end

        % 统一宽度，防止合成图过大
        scale = target_width / size(img, 2);
        new_h = max(1, round(size(img, 1) * scale));
        img = imresize(img, [new_h, target_width]);

        imgs{i} = img;
    end

    gap = 40;
    total_h = sum(cellfun(@(x) size(x,1), imgs)) + gap * (numel(imgs) - 1);
    canvas = uint8(255 * ones(total_h, target_width, 3));

    y = 1;
    for i = 1:numel(imgs)
        h = size(imgs{i}, 1);
        canvas(y:y+h-1, :, :) = imgs{i};
        y = y + h + gap;
    end

    % JPG 比大尺寸 PNG 更适合 Overleaf 快速编译
    imwrite(canvas, out_file, 'jpg', 'Quality', 92);
end
