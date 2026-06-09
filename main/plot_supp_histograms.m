clc; clear; close all;

% =========================================================
% Supplementary Figure S1a/S1b/S1c
% Histograms of original and encrypted images (fixed version)
% Outputs:
%   Figure_S1a_House_Histograms.png
%   Figure_S1b_Peppers_Histograms.png
%   Figure_S1c_Baboon_Histograms.png
% =========================================================

out_dir = 'supplementary_figures';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

image_names = {'House', 'Peppers', 'Baboon'};
output_names = {'Figure_S1a_House_Histograms.png', ...
                'Figure_S1b_Peppers_Histograms.png', ...
                'Figure_S1c_Baboon_Histograms.png'};
image_candidates = {
    {'photos\house.tiff', 'photos\house.tif', 'house.tiff', 'house.tif', 'House.tiff', 'House.tif'}, ...
    {'photos\peppers.png', 'peppers.png', 'Peppers.png'}, ...
    {'photos\baboon.jpg', 'photos\baboon.tiff', 'photos\baboon.png', 'baboon.jpg', 'baboon.tiff', 'baboon.tif', 'Baboon.jpg', 'Baboon.tiff', 'Baboon.tif'}
};

dna_complementary_principle = [
    84, 71, 67, 65;
    84, 67, 71, 65;
    71, 84, 65, 67;
    67, 84, 65, 71;
    84, 65, 71, 67;
    71, 65, 84, 67;
    67, 65, 84, 71;
    65, 84, 67, 71
];
dna_complementary_principle = uint8(dna_complementary_principle);

channel_names = {'R', 'G', 'B'};

for idx = 1:3
    img_path = find_existing_file(image_candidates{idx});
    if isempty(img_path)
        error('找不到图像 %s。请检查文件是否放在当前目录或 photos 文件夹下。', image_names{idx});
    end

    original_image = imread(img_path);
    if size(original_image, 3) == 1
        original_image = repmat(original_image, [1, 1, 3]);
    end
    if ~isa(original_image, 'uint8')
        original_image = im2uint8(original_image);
    end

    [key_stream_diffusion, key_stream_scrambling, key_stream_dna] = generate_chaotic_quence_ten(original_image);
    scrambled_image = scramble_img(original_image, key_stream_scrambling);
    diffused_image  = diffuse_img(scrambled_image, key_stream_diffusion);
    encoded_image   = encode_img(diffused_image, dna_complementary_principle, key_stream_dna);

    fig = figure('Color', 'w', 'Position', [100, 100, 1500, 850]);
    t = tiledlayout(2, 3, 'TileSpacing', 'loose', 'Padding', 'compact');

    % 原图 R/G/B 直方图
    for ch = 1:3
        nexttile;
        plot_hist_bar(original_image(:,:,ch));
        title(sprintf('%s - Original %s', image_names{idx}, channel_names{ch}), ...
            'FontSize', 12, 'FontWeight', 'normal');
        xlabel('Pixel value', 'FontSize', 11);
        ylabel('Frequency', 'FontSize', 11);
        set(gca, 'FontSize', 10);
        box on;
    end

    % 密文 R/G/B 直方图
    for ch = 1:3
        nexttile;
        plot_hist_bar(encoded_image(:,:,ch));
        title(sprintf('%s - Encrypted %s', image_names{idx}, channel_names{ch}), ...
            'FontSize', 12, 'FontWeight', 'normal');
        xlabel('Pixel value', 'FontSize', 11);
        ylabel('Frequency', 'FontSize', 11);
        set(gca, 'FontSize', 10);
        box on;
    end

    title(t, sprintf('Supplementary Figure S1%s. Histograms of original and encrypted %s image', ...
        char('a'+idx-1), image_names{idx}), 'FontWeight', 'bold', 'FontSize', 16);

    exportgraphics(fig, fullfile(out_dir, output_names{idx}), 'Resolution', 600);
    savefig(fig, fullfile(out_dir, strrep(output_names{idx}, '.png', '.fig')));
    fprintf('已保存: %s\n', fullfile(out_dir, output_names{idx}));
    close(fig);
end

function plot_hist_bar(channel_img)
    % 不用 imhist 直接画图，避免 tiledlayout 下 Position 警告和显示异常
    counts = imhist(channel_img, 256);
    x = 0:255;
    bar(x, counts, 1.0, 'EdgeColor', 'none');
    xlim([0 255]);
end

function filepath = find_existing_file(candidate_list)
    filepath = '';
    for i = 1:numel(candidate_list)
        if exist(candidate_list{i}, 'file') == 2
            filepath = candidate_list{i};
            return;
        end
    end
    for i = 1:numel(candidate_list)
        w = which(candidate_list{i});
        if ~isempty(w)
            filepath = w;
            return;
        end
    end
end

