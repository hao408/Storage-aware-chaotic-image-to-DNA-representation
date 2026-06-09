clc; clear; close all;

% =========================================================
% Supplementary Figure S2a/S2b/S2c
% Adjacent-pixel correlation scatter plots (split version)
% Outputs:
%   Figure_S2a_House_Correlation.png
%   Figure_S2b_Peppers_Correlation.png
%   Figure_S2c_Baboon_Correlation.png
% =========================================================

out_dir = 'supplementary_figures';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

image_names = {'House', 'Peppers', 'Baboon'};
output_names = {'Figure_S2a_House_Correlation.png', ...
                'Figure_S2b_Peppers_Correlation.png', ...
                'Figure_S2c_Baboon_Correlation.png'};
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

direction_names = {'H', 'V', 'D'};
rng(1);

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

    orig_ch = original_image(:,:,1);
    enc_ch  = encoded_image(:,:,1);

    fig = figure('Color', 'w', 'Position', [100, 100, 1300, 700]);
    t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    for d = 1:3
        nexttile;
        [x, y, r] = get_adjacent_pairs_and_corr(orig_ch, direction_names{d});
        [xs, ys] = sample_points(x, y, 5000);
        plot(xs, ys, '.', 'MarkerSize', 2);
        xlim([0 255]); ylim([0 255]); axis square; box on;
        title(sprintf('%s - Original %s (r = %.4f)', image_names{idx}, direction_names{d}, r), ...
            'FontSize', 11, 'FontWeight', 'normal');
        xlabel('Pixel value'); ylabel('Adjacent value');
        set(gca, 'FontSize', 10);
    end

    for d = 1:3
        nexttile;
        [x, y, r] = get_adjacent_pairs_and_corr(enc_ch, direction_names{d});
        [xs, ys] = sample_points(x, y, 5000);
        plot(xs, ys, '.', 'MarkerSize', 2);
        xlim([0 255]); ylim([0 255]); axis square; box on;
        title(sprintf('%s - Encrypted %s (r = %.4f)', image_names{idx}, direction_names{d}, r), ...
            'FontSize', 11, 'FontWeight', 'normal');
        xlabel('Pixel value'); ylabel('Adjacent value');
        set(gca, 'FontSize', 10);
    end

    title(t, sprintf('Supplementary Figure %s. Adjacent-pixel correlation scatter plots for %s', ...
        char('a'+idx-1), image_names{idx}), 'FontWeight', 'bold', 'FontSize', 15);

    exportgraphics(fig, fullfile(out_dir, output_names{idx}), 'Resolution', 600);
    savefig(fig, fullfile(out_dir, strrep(output_names{idx}, '.png', '.fig')));
    fprintf('已保存: %s\n', fullfile(out_dir, output_names{idx}));
    close(fig);
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

function [x, y, r] = get_adjacent_pairs_and_corr(channel_img, direction_type)
    channel_img = double(channel_img);
    switch upper(direction_type)
        case 'H'
            x = channel_img(:, 1:end-1);
            y = channel_img(:, 2:end);
        case 'V'
            x = channel_img(1:end-1, :);
            y = channel_img(2:end, :);
        case 'D'
            x = channel_img(1:end-1, 1:end-1);
            y = channel_img(2:end, 2:end);
        otherwise
            error('direction_type 必须是 H、V 或 D。');
    end
    x = x(:); y = y(:);
    C = corrcoef(x, y);
    r = C(1, 2);
end

function [xs, ys] = sample_points(x, y, max_points)
    n = numel(x);
    if n <= max_points
        xs = x; ys = y;
    else
        idx = randperm(n, max_points);
        xs = x(idx); ys = y(idx);
    end
end
