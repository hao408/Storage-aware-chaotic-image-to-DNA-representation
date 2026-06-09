%% run_multi_image_validation.m
% Multi-image validation for CIS-T-based secure image representation
% This script evaluates image-domain security and DNA/oligo sequence quality
% on multiple standard images under the same pipeline.

clear; clc; close all;

%% =========================
% 1. Image list
% ==========================
image_files = {
    'photos/house.tiff'
    'photos/Sailboat on lake.tiff'
    'photos/peppers.tiff'
    'photos/baboon.png'
    'photos/airplane.tiff'
};

resize_to = [512, 512];   % keep comparable with the House setting
out_dir = 'multi_image_results';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% =========================
% 2. DNA encoding rules
% ==========================
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
% 3. Oligo parameters
% ==========================
oligo_params = struct();
oligo_params.payload_len  = 120;
oligo_params.index_len    = 8;
oligo_params.ecc_len      = 8;
oligo_params.primer_left  = 'ACGTACGTACGTACGTACGT';
oligo_params.primer_right = 'TGCATGCATGCATGCATGCA';

%% =========================
% 4. Result containers
% ==========================
result_rows = {};

fprintf('\n=============================================\n');
fprintf(' Multi-image validation starts\n');
fprintf('=============================================\n');

for idx = 1:numel(image_files)

    img_path = image_files{idx};

    if ~exist(img_path, 'file')
        fprintf('[Skip] File not found: %s\n', img_path);
        continue;
    end

    [~, img_name, ext] = fileparts(img_path);
    fprintf('\n\n===== Processing %s%s =====\n', img_name, ext);

    t_start = tic;

    %% Read and preprocess image
    original_image = local_read_rgb_uint8(img_path, resize_to);

    %% Generate CIS-T key streams
    [key_stream_diffusion, key_stream_scrambling, key_stream_dna] = ...
        generate_chaotic_quence_ten(original_image);

    %% Encryption pipeline
    scrambled_image = scramble_img(original_image, key_stream_scrambling);
    diffused_image  = diffuse_img(scrambled_image, key_stream_diffusion);
    encoded_image   = encode_img(diffused_image, dna_complementary_principle, key_stream_dna);

    %% Decryption check
    decoded_image     = decode_img(encoded_image, dna_complementary_principle, key_stream_dna);
    dediffused_image  = dediffuse_img(decoded_image, key_stream_diffusion);
    recovered_image   = descramble_img(dediffused_image, key_stream_scrambling);

    is_reversible = isequal(original_image, recovered_image);
    if is_reversible
        fprintf('[OK] Noise-free reconstruction successful.\n');
    else
        warning('[Warning] Noise-free reconstruction failed for %s.', img_name);
    end

    %% Differential attack test: change one RGB pixel as in the current main program
    plain2 = original_image;
    plain2(1,1,:) = uint8(mod(double(plain2(1,1,:)) + 1, 256));

    scrambled2 = scramble_img(plain2, key_stream_scrambling);
    diffused2  = diffuse_img(scrambled2, key_stream_diffusion);
    encoded2   = encode_img(diffused2, dna_complementary_principle, key_stream_dna);

    %% Image-domain security metrics
    entropy_val = local_entropy_uint8(encoded_image);
    [corr_h, corr_v, corr_d] = local_corr_3dir(encoded_image);
    [npcr_val, uaci_val] = local_npcr_uaci(encoded_image, encoded2);

%% Fast DNA-symbol mapping for multi-image validation
% To keep the same length scale as the current paper results, each element
% of encoded_image is mapped to one DNA symbol. This avoids the expensive
% iterative constraint-repair step and keeps the oligo count comparable
% with the House experiment.

dna_alphabet = 'ACGT';
dna_idx = mod(double(encoded_image), 4) + 1;
dna_constrained = dna_alphabet(dna_idx);
dna_constrained = upper(char(dna_constrained));

% Basic DNA statistics in fast mode
dna_stats = struct();
dna_stats.total_bases = numel(dna_constrained);
dna_stats.gc_percent = (sum(dna_constrained(:) == 'G') + sum(dna_constrained(:) == 'C')) ...
                       / numel(dna_constrained) * 100;
dna_stats.max_homopolymer = local_max_homopolymer(dna_constrained(:)');

    [oligos, oligo_meta] = build_oligos(dna_constrained, oligo_params);

    [mean_gc, min_gc, max_gc, mean_max_hp, max_hp] = local_oligo_stats(oligos);

    runtime_sec = toc(t_start);

    fprintf('Entropy      : %.4f\n', entropy_val);
    fprintf('Corr-H/V/D   : %.4f / %.4f / %.4f\n', corr_h, corr_v, corr_d);
    fprintf('NPCR / UACI  : %.4f / %.4f\n', npcr_val, uaci_val);
    fprintf('Mean GC      : %.2f%%\n', mean_gc);
    fprintf('GC range     : %.2f--%.2f%%\n', min_gc, max_gc);
    fprintf('Mean max HP  : %.2f\n', mean_max_hp);
    fprintf('Max HP       : %d\n', max_hp);
    fprintf('Num oligos   : %d\n', oligo_meta.num_oligos);
    fprintf('Runtime      : %.4f s\n', runtime_sec);

    %% Save per-image workspace
    save(fullfile(out_dir, [img_name '_workspace.mat']), ...
        'original_image', 'encoded_image', 'encoded2', ...
        'dna_constrained', 'dna_stats', 'oligos', 'oligo_meta', ...
        'entropy_val', 'corr_h', 'corr_v', 'corr_d', ...
        'npcr_val', 'uaci_val', ...
        'mean_gc', 'min_gc', 'max_gc', 'mean_max_hp', 'max_hp', ...
        'runtime_sec', 'is_reversible');

    %% Collect row
    result_rows(end+1,:) = {
        img_name, ...
        size(original_image,1), size(original_image,2), size(original_image,3), ...
        is_reversible, ...
        entropy_val, corr_h, corr_v, corr_d, ...
        npcr_val, uaci_val, ...
        mean_gc, min_gc, max_gc, mean_max_hp, max_hp, ...
        oligo_meta.num_oligos, runtime_sec
    };

end

%% =========================
% 5. Export results
% ==========================
var_names = {
    'Image', ...
    'Height', 'Width', 'Channels', ...
    'Reversible', ...
    'Entropy', 'Corr_H', 'Corr_V', 'Corr_D', ...
    'NPCR_percent', 'UACI_percent', ...
    'Mean_GC_percent', 'Min_GC_percent', 'Max_GC_percent', ...
    'Mean_Max_HP', 'Max_HP', ...
    'Num_Oligos', 'Runtime_sec'
};

results_table = cell2table(result_rows, 'VariableNames', var_names);

csv_all = fullfile(out_dir, 'multi_image_all_results.csv');
writetable(results_table, csv_all);

security_table = results_table(:, {
    'Image', 'Entropy', 'Corr_H', 'Corr_V', 'Corr_D', 'NPCR_percent', 'UACI_percent'
});
csv_security = fullfile(out_dir, 'multi_image_security_results.csv');
writetable(security_table, csv_security);

sequence_table = results_table(:, {
    'Image', 'Mean_GC_percent', 'Min_GC_percent', 'Max_GC_percent', ...
    'Mean_Max_HP', 'Max_HP', 'Num_Oligos', 'Runtime_sec'
});
csv_sequence = fullfile(out_dir, 'multi_image_sequence_results.csv');
writetable(sequence_table, csv_sequence);

fprintf('\n=============================================\n');
fprintf(' Multi-image validation finished.\n');
fprintf(' Results saved to:\n');
fprintf('   %s\n', csv_all);
fprintf('   %s\n', csv_security);
fprintf('   %s\n', csv_sequence);
fprintf('=============================================\n');

disp(results_table);


%% ============================================================
% Local helper functions
% ============================================================

function img = local_read_rgb_uint8(img_path, resize_to)
    img = imread(img_path);

    % Remove alpha channel if present
    if ndims(img) == 3 && size(img,3) > 3
        img = img(:,:,1:3);
    end

    % Convert grayscale to RGB
    if ndims(img) == 2
        img = repmat(img, [1, 1, 3]);
    end

    % Convert to uint8
    if ~isa(img, 'uint8')
        img = im2uint8(img);
    end

    % Resize for fair comparison
    if ~isempty(resize_to)
        img = imresize(img, resize_to);
    end
end

function H = local_entropy_uint8(img)
    img = uint8(img);
    vals = img(:);
    counts = accumarray(double(vals)+1, 1, [256, 1]);
    p = counts / sum(counts);
    p = p(p > 0);
    H = -sum(p .* log2(p));
end

function [rH, rV, rD] = local_corr_3dir(img)
    img = double(img);

    % Horizontal adjacent pairs
    xH = img(:,1:end-1,:);
    yH = img(:,2:end,:);
    rH = local_corr_vec(xH(:), yH(:));

    % Vertical adjacent pairs
    xV = img(1:end-1,:,:);
    yV = img(2:end,:,:);
    rV = local_corr_vec(xV(:), yV(:));

    % Diagonal adjacent pairs
    xD = img(1:end-1,1:end-1,:);
    yD = img(2:end,2:end,:);
    rD = local_corr_vec(xD(:), yD(:));
end

function r = local_corr_vec(x, y)
    x = double(x(:));
    y = double(y(:));

    if numel(x) > 50000
        rng(1);
        idx = randperm(numel(x), 50000);
        x = x(idx);
        y = y(idx);
    end

    x = x - mean(x);
    y = y - mean(y);

    denom = sqrt(sum(x.^2) * sum(y.^2));
    if denom == 0
        r = 0;
    else
        r = sum(x .* y) / denom;
    end
end

function [npcr, uaci] = local_npcr_uaci(C1, C2)
    C1 = uint8(C1);
    C2 = uint8(C2);

    diff_pixels = C1 ~= C2;
    npcr = sum(diff_pixels(:)) / numel(C1) * 100;

    uaci = mean(abs(double(C1(:)) - double(C2(:))) / 255) * 100;
end

function [mean_gc, min_gc, max_gc, mean_max_hp, max_hp] = local_oligo_stats(oligos)
    n = numel(oligos);
    gc_vals = zeros(n,1);
    hp_vals = zeros(n,1);

    for i = 1:n
        s = upper(char(oligos{i}));
        gc_vals(i) = (sum(s == 'G') + sum(s == 'C')) / length(s) * 100;
        hp_vals(i) = local_max_homopolymer(s);
    end

    mean_gc = mean(gc_vals);
    min_gc = min(gc_vals);
    max_gc = max(gc_vals);
    mean_max_hp = mean(hp_vals);
    max_hp = max(hp_vals);
end

function hp = local_max_homopolymer(s)
    if isempty(s)
        hp = 0;
        return;
    end

    hp = 1;
    cur = 1;
    for i = 2:length(s)
        if s(i) == s(i-1)
            cur = cur + 1;
        else
            if cur > hp
                hp = cur;
            end
            cur = 1;
        end
    end

    if cur > hp
        hp = cur;
    end
end