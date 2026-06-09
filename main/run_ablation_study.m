%% run_ablation_study_with_rule_dependence.m
% Ablation study for CIS-T-based secure image representation
%
% Variants:
%   1) Full model
%   2) w/o CIS-T coupling: replace CIS-T with single-branch Logistic stream
%   3) w/o DNA rule control: use a fixed DNA rule
%   4) w/o balanced padding: use repeated A padding / ECC placeholder
%   5) w/o oligo organization: output only a long DNA sequence
%
% Added module:
%   DNA_Rule_Key_Diff_percent
%   This metric measures whether the nucleotide representation remains
%   key-dependent after image-domain encryption. For the same diffused image,
%   two DNA representations are generated using two different DNA-rule streams.
%   Dynamic DNA-rule control should give a high difference rate, while fixed
%   DNA-rule control should give approximately 0%.
%
% Output:
%   ablation_results/ablation_study_results.csv

clear; clc; close all;

%% =========================
% 1. Input image and output directory
% ==========================
img_path = 'photos/house.tiff';
resize_to = [512, 512];

out_dir = 'ablation_results';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

original_image = local_read_rgb_uint8(img_path, resize_to);

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
% 4. Generate full CIS-T key streams
% ==========================
master_key = 'CIST_DNA_STORAGE_KEY_000';
salt       = 'DNA_STORAGE_SALT_000';

[key_diff_cist, key_scr_cist, key_dna_cist] = ...
    generate_chaotic_quence(original_image, master_key, salt);

% Generate a second DNA-rule stream using a one-bit-changed master key.
% It is used only to measure whether DNA representation is still key-dependent
% after image-domain encryption.
master_key_alt = 'CIST_DNA_STORAGE_KEY_001';
[~, ~, key_dna_cist_alt] = ...
    generate_chaotic_quence(original_image, master_key_alt, salt);

%% =========================
% 5. Prepare differential plaintext
% ==========================
plain2 = original_image;
plain2(1,1,:) = uint8(mod(double(plain2(1,1,:)) + 1, 256));

%% =========================
% 6. Define and run variants
% ==========================
rows = {};

% ---------- Variant 1: Full model ----------
rows(end+1,:) = local_run_variant( ...
    'Full', ...
    original_image, plain2, ...
    key_diff_cist, key_scr_cist, key_dna_cist, key_dna_cist_alt, ...
    dna_complementary_principle, oligo_params, ...
    'dynamic', 'balanced', true);

% ---------- Variant 2: w/o CIS-T coupling ----------
% Use a single-branch Logistic generator to replace the coupled CIS-T source.
[key_diff_log, key_scr_log, key_dna_log] = ...
    local_generate_logistic_streams(size(original_image), 3.99, 0.3141592653);

% A slightly changed initial value provides an alternative Logistic DNA-rule stream.
[~, ~, key_dna_log_alt] = ...
    local_generate_logistic_streams(size(original_image), 3.99, 0.3141592654);

rows(end+1,:) = local_run_variant( ...
    'w/o CIS-T coupling', ...
    original_image, plain2, ...
    key_diff_log, key_scr_log, key_dna_log, key_dna_log_alt, ...
    dna_complementary_principle, oligo_params, ...
    'dynamic', 'balanced', true);

% ---------- Variant 3: w/o DNA rule control ----------
% Fixed rule 1 is used for all encoded units. The alternative rule stream is
% also fixed, so the DNA-rule key-dependence should be approximately 0%.
key_dna_fixed = ones(size(key_dna_cist), 'like', key_dna_cist);
key_dna_fixed_alt = ones(size(key_dna_cist), 'like', key_dna_cist);

rows(end+1,:) = local_run_variant( ...
    'w/o DNA rule control', ...
    original_image, plain2, ...
    key_diff_cist, key_scr_cist, key_dna_fixed, key_dna_fixed_alt, ...
    dna_complementary_principle, oligo_params, ...
    'fixed', 'balanced', true);

% ---------- Variant 4: w/o balanced padding ----------
% Use repeated A padding and repeated A ECC placeholder.
rows(end+1,:) = local_run_variant( ...
    'w/o balanced padding', ...
    original_image, plain2, ...
    key_diff_cist, key_scr_cist, key_dna_cist, key_dna_cist_alt, ...
    dna_complementary_principle, oligo_params, ...
    'dynamic', 'plainA', true);

% ---------- Variant 5: w/o oligo organization ----------
% Only output a long DNA sequence; oligo-level evaluation is unavailable.
rows(end+1,:) = local_run_variant( ...
    'w/o oligo organization', ...
    original_image, plain2, ...
    key_diff_cist, key_scr_cist, key_dna_cist, key_dna_cist_alt, ...
    dna_complementary_principle, oligo_params, ...
    'dynamic', 'none', false);

%% =========================
% 7. Export results
% ==========================
var_names = { ...
    'Variant', ...
    'Reversible', ...
    'Entropy', ...
    'Corr_H', 'Corr_V', 'Corr_D', ...
    'NPCR_percent', 'UACI_percent', ...
    'Mean_GC_percent', ...
    'Min_GC_percent', 'Max_GC_percent', ...
    'Mean_Max_HP', 'Max_HP', ...
    'DNA_2mer_Entropy', ...
    'DNA_Rule_Key_Diff_percent', ...
    'Num_Oligos', ...
    'Runtime_sec' ...
};

ablation_table = cell2table(rows, 'VariableNames', var_names);

csv_path = fullfile(out_dir, 'ablation_study_results.csv');
writetable(ablation_table, csv_path);

fprintf('\n=============================================\n');
fprintf(' Ablation study finished.\n');
fprintf(' Results saved to:\n');
fprintf('   %s\n', csv_path);
fprintf('=============================================\n');

disp(ablation_table);


%% ============================================================
% Local functions
% ============================================================

function row = local_run_variant(variant_name, original_image, plain2, ...
    key_diff, key_scr, key_dna, key_dna_alt, dna_rules, oligo_params, ...
    rule_mode, padding_mode, use_oligo)

    fprintf('\n===== Running variant: %s =====\n', variant_name);
    t_start = tic;

    %% Main encryption
    scrambled = scramble_img(original_image, key_scr);
    diffused  = diffuse_img(scrambled, key_diff);
    cipher    = encode_img(diffused, dna_rules, key_dna);

    %% Reconstruction check
    decoded    = decode_img(cipher, dna_rules, key_dna);
    dediffused = dediffuse_img(decoded, key_diff);
    recovered  = descramble_img(dediffused, key_scr);
    reversible = isequal(original_image, recovered);

    %% Differential plaintext encryption
    scrambled2 = scramble_img(plain2, key_scr);
    diffused2  = diffuse_img(scrambled2, key_diff);
    cipher2    = encode_img(diffused2, dna_rules, key_dna);

    %% Image-domain metrics
    entropy_val = local_entropy_uint8(cipher);
    [corr_h, corr_v, corr_d] = local_corr_3dir(cipher);
    [npcr_val, uaci_val] = local_npcr_uaci(cipher, cipher2);

    %% DNA representation
    % Fast DNA-symbol mapping: one encoded element -> one DNA symbol.
    % This keeps the same length scale as the current paper table.
    dna_seq = local_fast_dna_symbol_mapping(cipher);

    dna_2mer_entropy = local_kmer_entropy(dna_seq, 2);

    %% DNA-rule key-dependence module
    % Keep the image-domain encrypted state fixed (diffused), and only change
    % the DNA-rule stream. If DNA rule control is dynamic, the resulting DNA
    % sequence should change significantly. If the rule is fixed, it should not.
    cipher_rule_alt = encode_img(diffused, dna_rules, key_dna_alt);
    dna_seq_rule_alt = local_fast_dna_symbol_mapping(cipher_rule_alt);
    dna_rule_key_diff = local_sequence_difference_percent(dna_seq, dna_seq_rule_alt);

    if use_oligo
        [oligos, oligo_meta] = local_build_oligos_ablation( ...
            dna_seq, oligo_params, padding_mode);

        [mean_gc, min_gc, max_gc, mean_max_hp, max_hp] = ...
            local_oligo_stats(oligos);

        num_oligos = oligo_meta.num_oligos;
    else
        % Without oligo organization, oligo-level metrics are not applicable.
        mean_gc = local_gc_percent(dna_seq);
        min_gc = NaN;
        max_gc = NaN;
        mean_max_hp = NaN;
        max_hp = NaN;
        num_oligos = NaN;
    end

    runtime_sec = toc(t_start);

    fprintf('Reversible     : %d\n', reversible);
    fprintf('Entropy        : %.4f\n', entropy_val);
    fprintf('Corr-H/V/D     : %.4f / %.4f / %.4f\n', corr_h, corr_v, corr_d);
    fprintf('NPCR / UACI    : %.4f / %.4f\n', npcr_val, uaci_val);
    fprintf('Mean GC        : %.2f%%\n', mean_gc);
    if use_oligo
        fprintf('GC range       : %.2f--%.2f%%\n', min_gc, max_gc);
        fprintf('Mean max HP    : %.2f\n', mean_max_hp);
        fprintf('Max HP         : %d\n', max_hp);
        fprintf('Num oligos     : %d\n', num_oligos);
    else
        fprintf('Oligo metrics  : N/A\n');
    end
    fprintf('DNA 2mer Ent.  : %.4f\n', dna_2mer_entropy);
    fprintf('DNA rule-key difference: %.4f%%\n', dna_rule_key_diff);
    fprintf('Runtime        : %.4f s\n', runtime_sec);

    row = { ...
        variant_name, ...
        reversible, ...
        entropy_val, ...
        corr_h, corr_v, corr_d, ...
        npcr_val, uaci_val, ...
        mean_gc, ...
        min_gc, max_gc, ...
        mean_max_hp, max_hp, ...
        dna_2mer_entropy, ...
        dna_rule_key_diff, ...
        num_oligos, ...
        runtime_sec ...
    };
end

function img = local_read_rgb_uint8(img_path, resize_to)
    img = imread(img_path);

    if ndims(img) == 3 && size(img,3) > 3
        img = img(:,:,1:3);
    end

    if ndims(img) == 2
        img = repmat(img, [1, 1, 3]);
    end

    if ~isa(img, 'uint8')
        img = im2uint8(img);
    end

    if ~isempty(resize_to)
        img = imresize(img, resize_to);
    end
end

function [key_diff, key_scr, key_dna] = local_generate_logistic_streams(img_size, mu, x0)
    H = img_size(1);
    W = img_size(2);
    C = img_size(3);

    total = H * W * C;
    warmup = 1000;
    need = total * 3 + warmup + 10;

    x = zeros(need,1);
    x(1) = x0;

    for i = 1:need-1
        x(i+1) = mu * x(i) * (1 - x(i));
        if x(i+1) <= 0 || x(i+1) >= 1 || isnan(x(i+1))
            x(i+1) = mod(abs(x(i)) + 0.123456789, 1);
        end
    end

    x = x(warmup+1:end);

    s1 = x(1:total);
    s2 = x(total+1:2*total);
    s3 = x(2*total+1:3*total);

    n = H * W;

    key_diff = uint8(floor(mod(s1 * 1e14, 256)));
    key_diff = reshape(key_diff, [H, W, C]);

    key_scr = floor(mod(s2 * 1e14, n)) + 1;
    key_scr = reshape(key_scr, [H, W, C]);

    key_dna = uint8(floor(mod(s3 * 1e14, 8)) + 1);
    key_dna = reshape(key_dna, [H, W, C]);
end

function dna_seq = local_fast_dna_symbol_mapping(encoded_img)
    dna_alphabet = 'ACGT';
    dna_idx = mod(double(encoded_img(:)), 4) + 1;
    dna_seq = dna_alphabet(dna_idx);
    dna_seq = upper(char(dna_seq(:)'));
end

function diff_percent = local_sequence_difference_percent(seq1, seq2)
    seq1 = char(seq1(:)');
    seq2 = char(seq2(:)');
    L = min(length(seq1), length(seq2));
    if L == 0
        diff_percent = NaN;
        return;
    end
    diff_percent = sum(seq1(1:L) ~= seq2(1:L)) / L * 100;
end

function [oligos, meta] = local_build_oligos_ablation(dna_seq, params, padding_mode)
    dna_seq = upper(char(dna_seq(:)'));
    payload_len = params.payload_len;

    total_bases = numel(dna_seq);
    num_oligos = ceil(total_bases / payload_len);

    oligos = cell(num_oligos, 1);

    for i = 1:num_oligos
        st = (i-1) * payload_len + 1;
        ed = min(i * payload_len, total_bases);
        payload = dna_seq(st:ed);

        % Payload padding
        if length(payload) < payload_len
            pad_len = payload_len - length(payload);
            payload = [payload, local_padding_pattern(pad_len, padding_mode)];
        end

        index_seq = local_index_to_dna(i-1, params.index_len);

        if strcmpi(padding_mode, 'plainA')
            ecc_seq = repmat('A', 1, params.ecc_len);
        else
            ecc_seq = local_padding_pattern(params.ecc_len, 'balanced');
        end

        oligos{i} = [params.primer_left, index_seq, payload, ecc_seq, params.primer_right];
    end

    meta = struct();
    meta.total_bases = total_bases;
    meta.num_oligos = num_oligos;
    meta.params = params;
end

function pad = local_padding_pattern(pad_len, padding_mode)
    if pad_len <= 0
        pad = '';
        return;
    end

    switch lower(padding_mode)
        case 'plaina'
            pad = repmat('A', 1, pad_len);
        otherwise
            base = 'ACGT';
            pad = repmat(base, 1, ceil(pad_len / 4));
            pad = pad(1:pad_len);
    end
end

function seq = local_index_to_dna(index_value, len)
    % Convert nonnegative integer to base-4 DNA index.
    alphabet = 'ACGT';
    seq = repmat('A', 1, len);

    value = index_value;
    for pos = len:-1:1
        digit = mod(value, 4);
        seq(pos) = alphabet(digit + 1);
        value = floor(value / 4);
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

    xH = img(:,1:end-1,:);
    yH = img(:,2:end,:);
    rH = local_corr_vec(xH(:), yH(:));

    xV = img(1:end-1,:,:);
    yV = img(2:end,:,:);
    rV = local_corr_vec(xV(:), yV(:));

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
        gc_vals(i) = local_gc_percent(s);
        hp_vals(i) = local_max_homopolymer(s);
    end

    mean_gc = mean(gc_vals);
    min_gc = min(gc_vals);
    max_gc = max(gc_vals);
    mean_max_hp = mean(hp_vals);
    max_hp = max(hp_vals);
end

function gc = local_gc_percent(s)
    s = upper(char(s(:)'));
    gc = (sum(s == 'G') + sum(s == 'C')) / length(s) * 100;
end

function hp = local_max_homopolymer(s)
    s = char(s(:)');
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
            hp = max(hp, cur);
            cur = 1;
        end
    end
    hp = max(hp, cur);
end

function H = local_kmer_entropy(seq, k)
    seq = upper(char(seq(:)'));
    if length(seq) < k
        H = NaN;
        return;
    end

    alphabet = 'ACGT';
    num_kmers = 4^k;
    counts = zeros(num_kmers, 1);

    for i = 1:(length(seq)-k+1)
        word = seq(i:i+k-1);
        idx = 0;
        valid = true;

        for j = 1:k
            d = find(alphabet == word(j), 1) - 1;
            if isempty(d)
                valid = false;
                break;
            end
            idx = idx * 4 + d;
        end

        if valid
            counts(idx+1) = counts(idx+1) + 1;
        end
    end

    p = counts / sum(counts);
    p = p(p > 0);
    H = -sum(p .* log2(p));
end
