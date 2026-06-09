function stats = analyse_dna_storage_density(original_img, encoded_img)
% analyse_dna_storage_density
%  计算 DNA 存储密度和相对开销（针对当前图像 + 编码方式）
%
% 用法：
%   stats = analyse_dna_storage_density(original_image, encoded_image);
%
% 输入：
%   original_img : 原始图像矩阵（uint8 或 double），大小 H x W x C
%   encoded_img  : 加密后/编码后的图像矩阵（即 encoded_image）
%
% 假设：
%   - 每个像素（一个字节）被 encode_img 打包成 4 个 2bit 的碱基
%     → 即每个像素对应 4 个碱基。
%
% 输出 stats 字段：
%   stats.orig_bytes                原始字节数
%   stats.orig_bits                 原始比特数 = orig_bytes * 8
%   stats.dna_bases                 实际使用的 DNA 碱基数
%   stats.bases_per_byte            每个原始字节需要的碱基数
%   stats.bits_per_base_effective   实际有效比特/碱基
%   stats.theoretical_min_bases     理论最少碱基数（2 bits/碱基）
%   stats.overhead_factor_vs_theory 相对理论的开销倍数
%   stats.efficiency_vs_theory      相对理论极限的效率（百分比）

    % ------- 1. 原始图像的大小 -------
    orig_bytes = numel(original_img);     % 一个像素一个字节（假定 uint8）
    orig_bits  = orig_bytes * 8;

    % ------- 2. DNA 碱基数 -------
    % 严格解码下：每个像素 = 4 个碱基
    % encoded_img 的元素个数 = 原图像像素个数（通常大小一致）
    dna_bases = 4 * numel(encoded_img);

    % ------- 3. 核心指标 -------
    bases_per_byte = dna_bases / orig_bytes;          % 每个原始字节对应多少碱基
    bits_per_base_effective = orig_bits / dna_bases;  % 有效比特/碱基

    % 理论最少碱基数（每个碱基 2 bit）
    theoretical_min_bases = orig_bits / 2;

    % 相对理论开销倍数：实际/理论 ≥ 1
    overhead_factor = dna_bases / theoretical_min_bases;

    % 相对理论效率：实际 bits/base / 2 bits/base
    efficiency_ratio = bits_per_base_effective / 2 * 100;  % 百分比

    % ------- 4. 打印结果 -------
    fprintf('\n[analyse_dna_storage_density]\n');
    fprintf('  Orig bytes                 : %d\n', orig_bytes);
    fprintf('  Orig bits                  : %d\n', orig_bits);
    fprintf('  DNA bases (actual)         : %d\n', dna_bases);
    fprintf('  Bases per byte             : %.3f\n', bases_per_byte);
    fprintf('  Effective bits/base        : %.3f (理想为 2 bits/base)\n', ...
        bits_per_base_effective);
    fprintf('  Theoretical min bases      : %.0f (2 bits/base)\n', ...
        theoretical_min_bases);
    fprintf('  Overhead factor vs theory  : %.3f x\n', overhead_factor);
    fprintf('  Efficiency vs theory       : %.2f %%\n', efficiency_ratio);

    % ------- 5. 输出结构体 -------
    stats = struct();
    stats.orig_bytes                = orig_bytes;
    stats.orig_bits                 = orig_bits;
    stats.dna_bases                 = dna_bases;
    stats.bases_per_byte            = bases_per_byte;
    stats.bits_per_base_effective   = bits_per_base_effective;
    stats.theoretical_min_bases     = theoretical_min_bases;
    stats.overhead_factor_vs_theory = overhead_factor;
    stats.efficiency_vs_theory      = efficiency_ratio;
end
