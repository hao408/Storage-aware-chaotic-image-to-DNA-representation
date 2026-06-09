function Mchar = to_dna_char_matrix_strict(X, targetSize)
% to_dna_char_matrix_strict  严格按位将 encoded 图像解码为 DNA 字符矩阵
%
% 用法：
%   Mchar = to_dna_char_matrix_strict(encoded_image, [H W4 C]);
%
% 输入：
%   X          : 数值矩阵（一般是 encoded_image），uint8/double 都可以
%   targetSize : 目标尺寸 [H W C]，要求 prod(targetSize) == 4 * numel(X)
%
% 规则：
%   - 将 X 看成一个字节流，每个像素 = 1 个字节(8bit)
%   - 每个字节拆成 4 个 2bit：
%         (b7 b6) -> d3
%         (b5 b4) -> d2
%         (b3 b2) -> d1
%         (b1 b0) -> d0
%   - 每个 d ∈ {0,1,2,3} 再映射到 A/C/G/T：
%         0->A, 1->C, 2->G, 3->T
%   - 最后把所有碱基 reshape 成 targetSize
%
% 输出：
%   Mchar : 大小为 targetSize 的 char 矩阵，每个元素是 'A','C','G','T'

    if nargin < 2
        error('必须提供 targetSize，例如 [H, W*4, C]。');
    end
    if numel(targetSize) ~= 3
        error('targetSize 必须是长度为3的向量 [H W C]。');
    end

    bytes = uint8(X(:));        % N 个字节
    N = numel(bytes);
    total_bases = 4 * N;

    if prod(targetSize) ~= total_bases
        error('targetSize 不匹配：prod(targetSize) = %d, 但 4*numel(X) = %d', ...
              prod(targetSize), total_bases);
    end

    % ---- 拆 bit：b7b6, b5b4, b3b2, b1b0 ----
    d3 = bitshift(bitand(bytes, uint8(192)), -6); % 192 = 1100 0000
    d2 = bitshift(bitand(bytes, uint8(48)),  -4); %  48 = 0011 0000
    d1 = bitshift(bitand(bytes, uint8(12)),  -2); %  12 = 0000 1100
    d0 = bitand(bytes, uint8(3));                 %   3 = 0000 0011

    digits = zeros(4*N, 1, 'uint8');
    digits(1:4:end) = d3;
    digits(2:4:end) = d2;
    digits(3:4:end) = d1;
    digits(4:4:end) = d0;

    % ---- 0/1/2/3 -> A/C/G/T ----
    map = 'ACGT';
    chars = map(double(digits) + 1);   % 4N x 1 char

    % 按列主序 reshape 成 targetSize
    Mchar = reshape(chars, targetSize);
end
