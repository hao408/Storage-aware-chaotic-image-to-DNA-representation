function diffused_img = diffuse_img(original_img, key_stream_diffusion)
% DIFFUSE_IMG  跨通道强扩散（串行 + XOR + 加法）
%
% 思路：
%   1) 将 RGB 三通道打平成一个序列：
%        V = [R1,G1,B1, R2,G2,B2, ...]
%   2) 对应密钥流也打平成同样顺序的 K
%   3) 扩散公式：
%        C(1) = (V(1) + K(1)) mod 256
%        对 k >= 2:
%           T   = bitxor(V(k), C(k-1));
%           C(k)= (T + K(k)) mod 256
%
%   这样：
%     - 任意一个通道的一个像素改变，都会影响后面所有通道/像素；
%     - 差值不再是简单的"+1 传递"，而是经过 XOR + 加法混在一起。

    original_img         = uint8(original_img);
    key_stream_diffusion = uint8(key_stream_diffusion);

    [M, N, C] = size(original_img);
    if C ~= 3
        error('当前实现假设为 RGB 图像 (C=3)。');
    end

    % 若密钥只有1通道，复制到3通道
    if size(key_stream_diffusion,3) == 1
        key_stream_diffusion = repmat(key_stream_diffusion, [1,1,3]);
    end

    % ---- 打平成 V = [R1,G1,B1, R2,G2,B2, ...] ----
    R = original_img(:,:,1); Rv = R(:);
    G = original_img(:,:,2); Gv = G(:);
    B = original_img(:,:,3); Bv = B(:);

    KvR = key_stream_diffusion(:,:,1); KvR = KvR(:);
    KvG = key_stream_diffusion(:,:,2); KvG = KvG(:);
    KvB = key_stream_diffusion(:,:,3); KvB = KvB(:);

    numPix = numel(Rv);
    L = 3 * numPix;

    V  = zeros(L,1,'uint8');
    Kv = zeros(L,1,'uint8');

    V(1:3:end)  = Rv;
    V(2:3:end)  = Gv;
    V(3:3:end)  = Bv;

    Kv(1:3:end) = KvR;
    Kv(2:3:end) = KvG;
    Kv(3:3:end) = KvB;

    % ---- 串行扩散 ----
    Cvec = zeros(L,1,'uint8');

    % 第一个元素
    Cvec(1) = uint8( mod(double(V(1)) + double(Kv(1)), 256) );

    % 后续元素：先 XOR 上一个密文，再加密
    for k = 2:L
        T = bitxor(V(k), Cvec(k-1));                 % 非线性掺入前一密文
        Cvec(k) = uint8( mod(double(T) + double(Kv(k)), 256) );
    end

    % ---- 还原到 RGB 图像 ----
    Rv2 = Cvec(1:3:end);
    Gv2 = Cvec(2:3:end);
    Bv2 = Cvec(3:3:end);

    diffused_img = zeros(M,N,3,'uint8');
    diffused_img(:,:,1) = reshape(Rv2, [M,N]);
    diffused_img(:,:,2) = reshape(Gv2, [M,N]);
    diffused_img(:,:,3) = reshape(Bv2, [M,N]);
end
