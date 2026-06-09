function dediffused_img = dediffuse_img(diffused_img, key_stream_diffusion)
% DEDIFFUSE_IMG  跨通道强扩散的逆过程
%
% 对应加密：
%   C(1) = (V(1) + K(1)) mod 256
%   T(k) = bitxor(V(k), C(k-1))
%   C(k) = (T(k) + K(k)) mod 256
%
% 解密：
%   V(1) = (C(1) - K(1)) mod 256
%   对 k >= 2:
%       T(k) = (C(k) - K(k)) mod 256
%       V(k) = bitxor(T(k), C(k-1))

    diffused_img         = uint8(diffused_img);
    key_stream_diffusion = uint8(key_stream_diffusion);

    [M, N, C] = size(diffused_img);
    if C ~= 3
        error('当前实现假设为 RGB 图像 (C=3)。');
    end

    if size(key_stream_diffusion,3) == 1
        key_stream_diffusion = repmat(key_stream_diffusion, [1,1,3]);
    end

    % ---- 打平成 Cvec、Kv ----
    R = diffused_img(:,:,1); Rv = R(:);
    G = diffused_img(:,:,2); Gv = G(:);
    B = diffused_img(:,:,3); Bv = B(:);

    KvR = key_stream_diffusion(:,:,1); KvR = KvR(:);
    KvG = key_stream_diffusion(:,:,2); KvG = KvG(:);
    KvB = key_stream_diffusion(:,:,3); KvB = KvB(:);

    numPix = numel(Rv);
    L = 3 * numPix;

    Cvec = zeros(L,1,'uint8');
    Kv   = zeros(L,1,'uint8');

    Cvec(1:3:end) = Rv;
    Cvec(2:3:end) = Gv;
    Cvec(3:3:end) = Bv;

    Kv(1:3:end) = KvR;
    Kv(2:3:end) = KvG;
    Kv(3:3:end) = KvB;

    % ---- 逆扩散 ----
    V = zeros(L,1,'uint8');

    % k = 1
    V(1) = uint8( mod(double(Cvec(1)) - double(Kv(1)), 256) );

    % k >= 2
    for k = 2:L
        Tk = uint8( mod(double(Cvec(k)) - double(Kv(k)), 256) );
        V(k) = bitxor(Tk, Cvec(k-1));
    end

    % ---- 还原成 RGB ----
    Rv2 = V(1:3:end);
    Gv2 = V(2:3:end);
    Bv2 = V(3:3:end);

    dediffused_img = zeros(M,N,3,'uint8');
    dediffused_img(:,:,1) = reshape(Rv2, [M,N]);
    dediffused_img(:,:,2) = reshape(Gv2, [M,N]);
    dediffused_img(:,:,3) = reshape(Bv2, [M,N]);
end
