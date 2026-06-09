function [key_stream_diffusion, key_stream_scrambling, key_stream_dna] = generate_chaotic_quence_logistic(img, master_key, salt)
% LOGISTIC-3TRACK chaotic key streams: diffusion, scrambling, DNA rule
% 与你现有的 generate_chaotic_quence_ten 用法一致，只是内部换成 Logistic 映射。
%
% Usage:
%   [Kdif, Kscr, Kdna] = generate_chaotic_quence_logistic(img)
%   [Kdif, Kscr, Kdna] = generate_chaotic_quence_logistic(img, master_key, salt)
%
% 输出：
%   key_stream_diffusion  [H,W,C] uint8 0..255
%   key_stream_scrambling [H,W,C] uint8 0..255
%   key_stream_dna        [H,W,C] uint8 0..7

    % ---------- 输入检查 ----------
    if ~isa(img,'uint8') || ndims(img)~=3
        error('img must be HxWxC uint8.');
    end
    [H,W,C] = size(img);
    N = H*W*C;

    % ---------- 默认 master_key / salt ----------
    if nargin < 2 || isempty(master_key)
        master_key = 'LOGISTIC-DEFAULT-KEY';
    end
    img_hash = sha256_bytes(uint8(img(:)));
    if nargin < 3 || isempty(salt)
        salt = img_hash;
    end

    % ---------- HKDF-like 派生 6 个 (0,1) 上的 u_i ----------
    seed_in   = [rowu8(to_bytes(master_key)), rowu8(img_hash), rowu8(to_bytes(salt))];
    needBytes = 6*16;           % 6 个参数，每个 16 字节
    kbytes    = kdf_expand(seed_in, 'LOGISTIC-3TRACK', needBytes);

    u  = zeros(1,6); off = 0;
    for i=1:6
        u(i) = bytes_to_unit(kbytes(off+(1:16)));
        off  = off + 16;
    end

    % ---------- 映射到 Logistic 参数 ----------
    % r in (3.9, 4.0)  — 混沌区
    r1 = 3.9 + 0.1*u(1);
    r2 = 3.9 + 0.1*u(2);
    r3 = 3.9 + 0.1*u(3);

    x01 = clamp01(u(4)*0.999999 + 1e-7);
    x02 = clamp01(u(5)*0.999999 + 1e-7);
    x03 = clamp01(u(6)*0.999999 + 1e-7);

    Lwarm = 1000;          % 预迭代
    Lneed = N + Lwarm;

    % ---------- 迭代 3 轨 Logistic ----------
    x1 = iterate_logistic(x01, r1, Lneed); x1 = x1(Lwarm+1:end);
    x2 = iterate_logistic(x02, r2, Lneed); x2 = x2(Lwarm+1:end);
    x3 = iterate_logistic(x03, r3, Lneed); x3 = x3(Lwarm+1:end);

    % ---------- 量化 ----------
    Kdif_vec = uint8( floor(mod(x1 * 1.0e14, 256)) );   % 0..255
    ranks    = rank_from_sequence(x2);                  % 1..N
    Kscr_vec = uint8( mod(double(ranks), 256) );        % 0..255
    Kdna_vec = uint8( mod(floor(x3 * 1.0e6), 8) );      % 0..7

    % ---------- reshape ----------
    key_stream_diffusion  = reshape(Kdif_vec, [H,W,C]);
    key_stream_scrambling = reshape(Kscr_vec, [H,W,C]);
    key_stream_dna        = reshape(Kdna_vec, [H,W,C]);
end

% =================== Logistic 迭代 ===================

function x = iterate_logistic(x0, r, L)
    x = zeros(L,1);
    x(1) = clamp01(x0);
    for n = 1:(L-1)
        x(n+1) = r * x(n) * (1 - x(n));
    end
    % 保证数值稳定在 (0,1)
    x = clamp01(x);
end

% =================== KDF & 工具函数 ===================

function u = bytes_to_unit(b16)
    if numel(b16) ~= 16
        error('bytes_to_unit: need 16 bytes');
    end
    v = double(b16) * (256.^(15:-1:0)).';
    u = v / (256^16 - 1);
end

function okm = kdf_expand(ikm, info, L)
    if nargin < 2 || isempty(info)
        info = 'KDF-EXPAND';
    end
    if nargin < 3
        L = 32;
    end
    T = uint8([]);
    out = uint8([]);
    ctr = uint8(1);
    while numel(out) < L
        blk = sha256_bytes([uint8(info), ikm, T, ctr]);
        out = [out, blk]; %#ok<AGROW>
        T = blk;
        ctr = uint8(mod(double(ctr)+1,256));
        if ctr==0, ctr=uint8(1); end
    end
    okm = out(1:L);
end

function h = sha256_bytes(msg_bytes)
    md = java.security.MessageDigest.getInstance('SHA-256');
    md.update(uint8(msg_bytes));
    h = typecast(md.digest(),'uint8'); 
    h = h(:).';
end

function y = rowu8(x)
    if ~isa(x,'uint8')
        x = uint8(x);
    end
    y = x(:).';
end

function b = to_bytes(x)
    if isempty(x)
        b = uint8([]);
        return;
    end
    if isa(x,'uint8')
        b = x(:).';
        return;
    end
    if ischar(x) || isstring(x)
        b = uint8(char(x));
        b = b(:).';
        return;
    end
    if isnumeric(x)
        b = uint8(mod(double(x(:).'), 256));
        return;
    end
    % fallback
    b = uint8(char(string(x)));
    b = b(:).';
end

function ranks = rank_from_sequence(x)
    [~, perm] = sort(x, 'ascend');
    ranks = zeros(size(perm));
    ranks(perm) = 1:numel(perm);
    ranks = ranks(:);
end

function y = clamp01(x)
    y = min(max(x, 0+eps), 1-eps);
end
