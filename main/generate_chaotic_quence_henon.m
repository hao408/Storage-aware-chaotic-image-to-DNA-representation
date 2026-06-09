function [key_stream_diffusion, key_stream_scrambling, key_stream_dna] = generate_chaotic_quence_henon(img, master_key, salt)
% HENON-based chaotic key streams: diffusion, scrambling, DNA rule
% 使用 1 组 Hénon 映射生成 (x,y,z) 三路序列：
%   x -> diffusion,  y -> scrambling,  z=frac(x+y+h) -> DNA 规则。
%
% Usage:
%   [Kdif, Kscr, Kdna] = generate_chaotic_quence_henon(img)
%   [Kdif, Kscr, Kdna] = generate_chaotic_quence_henon(img, master_key, salt)

    if ~isa(img,'uint8') || ndims(img)~=3
        error('img must be HxWxC uint8.');
    end
    [H,W,C] = size(img);
    N = H*W*C;

    if nargin < 2 || isempty(master_key)
        master_key = 'HENON-DEFAULT-KEY';
    end
    img_hash = sha256_bytes(uint8(img(:)));
    if nargin < 3 || isempty(salt)
        salt = img_hash;
    end

    % ---------- HKDF-like 派生参数 ----------
    seed_in   = [rowu8(to_bytes(master_key)), rowu8(img_hash), rowu8(to_bytes(salt))];
    needBytes = 6*16;                          % a,b,x0,y0,hshift 等
    kbytes    = kdf_expand(seed_in, 'HENON-1SYSTEM', needBytes);

    u  = zeros(1,6); off = 0;
    for i=1:6
        u(i) = bytes_to_unit(kbytes(off+(1:16)));
        off  = off + 16;
    end

    % a,b 取在典型混沌区附近
    a = 1.2 + 0.4*u(1);    % ≈ [1.2,1.6]，常用 1.4
    b = 0.2 + 0.2*u(2);    % ≈ [0.2,0.4]，常用 0.3

    x0 = 2*u(3) - 1;       % [-1,1]
    y0 = 2*u(4) - 1;
    hshift = u(5);         % [0,1) 用于 z = frac(x+y+h)

    Lwarm = 1000;
    Lneed = N + Lwarm;

    % ---------- 迭代 Hénon ----------
    [x, y] = iterate_henon(x0, y0, a, b, Lneed);
    x = x(Lwarm+1:end);
    y = y(Lwarm+1:end);
    z = frac(x + y + hshift);

    % ---------- 量化 ----------
    Kdif_vec = uint8( floor(mod(x * 1.0e14, 256)) );   % 0..255
    ranks    = rank_from_sequence(y);
    Kscr_vec = uint8( mod(double(ranks), 256) );       % 0..255
    Kdna_vec = uint8( mod(floor(z * 1.0e6), 8) );      % 0..7

    key_stream_diffusion  = reshape(Kdif_vec, [H,W,C]);
    key_stream_scrambling = reshape(Kscr_vec, [H,W,C]);
    key_stream_dna        = reshape(Kdna_vec, [H,W,C]);
end

% =================== Hénon 迭代 ===================

function [x,y] = iterate_henon(x0, y0, a, b, L)
    x = zeros(L,1);
    y = zeros(L,1);
    x(1) = x0;
    y(1) = y0;
    for n = 1:(L-1)
        x(n+1) = 1 - a*x(n)^2 + y(n);
        y(n+1) = b*x(n);
    end
    % 简单归一到 (0,1) 区间，保证后续量化
    x = frac(x*0.5 + 0.5);
    y = frac(y*0.5 + 0.5);
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
    b = uint8(char(string(x)));
    b = b(:).';
end

function ranks = rank_from_sequence(x)
    [~, perm] = sort(x, 'ascend');
    ranks = zeros(size(perm));
    ranks(perm) = 1:numel(perm);
    ranks = ranks(:);
end

function z = frac(x)
    z = x - floor(x);
end
