function [key_stream_diffusion, key_stream_scrambling, key_stream_dna] = generate_chaotic_quence_6d(img, master_key, salt)
% GENERATE_CHAOTIC_QUENCE_6D
% 6 维离散混沌系统产生三路密钥流：扩散 / 置乱 / DNA 规则。
%
% 接口与其它 generate_chaotic_quence_* 完全一致：
%   [Kdif, Kscr, Kdna] = generate_chaotic_quence_6d(img)
%   [Kdif, Kscr, Kdna] = generate_chaotic_quence_6d(img, master_key, salt)
%
% 如果你有自己原来的 6D 混沌方程，只需要在下面的 iterate_6d_map
% 里把更新公式换成你的那一套即可，输入输出保持不变。

    if ~isa(img,'uint8') || ndims(img)~=3
        error('img must be HxWxC uint8.');
    end
    [H,W,C] = size(img);
    N = H*W*C;

    % ---------- 默认密钥 ----------
    if nargin < 2 || isempty(master_key)
        master_key = 'SIXDIM-DEFAULT-KEY';
    end
    img_hash = sha256_bytes(uint8(img(:)));
    if nargin < 3 || isempty(salt)
        salt = img_hash;
    end

    % ---------- HKDF-like 派生参数 ----------
    seed_in   = [rowu8(to_bytes(master_key)), rowu8(img_hash), rowu8(to_bytes(salt))];
    % 需要 12 个参数（a1..a6, b1..b6）+ 6 个初值共 18 个 u，预留多点
    needBytes = 24*16;
    kbytes    = kdf_expand(seed_in, 'SIXD-MAP', needBytes);

    u  = zeros(1,24); off = 0;
    for i=1:24
        u(i) = bytes_to_unit(kbytes(off+(1:16)));
        off  = off + 16;
    end

    % 参数映射：a_i in (1.2,1.9), b_i in (0.5,1.0)
    par = struct();
    for i = 1:6
        par.(['a',num2str(i)]) = 1.2 + 0.7*u(i);
        par.(['b',num2str(i)]) = 0.5 + 0.5*u(6+i);
    end

    % 初值 x1..x6 in (0,1)
    x0 = zeros(6,1);
    for i = 1:6
        x0(i) = clamp01(u(12+i));
    end

    hshift = u(24);   % 用于第三路 z = frac(x1+x2+h)

    Lwarm = 1000;
    Lneed = N + Lwarm;

    % ---------- 迭代 6D 系统 ----------
    X = iterate_6d_map(x0, par, Lneed);   % 返回 L×6 矩阵
    X = X(Lwarm+1:end, :);
    x1 = X(:,1);
    x2 = X(:,2);
    x3 = X(:,3);

    % 第三路：x1 与 x2 组合后平移
    z  = frac(x1 + x2 + hshift);

    % ---------- 量化 ----------
    Kdif_vec = uint8( floor(mod(x1 * 1.0e14, 256)) );   % 0..255
    ranks    = rank_from_sequence(x2);
    Kscr_vec = uint8( mod(double(ranks), 256) );        % 0..255
    Kdna_vec = uint8( mod(floor(z * 1.0e6), 8) );       % 0..7

    key_stream_diffusion  = reshape(Kdif_vec, [H,W,C]);
    key_stream_scrambling = reshape(Kscr_vec, [H,W,C]);
    key_stream_dna        = reshape(Kdna_vec, [H,W,C]);
end

% =====================================================================
% 6D 耦合映射迭代
% 如果你有自己的 6D 系统，把这里改成你的原始更新方程即可。
% =====================================================================
function X = iterate_6d_map(x0, par, L)
    X = zeros(L, 6);
    X(1,:) = x0(:).';
    for n = 1:(L-1)
        x1 = X(n,1); x2 = X(n,2); x3 = X(n,3);
        x4 = X(n,4); x5 = X(n,5); x6 = X(n,6);

        % 一个示例 6 维耦合映射（全部在 (0,1) 上）
        y1 = par.a1*x1 + par.b1*sin(pi*(x2 + x6));
        y2 = par.a2*x2 + par.b2*sin(pi*(x3 + x1));
        y3 = par.a3*x3 + par.b3*sin(pi*(x4 + x2));
        y4 = par.a4*x4 + par.b4*sin(pi*(x5 + x3));
        y5 = par.a5*x5 + par.b5*sin(pi*(x6 + x4));
        y6 = par.a6*x6 + par.b6*sin(pi*(x1 + x5));

        X(n+1,1) = frac(y1);
        X(n+1,2) = frac(y2);
        X(n+1,3) = frac(y3);
        X(n+1,4) = frac(y4);
        X(n+1,5) = frac(y5);
        X(n+1,6) = frac(y6);
    end
end

% =====================================================================
% 通用工具函数（和其它 generate_chaotic_quence_* 保持一致）
% =====================================================================
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

function y = clamp01(x)
    y = min(max(x, 0+eps), 1-eps);
end
