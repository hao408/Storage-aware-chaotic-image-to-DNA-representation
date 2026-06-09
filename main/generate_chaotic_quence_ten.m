function [key_stream_diffusion, key_stream_scrambling, key_stream_dna] = generate_chaotic_quence_ten(img, master_key, salt)
% CIS-T: Coupled IS-map <-> Skew-Tent (2D, bi-directional diff-coupling)
% Outputs:
%   key_stream_diffusion  [H,W,C] uint8 0..255    (for diffusion, from x)
%   key_stream_scrambling [H,W,C] uint8 0..255    (for scrambling, rank mod 256 from y)
%   key_stream_dna        [H,W,C] uint8 0..7      (DNA rule index from frac(x+y))
%
% Usage:
%   [Kdif, Kscr, Kdna] = generate_chaotic_quence(img)
%   [Kdif, Kscr, Kdna] = generate_chaotic_quence(img, 'myKey', 'mySalt')

    % ---------- checks ----------
    if ~isa(img,'uint8') || ndims(img)~=3
        error('img must be HxWxC uint8.');
    end
    [H,W,C] = size(img); N = H*W*C;

    % ---------- defaults ----------
    if nargin < 2 || isempty(master_key), master_key = 'CIS-T-DEFAULT-KEY'; end

    % bind to image hash by default
    img_hash = sha256_bytes(uint8(img(:)));
    if nargin < 3 || isempty(salt), salt = img_hash; end

    % ---------- derive raw bytes (HKDF-like expander) ----------
    % NOTE: to_bytes() is defined below (now included)
    seed_in  = [rowu8(to_bytes(master_key)), rowu8(img_hash), rowu8(to_bytes(salt))];
    needBytes = 6*16; % a,m,k + x0,y0,h (each 16 bytes -> 0..1)
    kbytes    = kdf_expand(seed_in, 'CIS-T-ISmap-Tent', needBytes);

    u  = zeros(1,6); off = 0;
    for i=1:6
        u(i) = bytes_to_unit(kbytes(off+(1:16)));
        off  = off + 16;
    end
    % map to parameter ranges
    par.a = 0.6 + 3.3*u(1);     % (0.6, 3.9)
    par.m = 0.3 + 0.4*u(2);     % (0.3, 0.7)
    par.k = 0.05 + 0.30*u(3);   % (0.05, 0.35)
    x0     = clamp01(u(4)*0.999999 + 1e-7);
    y0     = clamp01(u(5)*0.999999 + 1e-7);
    hshift = u(6);              % for third track frac(x+y+hshift)

    Lwarm = 1000;      % warm-up steps
    Lneed = N + Lwarm; % total steps

    % ---------- iterate 2D coupled system ----------
    [x, y] = iterate_cis_t(x0, y0, par.a, par.m, par.k, Lneed);
    x = x(Lwarm+1:end);
    y = y(Lwarm+1:end);

    % third track
    z = frac(x + y + hshift);

    % ---------- quantization ----------
    Kdif_vec = uint8( floor(mod(x * 1.0e14, 256)) );    % 0..255
    ranks    = rank_from_sequence(y);                   % 1..N
    Kscr_vec = uint8( mod(double(ranks), 256) );        % 0..255 (mask)
    Kdna_vec = uint8( mod(floor(z * 1.0e6), 8) );       % 0..7

    % ---------- reshape ----------
    key_stream_diffusion  = reshape(Kdif_vec, [H,W,C]);
    key_stream_scrambling = reshape(Kscr_vec, [H,W,C]);
    key_stream_dna        = reshape(Kdna_vec, [H,W,C]);
end

% =====================  Helpers  =====================

function [x,y] = iterate_cis_t(x0, y0, a, m, k, L)
    x = zeros(L,1); y = zeros(L,1);
    x(1) = clamp01(x0); y(1) = clamp01(y0);
    for n = 1:(L-1)
        % IS(x)
        fis = a*sin(pi*x(n)) + (4-a)*x(n)*(1 - x(n));
        % skew-Tent(y)
        if y(n) < m
            Ty = y(n)/m;
        else
            Ty = (1 - y(n)) / (1 - m);
        end
        xn1 = fis + k*(y(n) - x(n));
        yn1 = Ty  + k*(x(n) - y(n));
        xn1 = frac(xn1);
        yn1 = frac(yn1);
        if xn1 <= 0 || xn1 >= 1, xn1 = frac(x(n) + 1e-12); end
        if yn1 <= 0 || yn1 >= 1, yn1 = frac(y(n) + 1e-12); end
        x(n+1) = xn1; y(n+1) = yn1;
    end
end

function r = rank_from_sequence(v)
    [~, perm] = sort(v, 'ascend');    % perm: indices in sorted order
    r = zeros(size(perm));
    r(perm) = 1:numel(perm);
    r = r(:);
end

function u = bytes_to_unit(b)
    b = rowu8(b);
    if numel(b) < 16
        b = [b, zeros(1,16-numel(b),'uint8')];
    else
        b = b(1:16);
    end
    v = 0.0;
    for i=1:16, v = v*256.0 + double(b(i)); end
    u = v / (256.0^16);
    if u==0, u = rand(); end
end

function okm = kdf_expand(ikm, info, L)
    % Lightweight expander using repeated SHA256 chaining (simple and portable)
    if nargin < 2 || isempty(info), info = 'CIS-T-KDF'; end
    if nargin < 3, L = 64; end
    out = uint8([]); T = uint8([]);
    ctr = uint8(1);
    while numel(out) < L
        blk = sha256_bytes([uint8(info), ikm, T, ctr]);
        out = [out, blk]; %#ok<AGROW>
        T = blk;
        ctr = uint8(mod(double(ctr)+1,256)); if ctr==0, ctr=uint8(1); end
    end
    okm = out(1:L);
end

function h = sha256_bytes(msg_bytes)
    md = java.security.MessageDigest.getInstance('SHA-256');
    md.update(uint8(msg_bytes));
    h = typecast(md.digest(),'uint8'); h = h(:).';
end

function y = rowu8(x)
    if ~isa(x,'uint8'), x = uint8(x); end
    y = x(:).';
end

function b = to_bytes(x)
% Robustly turn common types into uint8 row vector
% - char/string -> ASCII bytes
% - numeric -> uint8 with mod 256 (vectorized)
% - uint8 stays as-is
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

function y = clamp01(x)
    y = min(max(x, 0+eps), 1-eps);
end

function z = frac(x)
    z = x - floor(x);
end
