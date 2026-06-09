function [key_stream_diffusion, key_stream_scrambling, key_stream_dna] = generate_chaotic_quence(img, master_key, salt)
% IS-map chaotic key streams: diffusion, scrambling, DNA rule (3 tracks + warmup)
% Usage:
%   [Kdif, Kscr, Kdna] = generate_chaotic_quence(img)
%   [Kdif, Kscr, Kdna] = generate_chaotic_quence(img, master_key, salt)
%
% Outputs shape/type match img:
%   Kdif: [H,W,C] uint8 0..255
%   Kscr: [H,W,C] uint8 (store permutation ranks mod 256)
%   Kdna: [H,W,C] uint8 0..7

    % -------- Shapes --------
    if ~isa(img,'uint8') || ndims(img)~=3
        error('img must be HxWxC uint8.');
    end
    [H, W, C] = size(img);
    N  = H * W * C;

    % -------- Parameters --------
    Lw = 1000; % warm-up steps
    if nargin < 2 || isempty(master_key)
        master_key = 'ISMAP-DEFAULT-KEY';
    end

    % salt: default bind to image hash
    img_hash = sha256_bytes(uint8(img(:))); % 1x32 row
    if nargin < 3 || isempty(salt)
        salt = img_hash;
    end

    % -------- KDF (HKDF-like) --------
    seed_in   = [rowu8(to_bytes(master_key)), rowu8(img_hash), rowu8(to_bytes(salt))];
    needBytes = 3*(16+16); % per track: 16B params + 16B x0
    kbytes    = hkdf_sha256(seed_in, 'ISMAP-3TRACKS', needBytes);

    % split into 3 tracks
    off = 0;
    [theta1, x01] = derive_params(kbytes(off+(1:16)), kbytes(off+(17:32))); off = off + 32;
    [theta2, x02] = derive_params(kbytes(off+(1:16)), kbytes(off+(17:32))); off = off + 32;
    [theta3, x03] = derive_params(kbytes(off+(1:16)), kbytes(off+(17:32))); off = off + 32;

    % -------- Iterate IS-map and discard warm-up --------
    x1 = iterate_is_map(x01, theta1, N + Lw); x1 = x1(Lw+1:end);
    x2 = iterate_is_map(x02, theta2, N + Lw); x2 = x2(Lw+1:end);
    x3 = iterate_is_map(x03, theta3, N + Lw); x3 = x3(Lw+1:end);

    % -------- Quantization --------
    % diffusion bytes
    Kdif_vec = uint8( floor(mod(x1 * 1.0e14, 256)) );
    % scrambling: ranks 1..N -> uint8 mask (mod 256)
    ranks    = rank_from_sequence(x2);
    Kscr_vec = uint8( mod(double(ranks), 256) );
    % dna rule 0..7
    Kdna_vec = uint8( mod(floor(x3 * 1.0e6), 8) );

    % -------- Reshape --------
    key_stream_diffusion  = reshape(Kdif_vec, [H, W, C]);
    key_stream_scrambling = reshape(Kscr_vec, [H, W, C]);
    key_stream_dna        = reshape(Kdna_vec, [H, W, C]);
end

% ===================== Helpers =====================

function y = rowu8(x)
    if ~isa(x,'uint8')
        x = uint8(x);
    end
    y = x(:).';
end

function y = to_bytes(x)
    if isempty(x)
        y = uint8([]);
        return;
    end
    if isa(x,'uint8')
        y = x;
        return;
    end
    if isstring(x) || ischar(x)
        y = uint8(char(x));
        return;
    end
    y = uint8(x);
end

function h = sha256_bytes(msg_bytes)
    md = java.security.MessageDigest.getInstance('SHA-256');
    md.update(uint8(msg_bytes));
    h = typecast(md.digest(),'uint8');
    h = h(:).'; % force row
end

function okm = hkdf_sha256(ikm, info, L)
    if nargin < 2 || isempty(info)
        info = 'ISMAP-HKDF';
    end
    if nargin < 3
        L = 64;
    end
    out = uint8([]);
    T   = uint8([]);
    i   = uint8(1);
    while numel(out) < L
        data = [T, uint8(info), i, sha256_bytes([ikm, T, uint8(info), i])];
        T = sha256_bytes(data);
        out = [out, T]; %#ok<AGROW>
        i = uint8(mod(double(i)+1,256));
        if i==0
            i = uint8(1);
        end
    end
    okm = out(1:L);
end

function [theta, x0] = derive_params(b1, b2)
    u1 = bytes_to_unit(b1); % in [0,1)
    u2 = bytes_to_unit(b2);
    theta.a      = 0.5 + 3.49 * u1;           % (0.5,3.99)
    theta.dither = 1e-12 + 1e-10 * u2;        % small dither
    x0 = clamp01(u2 * 0.999999 + 1.0e-7);     % (0,1)
end

function u = bytes_to_unit(b)
    % map 16 bytes -> [0,1)
    b = rowu8(b);
    if numel(b) < 16
        b = [b, zeros(1,16-numel(b),'uint8')];
    else
        b = b(1:16);
    end
    % accumulate as double to avoid uint64 corner parsing/overflow issues
    v = 0.0;
    for k = 1:16
        v = v*256.0 + double(b(k));
    end
    u = v / (256.0^16); % in [0,1)
    if u == 0
        u = rand();
    end
end

function x = iterate_is_map(x0, theta, L)
    % x_{n+1} = mod( a*sin(pi*x_n) + (4-a)*x_n*(1-x_n), 1 )
    a  = theta.a;
    dz = theta.dither;
    x  = zeros(L,1);
    x(1) = clamp01(x0);
    for n = 1:(L-1)
        xn = x(n);
        xn1 = a*sin(pi*xn) + (4.0 - a)*xn*(1.0 - xn);
        % reduce to (0,1)
        xn1 = xn1 - floor(xn1);
        if xn1 <= 0 || xn1 >= 1
            xn1 = xn + dz;
            xn1 = xn1 - floor(xn1);
        end
        if xn1 <= 0
            xn1 = eps;
        elseif xn1 >= 1
            xn1 = 1 - eps;
        end
        x(n+1) = xn1;
    end
end

function ranks = rank_from_sequence(x)
    % stable ranking: return rank 1..N for each element
    [~, perm] = sort(x, 'ascend'); % perm: positions in sorted order
    ranks = zeros(size(perm));
    ranks(perm) = 1:numel(perm);
    ranks = ranks(:);
end

function y = clamp01(x)
    y = min(max(x, 0+eps), 1-eps);
end
