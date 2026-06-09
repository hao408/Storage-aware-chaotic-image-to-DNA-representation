function summary = export_oligos(dna_matrix, out_basename, varargin)
%EXPORT_OLIGOS  Export DNA-encoded image (A/C/G/T matrix) to oligo FASTA & CSV.
%
%   SUMMARY = EXPORT_OLIGOS(DNA_MATRIX, OUT_BASENAME, 'Name',Value,...)
%
%   Inputs:
%     dna_matrix   : HxW (or HxWxC) char array of 'A','C','G','T' produced by your encode_img.
%     out_basename : output path WITHOUT extension. e.g. 'out/lena_cipher_dna'
%
%   Name-Value options (all optional):
%     'OligoLength'      : target length per oligo (including headers & CRC). Default 220.
%     'PrefixSync'       : fixed sync motif (char), default 'ACGTAC' (6 nt).
%     'BatchID'          : any string/char to tag the pool; mapped to 8 nt via hash. Default 'POOL-001'.
%     'IndexFieldLen'    : length (nt) for index and total fields (each). Default 6 (base-4).
%     'CRCFieldLen'      : length (nt) of CRC16 field (base-4). Default 8 (encodes 16 bits -> 8 nt).
%     'Order'            : 'row-major' or 'channel-stacked' (for HxWx3). Default 'row-major'.
%     'LineWrap'         : FASTA line wrap (0=no wrap). Default 0.
%
%   Outputs:
%     summary : struct with counts and basic QC (GC%, homopolymer stats), and file paths.
%
%   Files written:
%     <out_basename>.fasta : with headers >OLI_<idx> containing batch, index, total, CRC
%     <out_basename>.csv   : columns: index,sequence,gc_percent,max_homopolymer,payload_start,payload_len
%
%   Example:
%     % Suppose encoded_image is your DNA char matrix from encode_img (HxW x 3 allowed)
%     S = export_oligos(encoded_image, 'exports/lena_cipher_dna', ...
%                       'OligoLength', 136, 'BatchID','EXP-2025-10-29', 'Order','row-major');
%     disp(S)
%
%   Notes:
%     1) 本版本不主动“改写”碱基以满足 GC/同碱基约束，而是完整记录 QC 指标；若需联合约束设计，
%        可在切片前对 payload 调用你现有的约束器或我再给你一个可逆的受约束规则选择模块。
%     2) CRC16 使用标准多项式 0x1021，初值 0xFFFF；编码为 base-4（2bit/nt）→ 8 nt。
%
%   Author: d’s assistant
%   Date  : 2025-10-29

% ---------- Parse options ----------
p = inputParser;
addParameter(p,'OligoLength',220,@(x)isnumeric(x)&&isscalar(x)&&x>=60);
addParameter(p,'PrefixSync','ACGTAC',@(s)ischar(s)||isstring(s));
addParameter(p,'BatchID','POOL-001',@(s)ischar(s)||isstring(s));
addParameter(p,'IndexFieldLen',6,@(x)isnumeric(x)&&isscalar(x)&&x>=4);
addParameter(p,'CRCFieldLen',8,@(x)isnumeric(x)&&isscalar(x)&&x>=6);
addParameter(p,'Order','row-major',@(s)ischar(s)||isstring(s));
addParameter(p,'LineWrap',0,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
parse(p,varargin{:});
opt = p.Results;

% ---------- Flatten DNA matrix to 1-D payload stream ----------
dna_vec = flatten_dna(dna_matrix, string(opt.Order));

% ---------- Build constant headers ----------
sync_nt   = upper(char(opt.PrefixSync(:).'));
batch8nt  = batch_to_8nt(char(opt.BatchID)); % 8 nt from SHA-256
idx_len   = opt.IndexFieldLen;   % nt for index (base-4)
tot_len   = opt.IndexFieldLen;   % nt for total (base-4)
crc_len   = opt.CRCFieldLen;     % 8 nt (16 bits)
fixed_hdr_len = length(sync_nt) + 8 + idx_len + tot_len; % total header nts
fixed_tail_len = crc_len;

oligo_len = opt.OligoLength;
payload_room = oligo_len - fixed_hdr_len - fixed_tail_len;
if payload_room <= 0
    error('OligoLength too small. Need > %d nt', fixed_hdr_len + fixed_tail_len + 1);
end

% ---------- Slice payload into oligos ----------
payload_slices = slice_payload(dna_vec, payload_room);
total_oligos = numel(payload_slices);

% ---------- Encode each oligo with addressing & CRC ----------
headers = cell(total_oligos,1);
seqs    = cell(total_oligos,1);
gc_list = zeros(total_oligos,1);
hp_list = zeros(total_oligos,1);
pstart  = zeros(total_oligos,1);
plen    = zeros(total_oligos,1);

total_nt_count = length(dna_vec);
cursor = 0;

for k = 1:total_oligos
    pay = payload_slices{k};
    pstart(k) = cursor + 1;
    plen(k)   = length(pay);
    cursor    = cursor + plen(k);

    idx_nt  = base4_encode_uint(k-1, idx_len);                % 0-based index
    tot_nt  = base4_encode_uint(total_oligos-1, tot_len);     % 0-based total
    head    = [sync_nt, batch8nt, idx_nt, tot_nt];

    % CRC over header+payload (in bits), then base4->nt (8 nt)
    crc_bits = crc16_ccitt_bits([nt_to_bits(head), nt_to_bits(pay)]);
    crc_nt   = bits_to_nt(crc_bits);  % 16 bits -> 8 nt

    seq      = [head, pay, crc_nt];

    % QC
    [gc, maxhp] = gc_homopolymer_stats(seq);
    gc_list(k)  = gc;
    hp_list(k)  = maxhp;

    headers{k} = sprintf('OLI_%05d|BATCH=%s|IDX=%d|TOT=%d|LEN=%d|GC=%.1f|HP=%d', ...
                         k-1, char(opt.BatchID), k-1, total_oligos-1, length(seq), gc, maxhp);
    seqs{k}    = char(seq);
end

% ---------- Write FASTA ----------
fasta_path = [char(out_basename), '.fasta'];
write_fasta(fasta_path, headers, seqs, opt.LineWrap);

% ---------- Write CSV ----------
csv_path = [char(out_basename), '.csv'];
write_csv(csv_path, seqs, gc_list, hp_list, pstart, plen);

% ---------- Summary ----------
summary = struct();
summary.num_oligos = total_oligos;
summary.oligo_length = oligo_len;
summary.payload_room = payload_room;
summary.total_payload_nt = total_nt_count;
summary.gc_mean = mean(gc_list);
summary.homopolymer_max = max(hp_list);
summary.fasta = fasta_path;
summary.csv   = csv_path;

end
% =========================== Helpers ===========================

function dna_vec = flatten_dna(dna_matrix, order)
    dna = upper(dna_matrix);
    if ndims(dna) == 3
        switch lower(order)
            case 'row-major'
                % HxWxC -> concatenate channels in RGB order (row-major inside each)
                [H,W,C] = size(dna);
                buf = cell(1,C);
                for c=1:C
                    buf{c} = reshape(dna(:,:,c).', 1, H*W); % transpose then linearize for row-major
                end
                dna_vec = [buf{:}];
            case 'channel-stacked'
                % Same as above (explicit), kept for clarity
                [H,W,C] = size(dna);
                buf = cell(1,C);
                for c=1:C
                    buf{c} = reshape(dna(:,:,c).', 1, H*W);
                end
                dna_vec = [buf{:}];
            otherwise
                error('Unknown Order option: %s', order);
        end
    else
        % HxW
        [H,W] = size(dna);
        dna_vec = reshape(dna.', 1, H*W); % row-major
    end
end

function slices = slice_payload(vec, chunk_len)
    N = length(vec);
    M = ceil(N / chunk_len);
    slices = cell(M,1);
    s = 1;
    for i=1:M
        e = min(s+chunk_len-1, N);
        slices{i} = vec(s:e);
        s = e + 1;
    end
end

function out8nt = batch_to_8nt(batch_str)
    % Map SHA-256(batch_str) 16 bytes -> 8 nt (base-4, 2 bits per nt)
    b = sha256_bytes(uint8(batch_str));
    b = b(1:16);                % 16B -> 128 bits -> 64 nt if fully used, we take 8 nt = 16 bits
    bits = bytes_to_bits(b(1:2)); % 16 bits only
    out8nt = bits_to_nt(bits);    % 8 nt
end

function nt = base4_encode_uint(u, nt_len)
    % Encode non-negative integer u into nt_len nucleotides (base-4, big-endian)
    % digits: 0->A, 1->C, 2->G, 3->T
    digs = zeros(1, nt_len);
    x = double(u);
    for i = nt_len:-1:1
        digs(i) = mod(x,4);
        x = floor(x/4);
    end
    map = 'ACGT';
    nt = map(digs+1);
end

function bits = nt_to_bits(nt)
    % A->00 C->01 G->10 T->11
    nt = upper(nt);
    map = containers.Map({'A','C','G','T'},{[0 0],[0 1],[1 0],[1 1]});
    bits = zeros(1, 2*numel(nt));
    p = 1;
    for i=1:numel(nt)
        v = map(nt(i));
        bits(p:p+1) = v;
        p = p + 2;
    end
end

function nt = bits_to_nt(bits)
    % bits length must be even
    if mod(numel(bits),2)~=0
        error('bits_to_nt: bits length must be even');
    end
    map = ['A','C','G','T']; % 00 01 10 11
    nt = repmat('A', 1, numel(bits)/2);
    p = 1;
    for i=1:2:numel(bits)
        d = bits(i)*2 + bits(i+1);
        nt(p) = map(d+1);
        p = p + 1;
    end
end

function bits = bytes_to_bits(by)
    % by: uint8 row vector
    if ~isa(by,'uint8'); by = uint8(by); end
    bits = zeros(1, numel(by)*8);
    k = 1;
    for i=1:numel(by)
        v = by(i);
        for b=7:-1:0
            bits(k) = bitget(v, b+1);
            k = k+1;
        end
    end
end

function b = bits_to_bytes(bits)
    if mod(numel(bits),8)~=0
        error('bits_to_bytes: length must be multiple of 8');
    end
    n = numel(bits)/8;
    b = zeros(1,n,'uint8');
    k = 1;
    for i=1:n
        v = uint8(0);
        for bitpos=7:-1:0
            v = bitor(bitshift(v,1), uint8(bits(k)));
            k = k+1;
        end
        b(i) = v;
    end
end

function [gc, maxhp] = gc_homopolymer_stats(seq_nt)
    % GC %
    gc = 100 * (sum(seq_nt=='G') + sum(seq_nt=='C')) / numel(seq_nt);
    % longest homopolymer
    maxhp = 1;
    runlen = 1;
    for i=2:numel(seq_nt)
        if seq_nt(i)==seq_nt(i-1)
            runlen = runlen + 1;
            if runlen > maxhp, maxhp = runlen; end
        else
            runlen = 1;
        end
    end
end

function crc_bits = crc16_ccitt_bits(bits)
    % Compute CRC-16/CCITT-FALSE over a bit vector (MSB-first per byte semantics)
    % Polynomial 0x1021, init 0xFFFF, no XORout. We interpret bits in big-endian byte order.
    % First convert bits->bytes with padding to whole bytes (pad zeros at end).
    if mod(numel(bits),8)~=0
        pad = 8 - mod(numel(bits),8);
        bits = [bits, zeros(1,pad)];
    end
    data = bits_to_bytes(bits);
    crc = uint16(hex2dec('FFFF'));
    poly = uint16(hex2dec('1021'));
    for i=1:numel(data)
        crc = bitxor(crc, bitshift(uint16(data(i)),8));
        for b=1:8
            if bitand(crc, uint16(hex2dec('8000'))) ~= 0
                crc = bitxor(bitshift(crc,1), poly);
            else
                crc = bitshift(crc,1);
            end
        end
    end
    % Convert 16-bit value to 16 bits MSB-first
    out = zeros(1,16);
    for k=1:16
        out(k) = bitget(crc, 17-k); % MSB first
    end
    crc_bits = out;
end

function write_fasta(path, headers, seqs, linewrap)
    fid = fopen(path,'w');
    if fid<0, error('Cannot open %s for writing', path); end
    cleaner = onCleanup(@()fclose(fid));
    for i=1:numel(seqs)
        fprintf(fid, '>%s\n', headers{i});
        if linewrap>0
            s = seqs{i};
            for k=1:linewrap:numel(s)
                fprintf(fid, '%s\n', s(k:min(k+linewrap-1, numel(s))));
            end
        else
            fprintf(fid, '%s\n', seqs{i});
        end
    end
end

function write_csv(path, seqs, gc_list, hp_list, pstart, plen)
    fid = fopen(path,'w');
    if fid<0, error('Cannot open %s for writing', path); end
    cleaner = onCleanup(@()fclose(fid));
    fprintf(fid, 'index,sequence,gc_percent,max_homopolymer,payload_start,payload_length\n');
    for i=1:numel(seqs)
        fprintf(fid, '%d,%s,%.2f,%d,%d,%d\n', i-1, seqs{i}, gc_list(i), hp_list(i), pstart(i), plen(i));
    end
end

function h = sha256_bytes(msg_bytes)
    md = java.security.MessageDigest.getInstance('SHA-256');
    md.update(uint8(msg_bytes));
    h = typecast(md.digest(),'uint8');
end
