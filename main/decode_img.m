function decoded_img = decode_img(encoded_img, complementary_principle, key_stream_dna)
%DECODE_IMG  Inverse DNA rule-based substitution layer (decryption).
%
%   decoded_img = DECODE_IMG(encoded_img, complementary_principle, key_stream_dna)
%
%   encoded_img            : H×W or H×W×C uint8，DNA规则加密后的像素域图像
%   complementary_principle: 8×4 uint8, 与加密端相同的 DNA 规则表 (ASCII A/C/G/T)
%   key_stream_dna         : 与加密端一致的规则密钥流
%
%   decoded_img            : 还原到 DNA 层之前（扩散层之后）的像素域图像

    [M, N, C] = size(encoded_img);
    decoded_img = uint8(zeros(M, N, C));

    %-------- 1. ASCII 规则 -> 数字规则 0..3 --------
    rules_num = zeros(size(complementary_principle), 'uint8');  % 8×4

    for r = 1:size(complementary_principle,1)
        for d = 0:3
            base = complementary_principle(r, d+1);
            switch base
                case 65  % 'A'
                    rules_num(r, d+1) = uint8(0);
                case 67  % 'C'
                    rules_num(r, d+1) = uint8(1);
                case 71  % 'G'
                    rules_num(r, d+1) = uint8(2);
                case 84  % 'T'
                    rules_num(r, d+1) = uint8(3);
                otherwise
                    error('DNA规则中出现非法字符 (不是 A/C/G/T 的 ASCII)。');
            end
        end
    end

    %-------- 2. 预计算逆规则表 inv_rules --------
    % inv_rules(r, y+1) = x  使得 rules_num(r, x+1) = y
    inv_rules = zeros(size(rules_num), 'uint8');
    for r = 1:size(rules_num,1)
        for x = 0:3
            y = rules_num(r, x+1);      % y = σ_r(x)
            inv_rules(r, y+1) = uint8(x); % σ_r⁻¹(y) = x
        end
    end

    %-------- 3. 对每个像素的 4 个2-bit 组做 σ⁻¹ --------
    for ch = 1:C
        channel_img = encoded_img(:,:,ch);
        key_channel = key_stream_dna(:,:,ch);
        decoded_channel = uint8(zeros(M, N));
        
        for i = 1:M
            for j = 1:N
                % 选择规则索引 1..8
                rule_idx = mod(double(key_channel(i, j)), 8) + 1;
                inv_rule = inv_rules(rule_idx, :);  % 1×4, 元素 0..3
                
                pixel = channel_img(i, j);
                decoded_pixel = uint8(0);
                
                for k = 1:4
                    % 提取编码后的2位
                    bit_pos1 = 2*k-1;
                    bit_pos2 = 2*k;
                    bit1 = bitget(pixel, bit_pos1);
                    bit2 = bitget(pixel, bit_pos2);
                    encoded_num = bit1 + 2*bit2; % 当前 σ_r(x) 的值 0-3
                    
                    % 应用逆映射 σ_r⁻¹: 0..3 -> 0..3
                    original_num = inv_rule(encoded_num + 1);
                    
                    % 还原原始 2 位
                    orig_bit1 = bitget(original_num, 1);
                    orig_bit2 = bitget(original_num, 2);
                    decoded_pixel = bitset(decoded_pixel, bit_pos1, orig_bit1);
                    decoded_pixel = bitset(decoded_pixel, bit_pos2, orig_bit2);
                end
                
                decoded_channel(i, j) = decoded_pixel;
            end
        end
        decoded_img(:,:,ch) = decoded_channel;
    end
end
