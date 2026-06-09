function encoded_img = encode_img(original_img, complementary_principle, key_stream_dna)
%ENCODE_IMG  DNA rule-based substitution layer (encryption).
%
%   encoded_img = ENCODE_IMG(original_img, complementary_principle, key_stream_dna)
%
%   original_img          : H×W or H×W×C uint8 图像（已经做完置乱+扩散）
%   complementary_principle : 8×4 uint8, 每行是一条 DNA 互补规则 (ASCII A/C/G/T)
%   key_stream_dna        : H×W×C 或 H×W 的密钥流，每个位置 0..255，最终 mod 8 选规则
%
%   encoded_img           : 同尺寸 uint8 图像，完成 DNA 规则替换后的像素域密文

    [M, N, C] = size(original_img);
    encoded_img = uint8(zeros(M, N, C));

    %-------- 1. 将 ASCII 碱基规则转成数字规则 0..3 --------
    % 0=A, 1=C, 2=G, 3=T
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

    %-------- 2. 对每个像素的 4 个 2-bit 组做 σ 映射 --------
    for ch = 1:C
        channel_img = original_img(:,:,ch);
        key_channel = key_stream_dna(:,:,ch);
        encoded_channel = uint8(zeros(M, N));
        
        for i = 1:M
            for j = 1:N
                % 选择规则索引 1..8
                rule_idx = mod(double(key_channel(i, j)), 8) + 1;
                rule_num = rules_num(rule_idx, :);   % 1×4, 元素 0..3
                
                pixel = channel_img(i, j);
                encoded_pixel = uint8(0);
                
                % 处理每个2位组 (4个 2-bit)
                for k = 1:4
                    % 提取原始 2 位：bit1 为 LSB，bit2 为次低位
                    bit_pos1 = 2*k-1;
                    bit_pos2 = 2*k;
                    bit1 = bitget(pixel, bit_pos1);
                    bit2 = bitget(pixel, bit_pos2);
                    original_num = bit1 + 2*bit2; % 0-3
                    
                    % σ_r: 0..3 -> 0..3
                    mapped_num = rule_num(original_num + 1);  % 0..3
                    
                    % 写回到编码像素的 2 位
                    mapped_bit1 = bitget(mapped_num, 1);
                    mapped_bit2 = bitget(mapped_num, 2);
                    encoded_pixel = bitset(encoded_pixel, bit_pos1, mapped_bit1);
                    encoded_pixel = bitset(encoded_pixel, bit_pos2, mapped_bit2);
                end
                
                encoded_channel(i, j) = encoded_pixel;
            end
        end
        encoded_img(:,:,ch) = encoded_channel;
    end
end
