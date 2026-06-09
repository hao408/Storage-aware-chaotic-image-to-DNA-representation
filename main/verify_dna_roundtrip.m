function verify_dna_roundtrip(original_image, key_stream_dna, dna_complementary_principle)
%VERIFY_DNA_ROUNDTRIP  在小区域上验证 encode_img / decode_img 的可逆性
%
% 用法：
%   verify_dna_roundtrip(original_image, key_stream_dna, dna_complementary_principle);

    fprintf('Verifying DNA encoding/decoding on small region...\n');

    % 1. 取一个小 patch 做测试
    patch_H = min(5, size(original_image,1));
    patch_W = min(5, size(original_image,2));
    test_region = original_image(1:patch_H, 1:patch_W, :);
    test_key_dna = key_stream_dna(1:patch_H, 1:patch_W, :);

    % 2. 编码 + 解码
    test_encoded = encode_img(test_region, dna_complementary_principle, test_key_dna);
    test_decoded = decode_img(test_encoded, dna_complementary_principle, test_key_dna);

    % 3. 比较差异
    diff_dna = sum(abs(double(test_region(:)) - double(test_decoded(:))));
    if diff_dna == 0
        fprintf('✅ DNA encoding/decoding test: PASS\n');
    else
        fprintf('❌ DNA encoding/decoding test: FAIL (sum diff = %d)\n', diff_dna);
        diff_mask = test_region ~= test_decoded;
        [row, col, ch] = ind2sub(size(test_region), find(diff_mask, 1));
        if ~isempty(row)
            fprintf(['First difference: original=%d, decoded=%d ' ...
                     'at (row=%d, col=%d, ch=%d)\n'], ...
                    test_region(row, col, ch), ...
                    test_decoded(row, col, ch), ...
                    row, col, ch);
        end
    end
end
