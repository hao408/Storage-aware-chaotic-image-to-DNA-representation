function verify_dna_rules(dna_complementary_principle)
%VERIFY_DNA_RULES  验证 DNA 编码/解码规则是否自洽
%
% 用法：
%   verify_dna_rules(dna_complementary_principle);

    fprintf('Testing DNA encoding/decoding rules...\n');
    test_dna_encoding();

    function test_dna_encoding()
    % 测试DNA编码解码的正确性
    complementary_principle = dna_complementary_principle;
    
    % 测试所有可能的2位组合
    test_values = [0, 85, 170, 255]; % 00, 01010101, 10101010, 11111111
    
    all_pass = true;
    
    for rule_idx = 1:8
        fprintf('\n=== Testing rule %d ===\n', rule_idx);
        rule = complementary_principle(rule_idx, :);
        fprintf('Rule: A(%d)->%d, C(%d)->%d, G(%d)->%d, T(%d)->%d\n', ...
                65, rule(1), 67, rule(2), 71, rule(3), 84, rule(4));
        
        key_test = (rule_idx - 1) * ones(1, 1, 3, 'uint8');
        
        rule_pass = true;
        for test_val = test_values
            test_img = repmat(uint8(test_val), [1, 1, 3]);
            
            encoded = encode_img(test_img, complementary_principle, key_test);
            decoded = decode_img(encoded, complementary_principle, key_test);
            
            if decoded(1,1,1) ~= test_val
                fprintf('  ❌ FAIL: %d -> %d -> %d\n', test_val, encoded(1,1,1), decoded(1,1,1));
                rule_pass = false;
                all_pass = false;
            else
                fprintf('  ✅ PASS: %d -> %d -> %d\n', test_val, encoded(1,1,1), decoded(1,1,1));
            end
        end
        
        if rule_pass
            fprintf('Rule %d: ✅ ALL PASS\n', rule_idx);
        else
            fprintf('Rule %d: ❌ SOME FAIL\n', rule_idx);
        end
    end
    
    if all_pass
        fprintf('\n🎉 ALL RULES PASSED! DNA encoding/decoding is working correctly.\n');
    else
        fprintf('\n❌ SOME RULES FAILED. Please check the implementation.\n');
    end
end

end



