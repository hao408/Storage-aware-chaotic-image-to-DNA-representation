clc; clear; close all;

%% =========================================================
%  Robustness experiment on encoded-image representation
%  说明：
%  1. 本脚本直接对 encoded_image 做 substitution-like 扰动
%  2. 然后按正确逆流程恢复：
%     decode_img -> dediffuse_img -> descramble_img
%% =========================================================

%% =========================
%  输出目录
%% =========================
out_dir = 'encoded_channel_results';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% =========================
%  读取原始图像（与 example_file.m 保持一致）
%% =========================
% 如果你之后换图，这里一起改
original_image = imread('photos\house.tiff');

if ~isa(original_image, 'uint8')
    original_image = uint8(original_image);
end

%% =========================
%  DNA互补规则（与 example_file.m 保持一致）
%% =========================
dna_complementary_principle = [
    84, 71, 67, 65;
    84, 67, 71, 65;
    71, 84, 65, 67;
    67, 84, 65, 71;
    84, 65, 71, 67;
    71, 65, 84, 67;
    67, 65, 84, 71;
    65, 84, 67, 71;
];
dna_complementary_principle = uint8(dna_complementary_principle);

%% =========================
%  重新生成密钥流
%% =========================
[key_stream_diffusion, key_stream_scrambling, key_stream_dna] = generate_chaotic_quence_ten(original_image);

%% =========================
%  重新走一遍正确加密主流程
%% =========================
scrambled_image = scramble_img(original_image, key_stream_scrambling);
diffused_image  = diffuse_img(scrambled_image, key_stream_diffusion);
encoded_image   = encode_img(diffused_image, dna_complementary_principle, key_stream_dna);

%% =========================
%  先验证：无扰动时恢复是否正确
%% =========================
decoded_0    = decode_img(encoded_image, dna_complementary_principle, key_stream_dna);
dediffused_0 = dediffuse_img(decoded_0, key_stream_diffusion);
recovered_0  = descramble_img(dediffused_0, key_stream_scrambling);
recovered_0  = uint8(recovered_0);

mse_0  = immse(recovered_0, original_image);
psnr_0 = psnr(recovered_0, original_image);
ssim_0 = ssim(recovered_0, original_image);

fprintf('\n========== Sanity Check (No Perturbation) ==========\n');
fprintf('MSE  = %.10f\n', mse_0);
fprintf('PSNR = %.10f dB\n', psnr_0);
fprintf('SSIM = %.10f\n', ssim_0);
fprintf('====================================================\n');

% 如果这里都不对，说明主链路还有别的问题
if psnr_0 < 40
    warning(['无扰动恢复的 PSNR 小于 40 dB，说明主链路恢复可能仍有问题。', ...
             '如果 example_file.m 中原图和解密图看起来一致，请优先检查这里是否与 example_file.m 完全同图、同参数。']);
end

%% =========================
%  扰动率设置
%% =========================
noise_rates  = [0, 0.001, 0.005, 0.01];   % 0%, 0.1%, 0.5%, 1%
noise_labels = {'0%', '0.1%', '0.5%', '1%'};
num_cases    = numel(noise_rates);

%% =========================
%  结果存储
%% =========================
MSE_list  = zeros(num_cases, 1);
PSNR_list = zeros(num_cases, 1);
SSIM_list = zeros(num_cases, 1);

%% =========================
%  图像展示
%% =========================
fig1 = figure('Position', [100, 80, 1500, 700]);

subplot(2,3,1);
imshow(original_image);
title('Original image');

%% =========================
%  主循环：对 encoded_image 做 substitution-like 扰动
%% =========================
for k = 1:num_cases
    rate = noise_rates(k);

    fprintf('\n==============================\n');
    fprintf('Running perturbation rate = %.4f (%s)\n', rate, noise_labels{k});
    fprintf('==============================\n');

    % 1) 对 encoded_image 做随机替换扰动
    encoded_image_noisy = apply_uint8_substitution_noise(encoded_image, rate);

    % 2) 正确逆流程恢复
    decoded_image_noisy    = decode_img(encoded_image_noisy, dna_complementary_principle, key_stream_dna);
    dediffused_image_noisy = dediffuse_img(decoded_image_noisy, key_stream_diffusion);
    recovered_image        = descramble_img(dediffused_image_noisy, key_stream_scrambling);
    recovered_image        = uint8(recovered_image);

    % 3) 指标
    mse_val  = immse(recovered_image, original_image);
    psnr_val = psnr(recovered_image, original_image);
    ssim_val = ssim(recovered_image, original_image);

    MSE_list(k)  = mse_val;
    PSNR_list(k) = psnr_val;
    SSIM_list(k) = ssim_val;

    fprintf('MSE  = %.10f\n', mse_val);
    fprintf('PSNR = %.10f dB\n', psnr_val);
    fprintf('SSIM = %.10f\n', ssim_val);

    % 4) 保存单张图
    img_name = sprintf('recovered_encoded_perturb_%s.png', strrep(noise_labels{k}, '%', 'pct'));
    imwrite(recovered_image, fullfile(out_dir, img_name));

    % 5) 展示
    subplot(2,3,k+1);
    imshow(recovered_image);
    title(sprintf('Perturbation %s', noise_labels{k}));
end

saveas(fig1, fullfile(out_dir, 'encoded_perturbation_recovery_overview.png'));

%% =========================
%  保存结果表
%% =========================
ResultTable = table( ...
    string(noise_labels(:)), ...
    noise_rates(:), ...
    MSE_list, ...
    PSNR_list, ...
    SSIM_list, ...
    'VariableNames', {'NoiseLevel', 'NoiseRate', 'MSE', 'PSNR_dB', 'SSIM'});

disp(' ');
disp('================ Robustness Results Table ================');
disp(ResultTable);
disp('==========================================================');

writetable(ResultTable, fullfile(out_dir, 'encoded_robustness_results.csv'));

%% =========================
%  保存 TXT
%% =========================
txt_file = fullfile(out_dir, 'encoded_robustness_results.txt');
fid = fopen(txt_file, 'w');

fprintf(fid, '================ Robustness Results Table ================\n');
fprintf(fid, 'NoiseLevel\tNoiseRate\tMSE\tPSNR_dB\tSSIM\n');
for i = 1:height(ResultTable)
    fprintf(fid, '%s\t%.6f\t%.6f\t%.6f\t%.6f\n', ...
        ResultTable.NoiseLevel(i), ...
        ResultTable.NoiseRate(i), ...
        ResultTable.MSE(i), ...
        ResultTable.PSNR_dB(i), ...
        ResultTable.SSIM(i));
end
fprintf(fid, '==========================================================\n');

fclose(fid);

%% =========================
%  画曲线
%% =========================
fig2 = figure('Position', [180, 180, 1200, 500]);

subplot(1,2,1);
plot(noise_rates * 100, PSNR_list, '-o', 'LineWidth', 1.8, 'MarkerSize', 8);
xlabel('Perturbation rate (%)');
ylabel('PSNR (dB)');
title('PSNR under encoded-image perturbation');
grid on;

subplot(1,2,2);
plot(noise_rates * 100, SSIM_list, '-o', 'LineWidth', 1.8, 'MarkerSize', 8);
xlabel('Perturbation rate (%)');
ylabel('SSIM');
title('SSIM under encoded-image perturbation');
grid on;

saveas(fig2, fullfile(out_dir, 'encoded_perturbation_psnr_ssim_curves.png'));

fprintf('\n已完成并保存：\n');
fprintf('1. 恢复总览图: %s\n', fullfile(out_dir, 'encoded_perturbation_recovery_overview.png'));
fprintf('2. 曲线图    : %s\n', fullfile(out_dir, 'encoded_perturbation_psnr_ssim_curves.png'));
fprintf('3. 结果 CSV  : %s\n', fullfile(out_dir, 'encoded_robustness_results.csv'));
fprintf('4. 结果 TXT  : %s\n', txt_file);

%% =========================================================
%  辅助函数：对 uint8 图像做 substitution-like 随机替换扰动
%% =========================================================
function noisy_img = apply_uint8_substitution_noise(img, rate)
    noisy_img = img;
    total_num = numel(img);
    num_change = round(total_num * rate);

    if num_change == 0
        return;
    end

    idx = randperm(total_num, num_change);

    % 对每个被选中的位置，替换成一个不同的随机 uint8 值
    for t = 1:length(idx)
        old_val = img(idx(t));
        new_val = uint8(randi([0,255]));
        while new_val == old_val
            new_val = uint8(randi([0,255]));
        end
        noisy_img(idx(t)) = new_val;
    end
end