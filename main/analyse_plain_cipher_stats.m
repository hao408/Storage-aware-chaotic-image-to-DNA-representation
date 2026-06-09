function stats = analyse_plain_cipher_stats(original_image, cipher_image_like)
%ANALYSE_PLAIN_CIPHER_STATS  比较原图与密文的熵 / 相关性（支持灰度或彩色）
%
% 用法：
%   stats = analyse_plain_cipher_stats(original_image, cipher_image_like);
%
% 输入：
%   original_image    : 原始图像 (HxWx1 或 HxWx3, uint8)
%   cipher_image_like : 像素域密文 或 由 encoded_image 映射得到的伪灰度密文
%
% 输出 stats 结构体示例：
%   stats.plain(i).entropy, stats.plain(i).rH, rV, rD
%   stats.cipher(i).entropy, stats.cipher(i).rH, rV, rD

    P = uint8(original_image);
    C = uint8(cipher_image_like);

    compute_entropy = @(X) ...
        - sum( ...
            (histcounts(X(:), 0:256, 'Normalization','probability') + eps) ...
            .* log2(histcounts(X(:), 0:256, 'Normalization','probability') + eps) ...
          );

    stats = struct('plain',[],'cipher',[]);
    chn = {'R','G','B'};
    if size(P,3)==1
        chn = {'Gray'};
    end

    % ---- 原图 ----
    for k = 1:size(P,3)
        X = P(:,:,k);
        H_plain = compute_entropy(X);

        X1 = double(X(:,1:end-1)); X2 = double(X(:,2:end));
        rH = corrcoef(X1(:), X2(:)); rH = rH(1,2);

        X1 = double(X(1:end-1,:)); X2 = double(X(2:end,:));
        rV = corrcoef(X1(:), X2(:)); rV = rV(1,2);

        X1 = double(X(1:end-1,1:end-1)); X2 = double(X(2:end,2:end));
        rD = corrcoef(X1(:), X2(:)); rD = rD(1,2);

        stats.plain(k).channel = chn{k};
        stats.plain(k).entropy = H_plain;
        stats.plain(k).rH = rH;
        stats.plain(k).rV = rV;
        stats.plain(k).rD = rD;

        fprintf('原图-%s: 熵=%.4f, 相关性[H V D]=[%.5f %.5f %.5f]\n', ...
            chn{k}, H_plain, rH, rV, rD);
    end

    % ---- 密文 ----
    for k = 1:size(C,3)
        X = C(:,:,k);
        H_cipher = compute_entropy(X);

        X1 = double(X(:,1:end-1)); X2 = double(X(:,2:end));
        rH = corrcoef(X1(:), X2(:)); rH = rH(1,2);

        X1 = double(X(1:end-1,:)); X2 = double(X(2:end,:));
        rV = corrcoef(X1(:), X2(:)); rV = rV(1,2);

        X1 = double(X(1:end-1,1:end-1)); X2 = double(X(2:end,2:end));
        rD = corrcoef(X1(:), X2(:)); rD = rD(1,2);

        stats.cipher(k).channel = chn{min(k,numel(chn))};
        stats.cipher(k).entropy = H_cipher;
        stats.cipher(k).rH = rH;
        stats.cipher(k).rV = rV;
        stats.cipher(k).rD = rD;

        fprintf('密文-%s: 熵=%.4f, 相关性[H V D]=[%.5f %.5f %.5f]\n', ...
            stats.cipher(k).channel, H_cipher, rH, rV, rD);
    end
end
