function results = run_experiments_table1(original_image, dna_complementary_principle, n_runs)
% RUN_EXPERIMENTS_TABLE1  —— 生成 Table 1 所需的实验数据（v2.1，含 6D）
%
% 特点：
%   1) 同一加密框架，切换不同混沌系统进行对比；
%   2) 相关系数采用“全图像邻域像素”计算，结果确定性；
%   3) 加密时间采用 n_runs 次实验的平均值和标准差；
%   4) 打印时熵、NPCR、UACI、相关系数、时间均保留 6 位小数。

    if nargin < 2
        error('需要 original_image 和 dna_complementary_principle 两个输入。');
    end
    if nargin < 3 || isempty(n_runs)
        n_runs = 10;
    end

    img = uint8(original_image);

    % ====== 要比较的混沌系统（你可以按需要增删顺序） ======
    schemes = {
        struct('tag','CIST',  'name','CIS-T IS+Tent',      'fun',@generate_chaotic_quence_ten), ...
        struct('tag','IS3',   'name','3-track IS-map',     'fun',@generate_chaotic_quence), ...
        struct('tag','LOG',   'name','Logistic-3track',    'fun',@generate_chaotic_quence_logistic), ...
        struct('tag','HENON', 'name','Henon-1system',      'fun',@generate_chaotic_quence_henon), ...
        struct('tag','SYS6D', 'name','SixD-chaotic',       'fun',@generate_chaotic_quence_6d)
    };

    num_schemes = numel(schemes);
    results = struct([]);

    % ====== 构造“单像素扰动”的原图，用来计算 NPCR/UACI ======
    img2 = img;
    img2(1,1,1) = uint8(mod(double(img2(1,1,1)) + 1, 256));

    fprintf('==== Table 1: 同一框架，不同混沌系统实验（v2.1）====\n');
    fprintf('Tag\tChaos\t\tH_mean\t\tNPCR\t\tUACI\t\tCorrH\t\tCorrV\t\tCorrD\t\tTime_mean(s)\tTime_std(s)\n');

    for k = 1:num_schemes
        fgen = schemes{k}.fun;
        tag  = schemes{k}.tag;
        name = schemes{k}.name;

        % ---------- 1) 指标：熵 / NPCR / UACI / 相关系数 ----------
        [Kdif1, Kscr1, ~] = fgen(img);
        cipher1 = diffuse_img(scramble_img(img, Kscr1), Kdif1);

        [Kdif2, Kscr2, ~] = fgen(img2);
        cipher2 = diffuse_img(scramble_img(img2, Kscr2), Kdif2);

        ent_arr = analyse_entropy(cipher1);
        H_mean  = mean(ent_arr);

        uaci_npcr = analyse_uaci_nprc(cipher1, cipher2, false, 255);
        NPCR = uaci_npcr.npcr_score;
        UACI = uaci_npcr.uaci_score;

        ch1 = cipher1(:,:,1);
        [rH, rV, rD] = compute_neighbor_correlation_full(ch1);

        % ---------- 2) 加密时间：重复 n_runs 次 ----------
        times = zeros(1, n_runs);
        for t = 1:n_runs
            tic;
            [Kdif_t, Kscr_t, ~] = fgen(img);
            cipher_t = diffuse_img(scramble_img(img, Kscr_t), Kdif_t); %#ok<NASGU>
            times(t) = toc;
        end
        time_mean = mean(times);
        time_std  = std(times);

        % ---------- 3) 保存 ----------
        results(k).tag       = tag;
        results(k).name      = name;
        results(k).H_mean    = H_mean;
        results(k).NPCR      = NPCR;
        results(k).UACI      = UACI;
        results(k).corrH     = rH;
        results(k).corrV     = rV;
        results(k).corrD     = rD;
        results(k).time_mean = time_mean;
        results(k).time_std  = time_std;
        results(k).time_all  = times;
        results(k).n_runs    = n_runs;

        % ---------- 4) 打印 ----------
        fprintf('%s\t%-10s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n', ...
            tag, name, H_mean, NPCR, UACI, rH, rV, rD, time_mean, time_std);
    end

    fprintf('=======================================================================\n');
end

function [rH, rV, rD] = compute_neighbor_correlation_full(ch)
    ch = double(ch);

    X = ch(:, 1:end-1);
    Y = ch(:, 2:end);
    r = corrcoef(X(:), Y(:));
    rH = r(1,2);

    X = ch(1:end-1, :);
    Y = ch(2:end,   :);
    r = corrcoef(X(:), Y(:));
    rV = r(1,2);

    X = ch(1:end-1, 1:end-1);
    Y = ch(2:end,   2:end);
    r = corrcoef(X(:), Y(:));
    rD = r(1,2);
end
