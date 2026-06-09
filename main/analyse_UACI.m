function UACI = analyse_UACI(original_img, encrypted_img)
    [M, N, C] = size(original_img);
    accumulate = 0;
    
    for ch = 1:C
        for i = 1:M
            for j = 1:N
                accumulate = accumulate + abs(double(original_img(i,j,ch)) - double(encrypted_img(i,j,ch))) / 255.0;
            end
        end
    end
    
    UACI = accumulate / (M * N * C);
end