function descrambled_img = descramble_img(scrambled_img, key_stream_scrambling)
    [M, N, dimension] = size(scrambled_img);
    n = M*N;
    descrambled_img = scrambled_img;
    
    % Reshape key stream and ensure valid indices
    key_reshaped = reshape(key_stream_scrambling, [M, N, 3]);
    
    red_stream = reshape(key_reshaped(:,:,1), [n, 1]);
    green_stream = reshape(key_reshaped(:,:,2), [n, 1]);
    blue_stream = reshape(key_reshaped(:,:,3), [n, 1]);
    
    % Ensure indices are within valid range [1, n]
    red_stream = mod(red_stream - 1, n) + 1;
    green_stream = mod(green_stream - 1, n) + 1;
    blue_stream = mod(blue_stream - 1, n) + 1;
    
    red_img = descrambled_img(:,:,1);
    green_img = descrambled_img(:,:,2);
    blue_img = descrambled_img(:,:,3);
    
    % Convert to linear indexing
    red_img_linear = red_img(:);
    green_img_linear = green_img(:);
    blue_img_linear = blue_img(:);
    
    % Perform inverse scrambling using temporary array
    temp_red = red_img_linear;
    temp_green = green_img_linear;
    temp_blue = blue_img_linear;
    
    for i = n:-1:1
        swap_idx = red_stream(i);
        temp = temp_red(i);
        temp_red(i) = temp_red(swap_idx);
        temp_red(swap_idx) = temp;
        
        swap_idx = green_stream(i);
        temp = temp_green(i);
        temp_green(i) = temp_green(swap_idx);
        temp_green(swap_idx) = temp;
        
        swap_idx = blue_stream(i);
        temp = temp_blue(i);
        temp_blue(i) = temp_blue(swap_idx);
        temp_blue(swap_idx) = temp;
    end
    
    % Reshape back to 2D
    descrambled_img(:,:,1) = reshape(temp_red, [M, N]);
    descrambled_img(:,:,2) = reshape(temp_green, [M, N]);
    descrambled_img(:,:,3) = reshape(temp_blue, [M, N]);
end