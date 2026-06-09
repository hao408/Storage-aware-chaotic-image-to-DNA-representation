function analyse_gc_content(encoded_img, complementary_principle)
%{
    Analyze GC content of DNA-encoded image
%}
    [M, N, C] = size(encoded_img);
    
    figure;
    sgtitle('GC Content Analysis');
    
    for ch = 1:3
        channel_img = encoded_img(:,:,ch);
        gc_content_map = zeros(M, N);
        
        for i = 1:M
            for j = 1:N
                pixel = channel_img(i, j);
                gc_count = 0;
                total_bases = 0;
                
                for k = 1:4
                    two_bits = bitget(pixel, [2*k-1, 2*k]);
                    num = two_bits(1) + 2 * two_bits(2);
                    base_value = complementary_principle(1, num + 1);
                    
                    % Count G and C bases (ASCII: G=71, C=67)
                    if base_value == 71 || base_value == 67
                        gc_count = gc_count + 1;
                    end
                    total_bases = total_bases + 1;
                end
                
                gc_content_map(i, j) = gc_count / total_bases;
            end
        end
        
        subplot(2,2,ch);
        histogram(gc_content_map(:), 20);
        title(sprintf('Channel %d GC Content Distribution', ch));
        xlabel('GC Ratio');
        ylabel('Frequency');
        grid on;
        
        mean_gc = mean(gc_content_map(:));
        std_gc = std(gc_content_map(:));
        fprintf('Channel %d - Mean GC: %.4f, Std: %.4f\n', ch, mean_gc, std_gc);
    end
    
    % Overall GC content
    overall_gc = [];
    for ch = 1:3
        channel_img = encoded_img(:,:,ch);
        for i = 1:M
            for j = 1:N
                pixel = channel_img(i, j);
                for k = 1:4
                    two_bits = bitget(pixel, [2*k-1, 2*k]);
                    num = two_bits(1) + 2 * two_bits(2);
                    base_value = complementary_principle(1, num + 1);
                    if base_value == 71 || base_value == 67
                        overall_gc = [overall_gc, 1];
                    else
                        overall_gc = [overall_gc, 0];
                    end
                end
            end
        end
    end
    
    overall_gc_ratio = mean(overall_gc);
    fprintf('Overall GC Content: %.4f (%.2f%%)\n', overall_gc_ratio, overall_gc_ratio * 100);
end