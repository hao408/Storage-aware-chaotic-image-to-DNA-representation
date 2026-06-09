function analyse_histogram(img)
% analyse_histogram — draw histogram of encrypted image without Image Processing Toolbox
%
% Usage:
%     analyse_histogram(img)
%
% Input:
%     img : encrypted or original RGB image matrix (uint8 or double)
%
% Description:
%     This version does not require imhist (Image Processing Toolbox).
%     Uses MATLAB built-in histogram() which works in all standard editions.

    % Ensure input exists
    if nargin < 1
        error('analyse_histogram: input image missing!');
    end

    % Convert to uint8 if needed (range 0–255)
    if ~isa(img,'uint8')
        img = uint8(img);
    end

    % If grayscale, convert to pseudo RGB for analysis
    if size(img,3) == 1
        img = cat(3, img, img, img);
    end

    % Separate three channels
    R = img(:,:,1);
    G = img(:,:,2);
    B = img(:,:,3);

    % Plot R histogram
    figure;
    histogram(R(:), 0:255);
    title('Red Channel Histogram');
    xlabel('Pixel Value'); ylabel('Frequency');
    grid on;

    % Plot G histogram
    figure;
    histogram(G(:), 0:255);
    title('Green Channel Histogram');
    xlabel('Pixel Value'); ylabel('Frequency');
    grid on;

    % Plot B histogram
    figure;
    histogram(B(:), 0:255);
    title('Blue Channel Histogram');
    xlabel('Pixel Value'); ylabel('Frequency');
    grid on;

    % Optional combined histogram
    figure;
    histogram(img(:), 0:255);
    title('Overall Pixel Histogram');
    xlabel('Pixel Value'); ylabel('Frequency');
    grid on;
end
