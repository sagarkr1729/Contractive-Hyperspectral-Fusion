clc;
clear;
close all;

% Parameters
blurSize = 25;            % Size of the Gaussian blur kernel
sigmaBlur = 1.6;          % Standard deviation for the Gaussian blur
maxIter = 400;             % Number of ISTA iterations
lambda = 0; 
l1=1;
lm=1;% Step size for ISTA
iptsetpref('ImshowBorder','tight');  % Adjust imshow preferences

% Initialize DnCNN denoiser network
net = denoisingNetwork('dncnn');  % Load pre-trained DnCNN

% Read the input image
data = imread('Set12/01.png');  % Use your image path here
inputImage = double(data) / 255.0;  % Normalize image to [0, 1]

% Convert to grayscale if the image is RGB
if size(inputImage, 3) == 3
    inputImage = rgb2gray(inputImage);
end

% Apply Gaussian blur to the image
h = fspecial('gaussian', blurSize, sigmaBlur);
blurredImage = imfilter(inputImage, h, 'circular');  % Apply the Gaussian blur

% Add Gaussian noise to the blurred image
blurredImage = imnoise(blurredImage, 'gaussian', 0, 0.000006);  % Adjust noise variance as needed
Noisyim=imnoise(inputImage, 'gaussian', 0, 0.03);
% Initialize the fused image (start with noisy image or blurred image)
%fusedImage = Noisyim;  % Or start with blurredImage as the initial guess
fusedImage = blurredImage;
% Initialize ISTA variables
x = fusedImage;  % Initial estimate (fused/noisy image)
pre = fusedImage;  % Store previous image for convergence check
sigma=1/255;
% Initialize variables for the best result and PSNR calculation
highest_psnr = -Inf;  % Start with a low PSNR
best_x = x;           % Store the best denoised/deblurred image
psnr_iter = zeros(maxIter, 1);  % To store PSNR at each iteration

for iter = 1:maxIter
    pre = x;  % Store the previous image
    
    % Apply blur to the current estimate (x) to simulate forward model
    x_blurred = imfilter(x, h, 'circular');
    
    % Calculate residual (Ax - b), where x is the current image and b is the noisy image
    residual = x_blurred - blurredImage;  % Residual between blurred and noisy image
    
    residual = l1*imfilter(residual,h,'circular') + lm * (x - Noisyim);
    % Proximal operator (shrinkage thresholding)
    % Apply DnCNN to the residual (remove noise in the residual)
    x = BM3D((pre - lambda * residual), sigma);
    
    % Ensure the image is within the valid range [0, 1]
   % x = max(0, min(1, x));
     % Compute PSNR and SSIM for the current iteration
    psnr_iter(iter) = psnr(double(x), double(inputImage));  % Ensure both are of the same class
    ssim_iter = ssim(x, inputImage);  % Compute SSIM for the current iteration
    
    % Print PSNR and SSIM for the current iteration
    fprintf('Iteration %d: PSNR = %.2f dB, SSIM = %.4f\n', iter, psnr_iter(iter), ssim_iter);
    % Compute PSNR for the current iteration
    psnr_iter(iter) = psnr(double(x), double(inputImage));  % Ensure both are of the same class
    
    % Update the best result if the current PSNR is the highest
    if psnr_iter(iter) > highest_psnr
        highest_psnr = psnr_iter(iter);
        best_x = x;  % Save the best image
    end
end

% Display results
figure;
subplot(1, 4, 1); imshow(inputImage); title('Original Image');
subplot(1, 4, 2); imshow(blurredImage); title('Blurred ');
subplot(1, 4, 3); imshow(Noisyim); title(' Noisy Image');
subplot(1, 4, 4); imshow(x); title('fused');

% Calculate PSNR and SSIM for the noisy image
psnr_noisy = psnr(double(Noisyim), double(inputImage));
ssim_noisy = ssim(Noisyim, inputImage);

% Calculate PSNR and SSIM for the reconstructed image with highest PSNR
psnr_recon = psnr(double(x), double(inputImage));
ssim_recon = ssim(x, inputImage);

% Output PSNR and SSIM results
fprintf('PSNR of noisy image: %.2f dB\n', psnr_noisy);
fprintf('SSIM of noisy image: %.4f\n', ssim_noisy);
fprintf('PSNR of deblurred image (Highest PSNR): %.2f dB\n', psnr_recon);
fprintf('SSIM of deblurred image (Highest PSNR): %.4f\n', ssim_recon);

% Plot PSNR over iterations
figure;
plot(1:maxIter, psnr_iter);
xlabel('Iteration');
ylabel('PSNR value');
title('PSNR Values over ISTA Iterations');

