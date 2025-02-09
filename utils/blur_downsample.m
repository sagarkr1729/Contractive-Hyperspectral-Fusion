function [Y] = blur_downsample(X,B,down_fact)
%BLUR_DOWNSAMPLE Apply blurring and downsampling operations (in order) to a
% multi-channel image
% Xmat = m-by-n-by-c input image, where c is the no. of channels
% B = Blurring kernel having the same size as the image (non-vectorized)
% down_fact = Downsampling factor

nrows = size(X,1);
Xmat = im2mat(X);
FB = fft2(B);       % Fourier transform of the filter
Y = ConvC(Xmat, FB, nrows);
Y = mat2im(Y, nrows);   % Convert to image
shift = 1;
Y = downsamp_HS(Y, down_fact, shift);

end

