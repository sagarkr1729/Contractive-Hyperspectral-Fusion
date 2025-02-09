function [Y] = upsample_blur(X,B,up_fact)
%UPSAMPLE_BLUR Apply upsampling and blurring operations (in order) to a
% multi-channel image
% X = m-by-n-by-c input image, where c is the no. of channels
% B = Time-reversed blurring kernel having the same size as the upsampled
% image (non-vectorized)
% up_fact = Downsampling factor

shift = 1;
Y = upsamp_HS(X, up_fact, size(X,1)*up_fact, size(X,2)*up_fact, shift);
nrows = size(Y,1);
Ymat = im2mat(Y);
FB = fft2(B);       % Fourier transform of the input filter
FB = conj(FB);      % Fourier transform of time-reversed filter
Y = ConvC(Ymat, FB, nrows);
Y = mat2im(Y, nrows);   % Convert to image

end

