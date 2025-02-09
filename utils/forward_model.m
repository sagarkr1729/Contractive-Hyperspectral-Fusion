function [Yhim,Ymim,Ypim,R,ms_bands,Z,B,Yhim_up] = forward_model(X,down_fact,SNRh,SNRm,SNRp,ikonos_sp)
%FORWARD_MODEL Generate HS, MS and PAN images from ground-truth data
% This function uses the IKONOS response to simulate the spectral
% degradation operator.
% Inputs:
% X = Ground-truth image with high spatial & spectral resolutions (unnormalized)
% down_fact = Spatial downsampling factor (to generate HS image)
% SNRh = Desired SNR of HS image
% SNRm = Desired SNR of MS image
% SNRp = Desired SNR of PAN image
% ikonos_data = Matrix containing IKONOS response
% Outputs:
% Yhim = HS image
% Ymim = MS image
% Ypim = PAN image
% R = Spectral degradation matrix
% 

[nl, nc, ~] = size(X);
Z = im2mat(X);      % Convert to matrix
Z = Z/max(Z(:));    % Normalize
Z = Z(11:end,:);
L = size(Z,1);

% Create HS image (low spatial & high spectral resolution)
middlel = round((nl+1)/2);
middlec = round((nc+1)/2);
B = zeros(nl,nc);
B(middlel-2:middlel+2, middlec-2:middlec+2) = ...
    [1 4 6 4 1;...
     4 16 24 16 4;...
     6 24 36 24 6;...
     4 16 24 16 4;...
     1 4 6 4 1];    % Starck-Murtagh filter
B = ifftshift(B);   % Circularly center B
B = B/sum(sum(B));  % Normalize
FB = fft2(B);       % Fourier transform of the filter
Yh = ConvC(Z, FB, nl);
Yhim_up = mat2im(Yh, nl);   % Convert to image
sigmah = sqrt(sum(Yhim_up(:).^2)/(10^(SNRh/10))/numel(Yhim_up));
Yhim_up = Yhim_up + sigmah*randn(size(Yhim_up));
Yhim = downsamp_HS(Yhim_up, down_fact, 1);

% Create MS & PAN images (high spatial & low spectral resolution)
ms_bands = 2:5;
pan_bands = 1;      % These parameters are hard-coded based on the IKONOS matrix structure
valid_wavelengths = 430:860;
[~, valid_ik_bands] = intersect(ikonos_sp(:,1), valid_wavelengths);   % Map IKONOS to ROSIS wavelengths (430-860nm)
no_wa = length(valid_ik_bands);         % No. of valid wavelengths
xx  = linspace(1, no_wa, L);
x = 1:no_wa;
R = zeros(5, L);
for i = 1:5 % 1 - pan; 2 - blue; 3 - green; 4 - red; 5 - NIR
    R(i,:) = spline(x, ikonos_sp(valid_ik_bands,i+1), xx);  % Spline interpolation
end

Ymp = R*Z;               % Spectral degradation
c = zeros(5);
for i=1:5    % Normalize all channels to 1 (NOTE: this changes the SNRm)
    c(i) = max(Ymp(i,:));
    Ymp(i,:) = Ymp(i,:)/c(i);
    R(i,:) =  R(i,:)/c(i);
end

% MS image
Ym = Ymp(ms_bands,:);
R = R(ms_bands,:);
sigmam = sqrt(sum(Ym(:).^2)/(10^(SNRm/10))/numel(Ym));
Ym = Ym + sigmam*randn(size(Ym));   % Add noise
Ymim = mat2im(Ym, nl);

% PAN image
if(~isempty(SNRp))
    Yp = Ymp(pan_bands,:);
    sigmap = sqrt(sum(Yp(:).^2)/(10^(SNRp/10))/numel(Yp));
    Yp = Yp + sigmap*randn(size(Yp));
    Ypim = mat2im(Yp, nl);
else
    Ypim = [];
end

end

