clear all; close all;clc;
addpath('../src', '../src/utils');
% % % % % % % % % % % % % 
%
% This script has four steps. 
% I. It starts by generating the observed hyperspectral and
% multispectral/panchromatic images. The following parameters can be
% modified to change the data generation:
rng(1);
downsamp_factor = 4;
down_fact=4;% Downsampling factor
SNRh = 20; % SNR (in dB) for the hyperspectral image
SNRm = 20; % SNR (in dB) for the multispectral/panchromatic image
% 
% II. Next, it estimates the spectral and spatial response of the sensors.
% The regularization parameters can be adjusted here:
lambda_R = 1e1;
lambda_B = 1e1;
% For the denoising with SVD, we need to specify the number of bands we
% want to keep
p = 10; % Corresponds to variable L_s in [1]; number of endmembers in VCA /
% number of non-truncated singular vectors
%
% III. The data fusion algorithm is then called using the estimated responses
% and the observed data. The following parameters can be
% modified:
% 
basis_type = 'SVD';
lambda_phi = 5e-4;
lambda_m = 1e0;
% 
% IV. In the end, three quality indices are computed using the ground
% truth image. These indices are ERGAS, SAM, and UIQI.
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                       %
% I. Observed data (simulation)                                         %
% -----------------------------                                         %
%                                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
iptsetpref('ImshowBorder','tight');
ms_bands = 2:5; % 5 bands: 1 - pan; 2 - blue; 3 - green; 4 - red; 5 - NIR

% Load ROSIS image (without the first ten bands and without the last
% collumn)
load  '/MATLAB Drive/ISTA_Fusion/ISTA FUSION/Fusion/HySure-master/data/rosis_small_Datacube'
[nl, nc, ~] = size(X);
% In matrix form
Z = im2mat(X);
Z = Z/max(Z(:));
Z=Z(11:end,:);
L= size(Z,1);
clear X;
        
% % % % % % % % % % % % % 
% Blur kernel
middlel = round((nl+1)/2);
middlec = round((nc+1)/2);
% Blur matrix
B = zeros(nl,nc);
% Starck-Murtagh filter
B(middlel-2:middlel+2, middlec-2:middlec+2) = [1 4 6 4 1; 4 16 24 16 4; 6 24 36 24 6; 4 16 24 16 4; 1 4 6 4 1];
% Circularly center B
B = ifftshift(B);
% Normalize
B = B/sum(sum(B));
% Fourier transform of the filters
FB = fft2(B);

% % % % % % % % % % % % % 
% Simulate the HS data
% Spatial degradation (blur)
Yh = ConvC(Z, FB, nl);
Yhim_up = mat2im(Yh, nl);
% Add noise
sigmah = sqrt(sum(Yhim_up(:).^2)/(10^(SNRh/10))/numel(Yhim_up));
Yhim_up = Yhim_up + sigmah*randn(size(Yhim_up));
% % % % % % % % % % % Downsampling (version with reduced size)
% Downsampling
Yhim = downsamp_HS(Yhim_up, downsamp_factor, 1);
    
% % % % % % % % % % % % % 
% Simulate the MS/PAN data
% Use IKONOS's spectral response 
% (wavelenghths, pan, blue, green, red, NIR, in nanometers)
load '/MATLAB Drive/ISTA_Fusion/ISTA FUSION/Fusion/HySure-master/data/ikonos_spec_resp.mat'
% Map IKONOS' wavelengths into ROSIS (430 - 860 nm) bands
% Find valid interval ikonos \subset rosis
[~, valid_ik_bands] = intersect(ikonos_sp(:,1), 430:860);
no_wa = length(valid_ik_bands);
% Spline interpolation
xx  = linspace(1, no_wa, L);
x = 1:no_wa;
R = zeros(5, L);
for i = 1:5 % 1 - pan; 2 - blue; 3 - green; 4 - red; 5 - NIR
    R(i,:) = spline(x, ikonos_sp(valid_ik_bands,i+1), xx);
end
% Use just the predefined bands
R = R(ms_bands,:);
% Spectral degradation
Ym = R*Z;
% Normalize all channels to 1 (NOTE: this changes the SNRm)
c = zeros(length(ms_bands));
for i=1:length(ms_bands);
    c(i) = max(Ym(i,:));
    Ym(i,:) = Ym(i,:)/c(i);
    R(i,:) =  R(i,:)/c(i);
end
% Add noise
sigmam = sqrt(sum(Ym(:).^2)/(10^(SNRm/10))/numel(Ym));
Ym = Ym + sigmam*randn(size(Ym));
Ymim = mat2im(Ym, nl);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                       %
% II. Spectral and spatial responses estimation                         %
% ---------------------------------------------                         %
%                                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


intersection = cell(1,length(ms_bands));
intersection{1} = 1:30;
intersection{2} = 10:44;
intersection{3} = 36:66;
intersection{4} = 60:93;
contiguous = intersection;

% Blur's support: [hsize_h hsize_w]
hsize_h = 10;
hsize_w = 10;
shift = 1; % the 'phase' parameter in MATLAB's 'upsample' function
blur_center = 0; % to center the blur kernel according to the simluated data
[V, R_est, B_est] = sen_resp_est(Yhim, Ymim, downsamp_factor, intersection, contiguous, p, lambda_R, lambda_B, hsize_h, hsize_w, shift, blur_center);

% Denoises the original image, since it is quite noisy as well
Z = (V*V')*Z;
% In image form
Zim = mat2im(Z, nl);

E = V;

lambda_phi = 1;
lambda_m = 6;
h = 20/255 ;
maxiters = 1500;
delta = 4;
scale_fact=1;
mu=.50;
sigma= 0.3;
PatchSizeHalf=3;
 WindowSizeHalf=5;
 Sigma=10/255;

SB = @(X) blur_downsample(X,B_est,down_fact);
SBadj = @(X) upsample_blur(X,B_est,down_fact);
fwdoper= @(X) compute_forward(X,SB,SBadj,E,R_est,lambda_m,lambda_phi);
%C=@(X) compute_C(X,E,Yh,Ym,lambda_m,lambda_phi,mu,sigma);
grad = @(X) compute_gradient(X,SB,SBadj,E,R_est,Yhim,Ymim,lambda_m,lambda_phi);

[rr,cc,~] = size(Ymim);

psnr_vals = nan(1,maxiters);
iters = 1;

bicubic_fusion = imresize(Yhim,down_fact);

tol = 1e-15;
max_iter=1000;
%W = compute_weights(Ymim,search_rad,patch_rad,h,1);
%load('W.mat');
%denoiser = @(x) perform_denoising(W,x);
% bb=[8 13 45];
% temp_fig = Zim;
% temp_rgb = imadjust(temp_fig(:,:,bb),stretchlim(temp_fig(:,:,bb)),[]);
% figure;
% imshow(temp_rgb)

X_curr = var_dim(ima_interp_spline(Yhim,downsamp_factor),pinv(E));
X_curr1 = var_dim(ima_interp_spline(Yhim,downsamp_factor),pinv(E));
x_Entry= cell(1,maxiters);
for itera=1:1:3
    V_curr = X_curr - delta * grad(X_curr);
denoiser = @(x) wrapper_FASTDSGNLM(x,h,X_curr1);
    V2= scale_fact * denoiser(V_curr);
    X_curr=V2;
 V_curr = X_curr ;
     B2=V_curr;
    for jj=1:size(B2,3)
    eigen_im=(  B2(:,:,jj));
     input=eigen_im;
    
     eigen_im_1=(  X_curr1(:,:,jj));
   input_1=eigen_im_1;

   res=Fast_DSG_NLM(input, input_1, PatchSizeHalf, WindowSizeHalf, Sigma);
   
   BB = gather(res);  
    
     V2(:,:,jj) = double(BB); 
    end
     X_next=V2;
end
X_curr = var_dim(ima_interp_spline(Yhim,downsamp_factor),pinv(E));
X_curr2 = zeros(rr,cc,p);
%X_curr=ones(rr,cc,p);
X_curr1=X_next;
X_curr2=X_next;
X_next=zeros(rr,cc,p);
 %C= compute_C(X_curr2,SBadj,E,R_est,Yhim, Ymim,lambda_m,lambda_phi,mu,sigma);
while(true)
    
   C= compute_C(X_next,SBadj,E,R_est,Yhim, Ymim,lambda_m,lambda_phi,mu,sigma);
    [X_k,~]= solve_AX_plus_XB_with_K(SB,SBadj,E,R_est,lambda_m,lambda_phi,mu,sigma, C, tol, max_iter,rr,cc,p);
   
   X_curr = reshape(X_k,rr,cc,p);
   X_next=X_curr;
   V_curr = X_curr;
   denoiser = @(x) wrapper_FASTDSGNLM(x,h,X_curr1);
    V2= scale_fact * denoiser(V_curr);
    X_curr=V2;
 V_curr = X_curr ;
     B2=V_curr;
    for jj=1:size(B2,3)
    eigen_im=(  B2(:,:,jj));
     input=eigen_im;
    
     eigen_im_1=(  X_curr1(:,:,jj));
   input_1=eigen_im_1;

   res=Fast_DSG_NLM(input, input_1, PatchSizeHalf, WindowSizeHalf, Sigma);
   
   BB = gather(res);  
    
     V2(:,:,jj) = double(BB); 
    end
     X_next=V2;
    Z_next = E*im2mat(X_next);
    Z_next = mat2im(Z_next,rr);
    [rmse_total_next, ergas_next, sam_next, uiqi_next] = quality_assessment(Zim, Z_next, 0, 1/down_fact);
    psnr_next = psnr(Z_next,Zim,1);
    a=norm(im2mat(X_next-X_curr2),'fro');
    %fprintf('Iteration = %d, PSNR = %f,iterationvalue=%f\n',iters,psnr_next,a);
    fprintf('Iteration = %d,PSNR=%f, RMSE = %f,ERGAS=%f,SAM=%f,UIQI=%f,iterationvalue=%f\n',iters,psnr_next,rmse_total_next,ergas_next,sam_next,uiqi_next,a);
    psnr_vals(iters) = psnr_next;
    x_Entry{iters}=X_curr2;
    
    if(iters==maxiters)
        break;
    end
    
    iters = iters+1;
    X_curr2 = X_next;
   % [rmse_total, ergas, sam, uiqi] = quality_assessment(Zim, Z_next, 0, 1/down_fact)

end

Z_hat = Z_next;
[rmse_total, ergas, sam, uiqi] = quality_assessment(Zim, Z_hat, 0, 1/down_fact)

bb=[3 2 1];
%bb=[20 80 100];
%bb=[15 50 100];

Z_hat_rgb = imadjust(Z_hat(:,:,bb),stretchlim(Z_hat(:,:,bb)),[]);
Zim_rgb = imadjust(Zim(:,:,bb),stretchlim(Zim(:,:,bb)),[]);
bicubic_rgb = imadjust(bicubic_fusion(:,:,bb),stretchlim(bicubic_fusion(:,:,bb)),[]);

figure;
subplot(1,3,1); imshow(Zim_rgb); title('Ground-truth');
subplot(1,3,2); imshow(Z_hat_rgb); title('Recovered');
subplot(1,3,3); imshow(bicubic_rgb); title('Bicubic');
iteration_Cou=1:1:maxiters;
% bb=[8 13 45];
% temp_fig = Zim;
% temp_rgb = imadjust(temp_fig(:,:,bb),stretchlim(temp_fig(:,:,bb)),[]);
% figure;
% imshow(temp_rgb)
X_conv=nan(1,maxiters);
for i=1:1:maxiters
    X_conv(i)=log(norm(im2mat(x_Entry{i}-X_curr2),'fro'));
end

figure;
plot(iteration_Cou,psnr_vals);xlabel('iteration');ylabel('PSNR')
figure;
plot(iteration_Cou,X_conv);xlabel('iteration');ylabel('$\log\|X_k-X_*\|$','interpreter','latex')

%%
function terma=compute_forward(X,SB,SBadj,E,R,lambda_m,lambda_phi,mu,sigma)
rows = size(X,1);

term1 = SB(X);
term1 = SBadj(term1);
term1 = lambda_phi*E'*E*im2mat(term1);
term2 = lambda_m * (R*E)' * (R*E) * im2mat(X);
term3=2*mu*sigma*sigma*im2mat(X);
 terma=(term1+term2+term3)';
 % terma=mat2im(terma,rows);
end

function C= compute_C(Y,SBadj,E,R,Yh,Ym,lambda_m,lambda_phi,mu,sigma)
%rows = size(Y,1);
term4 = SBadj(Yh);
term5 = lambda_phi*(E')*im2mat(term4);

term6 = lambda_m * (R*E)' * im2mat(Ym);
term7=2*mu*sigma*sigma*im2mat(Y);

C=(term5+term6+term7)';
%C=mat2im(C,rows);
end
%%
%% Function to solve A*X + X*B = C using Conjugate Gradient
function [X, iter] = solve_AX_plus_XB_with_K(SB,SBadj,E,R,lambda_m,lambda_phi,mu,sigma, C, tol, max_iter,n,m,p)
    % Solve A X + X B = C using Conjugate Gradient
    
    % Get the dimensions of C
    %[n, m] = size(C);  % n is the number of rows, m is the number of columns
    
    % Vectorize the right-hand side matrix (C)
    b = C(:);  % Convert C to a vector
    
   
    
    % Perform Conjugate Gradient
    [vec_X, iter] = conjugate_gradient_operator(SB,SBadj,E,R,lambda_m,lambda_phi,mu,sigma, b, n, m,p, tol, max_iter);
    
    % Reshape the solution vector back to matrix form
    X = reshape(vec_X, n*m, p);  % Convert vector back to matrix form
end

%% Conjugate Gradient Method with Operator-based Approach
function [x, iter] = conjugate_gradient_operator(SB,SBadj,E,R,lambda_m,lambda_phi,mu,sigma, b, n, m,pp, tol, max_iter)
    % Conjugate Gradient method to solve Ax = b, where the operator is applied directly
    
    x = zeros(n * m* pp , 1);  % Initial guess (vectorized form)
    X=zeros(n,m,pp);
    com=compute_forward(X,SB,SBadj,E,R,lambda_m,lambda_phi,mu,sigma);
    com=reshape(com,n*m,pp);
    com=com(:);
    r = b - com;  % Initial residual
    p = r;  % Initial search direction
    rs_old = r' * r;  % Initial residual squared
    
    for iter = 1:max_iter
        P=reshape(p,n,m,pp);
        com=compute_forward(P,SB,SBadj,E,R,lambda_m,lambda_phi,mu,sigma);
        com=reshape(com,n*m,pp);
        com=com(:);
        Ap = com;  % Apply the operator to p
        
        % Step size (alpha)
        alpha = rs_old / (p' * Ap);
        
        % Update the solution
        x = x + alpha * p;
        
        % Update the residual
        r = r - alpha * Ap;
        
        % Compute new residual squared
        rs_new = r' * r;
        
        % Check for convergence
        if sqrt(rs_new) < tol
            fprintf('Converged after %d iterations\n', iter);
            break;
        end
        
        % Update the search direction
        p = r + (rs_new / rs_old) * p;
        
        % Update old residual
        rs_old = rs_new;
    end
    
    % If the method doesn't converge within max_iter, print a message
    if iter == max_iter
        fprintf('Reached maximum number of iterations (%d)\n', max_iter);
    end
end
