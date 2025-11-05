clear all; close all;clc;
addpath('../Contractive-Hyperspectral-Fusion-main/utils', '../Contractive-Hyperspectral-Fusion-main/data');
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
load  'Contractive-Hyperspectral-Fusion-main/data/rosis_small_Datacube'
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
load 'Contractive-Hyperspectral-Fusion-main/data/ikonos_spec_resp.mat'
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
Zim = mat2im(Z,size(Ymim,1));

SB = @(X) blur_downsample(X,B_est,down_fact);
SBadj = @(X) upsample_blur(X,B_est,down_fact);
grad = @(X) compute_gradient(X,SB,SBadj,E,R_est,Yhim,Ymim,lambda_m,lambda_phi);
[rr,cc,~] = size(Ymim);
Zim = mat2im(Z,size(Ymim,1));

SB = @(X) blur_downsample(X,B_est,down_fact);
SBadj = @(X) upsample_blur(X,B_est,down_fact);
grad = @(X) compute_gradient(X,SB,SBadj,E,R_est,Yhim,Ymim,lambda_m,lambda_phi);
[rr,cc,~] = size(Ymim);
Zim = mat2im(Z,size(Ymim,1));

SB = @(X) blur_downsample(X,B_est,down_fact);
SBadj = @(X) upsample_blur(X,B_est,down_fact);
grad = @(X) compute_gradient(X,SB,SBadj,E,R_est,Yhim,Ymim,lambda_m,lambda_phi);
[rr,cc,~] = size(Ymim);


lambda_phi = 1;
lambda_m = 6;
h =25/255;
maxiters = 1000;
delta =4;
scale_fact=1;
%B_est = B;
%R_est = ones(1,4);
PatchSizeHalf=3;
 WindowSizeHalf=5;
 Sigma=10/255;
bicubic_fusion = imresize(Yhim,down_fact);
%X_curr = zeros(rr,cc,p);
psnr_vals = nan(1,maxiters);
iters = 1;
X_curr = var_dim(ima_interp_spline(Yhim,downsamp_factor),pinv(E));
%% Contractive factor calcualtion

 bicubic_fusion;
max_iter=50;
tol_norm=1e-20;
tol_dot=1e-20;
verbose=true;
Sample=rand(nl*nc,p);
stepsize=1.8;



X_curr = var_dim(ima_interp_spline(Yhim,downsamp_factor),pinv(E));
X_curr1 = var_dim(ima_interp_spline(Yhim,downsamp_factor),pinv(E));
%X_curr  = reshape(V,nl*nc,p)';
%X_curr =
beta=power_method_for_beta(nl,SB,SBadj,E,R_est,Yhim,Ymim,lambda_m,lambda_phi,stepsize, Sample, max_iter, tol_norm, tol_dot, verbose,p);
fprintf('The Lipschitz L for gradient_operator calculated is %d\n',beta);
%stepsize=1/beta;
%contrac=power_method_for_denoiser1(nl, SB, SBadj, E, R_est, Yhim, Ymim, lambda_m, lambda_phi, stepsize, Sample, max_iter, tol_norm, tol_dot, verbose, p, h, X_curr1)




lambda_m_values = 0:1:15;
stepsize_values = 0.1:0.1:1.95;
h_values = [5/255,10/255,30/255];
s_values=[5/255,10/255,30/255];

lambda_m_values = 0:1:5;
stepsize_values = 0.5;
h_values = [5/255,10/255];
s_values=[5/255,10/255];
%% Save ss_matrix and contra_matrix for all combinations of s_values and h_values

% Initialize matrices to store contractive factors
ss_matrix = zeros(length(s_values) * length(h_values), length(stepsize_values));
contra_matrix = zeros(length(s_values) * length(h_values), length(lambda_m_values));

% Initialize index counter
counter = 1;

% Loop through s_values first, then h_values
for s_idx = 1:length(s_values)
    Sigma = s_values(s_idx);
    
    for h_idx = 1:length(h_values)
        h = h_values(h_idx);
        ss_values = zeros(size(stepsize_values));
        contra_values = zeros(size(lambda_m_values));

        % Calculate step size vs contractive factor
        for idx = 1:length(stepsize_values)
            stepsize = stepsize_values(idx) / beta;
            ss_values(idx) = power_method_for_denoiser1(nl, SB, SBadj, E, R_est, Yhim, Ymim, lambda_m, lambda_phi, stepsize, Sample, max_iter, tol_norm, tol_dot, verbose, p, h, X_curr1, Sigma);
        end
        
        % Calculate lambda vs contractive factor
        for idx = 1:length(lambda_m_values)
            lambda_m = lambda_m_values(idx);
            beta = power_method_for_beta(nl, SB, SBadj, E, R_est, Yhim, Ymim, lambda_m, lambda_phi, stepsize, Sample, max_iter, tol_norm, tol_dot, verbose, p);
            stepsize = 1.9 / beta;
            contra_values(idx) = power_method_for_denoiser1(nl, SB, SBadj, E, R_est, Yhim, Ymim, lambda_m, lambda_phi, stepsize, Sample, max_iter, tol_norm, tol_dot, verbose, p, h, X_curr1, Sigma);
        end
        
        % Store results in matrices as rows
        ss_matrix(counter, :) = ss_values;
        contra_matrix(counter, :) = contra_values;

        % Save individual results
        ss_filename = sprintf('ss_matrix_s%.3f_h%.3f.mat', Sigma, h);
        contra_filename = sprintf('contra_matrix_s%.3f_h%.3f.mat', Sigma, h);
        save(ss_filename, 'ss_values', 'stepsize_values', 'h', 'Sigma');
        save(contra_filename, 'contra_values', 'lambda_m_values', 'h', 'Sigma');
        
        fprintf('Saved ss_matrix and contra_matrix for s = %.3f, h = %.3f\n', Sigma, h);

        % Update counter
        counter = counter + 1;
    end
end

% Save the full matrices
save('ss_matrix_all_combinations.mat', 'ss_matrix', 's_values', 'h_values', 'stepsize_values');
save('contra_matrix_all_combinations.mat', 'contra_matrix', 's_values', 'h_values', 'lambda_m_values');

fprintf('Saved full ss_matrix and contra_matrix for all (s, h) combinations as 2D matrices!\n');

% Let me know if you’d like me to refine or add anything! 🚀





%% Function for contractive functiom






%  matlab function for power method
function [eigenvalue] = power_method_for_beta(b,SB,SBadj,E,R_est,Yhim,Ymim,lambda_m,lambda_phi,stepsize, input_image, max_iter, tol_norm, tol_dot, verbose,p)
  % Get image dimensions
  [m, n] = size(input_image);
  % Initialize storage for convergence metrics 
  

  % Initialize variables
  x = rand(m, n);
  x = double(x);  % Convert to double precision
  x= x / norm(x, 'fro');  % Normalize
  x1=mat2im(x',b);
  grad= compute_gradient1(x1,SB,SBadj,E,R_est,Yhim,Ymim,lambda_m,lambda_phi);
  y=grad;

  for i = 1:max_iter
    % Update with new eigenvector
    y_old = x;  % Store previous state for distance calculation
    x=mat2im(x',b);
    grad= compute_gradient1(x,SB,SBadj,E,R_est,Yhim,Ymim,lambda_m,lambda_phi);
     gamma= grad;
     y=gamma;
     y=im2mat(gamma)';
    y_norm = norm(y, 'fro');
  
    x = y / y_norm;

  

    % Check for convergence based on norm difference
    if abs(norm(x-y_old, 'fro') ) <= tol_norm
      if verbose
        disp('Converged: Norm difference less than tolerance');
      end
      break;
    end
    
    
  end
 

  % Print final information if verbose
  if verbose
    if i == max_iter
      disp('Not converged');
    end
    fprintf('Iteration: %d\n', i);
  end
  
  % Return the final eigenvalue estimate (norm of y)
  eigenvalue = norm(y, 'fro');
end


function grad = compute_gradient1(X,SB,SBadj,E,R,Yh,Ym,lambda_m,lambda_phi)


rows = size(X,1);

term1 = SB(X);
term1 = SBadj(term1);
term1 = lambda_phi*E'*E*im2mat(term1);

term2 = lambda_m * (R*E)' * (R*E) * im2mat(X);

term3 = SBadj(Yh);
term3 = lambda_phi*(E')*im2mat(term3);

term4 = lambda_m * (R*E)' * im2mat(Ym);

grad = term1 + term2 ;
grad = mat2im(grad,rows);

end



function [eigenvalue] = power_method_for_denoiser1(b,SB,SBadj,E,R_est,Yhim,Ymim,lambda_m,lambda_phi,stepsize, input_image, max_iter, tol_norm, tol_dot, verbose,p,h,X_curr1,Sigma)
  % Get image dimensions
  [m, n] = size(input_image);
  % Initialize storage for convergence metrics (optional)
  denoiser = @(a) wrapper_FASTDSGNLM(a,h,X_curr1);
PatchSizeHalf=3;
 WindowSizeHalf=5;
 
  % Initialize variables
  x = rand(m, n);
  x = double(x);  % Convert to double precision
  x= x / norm(x, 'fro');  % Normalize
  x1=mat2im(x',b);
  grad= compute_gradient1(x1,SB,SBadj,E,R_est,Yhim,Ymim,lambda_m,lambda_phi);
  gamma= x1-stepsize*grad;
  
   y=im2double(gamma);
   
  y=denoiser(y);
B2=y;
    for jj=1:size(B2,3)
    eigen_im=(  B2(:,:,jj));
    
    input =eigen_im;
     %input = gpuArray(input);
     eigen_im_1=(  X_curr1(:,:,jj));
   
    input_1 =eigen_im_1;
  
   res=Fast_DSG_NLM(input, input_1, PatchSizeHalf, WindowSizeHalf, Sigma);
   
   BB = gather(res);  % Removed .x as it's unnecessary here
    V2(:,:,jj) = double(BB) ;
    end
    
    y=V2;
  for i = 1:max_iter
    % Update with new eigenvector
    y_old = x;  % Store previous state for distance calculation
    x=mat2im(x',b);
    grad= compute_gradient1(x,SB,SBadj,E,R_est,Yhim,Ymim,lambda_m,lambda_phi);
     gamma= x-stepsize*grad;
    y=gamma;
     y=denoiser(y);
      B2=y;
    for jj=1:size(B2,3)
    eigen_im=(  B2(:,:,jj));
    
    input =eigen_im;
     %input = gpuArray(input);
     eigen_im_1=(  X_curr1(:,:,jj));
   
    input_1 =eigen_im_1;
  
   res=Fast_DSG_NLM(input, input_1, PatchSizeHalf, WindowSizeHalf, Sigma);
   
   BB = gather(res);  % Removed .x as it's unnecessary here
    V2(:,:,jj) = double(BB) ;
    end
    
    y=V2;
      y=im2mat(y);
    y_norm = norm(y, 'fro');
 
    x = y' / y_norm;

  

    % Check for convergence based on norm difference
    if abs(norm(x-y_old, 'fro') ) <= tol_norm
      if verbose
        disp('Converged: Norm difference less than tolerance');
      end
      break;
    end
   
    
  end
 

  % Print final information if verbose
  if verbose
    if i == max_iter
      disp('Not converged');
    end
    fprintf('Iteration: %d\n', i);
  end
  
  % Return the final eigenvalue estimate (norm of y)
  eigenvalue = norm(y, 'fro');
end

%%

function [eigenvalue] = power_method_for_denoiser_only(b,input_image, max_iter, tol_norm, verbose,h,X_curr1)
  % Get image dimensions
  [m, n] = size(input_image);
  % Initialize storage for convergence metrics (optional)
  denoiser = @(a) wrapper_FASTDSGNLM(a,h,X_curr1);

  % Initialize variables
  x = rand(m, n);
  x = double(x);  % Convert to double precision
  x= x / norm(x, 'fro');  % Normalize
  x1=mat2im(x',b);
  y=denoiser(x1);

  for i = 1:max_iter
    % Update with new eigenvector
    y_old = x;  % Store previous state for distance calculation
     x1=mat2im(x',b);
     y=x1;
     y=denoiser(y)-mat2im((1/b).*ones(1,b)*ones(1,b)'*im2mat(y),b);
    % y=denoiser(y);
      y=im2mat(y);
      
    y_norm = norm(y, 'fro');
 
    x = y' / y_norm;

  

    % Check for convergence based on norm difference
    if abs(norm(x-y_old, 'fro') ) <= tol_norm
      if verbose
        disp('Converged: Norm difference less than tolerance');
      end
      break;
    end
   
    
  end
 

  % Print final information if verbose
  if verbose
    if i == max_iter
      disp('Not converged');
    end
    fprintf('Iteration: %d\n', i);
  end
  
  % Return the final eigenvalue estimate (norm of y)
  eigenvalue = norm(y, 'fro');
end
