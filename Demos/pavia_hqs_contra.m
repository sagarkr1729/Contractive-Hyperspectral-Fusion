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

[rr,cc,~] = size(Ymim);

lambda_phi = 1;
lambda_m = 6;
h = 10/255 ;
maxiters = 50;
delta = 4;
scale_fact=1;
mu=2;
sigma= 1;
PatchSizeHalf=3;
 WindowSizeHalf=5;
 Sigma=10/255;

bicubic_fusion = imresize(Yhim,down_fact);

tol = 1e-10;

max_iter=1;
tol_norm=1e-5;
tol_dot=1e-5;
verbose=true;
Sample=rand(nl*nc,p);
stepsize=1.8;
psnr_vals = nan(1,maxiters);
iters = 1;
X_curr = var_dim(ima_interp_spline(Yhim,downsamp_factor),pinv(E));
X_curr1 = var_dim(ima_interp_spline(Yhim,downsamp_factor),pinv(E));
%X_curr  = reshape(V,nl*nc,p)';
%X_curr =

%stepsize=1/beta;
%contrac=power_method_for_denoiser1(nl, SB, SBadj, E, R_est, Yhim, Ymim, lambda_m, lambda_phi, stepsize, Sample, max_iter, tol_norm, tol_dot, verbose, p, h, X_curr1)
lambda_m_values = 0:1:15;
mu_values = [0.5,1,2,4,5,10];
h_values = [5/255,10/255, 50/255];
% Initialize matrix to store ss_values (contractive factors) for each h (Step Size vs Contractive Factor)
ss_matrix = zeros(length(h_values), length(mu_values));

% Initialize matrix to store contra_values (contractive factors) for each h (Lambda vs Contractive Factor)
contra_matrix = zeros(length(h_values), length(lambda_m_values));

% Contractive Factor Calculation - Step Size vs Contractive Factor for different h values
 % Adjust the range as needed
figure;
hold on; % Hold the figure to plot multiple lines for different h values

for h_idx = 1:length(h_values)
    h = h_values(h_idx); % Update h value
    ss_values = zeros(size(mu_values)); % Reset ss_values array for each h

   for idx = 1:length(mu_values)
    mu = mu_values(idx); % Now we use mu for each individual value
    ss_values(idx) = power_method_for_denoiser1(nl, SB, SBadj, E, R_est, Yhim, Ymim, lambda_m, lambda_phi, mu, sigma, Sample, max_iter, tol_norm, verbose, p, h, X_curr1, tol, rr, cc);
end

    % Store the ss_values as a row in ss_matrix
    ss_matrix(h_idx, :) = ss_values;

    % Plot the results for each h value
   % plot(stepsize_values, ss_values, 'DisplayName', sprintf('h = %.3f', h));
end

xlabel('Step Size ($\frac{i}{\beta}$)','Interpreter','latex');
ylabel('Contractive factor');
legend show; % Show legend for h values
hold off;

% Contractive Factor Calculation - Lambda vs Contractive Factor for different h values
  % Define the range of lambda_m values
figure;
hold on;

for h_idx = 1:length(h_values)
    h = h_values(h_idx); % Update h value
    contra_values = zeros(size(lambda_m_values));
    step_values = zeros(size(lambda_m_values));

    for idx = 1:length(lambda_m_values)
        lambda_m = lambda_m_values(idx);
        mu=1;
        mu=0.5;
    
   contra_values(idx)= power_method_for_denoiser1(nl, SB, SBadj, E, R_est, Yhim, Ymim, lambda_m, lambda_phi,mu,sigma, Sample, max_iter, tol_norm, verbose, p, h, X_curr1,tol,rr,cc);
    end

    % Store the contra_values as a row in contra_matrix
    contra_matrix(h_idx, :) = contra_values;

    % Plot the results for each h value
    plot(lambda_m_values, contra_values, 'DisplayName', sprintf('h = %.3f', h));
end

xlabel('$\lambda_m$','Interpreter','latex');
ylabel('Contractive factor');
legend show; % Show legend for h values
hold off;



%%
function [eigenvalue] = power_method_for_denoiser1(b,SB,SBadj,E,R_est,Yhim,Ymim,lambda_m,lambda_phi,mu,sigma, input_image, max_iter, tol_norm,  verbose,p,h,X_curr1,tol,rr,cc)
  % Get image dimensions
  [m, n] = size(input_image);
  % Initialize storage for convergence metrics (optional)
  denoiser = @(a) wrapper_FASTDSGNLM(a,h,X_curr1);
PatchSizeHalf=3;
 WindowSizeHalf=5;
 Sigma=10/255;
  % Initialize variables
  x = rand(m, n);
  x = double(x);  % Convert to double precision
  x= x / norm(x, 'fro');  % Normalize
  x1=mat2im(x',b);
  C= compute_C(x1,SBadj,E,R_est,Yhim,Ymim,lambda_m,lambda_phi,mu,sigma);
  gamma = solve_AX_plus_XB_with_K(SB,SBadj,E,R_est,lambda_m,lambda_phi,mu,sigma, C, tol, max_iter,rr,cc,p);
  y=gamma;
  y=mat2im(y',b);
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
    C= compute_C(x,SBadj,E,R_est,Yhim,Ymim,lambda_m,lambda_phi,mu,sigma);
     gamma = solve_AX_plus_XB_with_K(SB,SBadj,E,R_est,lambda_m,lambda_phi,mu,sigma, C, tol, max_iter,rr,cc,p);
  y=gamma;
   y=mat2im(y',b);
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

C=(term7)';
%C=mat2im(C,rows);
end

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
