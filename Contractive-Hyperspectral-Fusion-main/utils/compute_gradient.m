function grad = compute_gradient(X,SB,SBadj,E,R,Yh,Ym,lambda_m,lambda_phi)
%COMPUTE_GRADIENT Gradient of data fidelity term for HS-MS image fusion
% X = Input image at which gradient is to be computed
% SB = Handle to function that performs spatial degradation (blurring and
% downsampling)
% SBadj = Handle to function that performs the adjoint operation
% (upsampling and blurring)
% E = Matrix of PCA/VCA basis vectors
% R = Spectral degradation matrix
% Yh = Observed hyperspectral image
% Ym = Observed multispectral image
% lambda = Regularization parameter which decides relative weightage of
% spectral degradation objective

rows = size(X,1);

term1 = SB(X);
term1 = SBadj(term1);
term1 = lambda_phi*E'*E*im2mat(term1);

term2 = lambda_m * (R*E)' * (R*E) * im2mat(X);

term3 = SBadj(Yh);
term3 = lambda_phi*(E')*im2mat(term3);

term4 = lambda_m * (R*E)' * im2mat(Ym);

grad = term1 + term2 - (term3 + term4);
grad = mat2im(grad,rows);

end

