function out = wrapper_FASTDSGNLM(in,sigma,inguide)
if ~exist('inguide','var')
     % guideimage is same as inputimage
     inguide=in;
end

Apca=inguide;

%% Clustering - For reducing computational overhead, we can cluster
%% the guide image as well
% pcadim=6;
% Apca=compute_pca(inguide, 3, pcadim);

out = FASTDSGNLM(in, 5, sigma, 40, Apca);
end