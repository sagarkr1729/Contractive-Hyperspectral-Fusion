function out=FASTDSGNLM(A,S,h,Cluster,Aguide)
% Main filtering code lies here
if ~exist('Aguide','var')
     % guideimage is same as inputimage
     Aguide=A;
end
[m,n,guided]=size(Aguide);

%% Get cluster centres and normalization weights
Ares=reshape(Aguide(1:2:m,1:2:n,:),m*ceil(n)/4,guided);
% Ares=reshape(Aguide(1:m,1:n,:),m*n,guided);

Centre=kmeans_recursive(Ares,Cluster);
B=zeros(size(A));
Cluster=size(Centre,1);
for i=1:Cluster
    W(:,:,i)=sum((Aguide-reshape(Centre(i,:),1,1,guided)).^2,3);   
end
W=exp(-bsxfun(@rdivide,W,(2*(h^2))));
% Wc=W(W~=0);
% W(W==0)=min(min(min(Wc)));
Wb=zeros(m,n);
C = zeros(m,n);
eps=10^-3;

%% Filtering using matlab command for convolutions
    filt     = ones(ceil((S+1)/2),ceil((S+1)/2)+1);  
    for i=1:Cluster
        Wb=Wb+bsxfun(@times,W(:,:,i),imfilter(imfilter(W(:,:,i),filt),filt));
    end
    Wb=real(sqrt(Wb));
    Wb(Wb==0)=eps;
    Achan=A./Wb;
    for i=1:Cluster  
        B=B+bsxfun(@times,W(:,:,i),imfilter(imfilter(bsxfun(@times,W(:,:,i),Achan),filt),filt));   
        C=C+bsxfun(@times,W(:,:,i),imfilter(imfilter(bsxfun(@times,W(:,:,i),1./Wb),filt),filt));                     
    end       
    B=B./Wb;
    C=C./Wb;
    alpha=1/max(max(C));
    B=alpha*B;
    out=B+A-alpha*(bsxfun(@times,C,A));    
end
