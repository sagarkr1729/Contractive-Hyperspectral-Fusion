%% Optimized version

function [img, den]=fastKmeansfiltapproxnystromtriangle(img,S,h,Centre,imgguide)
if ~exist('imgguide','var')
     % guideimage is same as inputimage
     imgguide=img;
end
[m,n,~]=size(img);
num=zeros(size(img));
guided=size(imgguide,3);
imgguide=reshape(imgguide,m*n,guided);
Cluster=size(Centre,1);
A=zeros(Cluster,Cluster);
for i=1:Cluster
    A(i,i)=1;
    for j=i+1:Cluster
        A(i,j)=exp(-sum((Centre(i,:)-Centre(j,:)).^2,2)/(h*h));
        A(j,i)=A(i,j);
    end
end
[Vhat,Dhat]=eig(A);
Vhat = flip(Vhat, 2);
Dhat = flip(flip(Dhat,2),1);
Dhat=diag(Dhat);
W=zeros(m*n,Cluster);
for i=1:Cluster
    W(:,i)=sum((imgguide-Centre(i,:)).^2,2);   
end
W=exp(-bsxfun(@rdivide,W,((h^2))));
Eigmat=reshape(W*Vhat,m,n,Cluster);
%% Interpolating numerator and denominator seperately
den=zeros(m,n);

% Uncomment for matlab implementation
filt     = ones(ceil((S+1)/2),ceil((S+1)/2)); 
%% Calculating Bilateral image for each cluster centre as index pixel
    for i=1:Cluster
         den=den+bsxfun(@times,(1/Dhat(i,1))*Eigmat(:,:,i),imfilter(imfilter(Eigmat(:,:,i),filt),filt));
         num=num+bsxfun(@times,(1/Dhat(i,1))*Eigmat(:,:,i),imfilter(imfilter(bsxfun(@times,Eigmat(:,:,i),img),filt),filt)); 
%          den=den+bsxfun(@times,(1/Dhat(i,1))*Eigmat(:,:,i),imfilter(Eigmat(:,:,i),filt));
%          num=num+bsxfun(@times,(1/Dhat(i,1))*Eigmat(:,:,i),imfilter(bsxfun(@times,Eigmat(:,:,i),img),filt));
    end
img=bsxfun(@rdivide,num,den);
end    

