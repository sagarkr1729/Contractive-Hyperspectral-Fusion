function [I] = patches2im(patches,rowinds,colinds,patchsize,imsize)
%PATCHES2IM Create image from a set of (possibly overlapping) patches

R = imsize(1);
C = imsize(2);
I = zeros(R+patchsize(1),C+patchsize(2));
counts = zeros(R+patchsize(1),C+patchsize(2));

for k = 1:size(patches,2)
    P = reshape(patches(:,k),patchsize);
    I(rowinds(k):rowinds(k)+patchsize(1)-1,...
      colinds(k):colinds(k)+patchsize(2)-1) = ...
    I(rowinds(k):rowinds(k)+patchsize(1)-1,...
      colinds(k):colinds(k)+patchsize(2)-1) + P;
    counts(rowinds(k):rowinds(k)+patchsize(1)-1,...
      colinds(k):colinds(k)+patchsize(2)-1) = ...
    counts(rowinds(k):rowinds(k)+patchsize(1)-1,...
      colinds(k):colinds(k)+patchsize(2)-1) + 1;
end
I = I(1:R,1:C);
counts = counts(1:R,1:C);
I = I./counts;

end

