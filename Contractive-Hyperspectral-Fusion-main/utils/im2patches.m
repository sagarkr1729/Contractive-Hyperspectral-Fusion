function [patches,rowinds,colinds] = im2patches(I,patchsize,strides)
%IM2PATCHES Summary of this function goes here
%   Detailed explanation goes here

[R,C] = size(I);
rowstride = strides(1);
colstride = strides(2);

rowpad = max(patchsize(1),rowstride);
colpad = max(patchsize(2),colstride);
I = padarray(I,[rowpad,colpad],'replicate','post');

valid_colinds = 1:strides(2):C;
valid_rowinds = 1:strides(1):R;
Clen = length(valid_colinds);
Rlen = length(valid_rowinds);

patches = zeros(patchsize(1)*patchsize(2),Rlen*Clen);

colinds = repmat(valid_colinds,1,Rlen);
rowinds = reshape(repmat(valid_rowinds,Clen,1),1,[]);

for k1 = 1:Rlen
    for k2 = 1:Clen
        P = I(valid_rowinds(k1):valid_rowinds(k1)+patchsize(1)-1,...
              valid_colinds(k2):valid_colinds(k2)+patchsize(2)-1);
        patches(:,(k1-1)*Clen+k2) = P(:);
    end
end

end
