function [DenoisedImg] = Fast_DSG_NLM(NoisyImg, GuideImg, PatchSizeHalf, WindowSizeHalf, Sigma)

% Date: May 14, 2018.

% % Details of the inputs
% NoisyImg: input image to be denoised at each iteration
% GuiedeImg: image from which the FIXED weight can be computed. (for default case, set, NoisyImg = GuideImg)
% Rest all are as usual ...


% % 
[Height,Width] = size(NoisyImg);
[M1,N] = size(NoisyImg);
S = WindowSizeHalf;
K = PatchSizeHalf;

h = 10*Sigma;

u = zeros(Height,Width);  % denoised 
M = u; % Initialize the weight max
W0 = M; % Initialize the accumlated weights

% PaddedImg = padarray(NoisyImg, [K, K],'symmetric','both');
PaddedImg = padarray(GuideImg, [K, K],'symmetric','both');
PaddedV = padarray(NoisyImg,[S, S],'symmetric','both');

% 0-th loop

for dx =  -S:S  
    for dy = -S:S  
        %if dx ~= 0 || dy ~= 0
        % 
        [Sd, diff, t] = integralImgSqDiff(PaddedImg,dx,dy);
        Hat = triangle(dx, dy, S);

        temp1 = img2DShift(Sd, K, K);
        temp2 = img2DShift(Sd, -K-1,- K-1);     % integralImgSqDiff(Sd, K, K);
        temp3 = img2DShift(Sd, -(K+1), K) ;    %integralImgSqDiff(Sd, 0, K);
        temp4 = img2DShift(Sd, K, -(K+1))  ;  %integralImgSqDiff(Sd, K, 0);
        Res = temp1 + temp2 -temp3 -temp4;

        SqDist1 = Res(K+1 : M1+K, K+1 : N+K);

      %     w = exp(-SqDist1/(h^2));
        w = Hat.*exp(-SqDist1/(h^2));    
        W0 = W0 + w; % image containing as sum weights
    end
end

% %
W1 = zeros(size(NoisyImg));
   
% 1st loop

for dx =  -S:S  
    for dy = -S:S  
        %if dx ~= 0 || dy ~= 0
        % 
        [Sd, diff, t] = integralImgSqDiff(PaddedImg,dx,dy);
        Hat = triangle(dx, dy, S);

        temp1 = img2DShift(Sd, K, K);
        temp2 = img2DShift(Sd, -K-1,- K-1);     % integralImgSqDiff(Sd, K, K);
        temp3 = img2DShift(Sd, -(K+1), K) ;    %integralImgSqDiff(Sd, 0, K);
        temp4 = img2DShift(Sd, K, -(K+1))  ;  %integralImgSqDiff(Sd, K, 0);
        Res = temp1 + temp2 -temp3 -temp4;

        SqDist1 = Res(K+1 : M1+K, K+1 : N+K);

      %     w = exp(-SqDist1/(h^2));
        w = Hat.*exp(-SqDist1/(h^2));        
        W0_pad = padarray(W0, [S, S],'symmetric','both');        
        W0_shift = img2DShift(W0_pad, dx, dy);        
        W0_temp = W0_shift(S+1:S+Height, S+1:Width+S);        
        w1 = w./(sqrt(W0).*sqrt(W0_temp));
        
        W1 = W1 + w1; % image containing as sum weights
    end
end

alpha = 1/max(W1(:));
W2 = zeros(size(NoisyImg));
% 2nd loop : Filtering step

for dx =  -S:S  %-WindowSizeHalf :WindowSizeHalf
    for dy = -S:S   %-WindowSizeHalf :WindowSizeHalf
        %
        if dx ~= 0 || dy ~= 0
        % Compute the Integral Image 
        [Sd, diff, t] = integralImgSqDiff(PaddedImg,dx,dy);
        Hat = triangle(dx, dy, S);

        temp1 = img2DShift(Sd, K, K);
        temp2 = img2DShift(Sd, -K-1,- K-1);     % integralImgSqDiff(Sd, K, K);
        temp3 = img2DShift(Sd, -(K+1), K) ;    %integralImgSqDiff(Sd, 0, K);
        temp4 = img2DShift(Sd, K, -(K+1))  ;  %integralImgSqDiff(Sd, K, 0);
        Res = temp1 + temp2 -temp3 -temp4;

        SqDist1 = Res(K+1 : M1+K, K+1 : N+K);

        % Compute the weights for every pixels
        % w = exp(-SqDist1/(h^2));
        w = Hat.*exp(-SqDist1/(h^2));  
        W0_pad = padarray(W0, [S, S],'symmetric','both');        
        W0_shift = img2DShift(W0_pad, dx, dy);        
        W0_temp = W0_shift(S+1:S+Height, S+1:Width+S);        
        w2 = (alpha*w)./(sqrt(W0).*sqrt(W0_temp));
             
       % Obtain the corresponding noisy pixels
        v = PaddedV((S+1+dx):(S+dx+Height),(S+1+dy):(S+dy+Width));
                
       % Compute and accumalate denoised pixels % v = NoisyImg;       
       
        u = u + w2.*v;
        W2 = W2 + w2;
       end      
       
    end
end

u = u + (1 - W2).*NoisyImg;

DenoisedImg = u;
%save W2


function [Sd, diff, t] = integralImgSqDiff(v,dx,dy)
% FUNCTION intergralImgDiff: Compute Integral Image of Squared Difference    
    t = img2DShift(v,dx,dy);    
    diff = (v-t).^2;    
    Sd = cumsum(diff,1);    
    Sd = cumsum(Sd,2);

function hat = triangle(dx, dy, Ns)
    r1 = abs(1 - abs(dx)/(Ns + 1));
    r2 = abs(1 - abs(dy)/(Ns + 1));    
    hat = r1 * r2;


function t = img2DShift(v,dx,dy)
% FUNCTION img2DShift: Shift Image with respect to x and y coordinates
    t = zeros(size(v));
    type = (dx>0)*2+(dy>0);
    switch type
        case 0 % dx<0,dy<0: move lower-right 
            t(-dx+1:end,-dy+1:end) = v(1:end+dx,1:end+dy);
        case 1 % dx<0,dy>0: move lower-left
            t(-dx+1:end,1:end-dy) = v(1:end+dx,dy+1:end);
        case 2 % dx>0,dy<0: move upper-right
            t(1:end-dx,-dy+1:end) = v(dx+1:end,1:end+dy);
        case 3 % dx>0,dy>0: move upper-left
            t(1:end-dx,1:end-dy) = v(dx+1:end,dy+1:end);
    end



% -------------------------------------------------------------------------
