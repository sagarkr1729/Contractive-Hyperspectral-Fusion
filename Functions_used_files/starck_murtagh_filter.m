function [B] = starck_murtagh_filter(rows,cols)
%STACK_MURTAGH_FILTER Summary of this function goes here
%   Detailed explanation goes here

middlel = round((rows+1)/2);
middlec = round((cols+1)/2);
B = zeros(rows,cols);
B(middlel-2:middlel+2, middlec-2:middlec+2) = ...
    [1 4 6 4 1;...
     4 16 24 16 4;...
     6 24 36 24 6;...
     4 16 24 16 4;...
     1 4 6 4 1];    % Starck-Murtagh filter
B = ifftshift(B);   % Circularly center B
B = B/sum(sum(B));  % Normalize

end

