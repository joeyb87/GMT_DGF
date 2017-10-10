function output = MRifft(input, dim)

% function output = MRifft(input, dim)
%
% MR fft convention which takes
% input:
% with DC Fourier coefficient at the center
%
% output:
% with the center pixel in the middle


% Ernie Yeh ENY 2/14/02

tmp = ifftshift(input);

for idim = dim,
    tmp = ifft(tmp,[],idim);
end

output = fftshift(tmp);