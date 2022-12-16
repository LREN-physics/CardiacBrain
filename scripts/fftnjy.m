function kdat = fftnjy ( kdat, DIMs )

% Orhtonormal N-Dimensional Fourier operator along the first and second
% dimension
% 
% dat = fftnjy ( kdat, DIMs )
%
% Input:
%       kdat  - Input Fourier coefficients
%       DIMs  - Array of dimensions along which to compute the Fourier
%               transform, e.g., [1 2 4] takes the Fourier
%               transform along the 1st, 2nd and 4th dimensions. If DIMs is
%               an empty arrya, then fftnjy takes the Fourier transform
%               along all dimensions of kdat;
%
% Output:
%       dat   - Output data
%
% (c) Jerome Yerly 2012


% NOTE: Taking the fftshift on larde complex dataset is extremely slow! It
%       is significantly faster to process real and imaginary data 
%       separately...
%
%         clear all;
%         clc;
% 
%         N=[320 320 10 12 10];
%         x0=rand(N)+1i*rand(N);
%         dim=2;
%         fprintf('fftshift\n');
%         tic; xr=fftshift(real(x0),dim); toc;
%         tic; xi=fftshift(imag(x0),dim); toc;
%         tic; x1 =complex(xr,xi); toc;
%         fprintf('fft\n');
%         tic; x1=fft(x1,[],dim); toc;
%         fprintf('fftshift\n');
%         tic; xr=fftshift(real(x1),dim); toc;
%         tic; xi=fftshift(imag(x1),dim); toc;
%         tic; x1 =complex(xr,xi); toc;
% 
%         clear xr xi;
% 
%         fprintf('\n');
%         fprintf('fftshift\n');
%         tic; x2=fftshift(x0,dim); toc;
%         fprintf('fft\n');
%         tic; x2=fft(x2,[],dim); toc;
%         fprintf('fftshift\n');
%         tic; x2=fftshift(x2,dim); toc;
% 
%         diff = x1 - x2;
%         max(abs(diff(:))),
% 

% kdat = dat; % Initialize output data
N    = 1;    % Number of Fourier coefficients

if nargin == 1
    DIMs = 1:ndims(kdat);
elseif nargin == 2
    DIMs = parseInputData(ndims(kdat), DIMs);
else
    error('The function fftnjy can only take 1 or 2 input arguments! ''kdat = fftnjy ( dat, DIMs )''')
end


for i = 1:length(DIMs)
    % FFT along the dimension dim
    dim  = DIMs(i);
    kdat = complex(ifftshift(real(kdat),dim),ifftshift(imag(kdat),dim));
    kdat = fft(kdat,[],dim);
    kdat = complex(fftshift(real(kdat),dim),fftshift(imag(kdat),dim));
    N    = N * size(kdat,dim); % Update number of Fourier coefficients
    
    %{
    % VERY SLOW FOR LARGE DATASET
    dim = DIMs(i);
    kdat = fftshift(fft(ifftshift(kdat,dim),[],dim),dim);
    N   = N * size(kdat,dim); % Update number of Fourier coefficients
    %}
end

% Normalize data
kdat = kdat / sqrt(N);



function DIMs = parseInputData(nDimensions, DIMs)

DIMs = floor(DIMs); % Only integer values

if isempty(DIMs)
    DIMs = 1:nDimensions;   % Take Fourier transform along all dimensions by default
end

if max(DIMs) > nDimensions
    warning('Warning in ifftnjy: DIMs cannot contain a value larger than the number of dimensions in kdat!')
    %error('Error in ifftnjy: DIMs cannot contain a value larger than the number of dimensions in kdat!')
elseif min(DIMs) < 0
    error('Error in ifftnjy: DIMs cannot contain a value smaller than 1!')    
end