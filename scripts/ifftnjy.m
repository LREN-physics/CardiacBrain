function dat = ifftnjy ( dat, DIMs )

% Orhtonormal N-Dimensional inverse Fourier operator along the first and second
% dimension
% 
% dat = ifftnjy ( kdat, DIMs )
%
% Input:
%       dat   - Input Fourier coefficients
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


%dat = kdat; % Initialize output data
N   = 1;    % Number of Fourier coefficients

if nargin == 1
    DIMs = 1:ndims(dat);
elseif nargin == 2
    DIMs = parseInputData(ndims(dat), DIMs);
else
    error('The function ifftnjy can only take 1 or 2 input arguments! ''kdat = fftnjy ( dat, DIMs )''')
end


for i = 1:length(DIMs)
    % iFFT along the dimension dim 
    dim = DIMs(i);
    dat = complex(ifftshift(real(dat),dim),ifftshift(imag(dat),dim));
    dat = ifft(dat,[],dim);
    dat = complex(fftshift(real(dat),dim),fftshift(imag(dat),dim));
    N   = N * size(dat,dim); % Update number of Fourier coefficients
    
    %{
    % VERY SLOW FOR LARGE DATASET
    dim = DIMs(i);
    dat = fftshift(ifft(ifftshift(dat,dim),[],dim),dim);
    N   = N * size(dat,dim); % Update number of Fourier coefficients
    %}
end

% Normalize data
dat = dat * sqrt(N);



function DIMs = parseInputData(nDimensions, DIMs)

DIMs = floor(DIMs); % Only integer values

if isempty(DIMs)
    DIMs = 1:nDimensions;   % Take Fourier transform along all dimensions by default
end

if max(DIMs) > nDimensions
    error('Error in ifftnjy: DIMs cannot contain a value larger than the number of dimensions in kdat!')
elseif min(DIMs) < 0
    error('Error in ifftnjy: DIMs cannot contain a value smaller than 1!')    
end




%{

%dat = kdat; % Initialize output data
N   = 1;    % Number of Fourier coefficients

parseInputData(ndims(dat), DIMs);

for i = 1:length(DIMs)
    % iFFT along the dimension dim
    dim = DIMs(i);
    dat = ifftshift(ifft(fftshift(dat,dim),[],dim),dim);
    N   = N * size(dat,dim); % Update number of Fourier coefficients
end

% Normalize data
dat = dat * sqrt(N);



function DIMs = parseInputData(nDimensions, DIMs)

DIMs = floor(DIMs); % Only integer values

if isempty(DIMs)
    DIMs = 1:nDimensions;   % Take Fourier transform along all dimensions by default
end

if max(DIMs) > nDimensions
    error('Error in ifftnjy: DIMs cannot contain a value larger than the number of dimensions in kdat!')
elseif min(DIMs) < 0
    error('Error in ifftnjy: DIMs cannot contain a value smaller than 1!')    
end



%}