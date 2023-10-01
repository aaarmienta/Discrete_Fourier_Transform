%% Reset
clc, clear all, close all;

%% Specifications
N = 4;                                  % N-Point DFT
L = 4;                                  % Sequence Length (Not currently used)
n = 0:N-1;                              % N Indexes
k = 0:N-1;                              % N Frequency Bins

%% Twiddle Factor
W_N = exp(-1*j*2*pi/N);                 % Twiddle Factor

% Input Sequence
xn = [1 2 3 4];                         % Real Sequence
xn_cplx = [1+2j 3+4j 2+2j 3-4j];        % Complex Sequence

%% -----------------------------------------------------------------------------
% Linearity Property -----------------------------------------------------------
xn_lin = 3*xn_cplx + 2*xn_cplx;

% Discrete Fourier Transform: xn -> Xk
Xk = zeros(1,N);
for k = 0:N-1
  for n = 0:N-1
    temp = xn_lin(n+1)*W_N^(n*k);
    Xk(k+1) = Xk(k+1) + temp;
  end
end

% Inverse Discrete Fourier Transform: Xk -> xn
xn_lin_i = zeros(1,N);
for n = 0:N-1
  for k = 0:N-1
    temp = Xk(k+1)*W_N^(-1*n*k);
    xn_lin_i(n+1) = xn_lin_i(n+1) + temp;
  end
  xn_lin_i(n+1) = (xn_lin_i(n+1)/N);
end
% ------------------------------------------------------------------------------

%% -----------------------------------------------------------------------------
% Complex Conjugate Property ---------------------------------------------------
xn_conj = conj(xn_cplx);

% Discrete Fourier Transform: xn_conj -> Xk_conj_rev
Xk_conj_rev = zeros(1,N);
for k = 0:N-1
  for n = 0:N-1
    temp = xn_conj(n+1)*W_N^(n*k);
    Xk_conj_rev(k+1) = Xk_conj_rev(k+1) + temp;
  end
end

% Inverse Discrete Fourier Transform: Xk_conj_rev -> xn_conj
xn_conj_i = zeros(1,N);
for n = 0:N-1
  for k = 0:N-1
    temp = Xk_conj_rev(k+1)*W_N^(-1*n*k);
    xn_conj_i(n+1) = xn_conj_i(n+1) + temp;
  end
  xn_conj_i(n+1) = (xn_conj_i(n+1)/N);
end
% ------------------------------------------------------------------------------

%% -----------------------------------------------------------------------------
% Circular Time Shift Property -------------------------------------------------
xn_circ_time = xn_cplx;
shift = 2;
l = shift;
if abs(l) > N
  l = mod(l,N);
end

if l < 0
  xn_circ_time = [xn_circ_time(l+1:N) xn_circ_time(1:l)];
elseif l > 0
  xn_circ_time = [xn_circ_time(N-l:N) xn_circ_time(1:N-l-1)];
end

% Discrete Fourier Transform: xn_circ_time -> Xk_circ_time
Xk_circ_time = zeros(1,N);
for k = 0:N-1
  for n = 0:N-1
    temp = xn_circ_time(n+1)*W_N^(n*k);
    Xk_circ_time(k+1) = Xk_circ_time(k+1) + temp;
  end
end

% Inverse Discrete Fourier Transform : Xk_circ_time -> xn_circ_time
xn_i_circ_time = zeros(1,N);
for n = 0:N-1
  for k = 0:N-1
    temp = Xk_circ_time(k+1)*W_N^(-1*n*k);
    xn_i_circ_time(n+1) = xn_i_circ_time(n+1) + temp;
  end
  xn_i_circ_time(n+1) = (xn_i_circ_time(n+1)/N);
end
% ------------------------------------------------------------------------------

%% -----------------------------------------------------------------------------
% Time Reversal Property -------------------------------------------------------
xn_rev = zeros(1,N);
for n = 0:N-1
  xn_rev(n+1) = xn_cplx(N-n);
end
xn_rev = [xn_rev(N) xn_rev(1:N-1)];

% Discrete Fourier Transform: xn_rev -> Xk_rev
Xk_rev = zeros(1,N);
for k = 0:N-1
  for n = 0:N-1
    temp = xn_rev(n+1)*W_N^(n*k);
    Xk_rev(k+1) = Xk_rev(k+1) + temp;
  end
end

% Inverse Discrete Fourier Transform: Xk_rev -> xn_rev
xn_i_rev = zeros(1,N);
for n = 0:N-1
  for k = 0:N-1
    temp = Xk_rev(k+1)*W_N^(-1*n*k);
    xn_i_rev(n+1) = xn_i_rev(n+1) + temp;
  end
  xn_i_rev(n+1) = (xn_i_rev(n+1)/N);
end
% ------------------------------------------------------------------------------

%% -----------------------------------------------------------------------------
% Frequency Shift Property -----------------------------------------------------
shift = 1;
l = shift;
n = 0:N-1;
fsp = exp(1j*2*pi*l.*n./N);
xn_circ_freq = xn_cplx.*fsp;

% Discrete Fourier Transform: xn_circ_freq -> Xk_l
Xk_l = zeros(1,N);
for k = 0:N-1
  for n = 0:N-1
    temp = xn_circ_freq(n+1)*W_N^(n*k);
    Xk_l(k+1) = Xk_l(k+1) + temp;
  end
end

% Inverse Discrete Fourier Transform: Xk_rev -> xn_rev
xn_i_circ_freq = zeros(1,N);
for n = 0:N-1
  for k = 0:N-1
    temp = Xk_l(k+1)*W_N^(-1*n*k);
    xn_i_circ_freq(n+1) = xn_i_circ_freq(n+1) + temp;
  end
  xn_i_circ_freq(n+1) = (xn_i_circ_freq(n+1)/N);
end
% ------------------------------------------------------------------------------

%% -----------------------------------------------------------------------------
%*******************************************************************************
% Issues
%   1) Taking transpose of complex vector causes it to be conjugated
%       https://www.mathworks.com/help/matlab/ref/ctranspose.html
%*******************************************************************************
