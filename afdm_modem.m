% This file is part of the AFDM distribution
% Author: Vincent Savaux, b<>com
% Copyright 2026 b<>com.
 
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
 
%       http://www.apache.org/licenses/LICENSE-2.0
 
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

clear all; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 512; % FFT/IFFT size
c1 = 1.5/(2*N); % post-chirp parameter
c2 = 0.0001; % pre-chirp parameter
cpp_length = N/8; % arbitrary value of the CPP length
snr_dB = 20; % SNR in dB

const_size = 2; % 2,4,6,8 -> QPSK,16,64,256-QAM
const_norm = sqrt([0,2,0,10,0,42,0,170]); % normalization coef for constellation 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modulation/demodulation variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vec_n = 0:N-1; % vector of index
vec_n_cpp = -cpp_length:-1; % index of CPP
vec_pre_chirp = exp(2i*pi*c2*vec_n.^2); 
vec_post_chirp = exp(2i*pi*c1*vec_n.^2); 
vec_cpp = exp(-2i*pi*c1*(N^2-vec_n_cpp).^2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AFDM Modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data generation
dataIn = randi([0 1],const_size*N,1); % binary data generation
data_tx = qammod(dataIn,2^const_size,'gray','InputType','bit'); % binary to complex 
data = data_tx/const_norm(const_size);

% IDAFT Modulation
afdm_signal = diag(vec_post_chirp)*ifft(diag(vec_pre_chirp)*data)*sqrt(N); 

% CCP 
afdm_signal_cpp = [diag(vec_cpp)*afdm_signal(end-cpp_length+1:end);afdm_signal]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel and noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AWGN Channel 
afdm_signal_rx_no_noise = afdm_signal_cpp; 

% Add noise
pow_signal = sum(abs(afdm_signal_rx_no_noise).^2)/(N+cpp_length); 
pow_noise = pow_signal*10^(-snr_dB/10); 
noise = sqrt(pow_noise/2)*(randn(N+cpp_length,1) + 1i*randn(N+cpp_length,1));
afdm_signal_rx = afdm_signal_rx_no_noise + noise; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demodulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CPP removal 
afdm_signal_rx_rem_cpp = afdm_signal_rx; 
afdm_signal_rx_rem_cpp(1:cpp_length) = []; 

% DAFT Demodulation
afdm_signal_demod = diag(vec_pre_chirp')*fft(diag(vec_post_chirp')*afdm_signal_rx_rem_cpp)/sqrt(N);
data_rx = afdm_signal_demod*const_norm(const_size); 
dataOut = qamdemod(data_rx,2^const_size,'gray','OutputType','bit');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BER = length(find(dataIn-dataOut))/(N*const_size); 

figure
plot(real(data_tx),imag(data_tx),'o')
hold
plot(real(data_rx),imag(data_rx),'x')
findfigs




