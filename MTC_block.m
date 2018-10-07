%% Main parameters
close all;

numberSyms = 2000; %Number of symbols
fc = 2e6;        %carrier frequency offset
T_SC = 1e-6;     %single carrier modulation symbol duration
t_obs = 2e-3;    %duration received signal is observed
fs = 4e6;        %sampling rate
Ts = 1/fs;     
n_samples = t_obs*fs;   %number of samples in y(n)

%OFDM parameters
N = 64;             %No. of carriers
delta_f = 15.625e3; %Sub-carrier spacing (Hz)
T_OFDM = 80e-6;     %OFDM symbol duration
cp = 16; %number of samples of cyclic prefix

%% Generate all the signals


%%%%%%%%%%%%%%%%%%%%%%%%
% 2. 4-QAM
M = 4;
%Generate 2000 random symbols
b_in = randi([0 1],log2(M)*numberSyms,1);
%QAM modulation
syms = qammod(b_in, M,'InputType','bit','UnitAveragePower',true);
%Upsample
syms_up = upsample(syms, 4);
size(syms_up)
%Pulse shaping filter
g = rcosdesign(0.5,8,4, 'normal');
%Symbols after upsampling and pulse shaping
x_conv = conv(syms_up, g);
%Normalize power of transmitted signal
x_4QAM = x_conv/(rms(x_conv));     
%x_4QAM = x_conv;


%%%%%%%%%%%%%%%%%%%%%%%%
% 3. 16-QAM
M = 16;
%Generate 2000 random symbols
b_in = randi([0 1],log2(M)*numberSyms,1);
%QAM modulation
syms = qammod(b_in, M,'InputType','bit','UnitAveragePower',true);
%Upsample by 4
syms_up = upsample(syms, 4);
%Pulse shaping filter
g = rcosdesign(0.5,8,4, 'normal');
%Symbols after upsampling and pulse shaping
x_conv2 = conv(syms_up, g);
%Normalize power of transmitted signal
x_16QAM = x_conv2./rms(x_conv2); 
%x_16QAM = x_conv;

%%%%%%%%%%%%%%%%%%%%%%%%
% 4. 2-PAM
M = 2;
%Generate 2000 random symbols
b_in = randi([0 M-1],log2(M)*numberSyms,1);
%QAM modulation
syms2p = pammod(b_in, M);
%Upsample by 4
syms_up = upsample(syms2p, 4);
%Pulse shaping filter
g = rcosdesign(0.5,8,4, 'normal');
%Symbols after upsampling and pulse shaping
x_conv3 = conv(syms_up, g);
%Normalize power of transmitted signal
x_2PAM = x_conv3/(rms(x_conv3)); 


%%%%%%%%%%%%%%%%%%%%%%%%
% 5. 4-PAM
M = 4;
%Generate 2000 random symbols
b_in = randi([0 M-1],log2(M)*numberSyms,1);
%QAM modulation
syms4p = pammod(b_in, M);
%Upsample by 4
syms_up = upsample(syms4p, 4);
%Pulse shaping filter
g = rcosdesign(0.5,8,4, 'normal');
%Symbols after upsampling and pulse shaping
x_conv4 = conv(syms_up, g);
%Normalize power of transmitted signal
x_4PAM = x_conv4/(rms(x_conv4)); 



%% Block 2: MTC

%%%%%%%%%%%
%%% PAM %%%
%%%%%%%%%%%
%Ideal case = no noise
N = length(x_2PAM);
n = 0:N-1;
y = transpose(x_2PAM).*exp(1i*2*pi*fc*Ts*n);


R_1 = 1/N * sum(y.*conj(y).* exp(-1i*2*pi*n*Ts/T_SC));
R_2 = 1/N * sum(y.^2 .* exp(-1i*2*pi*n*Ts*(2*fc-1/T_SC)));
R_3 = 1/N * sum(y.^2 .* exp(-1i*2*pi*n*Ts*(2*fc-1/(2*T_SC))));
R_4 = 1/N * sum(y.^2 .* exp(-1i*2*pi*n*Ts*2*fc));
R_5 = 1/N * sum(y.^2 .* exp(-1i*2*pi*n*Ts*(2*fc+1/(2*T_SC))));
R_6 = 1/N * sum(y.^2 .* exp(-1i*2*pi*n*Ts*(2*fc+1/T_SC)));

V_PAM = [abs(R_1), abs(R_2), abs(R_3), abs(R_4), abs(R_5), abs(R_6)];
    
%Normalize vector V_PAM
V_PAM = V_PAM/norm(V_PAM); 

%%%%%%%%%%%
%%% QAM %%%
%%%%%%%%%%%
%Ideal case = no noise
y = ch_out_SC;  %4-QAM
N = length(y);
n = 0:N-1;

R_1 = 1/N * sum(y.*conj(y).* exp(-1i*2*pi*n*Ts/T_SC));
R_2 = 1/N * sum(y.^2 .* exp(-1i*2*pi*n*Ts*(2*fc-1/T_SC)));
R_3 = 1/N * sum(y.^2 .* exp(-1i*2*pi*n*Ts*(2*fc-1/(2*T_SC))));
R_4 = 1/N * sum(y.^2 .* exp(-1i*2*pi*n*Ts*2*fc));
R_5 = 1/N * sum(y.^2 .* exp(-1i*2*pi*n*Ts*(2*fc+1/(2*T_SC))));
R_6 = 1/N * sum(y.^2 .* exp(-1i*2*pi*n*Ts*(2*fc+1/T_SC)));

V_QAM = [abs(R_1), abs(R_2), abs(R_3), abs(R_4), abs(R_5), abs(R_6)];

%Normalize vector V_QAM
V_QAM = V_QAM/norm(V_QAM); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SNR_vect = -10:2:20;

%Real case (with noise)
F = zeros(1, 6);
p_QAM = zeros(1, length(SNR_vect));
p_PAM = zeros(1, length(p_QAM));

%We are sending 16-QAM
N = length(x_16QAM);
n = 0:N-1;
ch_out_QAM = transpose(x_16QAM).*exp(1i*2*pi*fc*Ts*n);

%We are sending 4-PAM
N2 = length(x_4PAM);
n2 = 0:N2-1;
ch_out_PAM = transpose(x_4PAM).*exp(1i*2*pi*fc*Ts*n2);

l = 1; %index of p_QAM and p_PAM

for SNR=SNR_vect
    SNR_lin = 10^(SNR/10); %SNR linear value
    P_N = 1/SNR_lin; %Noise power
    
    count_QAM = 0;
    count_PAM = 0;
    
    for iter=1:1000
        
        %PAM
        w2 = sqrt(P_N/2)* (randn(size(ch_out_QAM)) + 1j* randn(size(ch_out_QAM))); 
        
        y2 = ch_out_QAM + w2;
            
        R_1 = 1/N * sum(y2.*conj(y2).* exp(-1i*2*pi*n*Ts/T_SC));
        R_2 = 1/N * sum(y2.^2 .* exp(-1i*2*pi*n*Ts*(2*fc-1/T_SC)));
        R_3 = 1/N * sum(y2.^2 .* exp(-1i*2*pi*n*Ts*(2*fc-1/(2*T_SC))));
        R_4 = 1/N * sum(y2.^2 .* exp(-1i*2*pi*n*2*fc*Ts));
        R_5 = 1/N * sum(y2.^2 .* exp(-1i*2*pi*n*Ts*(2*fc+1/(2*T_SC))));
        R_6 = 1/N * sum(y2.^2 .* exp(-1i*2*pi*n*Ts*(2*fc+1/T_SC)));

        F = [abs(R_1), abs(R_2), abs(R_3), abs(R_4), abs(R_5), abs(R_6)];
        
        F = F/norm(F); %we normalize F to compare it to V
            
        
        if(norm(F-V_QAM)^2 < norm(F-V_PAM)^2)
            count_QAM = count_QAM + 1;
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% PAM
        %Complex gaussian noise
        w3 = sqrt(P_N/2)* (randn(size(ch_out_PAM)) + 1j* randn(size(ch_out_PAM))); 
        
        y3 = ch_out_PAM + w3;
            
        R_1 = 1/N2 * sum(y3.*conj(y3).* exp(-1i*2*pi*n2*Ts/T_SC));
        R_2 = 1/N2 * sum(y3.^2 .* exp(-1i*2*pi*n2*Ts*(2*fc-1/T_SC)));
        R_3 = 1/N2 * sum(y3.^2 .* exp(-1i*2*pi*n2*Ts*(2*fc-1/(2*T_SC))));
        R_4 = 1/N2 * sum(y3.^2 .* exp(-1i*2*pi*n2*2*fc*Ts));
        R_5 = 1/N2 * sum(y3.^2 .* exp(-1i*2*pi*n2*Ts*(2*fc+1/(2*T_SC))));
        R_6 = 1/N2 * sum(y3.^2 .* exp(-1i*2*pi*n2*Ts*(2*fc+1/T_SC)));

        F = [abs(R_1), abs(R_2), abs(R_3), abs(R_4), abs(R_5), abs(R_6)];
        
        F = F/norm(F); %we normalize F to compare it to V
            
        
        if(norm(F-V_QAM)^2 > norm(F-V_PAM)^2)
            count_PAM = count_PAM + 1;
        end
        
    end
    p_QAM(l) = count_QAM/1000;
    p_PAM(l) = count_PAM/1000;
    l = l + 1;
end

figure
plot(SNR_vect, p_QAM, SNR_vect, p_PAM, 'LineWidth', 1.5)
title('Average probability of correct classification vs SNR')
legend('M-QAM', 'M-PAM')
xlabel('SNR (dB)')
ylabel('P_{detection}')
grid on


