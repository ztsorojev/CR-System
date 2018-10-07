%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ECE233 Project 4         %%%
%%% Author: Zaurbek Tsorojev %%%   
%%% Date: 06/03/2018         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Main parameters
close all;

numberSyms = 2000; %Number of symbols
fc = 2e6;        %carrier frequency offset
T_SC = 1e-6;     %single carrier modulation symbol duration
t_obs = 2e-3;    %duration received signal is observed
Ts = 2.5e-7;     %sampling rate
fs = 1/Ts;
n_samples = t_obs*fs;   %number of samples in y(n)

%OFDM parameters
N = 64;             %No. of carriers
delta_f = 15.625e3; %Sub-carrier spacing (Hz)
T_OFDM = 80e-6;     %OFDM symbol duration
cp = 16; %number of samples of cyclic prefix

%% Generate all the signals

% 1. OFDM
M = 4;
%Generate 25 random symbols
b_in = randi([0 1],log2(M)*25,64);
%QAM modulation
syms = qammod(b_in, M,'InputType','bit','UnitAveragePower',true);
%64-pt IFFT (argument 2 makes it do ifft on rows)
x_ifft = ifft(syms, N, 2);
%P to S
%x_serial = transpose(x_ifft);
x_serial = x_ifft;
%Add cyclic prefix
x_wCP = [x_serial(:,N-cp+1:N), x_serial];
%P to S
x_wCP = transpose(x_wCP);
%Upsample by 4
x_up = upsample(x_wCP, 4);

%Low Pass Filter
%w_c = f_cutoff/(fs/2);
[b, a] = butter(8, 0.25);
x_filter = filter(b, a, x_up);
x_filter = reshape(x_filter, 1, []);
x_OFDM = x_filter./(rms(x_filter));     %Normalize power of transmitted signal

x_OFDM = transpose(x_OFDM);



%%%%%%%%%%%%%%%%%%%%%%%%
% 2. 4-QAM
M = 4;
%Generate 1000 random symbols
b_in = randi([0 1],log2(M)*numberSyms,1);
%QAM modulation
syms = qammod(b_in, M,'InputType','bit','UnitAveragePower',true);
%Upsample
syms_up = upsample(syms, 4);
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
%Generate 1000 random symbols
b_in = randi([0 1],log2(M)*numberSyms,1);
%QAM modulation
syms = qammod(b_in, M,'InputType','bit','UnitAveragePower',true);
%Upsample by 4
syms_up = upsample(syms, 4);
%Pulse shaping filter
g = rcosdesign(0.5,8,4, 'normal');
%Symbols after upsampling and pulse shaping
x_conv = conv(syms_up, g);
%Normalize power of transmitted signal
x_16QAM = x_conv/rms(x_conv); 
%x_16QAM = x_conv;



%%%%%%%%%%%%%%%%%%%%%%%%
% 4. 2-PAM
M = 2;
%Generate 1000 random symbols
b_in = randi([0 M-1],log2(M)*numberSyms,1);
%QAM modulation
syms2p = pammod(b_in, M);
%Upsample by 4
syms_up = upsample(syms2p, 4);
%Pulse shaping filter
g = rcosdesign(0.5,8,4, 'normal');
%Symbols after upsampling and pulse shaping
x_conv2 = conv(syms_up, g);
%Normalize power of transmitted signal
x_2PAM = x_conv2/(rms(x_conv2));


%%%%%%%%%%%%%%%%%%%%%%%%
% 5. 4-PAM
M = 4;
%Generate 1000 random symbols
b_in = randi([0 M-1],log2(M)*numberSyms,1);
%QAM modulation
syms4p = pammod(b_in, M);
%Upsample by 4
syms_up = upsample(syms4p, 4);
%Pulse shaping filter
g = rcosdesign(0.5,8,4, 'normal');
%Symbols after upsampling and pulse shaping
x_conv3 = conv(syms_up, g);
%Normalize power of transmitted signal
x_4PAM = x_conv3/(rms(x_conv3)); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MTC: calculate theoritical V feature vectors


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
N = length(x_4QAM);
n = 0:N-1;
y = transpose(x_4QAM).*exp(1i*2*pi*fc*Ts*n);

R_1 = 1/N * sum(y.*conj(y).* exp(-1i*2*pi*n*Ts/T_SC));
R_2 = 1/N * sum(y.^2 .* exp(-1i*2*pi*n*Ts*(2*fc-1/T_SC)));
R_3 = 1/N * sum(y.^2 .* exp(-1i*2*pi*n*Ts*(2*fc-1/(2*T_SC))));
R_4 = 1/N * sum(y.^2 .* exp(-1i*2*pi*n*Ts*2*fc));
R_5 = 1/N * sum(y.^2 .* exp(-1i*2*pi*n*Ts*(2*fc+1/(2*T_SC))));
R_6 = 1/N * sum(y.^2 .* exp(-1i*2*pi*n*Ts*(2*fc+1/T_SC)));

V_QAM = [abs(R_1), abs(R_2), abs(R_3), abs(R_4), abs(R_5), abs(R_6)];

%Normalize vector V_QAM
V_QAM = V_QAM/norm(V_QAM); 



%************************************************************************
%%% MF: Matched filter and symbol extractor

% sample corresponding to the first local max
T_offset = length(g);

%Number of constellation points
M_4QAM = 4;
x_4 = 1/sqrt(2) * [1 + 1i, 1 - 1i, -1 + 1i, -1 - 1i];

M_16QAM = 16;
x_16 = 1/sqrt(10) * [1 + 1i, 1 - 1i, -1 + 1i, -1 - 1i, 3 + 3i, 3 - 3i, ...
    -3 + 3i, -3 - 3i, 3 + 1i, 3 - 1i, -3 + 1i, -3 - 1i, 1 + 3i, 1 - 3i, -1 + 3i, -1 - 3i];


%%%%%%%%%%%%%%%%%%%%%%%%
% 2-PAM with larger number of symbols to calculate the ECDF (that
% approximates the theoretical CDF)
M = 2;
%Generate 1000 random symbols
b_in = randi([0 M-1],log2(M)*numberSyms*5,1);
%QAM modulation
syms2p = pammod(b_in, M);
%Upsample by 4
syms_up = upsample(syms2p, 4);
%Pulse shaping filter
g = rcosdesign(0.5,8,4, 'normal');
%Symbols after upsampling and pulse shaping
x_conv2 = conv(syms_up, g);
%Normalize power of transmitted signal
x_2PAM_large = x_conv2/(rms(x_conv2)); 


%%%%%%%%%%%%%%%%%%%%%%%%
% 4-PAM with larger number of symbols to calculate the ECDF (that
% approximates the theoretical CDF)
M = 4;
%Generate 1000 random symbols
b_in = randi([0 M-1],log2(M)*numberSyms*5,1);
%QAM modulation
syms4p = pammod(b_in, M);
%Upsample by 4
syms_up = upsample(syms4p, 4);
%Pulse shaping filter
g = rcosdesign(0.5,8,4, 'normal');
%Symbols after upsampling and pulse shaping
x_conv = conv(syms_up, g);
%Normalize power of transmitted signal
x_4PAM_large = x_conv/(rms(x_conv)); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*************************************************************************%

%% Block 1: Mc-SC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of samples used for distinguishing MC and SC signals

%OFDM:
%2ms observation with symbol duration = 80us => 25 symbols = 50 samples for
%OFDM (4-QAM)
Nm = 8000;  

q = 1;          %Index in for loop SNR
probaPts = 100;  %Number of points used for the P_fa vs P_d plot

SNR_vect = 0:10:20;
%SNR_vect = 0;

%Preallocation of probability vectors
P_OFDM = zeros(5, length(SNR_vect));
P_fa = zeros(5, length(SNR_vect));

%Count the proba to detect or false alarm
p_4QAM = zeros(5, length(SNR_vect));
p_16QAM = zeros(5, length(SNR_vect));

p_2PAM = zeros(5, length(SNR_vect));
p_4PAM = zeros(5, length(SNR_vect));

%%%%%%%% MTC %%%%%%%%
F = zeros(1, 6);
p_QAM = zeros(1, length(SNR_vect));
p_PAM = zeros(1, length(p_QAM));

p_c = zeros(1, length(SNR_vect));

%figure; hold on
for SNR = SNR_vect  
    
    SNR_lin = 10^(SNR/10); %SNR linear value

    %OFDM
    n = 1:length(x_OFDM);
    ch_out_OFMD = transpose(x_OFDM).*exp(1i*2*pi*fc*Ts*n);
    P_N_OFDM = (rms(ch_out_OFMD).^2)/SNR_lin; %Noise power

    %SC
    n = 1:length(x_16QAM);
    ch_out_SC = transpose(x_16QAM).*exp(1i*2*pi*fc*Ts*n);
    P_N_SC = (rms(ch_out_SC).^2)/SNR_lin; %Noise power

    C_42 = zeros(1, 10000);
    C_42_SC = zeros(1, 10000);
    
    C_20 = 0;
    C_21 = 0;
    C_20_SC = 0;
    C_21_SC = 0;
    
    P_N = P_N_SC; %Noise power

    %*******************************************************************%
    % Generate theoretical CDF 
    %*******************************************************************%
    
    %Standard deviation of noise
    sigma = sqrt(P_N); 

    %QAM
    F_0_4 = zeros(1, length(x_4QAM));
    F_0_16 = zeros(1, length(x_16QAM));   

    z_CDF_4 = abs(x_4QAM);
    z_CDF_4 = sort(z_CDF_4);

    z_CDF_16 = abs(x_16QAM);
    z_CDF_16 = sort(z_CDF_16);

    for k = 1:length(F_0_4)   

        %QAM
        for l = 1:M_4QAM
            F_0_4(k) = F_0_4(k) + marcumq(sqrt(2)*abs(x_4(l))/sigma, sqrt(2)*z_CDF_16(k)/sigma);
        end
        F_0_4(k) = 1 - 1/M_4QAM * F_0_4(k);

        for l = 1:M_16QAM
            F_0_16(k) = F_0_16(k) + marcumq(sqrt(2)*abs(x_16(l))/sigma, sqrt(2)*z_CDF_16(k)/sigma);
        end
        F_0_16(k) = 1 - 1/M_16QAM * F_0_16(k);
    end

    D = F_0_4 - F_0_16;


    %Tolerance on finding indexes
    tol = 1e-9;
    
    [a, t_p1] = min(abs(D-max(D)));
    [b, t_p2] = min(abs(D-min(D)));
    
    %In case find() returns a vector, we only take the first element
    t_p1 = t_p1(1);
    t_p2 = t_p2(1);
    
    %2PAM
    for iter = 1:1000
        w = sqrt(P_N/2)* (randn(size(x_2PAM_large)) + 1j* randn(size(x_2PAM_large))); 
        y_tilde = x_2PAM_large + w;
        
        y_mf = conv(y_tilde, g);
        
        n = 0:numberSyms-1;
        z_2PAM = real(y_mf(T_offset + (n*4))) + 1i*imag(y_mf(T_offset + (n*4)));

        %Normalization of the constellation
        z_2PAM = z_2PAM/rms(z_2PAM);
        
        %Calculate ECDF
        z_2PAM = abs(z_2PAM);
        [f_2PAM,z_CDF_2PAM] = ecdf(z_2PAM);
    end
    
    
    %4PAM
    for iter = 1:1000
        w = sqrt(P_N/2)* (randn(size(x_4PAM_large)) + 1j* randn(size(x_4PAM_large))); 
        y_tilde = x_4PAM_large + w;
        
        y_mf = conv(y_tilde, g);
        
        n = 0:numberSyms-1;
        z_4PAM = real(y_mf(T_offset + (n*4))) + 1i*imag(y_mf(T_offset + (n*4)));

        %Normalization of the constellation
        z_4PAM = z_4PAM/rms(z_4PAM);
        
        %Calculate ECDF
        z_4PAM = abs(z_4PAM);
        [f_4PAM,z_CDF_4PAM] = ecdf(z_4PAM);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%
    %%% TEST POINTS for PAM
    ref1 = 0.5 + SNR/100;
    ref2 = 1.5 - SNR/60;
    %Index of test points
    [~, t_p1_4pam] = min(abs(z_CDF_4PAM-ref1));  
    [~, t_p2_4pam] = min(abs(z_CDF_4PAM-ref2)); 
    
    [~, t_p1_2pam] = min(abs(z_CDF_2PAM-ref1));  
    [~, t_p2_2pam] = min(abs(z_CDF_2PAM-ref2));
    
    %In case find() returns a vector, we only take the first element
    t_p1_4pam = t_p1_4pam(1);
    t_p2_4pam = t_p2_4pam(1);
    
%**************************************************************************
% CALCULATE C_42
%**************************************************************************
    for iter=1:10000

            %Complex gaussian noise
            w = sqrt(P_N_OFDM/2)* (randn(size(ch_out_OFMD)) + 1j* randn(size(ch_out_OFMD))); 

            %Received signal in the presence of AWGN
            y = ch_out_OFMD + w;

            C_20 = 1/Nm * sum(y(1:Nm).^2);
            C_21 = 1/Nm * sum(abs(y(1:Nm)).^2);
            
            C_42(iter) = 1/Nm * sum(abs(y(1:Nm)).^4) - abs(C_20).^2 - 2*(C_21).^2;
           
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %False alarm
            w = sqrt(P_N_SC/2)* (randn(size(ch_out_SC)) + 1j* randn(size(ch_out_SC))); 
            y_SC = ch_out_SC + w;
            
            C_20_SC = 1/Nm * sum(y_SC(1:Nm).^2);
            C_21_SC = 1/Nm * sum(abs(y_SC(1:Nm)).^2);

            C_42_SC(iter) = 1/Nm * sum(abs(y_SC(1:Nm)).^4) - abs(C_20_SC).^2 - 2*(C_21_SC).^2;
            
    end

    %Threshold
    gamma= (max(C_42_SC) + min(C_42))/2;

    
    for i=1:5
        count_d = 0;
        count_fa = 0;
        count_QAM = 0;
        count_PAM = 0;
    
    %**************************************************************************
    % SENDING 10000 SIGNALS
    % CHANGE THE SIGNAL HERE (in the two lines below and make sure to update 
    % the number of sample Nm is needed)
    if(i==1) 
        N = length(x_OFDM);
        n2 = 1:N;
        signal = x_OFDM;
    elseif(i==2)
        N = length(x_4QAM);
        n2 = 1:N;
        signal = x_4QAM;
    elseif(i==3)
        N = length(x_16QAM);
        n2 = 1:N;
        signal = x_16QAM;
    elseif(i==4)
        N = length(x_2PAM);
        n2 = 1:N;
        signal = x_2PAM;
    else
        N = length(x_4PAM);
        n2 = 1:N;
        signal = x_4PAM;
    end
    
    Nm = 8000; %number of samples
    P_N = (rms(signal)^2)/SNR_lin; %Noise power
    C_42_b = zeros(1, 10000);
    C_20_b = 0;
    C_21_b = 0;

    %**************************************************************************
        for iter=1:10000

            w2 = sqrt(P_N/2)* (randn(size(signal)) + 1j* randn(size(signal))); 
            y2 = signal.*transpose(exp(1i*2*pi*fc*Ts*n2)) + w2;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Calculate cumulants of signal
            C_20_b = 1/Nm * sum(y2(1:Nm).^2);
            C_21_b = 1/Nm * sum(abs(y2(1:Nm)).^2);

            C_42_b(iter) = 1/Nm * sum(abs(y2(1:Nm)).^4) - abs(C_20_b).^2 - 2*(C_21_b).^2;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %Multi-carrier
            if(C_42_b(iter) > gamma) 
                P_OFDM(i,q) = P_OFDM(i,q) + 1;
            else 
             %*****************************************************************
             %%% Block 2: MTC

                R_1 = 1/N * sum(y2.*conj(y2).* transpose(exp(-1i*2*pi*n2*Ts/T_SC)));
                R_2 = 1/N * sum(y2.^2 .* transpose(exp(-1i*2*pi*n2*Ts*(2*fc-1/T_SC))));
                R_3 = 1/N * sum(y2.^2 .* transpose(exp(-1i*2*pi*n2*Ts*(2*fc-1/(2*T_SC)))));
                R_4 = 1/N * sum(y2.^2 .* transpose(exp(-1i*2*pi*n2*2*fc*Ts)));
                R_5 = 1/N * sum(y2.^2 .* transpose(exp(-1i*2*pi*n2*Ts*(2*fc+1/(2*T_SC)))));
                R_6 = 1/N * sum(y2.^2 .* transpose(exp(-1i*2*pi*n2*Ts*(2*fc+1/T_SC))));

                F = [abs(R_1), abs(R_2), abs(R_3), abs(R_4), abs(R_5), abs(R_6)];

                F = F/norm(F); %we normalize F to compare it to V

                %If QAM
                if(norm(F-V_QAM)^2 < norm(F-V_PAM)^2)
                    count_QAM = count_QAM + 1;

                    %************************************************************************
                    %%% MF: Matched filter and symbol extractor

                    y_tilde = y2.*transpose(exp(-1i*2*pi*fc*Ts*n2));
                    %y_tilde = (transpose(x_4QAM).*exp(1i*2*pi*fc*Ts*n2) + w2).*exp(-1i*2*pi*fc*Ts*n2);
                    %y_tilde = x_16QAM + transpose(w2);
                    y_mf = conv(y_tilde, q);

                    if(length(y_mf)<8030)
                        n3 = 0:1990; %case of OFDM that has 8000 samples
                        z = real(y_mf(T_offset + (n3*4))) + 1i*imag(y_mf(T_offset + (n3*4)));
                    else
                        n3 = 0:numberSyms-1;
                        z = real(y_mf(T_offset + (n3*4))) + 1i*imag(y_mf(T_offset + (n3*4)));
                    end

                    %Normalization of the constellation
                    z = z/rms(z);

                    %******************************************%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% MLC: Modulation level classification %%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %Calculate ECDF
                    z = abs(z);
                    [f,x] = ecdf(z);

                    tol = 1e-2;

                    [a, k1_ecdf] = min(abs(x-z_CDF_16(t_p1)));
                    [b, k2_ecdf] = min(abs(x-z_CDF_16(t_p2)));

                    %In case find() returns vector, we only take 1st element
                    k1_ecdf = k1_ecdf(1);
                    k2_ecdf = k2_ecdf(1);

                    %4-QAM
                    D_4q = [ abs(F_0_4(t_p1) - f(k1_ecdf)), abs(F_0_4(t_p2) - f(k2_ecdf))];
                    %16-QAM
                    D_16q = [ abs(F_0_16(t_p1) - f(k1_ecdf)), abs(F_0_16(t_p2) - f(k2_ecdf))];

                    V_4q = norm(D_4q(1) + D_4q(2));
                    V_16q = norm(D_16q(1) + D_16q(2));

                    tolerance = 0.03;

                    if(max(D_4q) <= max(D_16q) || V_4q <= V_16q)
                        p_4QAM(i,q) = p_4QAM(i,q) + 1;
                    %if the CDFs are nearly the same, it's 50% chance to be right so we only add 1/2    
                    elseif(abs((max(D_4q) - max(D_16q))) < tolerance)
                        p_4QAM(i,q) = p_4QAM(i,q) + 1/2; 
                        p_16QAM(i,q) = p_16QAM(i,q) + 1/2;
                    else
                        p_16QAM(i,q) = p_16QAM(i,q) + 1;
                    end

                %If PAM    
                else
                    count_PAM = count_PAM + 1;

                    y_tilde2 = y2.*transpose(exp(-1i*2*pi*fc*Ts*n2));
                    y_mf2 = conv(y_tilde2, q);

                    if(length(y_mf2)<8030)
                        n4 = 0:1990; %case of OFDM that has 8000 samples
                        z2 = real(y_mf2(T_offset + (n4*4))) + 1i*imag(y_mf2(T_offset + (n4*4)));
                    else
                        n4 = 0:numberSyms-1;
                        z2 = real(y_mf2(T_offset + (n4*4))) + 1i*imag(y_mf2(T_offset + (n4*4)));
                    end

                    %Normalization of the constellation
                    z2 = z2/rms(z2);

                    %******************************************%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% MLC: Modulation level classification %%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %Calculate ECDF
                    z2 = abs(z2);
                    [f2,x2] = ecdf(z2);

                    tol = 1e-2;

                    [a, k1_ecdf2] = min(abs(x2-z_CDF_4PAM(t_p1_4pam)));
                    [b, k2_ecdf2] = min(abs(x2-z_CDF_4PAM(t_p2_4pam)));

                    %In case find() returns vector, we only take 1st element
                    k1_ecdf2 = k1_ecdf2(1);
                    k2_ecdf2 = k2_ecdf2(1);

                    %2PAM
                    D_2p2 = [ abs(f_2PAM(t_p1_2pam) - f2(k1_ecdf2)), abs(f_2PAM(t_p2_2pam) - f2(k2_ecdf2))];
                    %4PAM
                    D_4p2 = [ abs(f_4PAM(t_p1_4pam) - f2(k1_ecdf2)), abs(f_4PAM(t_p2_4pam) - f2(k2_ecdf2))];

                    V_2p2 = norm(D_2p2(1) + D_2p2(2));
                    V_4p2 = norm(D_4p2(1) + D_4p2(2));

                    tolerance = 0.02;

                    %if the CDFs are nearly the same, it's 50% chance to be right so we only add 1/2
                    if(abs((max(D_2p2) - max(D_4p2))) < tolerance)
                        p_4PAM(i,q) = p_4PAM(i,q) + 1/2; 
                        p_2PAM(i,q) = p_2PAM(i,q) + 1/2;                   
                    elseif(max(D_4p2) <= max(D_2p2) || V_4p2 <= V_2p2)
                        p_4PAM(i,q) = p_4PAM(i,q) + 1;
                    else
                        p_2PAM(i,q) = p_2PAM(i,q) + 1;
                    end

                end
            end

            %Sent SC but detected OFDM
            if(C_42_SC(iter) > gamma) 
                count_fa = count_fa + 1;
            end
        end
    
  
    P_OFDM(i,q) = P_OFDM(i,q)/10000;
    P_fa(i,q) = count_fa/10000;
    
    p_QAM(i,q) = count_QAM/10000;
    p_PAM(i,q) = count_PAM/10000;
    
    p_4QAM(i,q) = p_4QAM(i,q)/10000;
    p_16QAM(i,q) = p_16QAM(i,q)/10000;
    
    p_2PAM(i,q) = p_2PAM(i,q)/10000;
    p_4PAM(i,q) = p_4PAM(i,q)/10000;
    
    end
    
    p_c(q) = (P_OFDM(1,q) + p_4QAM(2,q) + p_16QAM(3,q) + p_2PAM(4,q)+ p_4PAM(5,q))/5;
    
    q = q + 1;
end



figure
plot(SNR_vect, p_c, 'LineWidth', 1.5);
title('Probability of correct classification P_c vs SNR for the entire classifier')
xlabel('SNR')
ylabel('P_c')
grid on;



