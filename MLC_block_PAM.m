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


%%%%%%%%%%%%%%%%%%%%%%%%
% 2-PAM
M = 2;
%Generate random symbols
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

figure; 
periodogram(x_2PAM); 
title('PSD of x(n) for 2-PAM')
grid on;

scatterplot(syms2p) 
title('Constellation diagram for 2-PAM')

%%%%%%%%%%%%%%%%%%%%%%%%
% 4-PAM
M = 4;
%Generate random symbols
b_in = randi([0 M-1],log2(M)*numberSyms,1);
%QAM modulation
syms4p = pammod(b_in, M);
%Upsample by 4
syms_up = upsample(syms4p, 4);
%Pulse shaping filter
g = rcosdesign(0.5,8,4, 'normal');
%Symbols after upsampling and pulse shaping
x_conv = conv(syms_up, g);
%Normalize power of transmitted signal
x_4PAM = x_conv/(rms(x_conv)); 


scatterplot(syms4p) 
title('Constellation diagram for 4-PAM')

figure; 
periodogram(x_4PAM); 
title('PSD of x(n) for 4-PAM')
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%
% 2-PAM with larger number of symbols to calculate the ECDF (that
% approximates the theoretical CDF)
M = 2;
%Generate random symbols
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
%Generate random symbols
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


%************************************************************************
%% MF: Matched filter and symbol extractor

close all;

%Pulse shaping filter
g = rcosdesign(0.5,8,4, 'normal');

% sample corresponding to the first local max
T_offset = length(g);


%Counter in for loop
g = 1;

SNR_vect = -4:2:20;
%SNR_vect = 0;

%Count the proba to detect or false alarm
count_2PAM = zeros(1, length(SNR_vect));
count_4PAM = zeros(1, length(SNR_vect));
p_detect = zeros(1, length(SNR_vect));

for SNR = SNR_vect
    SNR_lin = 10^(SNR/10); %SNR linear value
    P_N = (rms(x_4QAM)^2)/SNR_lin; %Noise power

    %*******************************************************************%
    % Generate theoretical CDF 
    %*******************************************************************%
    
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
     
    %*******************************************************************%
    %% 2PAM
    n2 = 1:length(x_2PAM);
    for iter = 1:1000
        w2 = sqrt(P_N/2)* (randn(size(x_2PAM)) + 1j* randn(size(x_2PAM))); 

        y_tilde = (x_2PAM.*transpose(exp(1i*2*pi*fc*Ts*n2)) + w2).*transpose(exp(-1i*2*pi*fc*Ts*n2));
        
        y_mf = conv(y_tilde, g);
        
        
        n = 0:numberSyms-1;
        z = real(y_mf(T_offset + (n*4))) + 1i*imag(y_mf(T_offset + (n*4)));

        %Normalization of the constellation
        z = z/rms(z);
        
        %{
        scatterplot(z)
        title('MF - Constellation diagram of received signal')
        
        figure
        plot(real(z), imag(z), '.', real(syms), imag(syms), 'o', 'MarkerSize',7, 'MarkerFaceColor', 'r')
        title('MF - Constellation diagram of received signal')
        legend('z', 'Ideal')
        %}
        
        %******************************************%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% MLC: Modulation level classification %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Calculate ECDF
        z = abs(z);
        [f,x] = ecdf(z);
        
        tol = 1e-2;
        
        [a, k1_ecdf] = min(abs(x-z_CDF_4PAM(t_p1_4pam)));
        [b, k2_ecdf] = min(abs(x-z_CDF_4PAM(t_p2_4pam)));
        
        %In case find() returns vector, we only take 1st element
        k1_ecdf = k1_ecdf(1);
        k2_ecdf = k2_ecdf(1);

        %{
        %4-QAM
        D_4q = [ abs(F_0_4(k1_4q) - f(k1_ecdf)), abs(F_0_4(k2_4q) - f(k2_ecdf))];
        %16-QAM
        D_16q = [ abs(F_0_16(k1_16q) - f(k1_ecdf)), abs(F_0_16(k2_16q) - f(k2_ecdf))];
        %}


        
        %2PAM
        D_2p_2 = [ abs(f_2PAM(t_p1_2pam) - f(k1_ecdf)), abs(f_2PAM(t_p2_2pam) - f(k2_ecdf))];
        %4PAM
        D_4p_2 = [ abs(f_4PAM(t_p1_4pam) - f(k1_ecdf)), abs(f_4PAM(t_p2_4pam) - f(k2_ecdf))];

        V_2p_2 = norm(D_2p_2(1) + D_2p_2(2));
        V_4p_2 = norm(D_4p_2(1) + D_4p_2(2));
        
        tolerance = 0.02;
        
        %if the CDFs are nearly the same, it's 50% chance to be right so we only add 1/2    
        if(abs((max(D_2p_2) - max(D_4p_2))) < tolerance) 
            count_2PAM(g) = count_2PAM(g) + 1/2;
        
        elseif(max(D_2p_2) <= max(D_4p_2) && V_2p_2 <= V_4p_2)
            count_2PAM(g) = count_2PAM(g) + 1;
        end
         
    end
    
    %% 4PAM
    n3 = 1:length(x_4PAM);
    for iter = 1:1000
        w3 = sqrt(P_N/2)* (randn(size(x_4PAM)) + 1j* randn(size(x_4PAM))); 

        y_tilde_4p = (x_4PAM.*transpose(exp(1i*2*pi*fc*Ts*n3)) + w3).*transpose(exp(-1i*2*pi*fc*Ts*n3));
        y_mf_4p = conv(y_tilde_4p, g);
        
        
        n = 0:numberSyms-1;
        z_4p = real(y_mf_4p(T_offset + (n*4))) + 1i*imag(y_mf_4p(T_offset + (n*4)));

        %Normalization of the constellation
        z_4p = z_4p/rms(z_4p);
        
        %Calculate ECDF
        z_4p = abs(z_4p);
        [f_4p,x_4p] = ecdf(z_4p);
        
        tol = 1e-2;
        
        [a, k1_ecdf_4p] = min(abs(x_4p-z_CDF_4PAM(t_p1_4pam)));
        [b, k2_ecdf_4p] = min(abs(x_4p-z_CDF_4PAM(t_p2_4pam)));
        
        %In case find() returns vector, we only take 1st element
        k1_ecdf_4p = k1_ecdf_4p(1);
        k2_ecdf_4p = k2_ecdf_4p(1);

        %{
        %4-QAM
        D_4q = [ abs(F_0_4(k1_4q) - f(k1_ecdf)), abs(F_0_4(k2_4q) - f(k2_ecdf))];
        %16-QAM
        D_16q = [ abs(F_0_16(k1_16q) - f(k1_ecdf)), abs(F_0_16(k2_16q) - f(k2_ecdf))];
        %}


        
        %2PAM
        D_2p = [ abs(f_2PAM(t_p1_2pam) - f_4p(k1_ecdf_4p)), abs(f_2PAM(t_p2_2pam) - f_4p(k2_ecdf_4p))];
        %4PAM
        D_4p = [ abs(f_4PAM(t_p1_4pam) - f_4p(k1_ecdf_4p)), abs(f_4PAM(t_p2_4pam) - f_4p(k2_ecdf_4p))];

        V_2p = norm(D_2p(1) + D_2p(2));
        V_4p = norm(D_4p(1) + D_4p(2));
        
        tolerance = 0.02;
        
        %if the CDFs are nearly the same, it's 50% chance to be right so we only add 1/2    
        if(abs((max(D_2p) - max(D_4p))) < tolerance)
            count_4PAM(g) = count_4PAM(g) + 1/2;    
        elseif(max(D_4p) <= max(D_2p) && V_4p <= V_2p)
            count_4PAM(g) = count_4PAM(g) + 1;
        end
         
    end
    
    figure
    plot(z_CDF_2PAM,f_2PAM, z_CDF_4PAM,f_4PAM, x, f, 'LineWidth', 1.5)
    title(['CDF for PAM for SNR = ', num2str(SNR), ' dB'])
    legend('2-PAM', '4-PAM', 'ECDF')
    grid on;
    xlabel('Amplitude')
    ylabel('CDF')
    
    count_2PAM(g) = count_2PAM(g)/1000;
    count_4PAM(g) = count_4PAM(g)/1000;
    
    p_detect(g) = (count_2PAM(g) + count_4PAM(g))/2;
    
    g = g + 1;
end

figure
plot(SNR_vect, count_4PAM, 'LineWidth', 1.5);
title('Average probability of correct level classification vs SNR for 4PAM')
xlabel('SNR')
ylabel('P_d')
grid on;

figure
plot(SNR_vect, count_2PAM, 'LineWidth', 1.5);
title('Average probability of correct level classification vs SNR for 2PAM')
xlabel('SNR')
ylabel('P_d')
grid on;

figure
plot(SNR_vect, p_detect, 'LineWidth', 1.5);
title('Average probability of correct level classification vs SNR for PAM')
xlabel('SNR')
ylabel('P_d')
grid on;
