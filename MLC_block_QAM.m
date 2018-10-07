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
%Generate random symbols
b_in = randi([0 1],log2(M)*numberSyms,1);
%QAM modulation
syms = qammod(b_in, M,'InputType','bit','UnitAveragePower',true);
%64-pt IFFT
x_ifft = ifft(syms, N);
%P to S
x_serial = transpose(x_ifft);
%Add cyclic prefix
x_wCP = [x_serial(N-cp+1:N), x_serial];
%Upsample by 4
x_up = upsample(x_wCP, 4);

%Low Pass Filter
%w_c = f_cutoff/(fs/2);
[b, a] = butter(8, 0.25);
x_filter = filter(b, a, x_up);
x_OFDM = x_filter/(rms(x_filter));     %Normalize power of transmitted signal


figure; 
periodogram(x_OFDM); 
title('PSD of x(n) for OFDM')
grid on;

scatterplot(syms) 
title('Constellation diagram for OFDM')


%%%%%%%%%%%%%%%%%%%%%%%%
% 2. 4-QAM
M = 4;
%Generate random symbols
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
%Generate random symbols
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

%************************************************************************
%% MF: Matched filter and symbol extractor

close all;

%Pulse shaping filter
g = rcosdesign(0.5,8,4, 'normal');

% sample corresponding to the first local max
T_offset = length(g);

%Number of constellation points
M_4QAM = 4;
x_4 = 1/sqrt(2) * [1 + 1i, 1 - 1i, -1 + 1i, -1 - 1i];

M_16QAM = 16;
x_16 = 1/sqrt(10) * [1 + 1i, 1 - 1i, -1 + 1i, -1 - 1i, 3 + 3i, 3 - 3i, ...
    -3 + 3i, -3 - 3i, 3 + 1i, 3 - 1i, -3 + 1i, -3 - 1i, 1 + 3i, 1 - 3i, -1 + 3i, -1 - 3i];

%Counter in for loop
q = 1;

SNR_vect = -4:2:20;
%SNR_vect = 0;

%Count the proba to detect or false alarm
count_16QAM = zeros(1, length(SNR_vect));

p_detect = zeros(1, length(SNR_vect));

for SNR = SNR_vect
    SNR_lin = 10^(SNR/10); %SNR linear value
    P_N = (rms(x_4QAM)^2)/SNR_lin; %Noise power

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
    
    z_CDF_4 = transpose(z_CDF_4);
    z_CDF_16 = transpose(z_CDF_16);
    

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

    %{
    figure
    plot(z_CDF_16,F_0_4, z_CDF_16,F_0_16, 'LineWidth', 1.5)
    title(['CDF for QAM for SNR = ', num2str(SNR), ' dB'])
    legend('4-QAM', '16-QAM')
    grid on;
    xlabel('Amplitude')
    ylabel('CDF')
    %}
    
    D = F_0_4 - F_0_16;
    %{
    figure
    plot(z_CDF_16, D , 'LineWidth', 1.5)
    title('F_0_4 - F_0_16')
    grid on;
    xlabel('Amplitude')
    ylabel('CDF')
    %}

    

    %Index of test points
    %t_p1 = find( abs(D-max(D))<tol );
    %t_p2 = find( abs(D-min(D))<tol );
    
    [a, t_p1] = min(abs(D-max(D)));
    [b, t_p2] = min(abs(D-min(D)));
 
    %*******************************************************************%
    n2 = 1:length(x_16QAM);
    for iter = 1:1000
        w2 = sqrt(1/(2*SNR_lin))* (randn(size(x_16QAM)) + 1j* randn(size(x_16QAM))); 

        y_tilde2 = (x_16QAM.*transpose(exp(1i*2*pi*fc*Ts*n2)) + w2).*transpose(exp(-1i*2*pi*fc*Ts*n2));
        
        y_mf2 = conv(y_tilde2, q);
        
        
        n = 0:numberSyms-1;
        z_16 = real(y_mf2(T_offset + (n*4))) + 1i*imag(y_mf2(T_offset + (n*4)));

        %Normalization of the constellation
        z_16 = z_16/rms(z_16);
        
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
        z_16 = abs(z_16);
        [f2,x2] = ecdf(z_16);
        
        tol = 1e-2;
        
        [a, k1_ecdf_16] = min(abs(x2-z_CDF_16(t_p1)));
        [b, k2_ecdf_16] = min(abs(x2-z_CDF_16(t_p2)));
        
        %In case find() returns vector, we only take 1st element
        k1_ecdf_16 = k1_ecdf_16(1);
        k2_ecdf_16 = k2_ecdf_16(1);

        %4-QAM
        D_4q_16 = [ abs(F_0_4(t_p1) - f2(k1_ecdf_16)), abs(F_0_4(t_p2) - f2(k2_ecdf_16))];
        %16-QAM
        D_16q_16 = [ abs(F_0_16(t_p1) - f2(k1_ecdf_16)), abs(F_0_16(t_p2) - f2(k2_ecdf_16))];

        V_4q_16 = norm(D_4q_16(1) + D_4q_16(2));
        V_16q_16 = norm(D_16q_16(1) + D_16q_16(2));
        
        tolerance = 0.02;
        
        %if the CDFs are nearly the same, it's 50% chance to be right so we only add 1/2
        if(abs((max(D_4q_16) - max(D_16q_16))) < tolerance)
            count_16QAM(q) = count_16QAM(q) + 1/2;  
        elseif(max(D_16q_16) <= max(D_4q_16) || V_16q_16 <= V_4q)
            count_16QAM(q) = count_16QAM(q) + 1;
        end
         
    end
    
    
    
    figure
    plot(z_CDF_16,F_0_4, z_CDF_16,F_0_16, x2, f2, 'LineWidth', 1.5)
    title(['CDF for QAM for SNR = ', num2str(SNR), ' dB'])
    legend('4-QAM', '16-QAM', 'ECDF')
    xlabel('Amplitude')
    ylabel('CDF')
    
    
    count_16QAM(q) = count_16QAM(q)/1000;
    
    p_detect(q) = count_16QAM(q);
    
    q = q + 1;
end

%{
figure
plot(SNR_vect, count_16QAM, 'LineWidth', 1.5);
title('Average probability of correct level classification vs SNR (16-QAM)')
xlabel('SNR')
ylabel('P_d')
grid on;


figure
plot(SNR_vect, count_4QAM, 'LineWidth', 1.5);
title('Average probability of correct level classification vs SNR (4-QAM)')
xlabel('SNR')
ylabel('P_d')
grid on;
%}

figure
plot(SNR_vect, p_detect, 'LineWidth', 1.5);
title('Average probability of correct level classification vs SNR for QAM')
xlabel('SNR')
ylabel('P_d')
grid on;

