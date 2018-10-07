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
%x_wCP = transpose(x_wCP);
%Upsample by 4
x_up = upsample(x_wCP, 4);

%Low Pass Filter
%w_c = f_cutoff/(fs/2);
[b, a] = butter(8, 0.25);
x_filter = filter(b, a, x_up);
x_filter = reshape(x_filter, 1, []);
x_OFDM = x_filter./(rms(x_filter));     %Normalize power of transmitted signal


figure; 
periodogram(x_OFDM); 
title('PSD of x(n) for OFDM')
grid on;

%%
%%%%%%%%%%%%%%%%%%%%%%%%
% 2. 4-QAM
M = 4;
%Generate 2000 random symbols
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

figure; 
periodogram(x_4QAM); 
title('PSD of x(n) for 4-QAM')
grid on;

scatterplot(syms) 
title('Constellation diagram for 4-QAM')

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

figure; 
periodogram(x_16QAM); 
title('PSD of x(n) for 16-QAM')
grid on;

scatterplot(syms) 
title('Constellation diagram for 16-QAM')

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
x_conv = conv(syms_up, g);
%Normalize power of transmitted signal
x_2PAM = x_conv/(rms(x_conv)); 

figure; 
periodogram(x_2PAM); 
title('PSD of x(n) for 2-PAM')
grid on;

scatterplot(syms2p) 
title('Constellation diagram for 2-PAM')

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
x_conv = conv(syms_up, g);
%Normalize power of transmitted signal
x_4PAM = x_conv/(rms(x_conv)); 


scatterplot(syms4p) 
title('Constellation diagram for 4-PAM')

figure; 
periodogram(x_4PAM); 
title('PSD of x(n) for 4-PAM')
grid on;


%% Block 1: Mc-SC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of samples used for distinguishing MC and SC signals

%2000 symbols upsampled by 4 = 8000 samples
Nm = 8000;    

q = 1;          %Index in for loop SNR
probaPts = 100;  %Number of points used for the P_fa vs P_d plot

SNR_vect = 0:10:20;

%Preallocation of probability vectors
P_d = zeros(length(SNR_vect), probaPts);
P_fa = zeros(length(SNR_vect), probaPts);

%figure; hold on
for SNR = SNR_vect  
    
    SNR_lin = 10^(SNR/10); %SNR linear value

    %OFDM
    n = 1:length(x_OFDM);
    ch_out_OFMD = x_OFDM.*exp(1i*2*pi*fc*Ts*n);
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
    
    for iter=1:10000

            %Complex gaussian noise
            w = sqrt(1/(2*SNR_lin))* (randn(size(ch_out_OFMD)) + 1j* randn(size(ch_out_OFMD))); 

            %Received signal in the presence of AWGN
            y = ch_out_OFMD + w;
            
            %{
            C_20 = 0;
            C_21 = 0;
            
            for m = 1:Nm
                C_20 =  C_20 + y(m).^2;
                C_21 = C_21 + abs(y(m)).^2;
                C_42(iter) = C_42(iter) + abs(y(m)).^4;
            end
            
            C_20 = C_20/Nm;
            C_21 = C_21/Nm;
            C_42(iter)= C_42(iter)/Nm;
            C_42(iter) = C_42(iter) - abs(C_20).^2 - 2*(C_21).^2;
            %}
            
            C_20 = 1/Nm * sum(y(1:Nm).^2);
            C_21 = 1/Nm * sum(abs(y(1:Nm)).^2);
            
            C_42(iter) = 1/Nm * sum(abs(y(1:Nm)).^4) - abs(C_20).^2 - 2*(C_21).^2;
            

            %False alarm
            w = sqrt(1/(2*SNR_lin))* (randn(size(ch_out_SC)) + 1j* randn(size(ch_out_SC))); 
            y_SC = ch_out_SC + w;

            C_20_SC = 1/Nm * sum(y_SC(1:Nm).^2);
            C_21_SC = 1/Nm * sum(abs(y_SC(1:Nm)).^2);

            C_42_SC(iter) = 1/Nm * sum(abs(y_SC(1:Nm)).^4) - abs(C_20_SC).^2 - 2*(C_21_SC).^2;
            
            
            %{
            C_20_SC = 0;
            C_21_SC = 0;
            
            for m = 1:Nm
                C_20_SC =  C_20_SC + y_SC(m).^2;
                C_21_SC = C_21_SC + abs(y_SC(m)).^2;
                C_42_SC(iter) = C_42_SC(iter) + 1/Nm * abs(y_SC(m)).^4;
            end

            C_20_SC = C_20_SC/Nm;
            C_21_SC = C_21_SC/Nm;
            C_42_SC(iter) = C_42_SC(iter) - abs(C_20_SC).^2 - 2*(C_21_SC).^2;
            %}
    end
    
    
    l = 1; %index of P_d and P_fa
    for gamma=linspace(min(C_42_SC),max(C_42),probaPts)
        for i=1:5
            %**************************************************************************
            % SENDING 10000 SIGNALS
            % CHANGE THE SIGNAL HERE (in the two lines below and make sure to update 
            % the number of sample Nm is needed)
            if(i==1) 
                n2 = 1:length(x_OFDM);
                signal = x_OFDM;
            elseif(i==2)
                n2 = 1:length(x_4QAM);
                signal = x_4QAM;
            elseif(i==3)
                n2 = 1:length(x_16QAM);
                signal = x_16QAM;
            elseif(i==4)
                n2 = 1:length(x_2PAM);
                signal = x_2PAM;
            else
                n2 = 1:length(x_4PAM);
                signal = x_4PAM;
            end
            count_d = 0;
            count_fa = 0;
            
            
            for iter=1:10000
                if(C_42(iter) > gamma) 
                    count_d = count_d + 1;
                end     

                if(C_42_SC(iter) > gamma) 
                       count_fa = count_fa + 1;
                end 
            end

            P_d(q, l) = count_d/10000;
            P_fa(q, l) = count_fa/10000;
            
        end
        
        l = l + 1;

    end

    
    
    figure 
    h1 = histogram(C_42);
    h1.FaceColor = 'c';
    h1.EdgeColor = 'b';
    hold on;
    h2 = histogram(C_42_SC);
    h2.FaceColor = [1,0.6,0.1];
    h2.EdgeColor = 'r';
    title(['C_{42} distribution with SNR = ', num2str(SNR), ' dB'])
    xlabel('C_{42}')
    ylabel('Frequency')
    legend('OFDM', '16-QAM')
    
    figure
    plot(P_fa(q,:), P_d(q,:), 'LineWidth', 1.5); hold on; 
    grid on;
    title('Probability of detection (P_d) vs probability of false alarm for OFDM (P_{fa})')
    xlabel('P_{FA}')
    ylabel('P_d')
    
    q = q + 1;

end


figure
plot(P_fa(1,:), P_d(1,:),'y*', P_fa(2,:), P_d(2,:), 'b--o', P_fa(3,:), P_d(3,:), 'rdiamond', 'LineWidth', 1.5);
grid on;
title('Probability of detection (P_d) vs probability of false alarm for OFDM (P_{fa})')
legend('SNR = 0 dB', 'SNR = 10 dB', 'SNR = 20 dB')
xlabel('P_{FA}')
ylabel('P_d')


%{
grid on;
title('Probability of detection (P_d) vs probability of false alarm for OFDM (P_{fa})')
legend('SNR = 0 dB', 'SNR = 10 dB', 'SNR = 20 dB')
xlabel('P_{FA}')
ylabel('P_d')
hold off;
%}