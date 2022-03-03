%% ECE 600 Final Project - RECEIVER Script
%   OFDM-based Wireless Communication Implementation on USRP Test
%   Student name: Adnan Quadri
%   ID: 5326490
%   Date: 

clear 
close all
clc

%% Recevier side

% Reading DAT file to load received signal

fileID1 = fopen('adnan_rx_ofdm_sig_nov_26_1.dat', 'rb');%fopen('ofdm_sig_q.dat', 'rb');%fopen('adnan_rx_ofdm_sig_nov_26_1.dat', 'rb');
frewind(fileID1);
fseek(fileID1, 1e6, 'bof');
ofdm_signal_float_rx = fread(fileID1, 5*2*1600, 'float');
fclose(fileID1);

rx_sig1 = transpose(ofdm_signal_float_rx(1:2:end) + 1i*ofdm_signal_float_rx(2:2:end));
rx_sig_cfo = rx_sig1;

%% Transmitter side preamble information
rx_sc = 64;
rx_n_preamb = 4;
data_pos = [2:3 5:17 19:27 39:45 47:60 62:64];
rx_preamb_pos = [2:27 39:64];
tx_ofdm_symb =20;
rx_cyc_block = 16;
pilot_pos = [4 18 46 61];

preamb_pos = [2:27 39:64]; %[7:32 34:59];

% Preamble Part

rx_preamb = [0,1,0,0,0,0,1,1,0,0,1,0,1,1,0,0,0,1,0,0,0,0,0,0,1,1,0,1,0,0,1,1,0,1,0,1,1,1,1,1,0,0,1,0,1,1,1,0,0,0,1,0]';

% BPSK MOdulation of Preamble========================================
for rx = 1:length(rx_preamb)
    if rx_preamb(rx) == 0
        rx_mod_preamb(rx) = -1;
    else
        rx_mod_preamb(rx) = 1;
    end
end

% 4 OFDM Symbols generated as Preambles=============================
for rx_pambs = 1:rx_n_preamb
    tx_preambs(:,rx_pambs) = rx_mod_preamb;
end

% Preambles inserted in the required positions for the OFDM WiFi frame
tx_preambles = zeros(rx_sc, rx_n_preamb);
tx_preambles(rx_preamb_pos,:) = tx_preambs;

preamb_ifft = ifft(tx_preambles);

preamb_stream_cp = [preamb_ifft(end - rx_cyc_block +1:end,:); preamb_ifft()];
ltfs = reshape(preamb_stream_cp, 1, []);


%% For time syncronization
cross_corr = zeros(1,length(rx_sig_cfo));

tx_ltfs = ltfs(1:320);%txSignal(1:320);

for corr = 1:length(rx_sig_cfo)-321
   cross_corr(1,corr) = tx_ltfs*rx_sig_cfo(corr:corr+(length(tx_ltfs)-1))';
end

[val1 idx1] = max(cross_corr);
[val2 idx2] = min(cross_corr);


figure(1);
plot(abs(cross_corr(:,1:length(rx_sig_cfo))),'b-');
xlabel('sample index');
ylabel('correlation');
title('Cross correlation');

%% Freq offset correction
% 
rx_signal = rx_sig_cfo(1,idx1:idx1-1+1600);

rxCFO = fft(rx_signal);
scatterplot(reshape(rxCFO,1,[]));
title('Before CFO')

% Frequency offset correcton
fr_off = rx_signal(1:160)*rx_signal(161:320)';
ang_off= angle(fr_off)/160;
fr_offset = exp(1i*ang_off*(1:length(rx_signal)));

rxSig_cfo = fr_offset.*rx_signal;

rxCFO_corr = fft(rxSig_cfo);
scatterplot(reshape(rxCFO_corr,1,[]));
title('After CFO')
%% Removing CP
rx_frame = reshape(rxSig_cfo, (rx_sc+rx_cyc_block), []);
rx_cyc_subtract = rx_frame(rx_cyc_block+1:end,:);

% Perfroming OFDM fft
rx_ofdm_fft = fft(rx_cyc_subtract);


%% Channel Estimation

rx_preamb1 = rx_ofdm_fft(rx_preamb_pos,1);

rx_preamb2 = rx_ofdm_fft(rx_preamb_pos,2);
rx_preamb3 = rx_ofdm_fft(rx_preamb_pos,3);
rx_preamb4 = rx_ofdm_fft(rx_preamb_pos,4);

rx_preambavg = rx_preamb1;%(rx_preamb1 + rx_preamb2 + rx_preamb3 + rx_preamb4)./4;

tx_preambleavg = tx_preambles(rx_preamb_pos,1);%(tx_preambles(rx_preamb_pos,1) + tx_preambles(rx_preamb_pos,2) + tx_preambles(rx_preamb_pos,3) + tx_preambles(rx_preamb_pos,4))./4;

Hhat = rx_preambavg./tx_preambleavg;

Hhat1 = repmat(Hhat,1,16);


H_mag = abs(Hhat([27:end 1:26]));

H_phase = angle(Hhat([27:end 1:26]))/pi;



figure(4)
plot(H_mag, '-.b*')
xlabel('Sub-carrier indexes');
ylabel('Magnitude');
title('Estimated channel')
axis([1 52 0 1.2*max(H_mag)]);
grid on

figure(5)
plot(H_phase, '-.b*')
xlabel('Sub-carrier indexes');
ylabel('Phase');
title('Estimated channel')
grid on


%% Channel Equalization

rx_data(rx_preamb_pos,:) = rx_ofdm_fft(rx_preamb_pos, 5:end)./Hhat1;
rx_pho = rx_data(data_pos,:);

scatterplot(reshape(rx_pho,1,[]))
title('Before phase correction');
%% Phase Correction

rx_pilots1 = rx_data(pilot_pos, 5);
rx_pilots2 = rx_data(pilot_pos, 6);

pho_est = angle(mean(rx_pilots1.*rx_pilots2));

rx_ph1 =rx_data.*exp(-1i*pho_est/(length(rx_data)));

rx_pho_corr = rx_ph1(data_pos,:);

scatterplot(reshape(rx_pho_corr,1,[]))
title('After phase correction');


rxData = reshape(rx_pho_corr,1,[]);
%% QPSK Demodulation 
demod_data = demod_qpsk(rxData);



%% BER calculation using TX side data

% Tranmitted data
data = csvread('adnan_tx_dataNov26.txt');
data;
ber2 = biterr(data, demod_data)/length(data);




