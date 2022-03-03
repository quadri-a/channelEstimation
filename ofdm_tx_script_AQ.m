%% ECE 600 Final Project - Transmitter Script
%   OFDM-based Wireless Communication Implementation on USRP Test
%   Student name: Adnan Quadri
%   ID: 5326490
%   Date: 
clear
close all
clc
%% STEP: 1 GENERATE OFDM SIGNAL
% Initializing Parameters

n_data_symb = 16;
sc = 64;
data_sc = 48;
n_pilots = 4;
n_preamb = 4;
pilot = [-1 1 -1 1];
pilots = repmat(pilot',1,n_data_symb);
nbits = 2*data_sc*n_data_symb;
cyc_block = 16;

M1 = 2;
M2 = 4;

data_pos = [2:3 5:17 19:27 39:45 47:60 62:64];
preamb_pos = [2:27 39:64]; 
pilot_pos = [4 18 46 61];

%% Preamble Part

preamb = [0,1,0,0,0,0,1,1,0,0,1,0,1,1,0,0,0,1,0,0,0,0,0,0,1,1,0,1,0,0,1,1,0,1,0,1,1,1,1,1,0,0,1,0,1,1,1,0,0,0,1,0]';

% BPSK MOdulation of Preamble========================================
for ii = 1:length(preamb)
    if preamb(ii) == 0
        mod_preamb(ii) = -1;
    else
        mod_preamb(ii) = 1;
    end
end

% 4 OFDM Symbols generated as Preambles=============================
for pambs = 1:n_preamb
    preambs(:,pambs) = mod_preamb;
end

% Preambles inserted in the required positions for the OFDM WiFi frame
preambles = zeros(sc, n_preamb);
preambles(preamb_pos,:) = preambs;

%% Data Part

data = randi([0 1], 1, nbits);
%csvwrite('adnan_tx_dataNov26.txt',data)

% QPSK modulation==================================================

d = zeros(1,length(data)/2);

for i = 1:(length(data)/2)
    in = data(2*i);
    quad = data(2*i -1);
    
    
    if (quad == 0) && (in == 0)
        d(i) = exp(j*pi/4);
    end
    if (quad == 1) && (in == 0)
        d(i) = exp(j*3*pi/4);
    end
    if (quad == 1) && (in == 1)
        d(i) = exp(j*5*pi/4);
    end
    if (quad == 0) && (in == 1)
        d(i) = exp(j*7*pi/4);
    end


end

scatterplot(d);
title('QPSK modulated Data');

modulated_data = d;
modulated_data = reshape(modulated_data, data_sc, n_data_symb);
data_stream = zeros(sc, n_data_symb);
data_stream(data_pos, :) = modulated_data;

% scatterplot(reshape(modulated_data,1,[]))
data_stream(pilot_pos,:) = pilots;

ofdm_stream = [preambles data_stream];

ofdm_stream_ifft = ifft(ofdm_stream);

ofdm_stream_cp = [ofdm_stream_ifft(end - cyc_block +1:end,:); ofdm_stream_ifft()];

plot(real(ofdm_stream_cp(:).'))


txSignal = reshape(ofdm_stream_cp, 1,[]);
figure;
plot(real(txSignal));

%% Creating DAT file for TX 


% % % fileID= fopen(['adnan_ofdm_signalnov_26_1.dat'], 'wb');
% % % 
% % % ofdm_signal_float = [real(txSignal); imag(txSignal)];
% % % 
% % % ofdm_signal_float = ofdm_signal_float(:);
% % % 
% % % fwrite(fileID, ofdm_signal_float, 'float');
% % % fclose(fileID);
