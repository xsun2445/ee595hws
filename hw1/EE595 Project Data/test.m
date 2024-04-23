% test.m

rf_freq = 915e6;    % carrier frequency 470 MHz, used to correct clock drift
sf = 8;             % spreading factor SF7
bw = 125e3;         % bandwidth 125 kHz
fs = 250e3;           % sampling rate 1 MHz

phy = LoRaEncoder(rf_freq, sf, bw, fs);
phy.has_header = 0;         % explicit header mode
phy.cr = 2;                 % code rate = 4/8 (1:4/5 2:4/6 3:4/7 4:4/8)
phy.crc = 0;                % enable payload CRC checksum
phy.preamble_len = 8;       % preamble: 8 basic upchirps

% Encode payload [1 2 3 4 5]
disp('payload:')
disp([[1 2 3 4]';double('PING')';(8:15)']')


symbols = phy.encode([[1 2 3 4]';double('PING')';(8:15)']);
fprintf("[encode] symbols:\n");
disp(symbols);
length(symbols)

% Baseband Modulation
sig = phy.modulate(symbols);
disp('sig size ')
size(sig)

upchirp1 = phy.chirp(true,sf,bw,fs,230,0,0);
downchirp1 = phy.chirp(false,sf,bw,fs,0,0,0);

length(upchirp1)
length(downchirp1)


% figure()
% plot(real(phy.chirp(true,sf,bw,fs,0,0,0)))
% hold on
% plot(imag(phy.chirp(true,sf,bw,fs,0,0,0)))

figure()
plot(abs(fft(upchirp1.*downchirp1)))

figure()
plot(abs(fftshift(fft(upchirp1.*downchirp1))))
