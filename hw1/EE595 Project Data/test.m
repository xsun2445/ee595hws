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


fprintf('ldro enabled? %d\n', phy.ldr)


original_data = [[1 2 3 4]';double('PING')';(8:15)'];
original_data = [[1 2 3 4 5]'];
% original_data = [1:500]';

% Encode payload [1 2 3 4 5]
disp('payload:')
disp(original_data)




symbols = phy.encode(original_data);
fprintf("[encode] symbols:\n");
disp(symbols);
disp(class(symbols))
length(symbols)
disp(dec2bin(symbols))


fprintf('\ngray coding in the decoding process\n')
disp(phy.gray_coding(symbols))


% fileName = 'dataSF8CR8packet1.bin';
% raw = phy.read(fileName);

% figure()
% plot(real(raw))

% figure()
% [s,f,t] = stft(raw,fs);

% nexttile
% mesh(t,f,abs(s).^2)
% title("stft")
% view(2), axis tight

% uc = LoRaEncoder.chirp(false, sf, bw, fs, 0, 0, 0);
% preamble = repmat(uc, 8, 1);

% figure()
% [s,f,t] = stft(xcorr(raw, preamble),fs);
% nexttile
% mesh(t,f,abs(s).^2)
% title("stft")
% view(2), axis tight

% size(t)
% size(f)
% size(s)




% % Baseband Modulation
% sig = phy.modulate(symbols);
% disp('sig size ')
% size(sig)

% upchirp1 = phy.chirp(true,sf,bw,fs,230,0,0);
% downchirp1 = phy.chirp(false,sf,bw,fs,0,0,0);

% length(upchirp1)
% length(downchirp1)


% % figure()
% % plot(real(phy.chirp(true,sf,bw,fs,0,0,0)))
% % hold on
% % plot(imag(phy.chirp(true,sf,bw,fs,0,0,0)))

% figure()
% plot(abs(fft(upchirp1.*downchirp1)))

% figure()
% plot(abs(fftshift(fft(upchirp1.*downchirp1))))
