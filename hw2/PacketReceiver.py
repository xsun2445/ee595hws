from pylab import *
from rtlsdr import *
from LoRaDecoder import LoRaDecoder

sdr = RtlSdr()
sdr.sample_rate = 250e3
sdr.center_freq = 915e6
# sdr.freq_correction = 0
sdr.gain = 'auto'

# res = sdr.read_samples(512*1024)
res = sdr.read_samples(512*1024)
# print(res)
# print(len(res), type(res), type(res[0]))
# psd(res, NFFT=1024, Fs=sdr.sample_rate/1e6, Fc=sdr.center_freq/1e6)
# show()

# 
# with open('Lora_packet_sf7_cp8_0.npy', 'wb') as f:
#     np.save(f, res)

phy = LoRaDecoder(sf=7, cp=8)
# res = phy.move_to_zero(res, 915e6)
res = res[::2]
symbols, loc_start = phy.demodulate(res)
# print(symbols)
codes = phy.decode(symbols)

print(codes[:16])
# print("".join([chr(x) for x in codes[:16]]))

# plt.figure()
# plt.plot(abs(res))
# plt.show()


