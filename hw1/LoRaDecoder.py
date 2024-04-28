import numpy as np
import sys
import commpy
import struct
import matplotlib.pyplot as plt
import scipy


class LoRaDecoder():
    def __init__(self, coding_para=None) -> None:
        self.coding_para = coding_para
        self.Fs = 250*1e3
        self.bandwidth = 125*1e3

    
    def decoding(self, fileName, len_payload):
        pass

    def read_complex(self, fileName):
        raw = open(fileName, mode='rb').read()
        data = np.array([struct.unpack('f',raw[4*x:4*x+4]) for x in range(len(raw)//4)]).astype(np.float32)
        # data = data[:len(data)//2] + 1j*data[len(data)//2:]
        data = data[0::2] + 1j*data[1::2]
        return data.T
    
    @staticmethod
    def chirp(is_up, sf, bw, fs, h, cfo=0, tdelta=0, tscale=1):
        ''' 
        ## Chirp generation:
            - is_up: True for up chirp, False for down chirp
            - sf: Spreading factor (# bits in one chirp)
            - bw: Bandwidth of chirp
            - fs: Sampling frequency
            - h: Starting freq, from 0 to 2^SF-1
            - cfo: Central Frequency Offset
            - tdelta: Time offset, from 0 to 1/fs
            - tscale: Scaling the sampling frequency
        ## return 
            generated LoRa symbol
        '''
        # print(sf)
        # number of bins
        N = 2**sf
        # Time for chirp
        T = N/bw
        # how many adc samples per T: T = fs * T = T / ts
        samp_per_sym = round(fs/bw*N)
        # print(samp_per_sym)
        # starting freq round up
        h = round(h)
        # cfo + starting freq error
        cfo = cfo + (np.round(h) - h) / N * bw
        if is_up:
            k = bw/T
            f0 = -bw/2+cfo
        else:
            k = -bw/T
            f0 = bw/2+cfo

        # f0 always is -half bw + cfo
        # true f0 is f0+k*T*h/N
        # retain last element to calculate phase

        t = np.arange(0, samp_per_sym*(N-h)//N + 1) / fs * tscale + tdelta
        c1 = np.exp(1j*2*np.pi*t*(f0+k*T*h/N+0.5*k*t))


        if len(t) == 0:
            phi = 0
        else:
            phi = np.angle(c1[-1])

        t = np.arange(0, samp_per_sym*h//N)/fs + tdelta
        c2 = np.exp(1j*(phi + 2*np.pi*(t * (f0 + 0.5 * k * t))))


        return np.concatenate((c1[:-1], c2))




        # t = (0:samp_per_sym*(N-h)/N)/fs*tscale + tdelta
        # snum = length(t)
        # c1 = exp(1j*2*pi*(t.*(f0+k*T*h/N+0.5*k*t)))

        # if snum == 0
        #     phi = 0;
        # else
        #     phi = angle(c1(snum));
        # end
        # % then start from -half bw + cfo
        # t = (0:samp_per_sym*h/N-1)/fs + tdelta;
        # c2 = exp(1j*(phi + 2*pi*(t.*(f0+0.5*k*t))));

        # y = cat(2, c1(1:snum-1), c2).';







if __name__ == '__main__':
    # print(sys.argv)
    if len(sys.argv) != 4:
        raise Exception('Incorrect number of inputs. (given {}, should be 3)'.format(len(sys.argv)-1))
    
    fileName = sys.argv[1]
    plen = int(sys.argv[2])
    coding_para = sys.argv[3]

    fileName = './EE595 Project Data/dataSF8CR8packet1.bin'
    decoder = LoRaDecoder(coding_para)
    data = decoder.read_complex(fileName)
    print(type(data))
    print(data.dtype)

    fftData = np.fft.fft(data)

    print(fftData[:10])
    print(fftData.shape)

    fftData = np.fft.fft(data.T, axis=0)
    print(fftData[:10])



    # plt.figure()
    # plt.plot(np.abs(data))
    # plt.show()

    # plt.figure()
    # plt.plot(data.real)
    # plt.show()

    # plt.figure()
    # plt.plot(np.fft.fft(data))
    # plt.title('fft')
    # plt.show()

    # import matplotlib.pyplot as plt
    # t = np.arange(256)
    # sp = np.fft.fft(np.sin(t))
    # freq = np.fft.fftfreq(t.shape[-1])
    # plt.plot(freq, sp.real, freq, sp.imag)
    # plt.show()

    # print(np.fft.fft(np.exp(2j * np.pi * np.arange(8) / 8)))


    

