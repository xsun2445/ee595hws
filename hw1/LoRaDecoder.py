import numpy as np
import sys
import commpy
import struct
import matplotlib.pyplot as plt
import scipy
from functools import reduce


class LoRaDecoder():
    def __init__(self, sf, plen=16, cp=6) -> None:

        self.fs = 250*1e3
        self.fs = self.fs/2
        self.bw = 125*1e3

        self.sf = sf
        self.cfo = 0
        self.cp = cp
        self.ldr = 0
        self.preamble_len = 8
        self.upsample_enable = False
        self.crc_enable = 1
        self.header_enable = 0

        self.samp_per_symbol = int(self.fs*2**self.sf//self.bw)
        self.num_bins = 2**self.sf
        # payload length is always 16 bytes
        self.payload_len = plen
        self.symbol_len = self.calc_symbol_num()
        # print(self.symbol_len)



        # self.whitening_seq = np.array([0xff, 0xfe, 0xfc, 0xf8, 0xf0, 0xe1, 0xc2, 0x85, 0xb, 0x17, 0x2f, 0x5e, 0xbc, 0x78, 0xf1, 0xe3, 0xc6, 0x8d, 0x1a, 0x34, 0x68, 0xd0, 0xa0, 0x40, 0x80, 0x1, 0x2, 0x4, 0x8, 0x11, 0x23, 0x47, 0x8e, 0x1c, 0x38, 0x71, 0xe2, 0xc4, 0x89, 0x12, 0x25, 0x4b, 0x97, 0x2e, 0x5c, 0xb8, 0x70, 0xe0, 0xc0, 0x81, 0x3, 0x6, 0xc, 0x19, 0x32, 0x64, 0xc9, 0x92, 0x24, 0x49, 0x93, 0x26, 0x4d, 0x9b, 0x37, 0x6e, 0xdc, 0xb9, 0x72, 0xe4, 0xc8, 0x90, 0x20, 0x41, 0x82, 0x5, 0xa, 0x15, 0x2b, 0x56, 0xad, 0x5b, 0xb6, 0x6d, 0xda, 0xb5, 0x6b, 0xd6, 0xac, 0x59, 0xb2, 0x65, 0xcb, 0x96, 0x2c, 0x58, 0xb0, 0x61, 0xc3, 0x87, 0xf, 0x1f, 0x3e, 0x7d, 0xfb, 0xf6, 0xed, 0xdb, 0xb7, 0x6f, 0xde, 0xbd, 0x7a, 0xf5, 0xeb, 0xd7, 0xae, 0x5d, 0xba, 0x74, 0xe8, 0xd1, 0xa2, 0x44, 0x88, 0x10, 0x21, 0x43, 0x86, 0xd, 0x1b, 0x36, 0x6c, 0xd8, 0xb1, 0x63, 0xc7, 0x8f, 0x1e, 0x3c, 0x79, 0xf3, 0xe7, 0xce, 0x9c, 0x39, 0x73, 0xe6, 0xcc, 0x98, 0x31, 0x62, 0xc5, 0x8b, 0x16, 0x2d, 0x5a, 0xb4, 0x69, 0xd2, 0xa4, 0x48, 0x91, 0x22, 0x45, 0x8a, 0x14, 0x29, 0x52, 0xa5, 0x4a, 0x95, 0x2a, 0x54, 0xa9, 0x53, 0xa7, 0x4e, 0x9d, 0x3b, 0x77, 0xee, 0xdd, 0xbb, 0x76, 0xec, 0xd9, 0xb3, 0x67, 0xcf, 0x9e, 0x3d, 0x7b, 0xf7, 0xef, 0xdf, 0xbf, 0x7e, 0xfd, 0xfa, 0xf4, 0xe9, 0xd3, 0xa6, 0x4c, 0x99, 0x33, 0x66, 0xcd, 0x9a, 0x35, 0x6a, 0xd4, 0xa8, 0x51, 0xa3, 0x46, 0x8c, 0x18, 0x30, 0x60, 0xc1, 0x83, 0x7, 0xe, 0x1d, 0x3a, 0x75, 0xea, 0xd5, 0xaa, 0x55, 0xab, 0x57, 0xaf, 0x5f, 0xbe, 0x7c, 0xf9, 0xf2, 0xe5, 0xca, 0x94, 0x28, 0x50, 0xa1, 0x42, 0x84, 0x9, 0x13, 0x27, 0x4f, 0x9f, 0x3f, 0x7f]).astype(np.uint8)
        self.whitening_seq = self.gen_whitening_seq()
        self.upchirp = self.chirp(True, sf, self.bw, self.fs, 0, 0, 0)
        self.downchirp = self.chirp(False, sf, self.bw, self.fs, 0, 0, 0)



    def demodulate(self, raw):
        # the input raw data should has fs 125e3
        loc_start, self.cfo = self.start_detection(raw)
        # always 16 bytes payload, then always self.symbol_len symbols
        chirp_list = raw[loc_start:loc_start+self.symbol_len*self.samp_per_symbol] 
        # correct the data with cfo
        chirp_list *= np.exp(-1j*2*np.pi*self.cfo*(1/self.fs)*np.arange(len(chirp_list)))

        # reshape s.t. each row is a chirp
        chirp_list = chirp_list.reshape(-1, self.samp_per_symbol)
        fftData = np.fft.fft(chirp_list*self.downchirp, axis=1)
        # print(fftData.shape)
        symbols = np.argmax(fftData, axis=1)
        # print(symbols)
        # print(symbols.shape)
        return symbols


    def start_detection(self, raw):
        preamble = np.matlib.repmat(self.downchirp, 1, self.preamble_len).flatten()
        # since the netid is unknown, so just avoid multiply it
        netid = np.zeros(2*self.samp_per_symbol)
        sync = np.concatenate((self.upchirp,self.upchirp,self.upchirp[:len(self.upchirp)//4]))
        preamble = np.concatenate((preamble, netid, sync))
        max_val = 0
        loc_time = 0
        loc_freq = 0
        for i in range(len(raw)-len(preamble)):
            temp_fft = abs(np.fft.fft(raw[i:i+len(preamble)]*preamble))
            # if self.enable_upsample:
                # temp_fft = temp_fft[:len(temp_fft)//2] + temp_fft[len(temp_fft)//2:]
            temp_freq = np.argmax(temp_fft)
            temp_max = temp_fft[temp_freq]
            if temp_max > max_val:
                max_val = temp_max
                loc_time = i
                loc_freq = temp_freq
        
        loc_freq %= self.num_bins
        loc_freq = loc_freq/len(preamble)*self.samp_per_symbol *self.fs/self.num_bins
        loc_time += len(preamble)

        return loc_time, loc_freq
            


    
    def decode(self, data_demod):
        # grey decoding
        data_interleaved = self.gray_coding(data_demod)

        # deinterleaving
        data_hamming = self.deinterleaving(data_interleaved)

        # hamming code
        data_nibble = self.hamming_decode(data_hamming)

        # nibble -> bytes
        data_whitened = np.array([(data_nibble[2*i+1]<<np.uint8(4))+data_nibble[2*i] for i in range(len(data_nibble)//2)], dtype=np.uint8)
        
        # dewhitening
        data_crc = self.whiten(data_whitened)

        # crc
        data_decoded = data_crc

        return data_decoded
    

    def gray_coding(self, data):
        # header
        data[:8] = np.floor(data[:8] / 4)
        if self.ldr:
            data[8:] = np.floor(data[8:]/4)
        else:
            data[8:] = (data[8:]-1) % (2**self.sf)
        data = data.astype(np.uint16)
        data = np.bitwise_xor(data, (data >> np.uint16(1)))
        
        return data


    def deinterleaving(self, interleaved):
        '''diagonol deinterleaving process
        input is uint16, output is uint8
        for each symbol:
        interleaving:   [cp, sf-2*ldr] -> [sf-2*ldr, cp] 
        deinterleaving: [sf-2*ldr, cp] -> [cp, sf-2*ldr] 
        '''
        def helper(symbol, ldr):
            # uint16 -> uint8 unpacked using big endian
            symbol = np.unpackbits(symbol.view(np.uint8), bitorder='little').reshape(len(symbol), -1)
            # print(symbol)
            # choose those valid bits containing information without ldro bits
            # symbol = symbol[:, -self.sf+2*ldr:]
            symbol = symbol[:,:self.sf-2*ldr]
            # symbol = symbol[:,::-1]
            # diagonol deinterleaving
            symbol = np.array([np.roll(x, idx) for idx, x in enumerate(symbol)], dtype=np.uint16)
            symbol = symbol.T
            symbol = np.packbits(symbol, bitorder='little', axis=1)
            return symbol.flatten()

        deinterleaved = []
        # header has cp=8, has shape [8, sf-2]
        # transposed back: [sf-2, 8]
        deinterleaved.append(helper(interleaved[:8], 1))

        # for every #cp coding parameter symbols, -> has shape [cp,sf-2*ldr]
        # transposed back: [sf-2*ldr, cp]
        for i in range(8, len(interleaved), self.cp):
            temp = interleaved[i:i+self.cp]
            # print(temp)
            deinterleaved.append(helper(temp, self.ldr))
        
        deinterleaved = np.concatenate([np.array(x) for x in deinterleaved])

        return deinterleaved


    def hamming_decode(self, data_nibble):
        ## TODO: hamming correction 
        # return np.array([nibble & np.uint8(0x0F) for nibble in data_nibble], dtype=np.uint8)
        decoded_header = np.array([self.hamming_decode_nibble(x, 8) for x in data_nibble[:self.sf-2]])
        decoded_payload = np.array([self.hamming_decode_nibble(x, self.cp) for x in data_nibble[self.sf-2:]])

        return np.concatenate((decoded_header, decoded_payload))
    
    def hamming_decode_nibble(self, nibble, cp):
        '''hamming decoding and correction process
        - nibble: uint8 hamming coded symbol
        - cp: coding parameter, cr = 4/cp
        output:
            decoded and corrected nibble (uint8)
        '''
        [d8, d7, d6, d5, d4, d3, d2, d1] = np.unpackbits(nibble)
        # p1 = d8^d4^d3^d1
        # p2 = d7^d4^d2^d1
        # p3 = d5^d3^d2^d1
        # p4 = d5^d4^d3^d2^d1
        # p5 = d6^d4^d3^d2
        p1 = d1^   d3^d4^d8
        p2 = d1^d2^   d4^d7
        p3 = d1^d2^d3^   d5
        p4 = d1^d2^d3^d4^d5
        p5 =    d2^d3^d4^d6

        def correction_47(parity):
            # p2p3p5
            if parity == 3:     #0011
                res = 4         #0100
            elif parity == 5:   #0101
                res = 8         #1000
            elif parity == 6:   #0110
                res = 1         #0001
            elif parity == 7:   #0111
                res = 2         #0010
            else:
                res = 0
            return np.uint8(res)
        
        # TODO: 4/5, 4/6, 4/8
        # finished: 4/7 hamming code
        if cp == 5:
            pass
        elif cp == 6:
            pass
        elif cp == 7:
            nibble = nibble ^ correction_47(p2*4+p3*2+p5)
        elif cp == 8:
            nibble = nibble ^ correction_47(p2*4+p3*2+p5)

        return nibble & np.uint8(0x0F)

    
    def hamming_encode(self, data_nibble):
        encoded_header = np.array([self.hamming_encode_nibble(x, 8) for x in data_nibble[:self.sf-2]])
        encoded_payload = np.array([self.hamming_encode_nibble(x, self.cp) for x in data_nibble[self.sf-2:]])

        return np.concatenate((encoded_header, encoded_payload))


    def hamming_encode_nibble(self, nibble, cp):
        '''hamming encode for uint8 nibble (only has lower 4-bit)
        - nibble: uint8 0x0-0xf, (upper 4 digits are 0s)
        - self.cp: coding parameter, coding rate = 4/cp
        output:
        hamming encodded nibble (uint8)
        '''
        # [d1,d2,d3,d4] = np.unpackbits(nibble)[-4:]
        # p1 = d1^d2^d4
        # p2 = d1^d3^d4
        # p3 = d2^d3^d4
        # p4 = d1^d2^d3^d4
        # p5 = d1^d2^d3

        [d4,d3,d2,d1] = np.unpackbits(nibble)[-4:]
        p1 = d1^   d3^d4
        p2 = d1^d2^   d4
        p3 = d1^d2^d3
        p4 = d1^d2^d3^d4
        p5 =    d2^d3^d4

        if cp == 5:
            codeword = (p4<<np.uint8(4)) + nibble
        elif cp == 6:
            codeword = (p5<<np.uint8(5)) + (p3<<np.uint8(4)) + nibble
        elif cp == 7:
            codeword = (p2<<np.uint8(6)) + (p5<<np.uint8(5)) + (p3<<np.uint8(4)) + nibble
        elif cp == 8:
            codeword = (p1<<np.uint8(7)) + (p2<<np.uint8(6)) + (p5<<np.uint8(5)) + (p3<<np.uint8(4)) + nibble
        else:
            raise Exception('wrong number of coding parameter (should be 5-8, {%d} given)'.format(cp))

        return codeword

    
    def gen_whitening_seq(self):
        '''generate the whitening sequence for randomize the signal by xor with it
        whitening: codeword xor whitening_seq = symbol
        inverse process: symbol xor whitening_seq = codeword
        '''
        reg = np.uint8(0xff)
        num = 255
        seq = []
        for i in range(num):
            seq.append(reg)
            reg = (reg<<np.uint8(1)) ^ (LoRaDecoder.getbit(reg,8) ^ (LoRaDecoder.getbit(reg,6) ^ (LoRaDecoder.getbit(reg,5) ^ LoRaDecoder.getbit(reg,4))))
        return np.array(seq).astype(np.uint8)
    

    def whiten(self, data_bytes):
        return np.array([x^y for x,y in zip(data_bytes, self.whitening_seq)]).astype(np.uint8)
    

    def calc_symbol_num(self):
        plen = self.payload_len
        res = 8 + max(self.cp*np.ceil((2*plen-self.sf+7+4*self.crc_enable-5*(1-self.header_enable))/(self.sf-2*self.ldr)), 0)
        return int(res)

    
    @staticmethod
    def getbit(data, n):
        '''get n-th bit of uint8 data
        '''
        return (data >> np.uint8(n-1)) & np.uint8(1)
    
    
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

    @staticmethod
    def read_complex(fileName):
        raw = open(fileName, mode='rb').read()
        data = np.array([struct.unpack('f',raw[4*x:4*x+4]) for x in range(len(raw)//4)]).astype(np.float32)
        # data = data[:len(data)//2] + 1j*data[len(data)//2:]
        data = data[0::2] + 1j*data[1::2]
        return data.flatten()



if __name__ == '__main__':
    # print(sys.argv)
    if len(sys.argv) != 4:
        raise Exception('Incorrect number of inputs. (given {}, should be 3)'.format(len(sys.argv)-1))
    
    fileName = sys.argv[1]
    plen = int(sys.argv[2])
    cp = sys.argv[3]

    fileName = './EE595 Project Data/dataSF8CR8packet1.bin'
    decoder = LoRaDecoder(sf, plen, cp)
    data = decoder.read_complex(fileName)
    print(type(data))
    print(data.dtype)

    fftData = np.fft.fft(data)

    print(fftData[:10])
    print(fftData.shape)

    fftData = np.fft.fft(data.T, axis=0)
    print(fftData[:10])






    

