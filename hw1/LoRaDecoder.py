import numpy as np
import sys
import commpy
import struct
import matplotlib.pyplot as plt

class LoRaDecoder():
    def __init__(self, coding_para) -> None:
        self.coding_para = coding_para
        self.Fs = 250*1e3
        self.bandwidth = 125*1e3

    
    def decoding(self, fileName, len_payload):
        pass

    def read_complex(self, fileName):
        raw = open(fileName, mode='rb').read()
        data = np.array([struct.unpack('f',raw[4*x:4*x+4]) for x in range(len(raw)//4)])
        data = data[:len(data)//2] + 1j*data[len(data)//2:]
        return data
    
    @staticmethod
    def chirp():
        pass






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
    print(data.shape)

    # plt.figure()
    # plt.plot(np.abs(data))
    # plt.show()

    plt.figure()
    plt.plot(np.abs(np.fft.fft(data)))
    plt.show()

    

