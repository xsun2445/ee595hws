import numpy as np
import sys
import commpy


class LoRaDecoder():
    def __init__(self, coding_para) -> None:
        self.coding_para = coding_para
        self.Fs = 250*1e3
        self.bandwidth = 125*1e3

    
    def decoding(fileName, len_payload):
        pass

    def 





if __name__ == '__main__':
    # print(sys.argv)
    if len(sys.argv) != 4:
        raise Exception('Incorrect number of inputs. (given {}, should be 3)'.format(len(sys.argv)-1))
    
    fileName = sys.argv[1]
    plen = int(sys.argv[2])
    coding_para = sys.argv[3]

