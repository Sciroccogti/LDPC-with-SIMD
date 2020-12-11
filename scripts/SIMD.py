'''
@file SIMD.py
@author Yifan Zhang (scirocco_gti@yeah.net)
@brief 
@date 2020-11-24 20:16:40
@modified: 2020-11-24 20:25:35
'''

import argparse
import time
from multiprocessing import Pool, Array

import bitarray as ba
import numpy as np
from bitarray.util import ba2int, int2ba

from alist import alist2sparse2


def loop(num: int):
    start = loop_num // threads_num * num
    print(start)
    for i in range(loop_num // threads_num):
        I = np.array(int2ba(start + i, 64).tolist(True))
        count[np.count_nonzero(np.matmul(I, G))] += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-H", "--dec-h-path",
                        help="path to H matrix in .alist")
    args = parser.parse_args()

    threads_num = 16
    loop_num = 4096000

    G = alist2sparse2(args.dec_h_path)
    G = np.flipud(G)
    # G_int = []
    # for col in G.T:
    #     G_int.append(ba2int(ba.bitarray(list(col))))

    iters = list(range(threads_num))
    count = Array('i', 129)
    pool = Pool()
    start = time.time()
    pool.map(loop, iters, 1)
    pool.close()
    pool.join()
    # for i in range(129):
    #     print("%d: %d" % (i, count[i]))
    print(time.time() - start)
