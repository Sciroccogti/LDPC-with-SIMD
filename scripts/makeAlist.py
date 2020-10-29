'''
File: makeAlist.py
File Created: Thursday, 29th October 2020 20:30:07
Author: Yifan Zhang (scirocco_gti@yeah.net)
Last Modified: Thursday, 29th October 2020 20:41:11
'''

import numpy as np

from alist import save_alist

nRow = int(input("nRow is: "))
nCol = int(input("nCol is: "))
mat = np.zeros((nRow, nCol))
for i in range(nRow):
    j = 0
    for n in input().split():
        mat[i][j] = n
        j += 1
fname = input("filename is: ")
save_alist(fname, mat)
