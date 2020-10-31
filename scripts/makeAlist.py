'''
@file makeAlist.py
@author Yifan Zhang (scirocco_gti@yeah.net)
@brief 
@date 2020-10-29 20:30:07
@modified: 2020-10-31 12:37:58
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
