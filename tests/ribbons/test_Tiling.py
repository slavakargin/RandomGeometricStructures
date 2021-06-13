'''
Created on Jun 13, 2021

@author: vladislavkargin
'''
#import numpy as np
#import matplotlib.pyplot as plt
#from itu.algs4.graphs import digraph as dgr


import rgs.ribbons.Ribbon as rb
#import rgs.ribbons.Utility as ut

from rgs.ribbons.Tiling import (Tiling, standardT, standardAztec2)

def test_Tiling():
    ribbon1 = rb.Ribbon(0, 0, [0, 1])
    ribbon2 = rb.Ribbon(0, 1, [1, 0])
    tiling = Tiling([ribbon1, ribbon2]) 
    assert tiling.__str__() == '\n' + str(ribbon1) + '\n' + str(ribbon2)

def test_standardT():
    n = 3
    M = 2
    N = 3
    t = standardT(n, M, N)
    assert str(t) == "\n(0,0), [1, 1]\n(1,0), [1, 1]\n(2,0), [1, 1]\n(0,3), [1, 1]\n(1,3), [1, 1]\n(2,3), [1, 1]"

test_Tiling()
test_standardT()