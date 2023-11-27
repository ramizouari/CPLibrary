import numpy as np
from numpy import fft

T=int(input())
if T>0:
    for k in range(T):
        line=input()
        a,b,c,d,e,alpha=line.split(" ")
        shape=[a,b,c,d,e]
        shape=list(map(int,shape))
        a,b,c,d,e=shape
        alpha=float(alpha)
        A=np.zeros(shape,dtype=np.complex64)
#        print()
        for i1 in range(a):
            for i2 in range(b):
                for i3 in range(c):
                    for i4 in range(d):
                        for i5 in range(e):
                            A[i1,i2,i3,i4,i5]=(i1^i2^i3^i4^i5)*np.exp((i1-i2+i3-i4+i5)*alpha*1.j)
#        print(f"A: {A.reshape(-1)}")
        B=fft.fftn(A,norm="ortho")
#        print(f"B: {B.reshape(-1)}")
        print(np.mean(np.abs(np.real(B))))
#        print()
