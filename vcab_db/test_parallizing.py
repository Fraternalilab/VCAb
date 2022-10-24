from joblib import Parallel, delayed
import multiprocessing
import numpy as np
import time

def test(rc,cc,x,y):
    result=0
    if cc>rc:
        result=x+y
    return 0

results = np.zeros((100, 100), float)

st1=time.time()
for rc,i in enumerate(range(100)):
    results[rc]=Parallel(n_jobs=1)(delayed(test)(rc,cc,i,j) for cc,j in enumerate(range(100)))
en1=time.time()
print (en1-st1)
#print (results)

results2=np.zeros((100, 100), float)
st2=time.time()
for rc,i in enumerate(range(100)):
    for cc,j in enumerate(range(100)):
        results2[rc][cc]=test(rc,cc,i,j)
en2=time.time()
print (en2-st2)
