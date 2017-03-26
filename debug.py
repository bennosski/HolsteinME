from numpy import *

t1 = load('test1/t.npy')
t2 = load('test2/t.npy')

print (t1-t2)/(t1+t2)*2
