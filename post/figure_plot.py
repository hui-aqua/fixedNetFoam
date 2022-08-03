import matplotlib.pyplot as plt
import numpy as np


f1=np.loadtxt('drag1.out',delimiter=',')
f2=np.loadtxt('drag2.out',delimiter=',')
time_end=115
time_list=[]
for i in range(time_end):
    time_list.append(i)
    for j in [0.2,0.4,0.6,0.8]:
        time_list.append(i+j)
print(time_list)

plt.figure()
plt.plot(time_list[10:],f1[10:,1],label='F1')
plt.plot(time_list[10:],f2[10:,1],label='F2')

print(np.mean(f1[-100:,0]),np.mean(f2[-100:,0]))

plt.legend()
plt.show()
