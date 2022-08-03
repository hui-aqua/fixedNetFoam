import numpy as np
import readof as rof

# processor
cage1=[0,4] 
cage2=[1,5]

time_end=115
time_list=[]
for i in range(time_end):
    time_list.append(i)
    for j in [0.2,0.4,0.6,0.8]:
        time_list.append(i+j)
print(time_list)



def get_drag(cage,t):
    drag=[]
    for item in t:
        part_drag=np.array([0,0,0.0])
        for c in cage:
            fh=rof.readvector('processor'+str(c),str(item),'Fh')    
            part_drag+=np.sum(fh,axis=1)       
        
        drag.append((part_drag).tolist())
    return np.array(drag)

drag1=get_drag(cage1,time_list)
drag2=get_drag(cage2,time_list)
print(drag1)

np.savetxt('drag1.out',  drag1, delimiter=',', fmt='%.6f')
np.savetxt('drag2.out',  drag2, delimiter=',', fmt='%.6f')




