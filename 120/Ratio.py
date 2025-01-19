import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def get_time(file_name):
     data = pd.read_csv(file_name, delimiter='\s+')
     Time = data['Time'].to_numpy()
     
     return Time

def get_ratio(file_name):
    data = pd.read_csv(file_name, delimiter='\s+')
    O = data['O'].to_numpy()
    OH = data['OH'].to_numpy()
    Al = data['Al'].to_numpy()
    AlXOH = data['AlXOH'].to_numpy()
    AlXOHH2O = data['AlXOHH2O'].to_numpy()

    ratio = (O+OH+AlXOH+AlXOHH2O)/(Al+AlXOH+AlXOHH2O)
 
    return ratio


Time = get_time('log.data')
ratio = get_ratio('log.data')
plt.plot(Time,ratio,color='red',linewidth = 3)
plt.axhline(y=1.5,xmin=0,xmax=8,ls='--')
plt.xticks(fontweight='bold',fontsize = 15)
plt.yticks(fontweight='bold',fontsize = 15)
plt.title('393K',fontweight='bold',size = 20)
plt.xlabel('Time',fontweight='bold',size = 20)
plt.ylabel('ratio',fontweight='bold',size = 20)
plt.show()