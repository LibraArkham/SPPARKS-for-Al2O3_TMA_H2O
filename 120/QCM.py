import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#读取log.spparks中表的行数
def read_lines(file_path, str):
  str_line = 0
  with open (file_path, 'r') as f:
       for line in f:
           str_line += 1
           if line.startswith(str):
               break
  
  return str_line

#保存为log.data
def save_table(file_path, start_line, end_line):
    with open (file_path, 'r') as f:
        lines = f.readlines()
            
    with open ('log.data', 'w') as d:
        d.writelines(lines[start_line - 1:end_line])

#作图
def draw_data(file_path, start_str, end_str, y, x, T):
    start_line = read_lines(file_path, start_str)
    end_line = read_lines(file_path, end_str)
    save_table(file_path, start_line, end_line)
    data = pd.read_csv('log.data', delimiter='\s+')
    y1 = data['QCM'].to_numpy()
    x1 = data['Time'].to_numpy()
    plt.plot(x1,(y1-y1[0])/(19.72*6.02*35.76)*10)
    plt.title(T)
    plt.xlabel(x)
    plt.ylabel(y)
    plt.savefig()



Temperature = '373K'
file_path = 'log.spparks'  
start_str = '      Time'
end_str = '         8'
y = 'Mass gain(ng/cm^2)'
x = 'Time(s)'
draw_data(file_path, start_str, end_str, y, x, Temperature)