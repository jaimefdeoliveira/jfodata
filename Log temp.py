# -*- coding: utf-8 -*-
"""
Created on Wed May 27 07:52:07 2020

@author: jaime
"""

import pyvisa
import time
import pandas as pd
import matplotlib.pyplot as plt
import csv
import os

rm = pyvisa.ResourceManager()
#it.close()
it=rm.open_resource(rm.list_resources()[0])

it.write_termination='\n'
time.sleep(1)
it.read_termination = '\n'
time.sleep(1)
it.timeout=2
time.sleep(1)
print(it.query('*IDN?',delay=1))
#%%
df=pd.DataFrame(columns= ['Time','T','H'])

filename = r"C:\Users\jaime\OneDrive\PHD\Dados\teste2.csv"
write_header = not os.path.exists(filename) or os.stat(filename).st_size == 0
for i in range(30):
    t=float(it.query(':MEASURE:TEMPrature?',delay=1))
    h=float(it.query(':MEASURE:HUMIdity?',delay=1))
    Time=time.time()
    df.loc[len(df)]=[Time,t,h]
    time.sleep(120)

    
    with open(filename,"a") as f:
        write_header = not os.path.exists(filename) or os.stat(filename).st_size == 0
        writer = csv.writer(f,delimiter="\t")  
        if write_header:
            writer.writerow(['Time','T','H' ])
        writer.writerow([Time,t,h])
plt.plot(df['Time'],df['T'])

it.close()

df.dropna()
