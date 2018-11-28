# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 14:05:00 2018

@author: jaime

pacote para plotes de dados com fuçoes rotineiras do pacote jfodata


"""
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy import interpolate
from scipy.optimize import curve_fit
import numpy as np
import math
import pandas as pd

plt.rcParams['axes.linewidth'] = 1.4
plt.rcParams['figure.figsize']= (12,9)

plt.rcParams['font.size'] = 32
plt.rcParams['legend.fontsize'] = 'large'
plt.rcParams['figure.titlesize'] = 'medium'
plt.rcParams['axes.labelsize']=32
plt.rcParams['savefig.bbox']='tight'

work=os.getcwd().split('PHD')[0]+ r'PHD\Medidas ADR\Dados py\functions'
os.chdir (r'C:\Users\jaime\Google Drive\PHD\Medidas ADR\Dados py\functions')

class plot:
        """
         seila


        """
        def __init__(self,x,y,dic={}):
             
            self.x = x
            self.y = y
            plt.plot(x,y,**dic)
            plt.xlabel(x.name)           
            plt.ylabel(y.name) 
        
        def nx(self,xn,dic={}):
            '''
            função faz plots normalizados dos dados em relação a um ponto em x
            '''
            yni=(np.abs(self.x-xn)).idxmin()
            plt.plot(self.x,self.y/self.y[yni],**dic)
        def plot(self,dic={}):
            plt.plot(self.x,self.y,**dic)
             
            
class urso:
    '''
        fazer coisas com o pandas
        Cria um data frame com as colunas importantes para se trabalhar 
    '''
    def __init__(self,df):
        self.df=df
        df
    
    def new(self,lista):
        dfn=pd.DataFrame(self.df[lista[0]])
        for i in np.arange(1,len(lista),1):
            dfn[lista[i]]=self.df[lista[i]]
        return dfn

 