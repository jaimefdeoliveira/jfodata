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
import matplotlib.gridspec as gridspec

plt.rcParams['axes.linewidth'] = 1.4
plt.rcParams['figure.figsize']= (12,9)

plt.rcParams['font.size'] = 32
plt.rcParams['legend.fontsize'] = 'large'
plt.rcParams['figure.titlesize'] = 'medium'
plt.rcParams['axes.labelsize']=32
plt.rcParams['savefig.bbox']='tight'

work=os.getcwd().split('PHD')[0]+ r'PHD\Medidas ADR\Dados py\functions'
os.chdir (work)


def pdinter(data1,dadox,dadoy,xmin,xmax,npoints=1000,smoothing=1,K=3):
    ''' 
    função que faz a interpolção em uA faixa 
    '''
    import numpy as np
    from scipy import interpolate
    import matplotlib.pyplot as plt
    data=data1.copy()
    data.sort_values(dadox,ascending=True,inplace=True)
    data.dropna(axis=1 ,how='all', inplace=True)
    data.drop_duplicates(dadox,inplace=True)
    df=data[(data[dadox]<xmax) & (data[dadox]>xmin)] 
    
    
    xs=np.linspace(xmin,xmax,npoints)
    
    spl = interpolate.UnivariateSpline(df[dadox], df[dadoy], k=K , s=smoothing,ext=3)
    return spl , xs


def foldershow(path,ext,Print=True):
    '''
    function to show all file to import files from an folder
    '''
    os.chdir (path)
    a= [name for name in os.listdir(".") if name.endswith("."+ ext)]
    if Print==True:
        for i in range(len(a)):
            print('%i %s'%(i,a[i]))
    return a, len(a)

def dataimport(path,ext,n,separador='\t',skip=1):
    '''
    function to show all file to import files from an folder
    '''
    os.chdir (path)
    a= [name for name in os.listdir(".") if name.endswith("."+ ext)]
    return pd.read_csv(a[n],sep=separador,skiprows=skip)

def saveRT(data,T,R):
    df=pd.DataFrame(columns=['T','R','rho','sigA'])
    df['T']=data[T]
    df['R']=data[R]
    df['rho']=(0.05/1.221)*data[R]
    df['sigA']=1/df['rho']
    return df

def saveMR(data,T,R,F,rate):
    df=pd.DataFrame(columns=['T','R','rho','sigA','Field'])
    df['T']=data[T]
    df['R']=data[R]
    df['rho']=(0.05/1.221)*data[R]
    df['sigA']=1/df['rho']
    df['Field']=rate*data[F]
    return df

def saveMREDC(data,T,V,F,I):
    df=pd.DataFrame(columns=['T(K)','R (Ohns)','Field (T)'])
    df['T(K)']=data[T]
    df['R (Ohns)']=data[V]/(I/1000)
    df['Field (T)']=data[F]/10000
    return df
#%%

def dex(df,columns,value):
    
    'isso returna o index do row com o valor que você quer'
    import numpy as np
    return (np.abs(df[columns]-value)).idxmin()
#%%
def separacurvaQV(df):
    
    try:
        df1=pd.DataFrame(data={'T(K)':df['Temperature (K)'],'B(T)':round(df['Field (Oe)']/10000,4),'Angle (deg)':df['Sample Position (deg)'],'V(V)':df['In Phase Voltage Ampl Ch1 (V)'],'4V(V)':df['Quadrature Voltage Ch1 (V)'],'R(ohms)':df['Resistance Ch1 (Ohms)'],"Phase(deg)":df["Phase Angle Ch1 (deg)"],"I(A)":df['AC Current Ch1 (mA)']/1000,"f(Hz)":df['Frequency Ch1 (Hz)'],"2 Harm":df['2nd Harmonic Ch1 (dB)'],"W(rad/s)":2*np.pi*df['Frequency Ch1 (Hz)']})
        df1.dropna(inplace=True)
        df1.reset_index(drop=True,inplace=True)
    
    except:
        df1=pd.DataFrame()
    try:
        df2=pd.DataFrame(data={'T(K)':df['Temperature (K)'],'B(T)':round(df['Field (Oe)']/10000,2),'Angle (deg)':df['Sample Position (deg)'],'V(V)':df['In Phase Voltage Ampl Ch2 (V)'],'4V(V)':df['Quadrature Voltage Ch2 (V)'],'R(ohms)':df['Resistance Ch2 (Ohms)'],"Phase(deg)":df["Phase Angle Ch2 (deg)"],"I(A)":df['AC Current Ch2 (mA)']/1000,"f(Hz)":df['Frequency Ch2 (Hz)'],"2 Harm":df['2nd Harmonic Ch2 (dB)'],"W(rad/s)":2*np.pi*df['Frequency Ch2 (Hz)']})
        df2.dropna(inplace=True)
        df2.reset_index(drop=True,inplace=True)
    except:
        df2=pd.DataFrame()
    return df1,df2

def MRlist(df_dict,T):
    listaMR=[]
    for i in df_dict.keys():
        md,ch1,ch2=data_class(df_dict[i])
        if md==1:
            if round(df_dict[i]['Temperature (K)'][0])==T:
                    listaMR.append(i)
    return listaMR

def data_class(df):
    """ Select type o data in MR, Ang scan ou RvsT
    
    
    Return:
     (mod,ch1,ch2)
     mod:
            0 RvsT
            1 MR
            2 Ang scan
    
    
    
    """
    ch1=ch2=0
    try:       
        if round(df['Temperature (K)'].var())!=0:
            mod=0
        elif round(df['Field (Oe)'].var())!=0:
            mod=1
        else: 
            mod=2
            
        if len(df['Resistance Ch1 (Ohms)'].dropna())!=0:
            ch1=1
        if len(df['Resistance Ch2 (Ohms)'].dropna())!=0:
            ch2=1
    except:
        mod=3
    return mod,ch1,ch2


def ASSYextractinter(dfinit,X='B(T)',Y='R(ohms)',Hmax=8,K=5,smoothing=0.000000002):
    import jfodata
    
    try:
        Hmax=dfinit['B(T)'].max()
        inter=jfodata.pdinter(dfinit,X,Y,-Hmax,Hmax,smoothing=smoothing,npoints=400,K=K)
        df=pd.DataFrame(columns=['T(K)','Angle (deg)',X,Y])
        for jj in range(len(inter[1])):
            df.loc[jj]=[dfinit['T(K)'][0],dfinit['Angle (deg)'][0],inter[1][jj],inter[0](inter[1][jj])]
        
        
        
        dfp=df[df[X]>0.01]
        dfp=dfp.sort_values(X)
        dfp.reset_index(drop=True,inplace=True)
        
        dfn=df[df[X]<-0.01]
        dfn=dfn.sort_values(X,ascending=False)
        dfn.reset_index(drop=True,inplace=True)
        
        dfnAS=(dfn[Y]-dfp[Y])/2
        dfpAS=(dfp[Y]-dfn[Y])/2
        
        
        
        dfnSY=(dfn[Y]+dfp[Y])/2
        dfpSY=(dfp[Y]+dfn[Y])/2
        
        
        
        dfn=dfn.assign(AS=dfnAS,SY=dfnSY)
        dfp=dfp.assign(AS=dfpAS,SY=dfpSY)
        dfn=dfn.dropna()
        dfp=dfp.dropna()
        
        dff=pd.concat([dfn,dfp]).sort_values(X).astype(float)
    except:
        dff=pd.DataFrame(columns=['T(K)','Angle (deg)',X,Y])
    return dff.reset_index(drop=True)


class DataDYMRT:
    ### Classe for plot adn analise Data for telurium projetc
    import jfodata as jfo 
    def __init__(self,df):
        ''' Convert tudo para dataframe legal
        
        '''
        import jfodata as jfo
        self.mod,self.ch1,self.ch2=jfo.data_class(df)
        self.CH1,self.CH2=jfo.separacurvaQV(df)
        self.CH1.drop_duplicates(subset='B(T)',inplace=True)
        self.CH2.drop_duplicates(subset='B(T)',inplace=True)
    def RT(self,CH=1):
        if CH==1:
            plt.plot(self.CH1["T(K)"],self.CH1["R(ohms)"])
            
        else:
            plt.plot(self.CH2["T(K)"],self.CH2["R(ohms)"])
            
            
    def MR(self,CH=1,descrip=True,label='f and I'):
        if CH==1:
            plt.plot(self.CH1["B(T)"],self.CH1["R(ohms)"],label=label)
            
        else:
            plt.plot(self.CH2["B(T)"],self.CH2["R(ohms)"],label=label)
        if descrip==True:
            plt.annotate(r'%s:%0.2f'%('Temperaure',self.CH1["T(K)"][0]),(0.05,0.9),xycoords='axes fraction',fontsize=18)
            plt.annotate(r'%s:%0.2f'%('Positions',self.CH1['Angle (deg)'][0]),(0.05,0.85),xycoords='axes fraction',fontsize=18)
            plt.annotate(r'%s:%0.2f'%('freq:',self.CH1["f(Hz)"][0]),(0.05,0.80),xycoords='axes fraction',fontsize=18)
            plt.annotate(r'%s:%0.2f'%('current:',self.CH1["I(A)"][0]),(0.05,0.75),xycoords='axes fraction',fontsize=18)
        plt.xlabel('H(T)')
        plt.ylabel('R')
    def QV(self,CH=1,descrip=True,label='f and I'):
        if CH==1:
            plt.plot(self.CH1["B(T)"],self.CH1["4V(V)"],label=label)
            
        else:
            plt.plot(self.CH2["B(T)"],self.CH2["4V(V)"],label=label)
        if descrip==True:
            plt.annotate(r'%s:%0.2f'%('Temperaure',self.CH1["T(K)"][0]),(0.05,0.9),xycoords='axes fraction',fontsize=18)
            plt.annotate(r'%s:%0.2f'%('Positions',self.CH1['Angle (deg)'][0]),(0.05,0.85),xycoords='axes fraction',fontsize=18)
            plt.annotate(r'%s:%0.2f'%('freq:',self.CH1["f(Hz)"][0]),(0.05,0.80),xycoords='axes fraction',fontsize=18)
            plt.annotate(r'%s:%0.2f'%('current:',self.CH1["I(A)"][0]),(0.05,0.75),xycoords='axes fraction',fontsize=18)
        plt.xlabel('H(T)')
        plt.ylabel('QV')

    def sArtplotchart(self,descrip=True):
        '''' 
            Criar um chart sArt do DF  com os valores de phase, In phase e QV e 4V e in phase
            
            Input: Dataframe 
            
            output fig com 4 graficos 
        '''
        
        if self.mod==0:
            Xaxis='T(K)'
        elif self.mod==1:
            Xaxis='B(T)'
        else:
            Xaxis='Angle (deg)'
            
            
        if self.ch1==1:
            data=self.CH1
        else:
            data=self.CH2
        
        fig, ax=plt.subplots(nrows=2,ncols=2,figsize=(36,36),constrained_layout=True)
        
        ax[0,0].plot(data[Xaxis],data['V(V)'])
        
        ax[0,0].set_xlabel(Xaxis)
        ax[0,0].set_ylabel("V(V)")
        
        
        ax[0,1].plot(data[Xaxis],data['4V(V)'])
        ax[0,1].set_xlabel(Xaxis)
        ax[0,1].set_ylabel("4V(V)")
        
        
        
        ax[1,0].plot(data[Xaxis],data['2 Harm'])
        ax[1,0].set_xlabel(Xaxis)
        ax[1,0].set_ylabel("2 Harm (dB)")
        
        
        
        ax[1,1].plot(data["V(V)"],data['4V(V)'])
        ax[1,1].set_xlabel("V(V)")
        ax[1,1].set_ylabel("4V(V)")
        
        if descrip==True and self.mod==1:
            ax[0,0].annotate(r'%s:%0.2f K'%('Temperature',self.CH1["T(K)"][0]),(0.4,0.9),xycoords='axes fraction',fontsize=28)
            ax[0,0].annotate(r'%s:%0.2f deg'%('Position',self.CH1['Angle (deg)'][0]),(0.4,0.85),xycoords='axes fraction',fontsize=28)
            ax[0,0].annotate(r'%s:%0.2f Hz'%('Freq',self.CH1["f(Hz)"][0]),(0.4,0.80),xycoords='axes fraction',fontsize=28)
            ax[0,0].annotate(r'%s:%0.2f A'%('Current',self.CH1["I(A)"][0]),(0.4,0.75),xycoords='axes fraction',fontsize=28)
        
        
        return

class DataChiral(DataDYMRT):
    
    
    def __init__(self,df):
        import jfodata as jfo 
        self.mod,self.ch1,self.ch2=jfo.data_class(df)
        self.CH1,self.CH2=jfo.separacurvaQV(df)
        self.CH1ASSY=jfo.ASSYextractinter(self.CH1)
        self.CH2ASSY=jfo.ASSYextractinter(self.CH2)
    
    def sArtplotchartSY(self,descrip=True):
        '''' 
            Criar um chart sArt do DF  com os valores de phase, In phase e QV e 4V e in phase
            
            Input: Dataframe 
            
            output fig com 4 graficos 
        '''
        

        if self.mod==1:
            Xaxis='B(T)'
        else:
            pass 
            
            
        if self.ch1==1:
            data=self.CH1ASSY
        else:
            data=self.CH2ASSY
        
        fig, ax=plt.subplots(nrows=2,ncols=2,figsize=(36,36),constrained_layout=True)
        
        ax[0,0].plot(data[Xaxis],data['R(ohms)'])
        
        ax[0,0].set_xlabel(Xaxis)
        ax[0,0].set_ylabel("R(ohms)")
        
        
        ax[0,1].plot(data[Xaxis],data['AS'])
        ax[0,1].set_xlabel(Xaxis)
        ax[0,1].set_ylabel("AS")
        
        
        
        ax[1,0].plot(data[Xaxis],data['SY'])
        ax[1,0].set_xlabel(Xaxis)
        ax[1,0].set_ylabel("SY")
        
        
        
        ax[1,1].plot(data["SY"],data["SY"])
        ax[1,1].set_xlabel("")
        ax[1,1].set_ylabel("")
        
        if descrip==True and self.mod==1:
            ax[0,0].annotate(r'%s:%0.2f K'%('Temperature',self.CH1["T(K)"][0]),(0.4,0.9),xycoords='axes fraction',fontsize=28)
            ax[0,0].annotate(r'%s:%0.2f deg'%('Position',self.CH1['Angle (deg)'][0]),(0.4,0.85),xycoords='axes fraction',fontsize=28)
            ax[0,0].annotate(r'%s:%0.2f Hz'%('Freq',self.CH1["f(Hz)"][0]),(0.4,0.80),xycoords='axes fraction',fontsize=28)
            ax[0,0].annotate(r'%s:%0.2f A'%('Current',self.CH1["I(A)"][0]),(0.4,0.75),xycoords='axes fraction',fontsize=28)
        
        
        return

class QVALL(DataChiral):
    
    def __init__(self,df):
        import jfodata as jfo 
        self.mod,self.ch1,self.ch2=jfo.data_class(df)
        self.CH1,self.CH2=jfo.separacurvaQV(df)
        self.CH1['Indtutance']=self.CH1['4V(V)']/(self.CH1['I(A)'].max())
        self.CH2['Indtutance']=self.CH2['4V(V)']/(self.CH2['I(A)'].max())
        self.CH1['X']=self.CH1['4V(V)']/(self.CH1['I(A)'].max())
        self.CH2['X']=self.CH2['4V(V)']/(self.CH2['I(A)'].max())
       
        self.dhall=1
        self.vdis=1
        self.croos=np.pi*(self.dhall/2)**2
        self.name='name'
        if len(self.CH1)!=0:
            self.T=round(self.CH1['T(K)'].max(),2)
        else:
            self.T=round(self.CH2['T(K)'].max(),2)
        try:
            self.angle=round(self.CH1['Angle (deg)'].max())
        except:
            self.angle=round(self.CH1['Angle (deg)'].max(),2)
        self.I=self.CH1['I(A)'].max()
        self.F=self.CH1['f(Hz)'].max()
        self.W=self.CH1['W(rad/s)'].max()
    def QI(self,CH=1,descrip=True,label='f and I'):
        if CH==1:
            plt.plot(self.CH1["B(T)"],self.CH1["Indtutance"],label=label)
            
        else:
            plt.plot(self.CH2["B(T)"],self.CH2["Indtutance"])
        if descrip==True:
            plt.annotate(r'%s:%0.2f'%('Temperaure',self.CH1["T(K)"][0]),(0.05,0.9),xycoords='axes fraction',fontsize=18)
            plt.annotate(r'%s:%0.2f'%('Positions',self.CH1['Angle (deg)'][0]),(0.05,0.85),xycoords='axes fraction',fontsize=18)
            plt.annotate(r'%s:%0.2f'%('freq:',self.CH1["f(Hz)"][0]),(0.05,0.80),xycoords='axes fraction',fontsize=18)
            plt.annotate(r'%s:%0.2f'%('current:',self.CH1["I(A)"][0]),(0.05,0.75),xycoords='axes fraction',fontsize=18)  

        plt.xlabel('B(T)')
        plt.ylabel('L($\Omega$)')
    def SyAydataV(self,s=0.000000002,K=5):

        from jfodata import ASSYextractinter
        self.CH1VASSY=ASSYextractinter(self.CH1,'B(T)','V(V)',8.9,K,s)
        self.CH2VASSY=ASSYextractinter(self.CH2,'B(T)','V(V)',8.9,K,s)
    def SyAydataR(self,s=0.000000002,K=5):
    
        from jfodata import ASSYextractinter
        self.CH1VASSY=ASSYextractinter(self.CH1,'B(T)','R(ohms)',8.9,K,s)
        self.CH2VASSY=ASSYextractinter(self.CH2,'B(T)','R(ohms)',8.9,K,s)
    def resistivit(self):
        from jfodata import dex
        if len(self.CH1)!=0:
            if self.mod==0: 
                self.CH1rho=(self.croos*self.CH1['R(ohms)'][self.CH1['T(K)'].idxmax()])/self.vdis
            elif self.mod==1: 
                self.CH1rho=(self.croos*self.CH1['R(ohms)'][dex(self.CH1,'B(T)',0)])/self.vdis
            elif self.mod==2:
                self.CH1rho=(self.croos*self.CH1['R(ohms)'][dex(self.CH1,'Angle (deg)',0)])/self.vdis
        else:
            self.CH1rho=0    

        if len(self.CH2)!=0:
            if self.mod==0: 
                self.CH2rho=(self.croos*self.CH2['R(ohms)'][self.CH2['T(K)'].idxmax()])/self.vdis
            elif self.mod==1: 
                self.CH2rho=(self.croos*self.CH2['R(ohms)'][dex(self.CH2,'B(T)',0)])/self.vdis
            elif self.mod==2:
                self.CH2rho=(self.croos*self.CH2['R(ohms)'][dex(self.CH2,'Angle (deg)',0)])/self.vdis
        else:
            self.CH2rho=0                  
                
                
    def SyAydataQV(self,s=0.000000002,K=5):
        from jfodata import ASSYextractinter
        self.CH1QVASSY=ASSYextractinter(self.CH1,'B(T)','4V(V)',8.8,K,s)
        self.CH2QVASSY=ASSYextractinter(self.CH2,'B(T)','4V(V)',8.8,K,s)
        
    def SyAydataX(self,s=0.000000002,K=5):
        from jfodata import ASSYextractinter
        self.CH1QVASSY=ASSYextractinter(self.CH1,'B(T)','X',8.8,K,s)
        self.CH2QVASSY=ASSYextractinter(self.CH2,'B(T)','X',8.8,K,s)   
        
    def SyAydataQVXLCH1(self):
        try:       
            self.CH1QVASSY['X']=self.CH1QVASSY['SY']/(self.CH1['I(A)'].max())
            self.CH1QVASSY['L']=self.CH1QVASSY['SY']/(self.CH1['f(Hz)'].max()*self.CH1['I(A)'].max())
        except:
            pass
        
    def SyAydataQVXLCH2(self,):    
        try:
            self.CH2QVASSY['X']=self.CH2QVASSY['SY']/(self.CH2['I(A)'].max())
            self.CH2QVASSY['L']=self.CH2QVASSY['SY']/(self.CH2['f(Hz)'].max()*self.CH2['I(A)'].max())
        except:
            pass

        
        
    def ASSYall(self,s=0.000000002,K=5):
        self.SyAydataV(s,K)
        self.SyAydataQV(s,K)
        self.SyAydataQVXLCH1()
        self.SyAydataQVXLCH2()
        
    def ASSYallX(self,s=0.000000002,K=5):
        self.SyAydataR(s,K)
        self.SyAydataX(s,K)
        
    def describ(self,ax,key):
        
        ax.annotate(r'$\bf{%s:%0.2f}$ $\bfK$'%('Temperaure',self.T),(0.05,0.9),xycoords='axes fraction',fontsize=25)
        ax.annotate(r'$\bf{%s:%0.2f}$ $\bfdeg$'%('Positions',self.angle),(0.05,0.85),xycoords='axes fraction',fontsize=25)
        ax.annotate(r'$\bf{%s:%0.2f}$ $\bfHz$'%('freq',self.F),(0.05,0.80),xycoords='axes fraction',fontsize=25)
        ax.annotate(r'$\bf{%s:%0.2f}$ $\bfmA$'%('current',self.I*1000),(0.05,0.75),xycoords='axes fraction',fontsize=25)     
        
        if round(CH[i].angle)==0:
            ax.annotate(r'$\bf Transverse$',(0.4,1.01),xycoords='axes fraction',fontsize=22)
        elif round(CH[i].angle)==90:
            ax.annotate(r'$\bfLongitudinal$',(0.4,1.01),xycoords='axes fraction',fontsize=22)     
            
    def plotRN(self,key,ax=1,desrib=1,label=''):
        import jfodata as jfo 
        if ax==1:
            fig, ax = plt.subplots()
        ax.plot(self.CH1VASSY['B(T)'],self.CH1VASSY['SY']/self.CH1VASSY['SY'][jfo.dex(self.CH1VASSY,'B(T)',0)],'o',label=label)            
        ax.set_xlabel('B(T)')
        ax.set_ylabel(r'$\frac{R}{R(B=0T)}$')
        if desrib==1:
            self.describ(ax,key)

    def plotXN(self,key,ax=1,desrib=1,label=''):
        import jfodata as jfo 
        if ax==1:
            fig, ax = plt.subplots()
        ax.plot(self.CH1QVASSY['B(T)'],self.CH1QVASSY['SY']/self.CH1QVASSY['SY'][jfo.dex(self.CH1QVASSY,'B(T)',0)],'o',label=label)            
        ax.set_xlabel('B(T)')
        ax.set_ylabel(r'$\frac{X}{X(B=0T)}$')
        if desrib==1:
            self.describ(ax,key)



"""
Funcoes de classes importnte para analises


"""

from scipy.constants import epsilon_0
from scipy.constants import mu_0
## Analises classe
def indutance(w,L,K,X_0):
    return w*L + K*np.log(w) + K*X_0


def indutanceCilindrical(w,tau,l,k3):
    k1=1/(2*np.pi*epsilon_0)
    k2=mu_0/(2*np.pi)
    return np.log(w*tau)*((k1/(w*l)) - k2*l*w) + k3

class  freqanasis:
    ### Classe for data analyse and plots
    import jfodata as jfo 
    def __init__(self,listdf,CH):
        self.dataXSY=pd.DataFrame(columns=['T(K)','Angle(deg)','Vcontac','croos','I(A)','f','w','B(T)','x','R','x_as','R_as'])
        for jj in np.arange(0,9.1,0.1):
            for i in listdf:
                self.dataXSY.loc[len(self.dataXSY)]=[round(CH[i].T),CH[i].angle,CH[i].vdis,CH[i].croos,CH[i].I,CH[i].F,CH[i].W,jj,CH[i].CH1QVASSY['SY'][jfo.dex(CH[i].CH1QVASSY,'B(T)',jj)],CH[i].CH1VASSY['SY'][jfo.dex(CH[i].CH1VASSY,'B(T)',jj)],CH[i].CH1QVASSY['AS'][jfo.dex(CH[i].CH1QVASSY,'B(T)',jj)],CH[i].CH1VASSY['AS'][jfo.dex(CH[i].CH1VASSY,'B(T)',jj)]]
    def dffilter(self,B=0,Ang=0,fcut=0,fundercut=200,I=0.1e-3,label='OI'):
        df=self.dataXSY[(self.dataXSY["B(T)"]==jj) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut) & (self.dataXSY["f"]<fundercut) & (self.dataXSY["I(A)"]==I)]
        return df
    def plotQVW(self,B=0,Ang=0,fundercut=200,fcut=0,I=0.1e-3,label='OI'):
        df=self.dataXSY[(self.dataXSY["B(T)"]==B) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut) & (self.dataXSY["f"]<fundercut) & (self.dataXSY["I(A)"]==I)]
        plt.plot(df["w"],df["x"]/df["w"],'-o',label=label)
        
        
    def plotZ(self,B=0,Ang=0,fundercut=200,fcut=0,T=2,I=0.1e-4,xaxis="w",label='OI'):
        df=self.dataXSY[(self.dataXSY["B(T)"]==B) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut) & (self.dataXSY["f"]<fundercut) & (self.dataXSY["I(A)"]==I) & (self.dataXSY["T(K)"]==T)]
        plt.plot(df[xaxis],df["x"],'-o',label=label,ms=10,lw=4)
        
    def plotZAS(self,B=0,Ang=0,fundercut=200,fcut=0,T=2,I=0.1e-4,xaxis="w",label='OI'):
        df=self.dataXSY[(self.dataXSY["B(T)"]==B) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut) & (self.dataXSY["f"]<fundercut) & (self.dataXSY["I(A)"]==I) & (self.dataXSY["T(K)"]==T)]
        plt.plot(df[xaxis],df["x_as"],'-o',label=label,ms=10,lw=4)

    def plotZnormal(self,B=0,Ang=0,fundercut=200,fcut=0,I=0.1e-4,xaxis="w",label='OI'):
        df=self.dataXSY[(self.dataXSY["B(T)"]==B) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut) & (self.dataXSY["f"]<fundercut) & (self.dataXSY["I(A)"]==I)]
        plt.plot(df[xaxis],df["x"]/df["x"][jfo.dex(df,"w",0)],'-o',label=label,ms=10,lw=4)
        
    def plotZR(self,B=0,Ang=0,fundercut=200,fcut=0,T=2,I=0.1e-4,xaxis="w",label='OI'):
        df=self.dataXSY[(self.dataXSY["B(T)"]==B) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut) & (self.dataXSY["f"]<fundercut) & (self.dataXSY["I(A)"]==I) & (self.dataXSY["T(K)"]==T)]
        plt.plot(df[xaxis],df["R"],'-o',label=label,ms=10,lw=4)
        
    def plotZRAS(self,B=0,Ang=0,fundercut=200,fcut=0,T=2,I=0.1e-4,xaxis="w",label='OI'):
        df=self.dataXSY[(self.dataXSY["B(T)"]==B) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut) & (self.dataXSY["f"]<fundercut) & (self.dataXSY["I(A)"]==I) & (self.dataXSY["T(K)"]==T)]
        plt.plot(df[xaxis],df["R_as"],'-o',label=label,ms=10,lw=4)
        
    def fitLiner(self,B=0,Ang=0,fcut1=0,fcut=0,label='OI'):
        df=self.dataXSY[(self.dataXSY["B(T)"]==B) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut) ]
        df1=self.dataXSY[(self.dataXSY["B(T)"]==B) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut1) ]

        plt.plot(df1["w"],df1["x"],'o',label=label)
        poly=np.polyfit(df["w"],df["x"],1)
        plt.plot(df["w"],np.polyval(poly,df["w"]),'-r')

        return poly
    
    
    def fitunderL(self,func,B=0,Ang=0,fundercut=200,fcut1=0,fcut=0,label='OI',I=0.1e-3):
        df=self.dataXSY[(self.dataXSY["B(T)"]==B) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut) & (self.dataXSY["f"]<fundercut) & (self.dataXSY["I(A)"]==I)]
        df1=self.dataXSY[(self.dataXSY["B(T)"]==B) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut1) ]

        plt.plot(df1["w"],df1["x"],'-o',label=label,ms=10,lw=4)
        popt,pcov=curve_fit(func,df["w"],df["x"])
        plt.plot(df["w"],func(df["w"],*popt),'r')

        return popt
    
    def fitunderL2(self,B=0,Ang=0,fundercut=200,fcut1=0,fcut=0,label='OI',I=0.1e-3):
        df=self.dataXSY[(self.dataXSY["B(T)"]==B) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut) & (self.dataXSY["f"]<fundercut) & (self.dataXSY["I(A)"]==I)]
        df1=self.dataXSY[(self.dataXSY["B(T)"]==B) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut1) ]

        plt.plot(df1["w"],df1["x"],'-o',label=label,ms=10,lw=4)
        popt,pcov=curve_fit(indutanceCilindrical,df["w"],df["x"])
        plt.plot(df["w"],indutanceCilindrical(df["w"],*popt),'r')

        return popt
    
    
    def fitunderLDF(self,Ang=0,fundercut=200,fcut1=0,fcut=0,label='OI',I=0.1e-3):
        coeffi=pd.DataFrame(columns=['B(T)','L','K','X_0'])
        for jj in np.arange(0,9.1,0.1):
            df=self.dataXSY[(self.dataXSY["B(T)"]==jj) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut) & (self.dataXSY["f"]<fundercut) & (self.dataXSY["I(A)"]==I)]
            df1=self.dataXSY[(self.dataXSY["B(T)"]==jj) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut1) ]
    
            
            popt,pcov=curve_fit(indutance,df["w"],df["x"])
            
            coeffi.loc[len(coeffi)]=[jj,popt[0],popt[1],popt[2]]
        self.coeffi=coeffi
        return coeffi
    
    def fitLinerDf(self,Ang=0,fcut1=0,fcut=0,label='OI'):
        
        X0LB=pd.DataFrame(columns=['T(K)','Angle(deg)','B(T)','Vcontac','croos','I(A)','f','w','X_0','L_B',])
        for jj in np.arange(0,9.1,0.1):
            
            df=self.dataXSY[(self.dataXSY["B(T)"]==jj) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut)]
            self.poly=np.polyfit(df["w"],df["x"],1)
            X0LB.loc[len(X0LB)]=[2,Ang,jj,CH[i].vdis,CH[i].croos,CH[i].I,CH[i].F,CH[i].W,self.poly[1],self.poly[0]]
            
        self.X0LB=X0LB
        return  X0LB

    def plotPuroC(self,B=0,Ang=0,fundercut=200,fcut=0,I=0.1e-3,xaxis="w",label='OI'):
        df=self.dataXSY[(self.dataXSY["B(T)"]==B) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut) & (self.dataXSY["f"]<fundercut) & (self.dataXSY["I(A)"]==I)]
        plt.plot(df[xaxis],df["x"]-self.X0LB['X_0'][jfo.dex(self.X0LB,"B(T)",B)]-(df[xaxis]*self.X0LB['L_B'][jfo.dex(self.X0LB,"B(T)",B)]),'-o',label=label,ms=10,lw=4)

        
    def lmfitodel0(self,Ang=0,fundercut=200,fcut1=0,fcut=0,label='OI',I=0.1e-3):   
        def indutancPower(w,alpha,beta,x_0):
            
            return alpha*(w) - 1/(beta*(w)) + x_0
        ind=Model(indutancPower)
        dfT=pd.DataFrame(columns=['L', 'K', 'X_0'])
        fitmodel={}
        for jj in np.arange(0,9.1,0.1):
            df=self.dataXSY[(self.dataXSY["B(T)"]==jj) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut) & (self.dataXSY["f"]<fundercut) & (self.dataXSY["I(A)"]==I)]
            params = ind.make_params(alpha=1,beta=1,x_0=1)
            fitmodel[jj]=ind.fit(df['x'],params=params,w=df['w'])
            dfT.loc[jj]=fitmodel.values
        return fitmodel,dfT        

        
    def lmfitodel(self,Ang=0,fundercut=200,fcut1=0,fcut=0,label='OI',I=0.1e-3):   
        def indutance(w,L,K,X_0):
            return w*L + K*np.log(w) + X_0
        ind=Model(indutance)
        dfT=pd.DataFrame(columns=['L', 'K', 'X_0'])
        fitmodel={}
      
        for jj in np.arange(0,9.1,0.1):
            df=self.dataXSY[(self.dataXSY["B(T)"]==jj) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut) & (self.dataXSY["f"]<fundercut) & (self.dataXSY["I(A)"]==I)]
            params = ind.make_params(L=-1.4,K=60,X_0=1.5)
            fitmodel[jj]=ind.fit(df['x'],params=params,w=df['w'])
            dfT.loc[jj]=fitmodel.values
        return fitmodel,dfT
    

  
    def lmfitodel2(self,Ang=0,fundercut=200,fcut1=0,fcut=0,label='OI',I=0.1e-3):   
        def indutanceCilindrical(w,tau,l):
            k1=1/(2*np.pi*epsilon_0)
            k2=mu_0/(2*np.pi)
            return np.log(w*tau)*((k1/(w*l)) - k2*l*w)
        ind=Model(indutanceCilindrical)
        dfT=pd.DataFrame(columns=['tau', 'l'])
        fitmodel={}
        for jj in np.arange(0,9.1,0.1):
            df=self.dataXSY[(self.dataXSY["B(T)"]==jj) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut) & (self.dataXSY["f"]<fundercut) & (self.dataXSY["I(A)"]==I)]
            params = ind.make_params(tau=1,l=0.001)
            fitmodel[jj]=ind.fit(df['x'],params=params,w=df['w'])
            dfT.loc[jj]=fitmodel.values
        return fitmodel,dfT
    
    def lmfitodel3(self,Ang=0,fundercut=200,fcut1=0,fcut=0,label='OI',I=0.1e-3):   
        def indutanceCilindrical(w,tau,l):
            k1=1/(2*np.pi*epsilon_0)
            k2=mu_0/(2*np.pi)
            return np.log(w*tau)*((k1/(w*l)) - k2*l*w)
        ind=Model(indutanceCilindrical)
        dfT=pd.DataFrame(columns=['tau', 'l'])
        fitmodel={}
        for jj in np.arange(0,9.1,0.1):
            df=self.dataXSY[(self.dataXSY["B(T)"]==jj) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut) & (self.dataXSY["f"]<fundercut) & (self.dataXSY["I(A)"]==I)]
            params = ind.make_params(tau=1,l=0.001)
            fitmodel[jj]=ind.fit(df['x'],params=params,w=df['w'])
            dfT.loc[jj]=fitmodel.values
        return fitmodel,dfT
    
  
    def lmfitodel4(self,Ang=0,fundercut=200,fcut1=0,fcut=0,label='OI',I=0.1e-3):   
        def indutancPower(w,alpha,beta,a,b):
            
            
            return alpha*(w**a) - beta*(w**b)
        ind=Model(indutancPower)
        dfT=pd.DataFrame(columns=['B(T)',"alpha","beta","a","b"])
        fitmodel={}
        for jj in np.arange(0,9.1,0.1):
            df=self.dataXSY[(self.dataXSY["B(T)"]==jj) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut) & (self.dataXSY["f"]<fundercut) & (self.dataXSY["I(A)"]==I)]
            #popt,pcov=curve_fit(indutancPower,df['w'],df['x'],maxfev = 50000)
            params = ind.make_params(alpha=1,beta=1,a=1/3,b=1)
            params['b'].vary=False
            params['a'].vary=False
            fitmodel[jj]=ind.fit(df['x'],params=params,w=df['w'])
            dic=fitmodel[jj].values 
            dic['B(T)']=jj
            dfT.loc[jj]=dic
            
        return fitmodel,dfT

    def lmfitodelunder(self,Ang=0,fundercut=200,fcut1=0,fcut=0,label='OI',I=0.1e-3):   
        def indutancPower(w,L,W_0):
            
            
            return w*L/(1+((w**2)/(W_0**2)))
        ind=Model(indutancPower)
        dfT=pd.DataFrame(columns=['B(T)',"L","W_0"])
        fitmodel={}
        for jj in np.arange(0,9.1,0.1):
            df=self.dataXSY[(self.dataXSY["B(T)"]==jj) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut) & (self.dataXSY["f"]<fundercut) & (self.dataXSY["I(A)"]==I)]
            #popt,pcov=curve_fit(indutancPower,df['w'],df['x'],maxfev = 50000)
            params = ind.make_params(L=1,W_0=1)

            fitmodel[jj]=ind.fit(df['x'],params=params,w=df['w'])
            dic=fitmodel[jj].values 
            dic['B(T)']=jj
            dfT.loc[jj]=dic
            
        return fitmodel,dfT
    
    
    def lmfitCarste1(self,Ang=0,fundercut=200,fcut1=0,fcut=0,label='OI',I=0.1e-3):   
        def indutancPower(w,L,W_0):
            
            
            return w*L/(1+((w**2)*(W_0**2)))
        ind=Model(indutancPower)
        dfT=pd.DataFrame(columns=['B(T)',"L","W_0"])
        fitmodel={}
        for jj in np.arange(0,9.1,0.1):
            df=self.dataXSY[(self.dataXSY["B(T)"]==jj) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut) & (self.dataXSY["f"]<fundercut) & (self.dataXSY["I(A)"]==I)]
            #popt,pcov=curve_fit(indutancPower,df['w'],df['x'],maxfev = 50000)
            params = ind.make_params(L=1,W_0=1)

            fitmodel[jj]=ind.fit(df['x'],params=params,w=df['w'])
            dic=fitmodel[jj].values 
            dic['B(T)']=jj
            dfT.loc[jj]=dic
            
        return fitmodel,dfT

    def lmfitmodes(self,Ang=0,fundercut=200,fcut1=0,fcut=0,label='OI',I=0.1e-3):   
        def indutancmodes(w,F_n,rho,alpha):
            
            
            return w*F_n/(1+((rho*F_n*w)**(2*alpha)))
        
        ind=Model(indutancmodes,prefix="1")
        ind2=Model(indutancmodes,prefix="2")
        dfT=pd.DataFrame(columns=['B(T)',"F_n","rho","alpha"])
        fitmodel={}
        for jj in np.arange(0,9.1,0.1):
            df=self.dataXSY[(self.dataXSY["B(T)"]==jj) & (self.dataXSY["Angle(deg)"]==Ang) & (self.dataXSY["f"]>fcut) & (self.dataXSY["f"]<fundercut) & (self.dataXSY["I(A)"]==I)]
            #popt,pcov=curve_fit(indutancPower,df['w'],df['x'],maxfev = 50000)
            params = ind.make_params(rho=0.04857153,F_n=0.92808498,alpha=0.8)

            fitmodel[jj]=ind.fit(df['x'],params=params,w=df['w'])
            dic=fitmodel[jj].values 
            dic['B(T)']=jj
            dfT.loc[jj]=dic
            
        return fitmodel,dfT
        
class plot2dnharmoci:

    def __init__(self,CH):

        self.CH=CH
        
    def describ(self,ax,key):
        
        ax.annotate(r'$\bf{%s:%0.2f}$ $\bfK$'%('Temperaure',self.CH[key].T),(0.05,0.9),xycoords='axes fraction',fontsize=25)
        ax.annotate(r'$\bf{%s:%0.2f}$ $\bfdeg$'%('Positions',self.CH[key].angle),(0.05,0.85),xycoords='axes fraction',fontsize=25)
        ax.annotate(r'$\bf{%s:%0.2f}$ $\bfHz$'%('freq',self.CH[key].F),(0.05,0.80),xycoords='axes fraction',fontsize=25)
        ax.annotate(r'$\bf{%s:%0.2f}$ $\bfmA$'%('current',self.CH[key].I*1000),(0.05,0.75),xycoords='axes fraction',fontsize=25)     
        
        if round(CH[i].angle)==0:
            ax.annotate(r'$\bf Transverse$',(0.4,1.01),xycoords='axes fraction',fontsize=22)
        elif round(CH[i].angle)==90:
            ax.annotate(r'$\bfLongitudinal$',(0.4,1.01),xycoords='axes fraction',fontsize=22)
        
    def realdataR(self,key,ax=1,desrib=1):
        if ax==1:
            fig, ax = plt.subplots()
            
        ax.plot(self.CH[key].CH1['B(T)'],self.CH[key].CH1['R(ohms)'],'ok',label='Data')
        ax.plot(self.CH[key].CH1VASSY['B(T)'],self.CH[key].CH1VASSY['R(ohms)'],'--r',lw=4,label='Data interpolate')
        ax.plot(self.CH[key].CH1VASSY['B(T)'],self.CH[key].CH1VASSY['SY'],'g--',lw=4,label='Data interpolate SY part')  
        #ax.plot(self.CH[key].CH1QVASSY['B(T)'],self.CH[key].CH1VASSY['AS'],'b--')  
        ax.set_xlabel('B(T)')
        ax.set_ylabel(r'$R(\Omega)$')
        if desrib==1:
            self.describ(ax,key)
        ax.legend(fontsize=18,loc='lower left')
        
    def realdataX(self,key,ax=1,desrib=1):
        if ax==1:
            fig, ax = plt.subplots()
            
        ax.plot(self.CH[key].CH1['B(T)'],self.CH[key].CH1['X'],'ok',label='Data')
        ax.plot(self.CH[key].CH1QVASSY['B(T)'],self.CH[key].CH1QVASSY['X'],'--r',lw=4,label='Data interpolate')
        ax.plot(self.CH[key].CH1QVASSY['B(T)'],self.CH[key].CH1QVASSY['SY'],'g--',lw=4,label='Data interpolate SY part')  
        ax.plot(self.CH[key].CH1QVASSY['B(T)'],self.CH[key].CH1QVASSY['AS'],'b--',label='Data interpolate AS part')  
        ax.set_xlabel('B(T)')
        ax.set_ylabel(r'$X(\Omega)$')
        if desrib==1:
            self.describ(ax,key)
        ax.legend(fontsize=18,loc='lower left')
        
        
    def plotX(self,key,ax=1,desrib=1,label=''):
        if ax==1:
            fig, ax = plt.subplots()
        ax.plot(self.CH[key].CH1QVASSY['B(T)'],self.CH[key].CH1QVASSY['SY'],'o',label=label)     
        ax.set_xlabel('B(T)')
        ax.set_ylabel(r'$X(\Omega)$')
        if desrib==1:
            self.describ(ax,key)
        
        

    def plotR(self,key,ax=1,desrib=1,label=''):
        if ax==1:
            fig, ax = plt.subplots()
        print(key)
        ax.plot(self.CH[key].CH1VASSY['B(T)'],self.CH[key].CH1VASSY['SY'],'o',label=label)            
        ax.set_xlabel('B(T)')
        ax.set_ylabel(r'R($\Omega$)')
        if desrib==1:
            self.describ(ax,key)
           
        
    def plot2nd(self,key,ax=1,desrib=1,label='',**kwargs):
        if ax==1:
            fig, ax = plt.subplots()
        ax.plot(self.CH[key].CH1['B(T)'],self.CH[key].CH1['2 Harm'],'o',label=label,**kwargs)        
        ax.set_xlabel('B(T)')
        ax.set_ylabel(r'2 Harmic (Db)')
        if desrib==1:
            self.describ(ax,key)
        
        
    def plotRdiff(self,key,ax=1,desrib=1,label='',**kwargs):
        if ax==1:
            fig, ax = plt.subplots()
        ax.plot(self.CH[key].CH1QVASSY['B(T)'][:-1],np.diff(self.CH[key].CH1VASSY['SY']),'o',label=label,**kwargs) 
        ax.set_xlabel('B(T)')
        ax.set_ylabel(r'$\frac{dR}{dB}$')
        if desrib==1:
            self.describ(ax,key)
            
       
        
    def plotXdiff(self,key,ax=1,desrib=1,label='',**kwargs):
        if ax==1:
            fig, ax = plt.subplots()
        ax.plot(self.CH[key].CH1QVASSY['B(T)'][:-1],np.diff(self.CH[key].CH1QVASSY['SY']),'o',label=label,**kwargs) 
        ax.set_xlabel('B(T)')
        ax.set_ylabel(r'$\frac{dX}{dB}$')
        if desrib==1:
            self.describ(ax,key)
        
    def plot2nddiff(self,key,ax=1,desrib=1,label=''):
        if ax==1:
            fig, ax = plt.subplots()
        ax.plot(self.CH[key].CH1['B(T)'][:-1],np.diff(self.CH[key].CH1['2 Harm']),'o',label=label) 
        ax.set_xlabel('B(T)')
        ax.set_ylabel(r'$\frac{d 2andH}{dB}$')
        if desrib==1:
            self.describ(ax,key)
        
    def tresfig3(self,key,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):
    
        fig = plt.figure(figsize=(36,24))
        gs0 = gridspec.GridSpec(2, 3, figure=fig,wspace=0.3,hspace=.30)
        
        ax0=fig.add_subplot(gs0[0,0])
        ax1=fig.add_subplot(gs0[0,1])
        ax2=fig.add_subplot(gs0[0,2])
        
        ax3=fig.add_subplot(gs0[1,0])
        ax4=fig.add_subplot(gs0[1,1])
        ax5=fig.add_subplot(gs0[1,2])

        
      
        
        
        self.plotR(key,ax0)
        self.plotX(key,ax1,desrib=0)
        self.plot2nd(key,ax2,desrib=0)
        
        self.plotRdiff(key,ax3,desrib=0)
        self.plotXdiff(key,ax4,desrib=0)
        self.plot2nddiff(key,ax5,desrib=0)
        
        
        
        plt.savefig(path + '\%s_T_%0.2f_K_ang_%0.2f_deg.jpg'%(key,self.CH[key].T,self.CH[key].angle))
        plt.show()
        
        
    def tresfig4(self,key,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):
    
        fig = plt.figure(figsize=(36,36))
        gs0 = gridspec.GridSpec(3, 3, figure=fig,wspace=0.3,hspace=.30)
        
        ax0=fig.add_subplot(gs0[0,0])
        ax1=fig.add_subplot(gs0[0,1])
        ax2=fig.add_subplot(gs0[0,2])
        
        ax3=fig.add_subplot(gs0[1,0])
        ax4=fig.add_subplot(gs0[1,1])
        ax5=fig.add_subplot(gs0[1,2])


        ax6=fig.add_subplot(gs0[2,0])
        ax7=fig.add_subplot(gs0[2,1])
        ax8=fig.add_subplot(gs0[2,2])
      
        
        self.realdataR(key,ax0)
        self.realdataX(key,ax1,desrib=0)
        
        
        self.plotR(key,ax3,desrib=0)
        self.plotX(key,ax4,desrib=0)
        self.plot2nd(key,ax5,desrib=0)
        
        self.plotRdiff(key,ax6,desrib=0)
        self.plotXdiff(key,ax7,desrib=0)
        self.plot2nddiff(key,ax8,desrib=0)
        

        
        plt.savefig(path + '\%s_T_%0.2f_K_ang_%0.2f_deg.jpg'%(key,self.CH[key].T,self.CH[key].angle))
        plt.show()
        
    def diffchart(self,listdf,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):
    
        fig = plt.figure(figsize=(36,32))
        gs0 = gridspec.GridSpec(2, 2, figure=fig,wspace=0.2,hspace=.20)
        
        ax0=fig.add_subplot(gs0[0,0])
        ax1=fig.add_subplot(gs0[0,1])
        ax2=fig.add_subplot(gs0[1,0])
        ax3=fig.add_subplot(gs0[1,1])
        
        for i in listdf:
            if CH[i].angle==0 and CH[i].I==0.0001:
                self.plotRdiff(i,ax0,desrib=0,label='%0.2f'%CH[i].F) 
                ax0.set_title('Transversal')
                ax0.legend()
        for i in listdf:
            if CH[i].angle==90 and CH[i].I==0.0001:
                self.plotRdiff(i,ax1,desrib=0) 
                ax1.set_title('logitudinal')
            
        for i in listdf:
            if CH[i].angle==0 and  CH[i].I==0.0001:
                self.plotXdiff(i,ax2,desrib=0) 
                
        for i in listdf:
            if CH[i].angle==90 and CH[i].I==0.0001:
                self.plotXdiff(i,ax3,desrib=0) 
                
    def RXchart(self,listdf,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):
    
        fig = plt.figure(figsize=(36,32))
        gs0 = gridspec.GridSpec(2, 2, figure=fig,wspace=0.2,hspace=.20)
        
        ax0=fig.add_subplot(gs0[0,0])
        ax1=fig.add_subplot(gs0[0,1])
        ax2=fig.add_subplot(gs0[1,0])
        ax3=fig.add_subplot(gs0[1,1])
        
        for i in listdf:
            if CH[i].angle==0 and  CH[i].I==0.0001:
                self.plotR(i,ax0,desrib=0,label='%0.2f'%CH[i].F) 
                ax0.set_title('Transversal')
                ax0.legend()
        for i in listdf:
            if CH[i].angle==90 and  CH[i].I==0.0001:
                self.plotR(i,ax1,desrib=0) 
                ax1.set_title('logitudinal')
            
        for i in listdf:
            if CH[i].angle==0 and  CH[i].I==0.0001:
                self.plotX(i,ax2,desrib=0) 
                
        for i in listdf:
            if CH[i].angle==90 and  CH[i].I==0.0001:
                self.plotX(i,ax3,desrib=0) 
                
                
    def Ichart(self,listdf,angle=0,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):
    
        fig = plt.figure(figsize=(36,24))
        gs0 = gridspec.GridSpec(2, 3, figure=fig,wspace=0.2,hspace=.20)
        
        ax0=fig.add_subplot(gs0[0,0])
        ax1=fig.add_subplot(gs0[0,1])
        ax2=fig.add_subplot(gs0[0,2])
        ax3=fig.add_subplot(gs0[1,0])
        ax4=fig.add_subplot(gs0[1,1])
        ax5=fig.add_subplot(gs0[1,2])
        
        Ilist=list(set([self.CH[i].I for i in listdf]))
        for i in listdf:
            if CH[i].angle==angle and Ilist[0]==CH[i].I:
                self.plotR(i,ax0,desrib=0,label='%0.2f'%CH[i].F) 
                if angle==0:
                    ax0.set_title('Transversel I=%0.2f mA'%(CH[i].I*1000))
                else:
                    ax0.set_title('Longitudinal I=%0.2f mA'%(CH[i].I*1000))
                ax0.legend()
        
        for i in listdf:
            if CH[i].angle==angle and Ilist[1]==CH[i].I:
                self.plotR(i,ax1,desrib=0,label='%0.2f'%CH[i].F) 
                if angle==0:
                    ax1.set_title('Transversel I=%0.2f mA'%(CH[i].I*1000))
                else:
                    ax1.set_title('Longitudinal I=%0.2f mA'%(CH[i].I*1000))
                ax1.legend()

                
        for i in listdf:
            if CH[i].angle==angle and Ilist[2]==CH[i].I:
                self.plotR(i,ax2,desrib=0,label='%0.2f'%CH[i].F) 
                if angle==0:
                    ax2.set_title('Transversel I=%0.2f mA'%(CH[i].I*1000))
                else:
                    ax2.set_title('Longitudinal I=%0.2f mA'%(CH[i].I*1000))
                ax2.legend()

        for i in listdf:
            if CH[i].angle==angle and Ilist[0]==CH[i].I:
                self.plotX(i,ax3,desrib=0,label='%0.2f'%CH[i].F) 
                if angle==0:
                    ax3.set_title('Transversel I=%0.2f mA'%(CH[i].I*1000))
                else:
                    ax3.set_title('Longitudinal I=%0.2f mA'%(CH[i].I*1000))
                ax0.legend()      
                
        for i in listdf:
            if CH[i].angle==angle and Ilist[1]==CH[i].I:
                self.plotX(i,ax4,desrib=0,label='%0.2f'%CH[i].F) 
                if angle==0:
                    ax4.set_title('Transversel I=%0.2f mA'%(CH[i].I*1000))
                else:
                    ax4.set_title('Longitudinal I=%0.2f mA'%(CH[i].I*1000))
                ax0.legend() 
                
        for i in listdf:
            if CH[i].angle==angle and Ilist[2]==CH[i].I:
                self.plotX(i,ax5,desrib=0,label='%0.2f'%CH[i].F) 
                if angle==0:
                    ax5.set_title('Transversel I=%0.2f mA'%(CH[i].I*1000))
                else:
                    ax5.set_title('Longitudinal I=%0.2f mA'%(CH[i].I*1000))
                ax0.legend() 
                
    def Idrchart(self,listdf,angle=0,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):
    
        fig = plt.figure(figsize=(36,24))
        gs0 = gridspec.GridSpec(2, 3, figure=fig,wspace=0.2,hspace=.20)
        
        ax0=fig.add_subplot(gs0[0,0])
        ax1=fig.add_subplot(gs0[0,1])
        ax2=fig.add_subplot(gs0[0,2])
        ax3=fig.add_subplot(gs0[1,0])
        ax4=fig.add_subplot(gs0[1,1])
        ax5=fig.add_subplot(gs0[1,2])
        
        Ilist=list(set([self.CH[i].I for i in listdf]))
        for i in listdf:
            if CH[i].angle==angle and Ilist[0]==CH[i].I:
                self.plotRdiff(i,ax0,desrib=0,label='%0.2f'%CH[i].F) 
                if angle==0:
                    ax0.set_title('Transversel I=%0.2f mA'%(CH[i].I*1000))
                else:
                    ax0.set_title('Longitudinal I=%0.2f mA'%(CH[i].I*1000))
                ax0.legend()
        
        for i in listdf:
            if CH[i].angle==angle and Ilist[1]==CH[i].I:
                self.plotRdiff(i,ax1,desrib=0,label='%0.2f'%CH[i].F) 
                if angle==0:
                    ax1.set_title('Transversel I=%0.2f mA'%(CH[i].I*1000))
                else:
                    ax1.set_title('Longitudinal I=%0.2f mA'%(CH[i].I*1000))
                ax1.legend()

                
        for i in listdf:
            if CH[i].angle==angle and Ilist[2]==CH[i].I:
                self.plotRdiff(i,ax2,desrib=0,label='%0.2f'%CH[i].F) 
                if angle==0:
                    ax2.set_title('Transversel I=%0.2f mA'%(CH[i].I*1000))
                else:
                    ax2.set_title('Longitudinal I=%0.2f mA'%(CH[i].I*1000))
                ax2.legend()

        for i in listdf:
            if CH[i].angle==angle and Ilist[0]==CH[i].I:
                self.plotXdiff(i,ax3,desrib=0,label='%0.2f'%CH[i].F) 
                if angle==0:
                    ax3.set_title('Transversel I=%0.2f mA'%(CH[i].I*1000))
                else:
                    ax3.set_title('Longitudinal I=%0.2f mA'%(CH[i].I*1000))
                ax0.legend()      
                
        for i in listdf:
            if CH[i].angle==angle and Ilist[1]==CH[i].I:
                self.plotXdiff(i,ax4,desrib=0,label='%0.2f'%CH[i].F) 
                if angle==0:
                    ax4.set_title('Transversel I=%0.2f mA'%(CH[i].I*1000))
                else:
                    ax4.set_title('Longitudinal I=%0.2f mA'%(CH[i].I*1000))
                ax0.legend() 
                
        for i in listdf:
            if CH[i].angle==angle and Ilist[2]==CH[i].I:
                self.plotXdiff(i,ax5,desrib=0,label='%0.2f'%CH[i].F) 
                if angle==0:
                    ax5.set_title('Transversel I=%0.2f mA'%(CH[i].I*1000))
                else:
                    ax5.set_title('Longitudinal I=%0.2f mA'%(CH[i].I*1000))
                ax0.legend() 
                
    def Fchart(self,listdf,angle=0,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):
    
        fig = plt.figure(figsize=(36,24))
        gs0 = gridspec.GridSpec(2, 3, figure=fig,wspace=0.2,hspace=.20)
        
        ax0=fig.add_subplot(gs0[0,0])
        ax1=fig.add_subplot(gs0[0,1])
        ax2=fig.add_subplot(gs0[0,2])
        ax3=fig.add_subplot(gs0[1,0])
        ax4=fig.add_subplot(gs0[1,1])
        ax5=fig.add_subplot(gs0[1,2])
        
        Ilist=list(set([self.CH[i].F for i in listdf]))
        Ilist=sorted(Ilist)
        
        for i in listdf:
            if CH[i].angle==angle and Ilist[0]==CH[i].F:
                self.plotR(i,ax0,desrib=0,label='%0.2f'%(CH[i].I*1000)) 
                if angle==0:
                    ax0.set_title('Transversel F=%0.2f Hz'%(CH[i].F))
                else:
                    ax0.set_title('Longitudinal F=%0.2f Hz'%(CH[i].F))
                ax0.legend()
        
        for i in listdf:
            if CH[i].angle==angle and Ilist[1]==CH[i].F:
                self.plotR(i,ax1,desrib=0,label='%0.2f'%(CH[i].I*1000)) 
                if angle==0:
                    ax1.set_title('Transversel F=%0.2f Hz'%(CH[i].F))
                else:
                    ax1.set_title('Longitudinal F=%0.2f Hz'%(CH[i].F))
                ax1.legend()

                
        for i in listdf:
            if CH[i].angle==angle and Ilist[2]==CH[i].F:
                self.plotR(i,ax2,desrib=0,label='%0.2f'%(CH[i].I*1000)) 
                if angle==0:
                    ax2.set_title('Transversel F=%0.2f Hz'%(CH[i].F))
                else:
                    ax2.set_title('Longitudinal F=%0.2f Hz'%(CH[i].F))
                ax2.legend()

        for i in listdf:
            if CH[i].angle==angle and Ilist[0]==CH[i].F:
                self.plotX(i,ax3,desrib=0,label='%0.2f'%(CH[i].I*1000)) 
                if angle==0:
                    ax3.set_title('Transversel F=%0.2f Hz'%(CH[i].F))
                else:
                    ax3.set_title('Longitudinal F=%0.2f Hz'%(CH[i].F))
                ax0.legend()      
                
        for i in listdf:
            if CH[i].angle==angle and Ilist[1]==CH[i].F:
                self.plotX(i,ax4,desrib=0,label='%0.2f'%(CH[i].I*1000)) 
                if angle==0:
                    ax4.set_title('Transversel I=%0.2f Hz'%(CH[i].F))
                else:
                    ax4.set_title('Longitudinal I=%0.2f Hz'%(CH[i].F))
                ax0.legend() 
                
        for i in listdf:
            if CH[i].angle==angle and Ilist[2]==CH[i].F:
                self.plotX(i,ax5,desrib=0,label='%0.2f'%(CH[i].I*1000)) 
                if angle==0:
                    ax5.set_title('Transversel F=%0.2f Hz'%(CH[i].F))
                else:
                    ax5.set_title('Longitudinal F=%0.2f Hz'%(CH[i].F))
                ax0.legend() 
                
                
    def Fdchart(self,listdf,angle=0,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):
    
        fig = plt.figure(figsize=(36,24))
        gs0 = gridspec.GridSpec(2, 3, figure=fig,wspace=0.2,hspace=.20)
        
        ax0=fig.add_subplot(gs0[0,0])
        ax1=fig.add_subplot(gs0[0,1])
        ax2=fig.add_subplot(gs0[0,2])
        ax3=fig.add_subplot(gs0[1,0])
        ax4=fig.add_subplot(gs0[1,1])
        ax5=fig.add_subplot(gs0[1,2])
        
        Ilist=list(set([self.CH[i].F for i in listdf]))
        Ilist=sorted(Ilist)
        
        for i in listdf:
            if CH[i].angle==angle and Ilist[0]==CH[i].F:
                self.plotRdiff(i,ax0,desrib=0,label='%0.2f'%(CH[i].I*1000)) 
                if angle==0:
                    ax0.set_title('Transversel I=%0.2f Hz'%(CH[i].F))
                else:
                    ax0.set_title('Longitudinal I=%0.2f Hz'%(CH[i].F))
                ax0.legend()
        
        for i in listdf:
            if CH[i].angle==angle and Ilist[1]==CH[i].F:
                self.plotRdiff(i,ax1,desrib=0,label='%0.2f'%(CH[i].I*1000)) 
                if angle==0:
                    ax1.set_title('Transversel I=%0.2f Hz'%(CH[i].F))
                else:
                    ax1.set_title('Longitudinal I=%0.2f Hz'%(CH[i].F))
                ax1.legend()

                
        for i in listdf:
            if CH[i].angle==angle and Ilist[2]==CH[i].F:
                self.plotRdiff(i,ax2,desrib=0,label='%0.2f'%(CH[i].I*1000)) 
                if angle==0:
                    ax2.set_title('Transversel I=%0.2f Hz'%(CH[i].F))
                else:
                    ax2.set_title('Longitudinal I=%0.2f Hz'%(CH[i].F))
                ax2.legend()

        for i in listdf:
            if CH[i].angle==angle and Ilist[0]==CH[i].F:
                self.plotXdiff(i,ax3,desrib=0,label='%0.2f'%(CH[i].I*1000)) 
                if angle==0:
                    ax3.set_title('Transversel I=%0.2f Hz'%(CH[i].F))
                else:
                    ax3.set_title('Longitudinal I=%0.2f Hz'%(CH[i].F))
                ax0.legend()      
                
        for i in listdf:
            if CH[i].angle==angle and Ilist[1]==CH[i].F:
                self.plotXdiff(i,ax4,desrib=0,label='%0.2f'%(CH[i].I*1000)) 
                if angle==0:
                    ax4.set_title('Transversel I=%0.2f Hz'%(CH[i].F))
                else:
                    ax4.set_title('Longitudinal I=%0.2f Hz'%(CH[i].F))
                ax0.legend() 
                
        for i in listdf:
            if CH[i].angle==angle and Ilist[2]==CH[i].F:
                self.plotXdiff(i,ax5,desrib=0,label='%0.2f'%(CH[i].I*1000)) 
                if angle==0:
                    ax5.set_title('Transversel I=%0.2f Hz'%(CH[i].F))
                else:
                    ax5.set_title('Longitudinal I=%0.2f Hz'%(CH[i].F))
                ax0.legend() 
                
                
    def Fanglechart(self,listdf,angle=0,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):    
        fig = plt.figure(figsize=(36,24))
        gs0 = gridspec.GridSpec(2, 3, figure=fig,wspace=0.2,hspace=.20)
        
        ax0=fig.add_subplot(gs0[0,0])
        ax1=fig.add_subplot(gs0[0,1])
        ax2=fig.add_subplot(gs0[0,2])
        ax3=fig.add_subplot(gs0[1,0])
        ax4=fig.add_subplot(gs0[1,1])
        ax5=fig.add_subplot(gs0[1,2])
        
        Ilist=list(set([self.CH[i].F for i in listdf]))
        Ilist=sorted(Ilist)
        
        
        for i in listdf:
            if CH[i].angle==0 and Ilist[0]==CH[i].F:
                self.plotRdiff(i,ax0,desrib=0,color='red',label='%0.2f'%(CH[i].I*1000)) 
            if CH[i].angle==90 and Ilist[0]==CH[i].F:
                self.plotRdiff(i,ax0,desrib=0,color='black',label='%0.2f'%(CH[i].I*1000)) 
                ax0.set_title('F=%0.2f Hz'%(CH[i].F))
                
        for i in listdf:
            if CH[i].angle==0 and Ilist[1]==CH[i].F:
                self.plotRdiff(i,ax1,desrib=0,color='red',label='%0.2f'%(CH[i].I*1000)) 
            if CH[i].angle==90 and Ilist[1]==CH[i].F:
                self.plotRdiff(i,ax1,desrib=0,color='black',label='%0.2f'%(CH[i].I*1000)) 
                ax0.set_title('F=%0.2f Hz'%(CH[i].F))
                
        for i in listdf:
            if CH[i].angle==0 and Ilist[2]==CH[i].F:
                self.plotRdiff(i,ax2,desrib=0,color='red',label='%0.2f'%(CH[i].I*1000)) 
            if CH[i].angle==90 and Ilist[2]==CH[i].F:
                self.plotRdiff(i,ax2,desrib=0,color='black',label='%0.2f'%(CH[i].I*1000)) 
                ax0.set_title('F=%0.2f Hz'%(CH[i].F))
                
                
        for i in listdf:
            if CH[i].angle==0 and Ilist[0]==CH[i].F:
                self.plotRdiff(i,ax0,desrib=0,color='red',label='%0.2f'%(CH[i].I*1000)) 
            if CH[i].angle==90 and Ilist[0]==CH[i].F:
                self.plotRdiff(i,ax0,desrib=0,color='black',label='%0.2f'%(CH[i].I*1000)) 
                ax0.set_title('F=%0.2f Hz'%(CH[i].F))
                
        for i in listdf:
            if CH[i].angle==0 and Ilist[1]==CH[i].F:
                self.plotRdiff(i,ax1,desrib=0,color='red',label='%0.2f'%(CH[i].I*1000)) 
            if CH[i].angle==90 and Ilist[1]==CH[i].F:
                self.plotRdiff(i,ax1,desrib=0,color='black',label='%0.2f'%(CH[i].I*1000)) 
                ax1.set_title('F=%0.2f Hz'%(CH[i].F))
                
        for i in listdf:
            if CH[i].angle==0 and Ilist[2]==CH[i].F:
                self.plotRdiff(i,ax2,desrib=0,color='red',label='%0.2f'%(CH[i].I*1000)) 
            if CH[i].angle==90 and Ilist[2]==CH[i].F:
                self.plotRdiff(i,ax2,desrib=0,color='black',label='%0.2f'%(CH[i].I*1000)) 
                ax3.set_title('F=%0.2f Hz'%(CH[i].F))
                
                
        for i in listdf:
            if CH[i].angle==0 and Ilist[0]==CH[i].F:
                self.plotXdiff(i,ax3,desrib=0,color='red',label='%0.2f'%(CH[i].I*1000)) 
            if CH[i].angle==90 and Ilist[0]==CH[i].F:
                self.plotXdiff(i,ax3,desrib=0,color='black',label='%0.2f'%(CH[i].I*1000)) 
                ax3.set_title('F=%0.2f Hz'%(CH[i].F))
                
        for i in listdf:
            if CH[i].angle==0 and Ilist[1]==CH[i].F:
                self.plotXdiff(i,ax4,desrib=0,color='red',label='%0.2f'%(CH[i].I*1000)) 
            if CH[i].angle==90 and Ilist[1]==CH[i].F:
                self.plotXdiff(i,ax5,desrib=0,color='black',label='%0.2f'%(CH[i].I*1000)) 
                ax4.set_title('F=%0.2f Hz'%(CH[i].F))
                
        for i in listdf:
            if CH[i].angle==0 and Ilist[2]==CH[i].F:
                self.plotXdiff(i,ax5,desrib=0,color='red',label='%0.2f'%(CH[i].I*1000)) 
            if CH[i].angle==90 and Ilist[2]==CH[i].F:
                self.plotXdiff(i,ax5,desrib=0,color='black',label='%0.2f'%(CH[i].I*1000)) 
                ax5.set_title('F=%0.2f Hz'%(CH[i].F))

class Datafit:
    ### Classe for data analyse and fit quadrature voltage and save the data.
    import jfodata as jfo 
    def __init__(self,listdf,CH):
        self.df={}
        self.intparam={}
        self.inter={}
        for i in listdf:
            self.df[i]=CH[i].CH1
        self.listdf=listdf
        for i in listdf:
            self.df[i]=CH[i].CH1
        for i in listdf:
            self.intparam[i]=[1,3]
    
    def plotX(self,key,ax=1):
        if ax==1:
            fig, ax = plt.subplots()
        
        ax.plot(self.df[key]['B(T)'],self.df[key]['X'],'ok')
        ax.set_xlabel('B(T)')
        ax.set_ylabel('X($\Omega$)')
        
    def plotR(self,key,ax=1):
        if ax==1:
            fig, ax = plt.subplots()
        
        ax.plot(self.df[key]['B(T)'],self.df[key]['R(ohms)'],'ok')
        ax.set_xlabel('B(T)')
        ax.set_ylabel('R($\Omega$)')
        
    def interpo(self,key,s=1,K=3,Y='X'):
        Hmax=self.df[key]['B(T)'].max()
        self.inter[key]=jfo.pdinter(self.df[key],'B(T)',Y,-Hmax,Hmax,npoints=200,smoothing=s,K=K)
        self.intparam[i]=[s,K]
    def plotinter(self,key,ax=1):
        if ax==1:
            fig, ax = plt.subplots()
        ax.plot(self.inter[key][1],self.inter[key][0](self.inter[key][1]),'-.r')
        
    def interfinal(self,key,s=1,K=3):
        fig, ax = plt.subplots()
        self.interpo(key,s,K)
        self.plotX(key,ax)
        self.plotinter(key,ax)
        plt.show()
    
    def datasave(self,path):
        np.save(path+'\intparam.npy',self.intparam)
        
    def dataload(self,path):
        self.intparam=np.load(path+'\intparam.npy',allow_pickle='TRUE').item()     

class Analiperpara(plot2dnharmoci):
    import matplotlib.gridspec as gridspec
    def __init__(self,CH):

        self.CH=CH
        
    def plotRN(self,key,ax=1,desrib=1,label=''):
        if ax==1:
            fig, ax = plt.subplots()
        ax.plot(self.CH[key].CH1VASSY['B(T)'],self.CH[key].CH1VASSY['SY']/self.CH[key].CH1VASSY['SY'][jfo.dex(self.CH[key].CH1VASSY,'B(T)',0)],'o',label=label)            
        ax.set_xlabel('B(T)')
        ax.set_ylabel(r'$\frac{R}{R(B=0T)}$')
        if desrib==1:
            self.describ(ax,key)

    def plotRNAS(self,key,ax=1,desrib=1,label=''):
        if ax==1:
            fig, ax = plt.subplots()
        ax.plot(self.CH[key].CH1VASSY['B(T)'],self.CH[key].CH1VASSY['AS']/self.CH[key].CH1VASSY['AS'][jfo.dex(self.CH[key].CH1VASSY,'B(T)',0)],'o',label=label)            
        ax.set_xlabel('B(T)')
        ax.set_ylabel(r'$\frac{R_AS}{R(B=0T)}$')
        if desrib==1:
            self.describ(ax,key)
    
    def plotRAS(self,key,ax=1,desrib=1,label=''):
        if ax==1:
            fig, ax = plt.subplots()
        ax.plot(self.CH[key].CH1VASSY['B(T)'],self.CH[key].CH1VASSY['AS'],'o',label=label)            
        ax.set_xlabel('B(T)')
        ax.set_ylabel(r'$\frac{R_AS}{R(B=0T)}$')
        if desrib==1:
            self.describ(ax,key)            
            
    def plotXN(self,key,ax=1,desrib=1,label=''):
        if ax==1:
            fig, ax = plt.subplots()
        ax.plot(self.CH[key].CH1QVASSY['B(T)'],self.CH[key].CH1QVASSY['SY']/self.CH[key].CH1QVASSY['SY'][jfo.dex(self.CH[key].CH1QVASSY,'B(T)',0)],'o',label=label)            
        ax.set_xlabel('B(T)')
        ax.set_ylabel(r'$\frac{X}{X(B=0T)}$')
        if desrib==1:
            self.describ(ax,key)
            
            
    def plotXNAS(self,key,ax=1,desrib=1,label=''):
        if ax==1:
            fig, ax = plt.subplots()
        ax.plot(self.CH[key].CH1QVASSY['B(T)'],self.CH[key].CH1QVASSY['AS']/self.CH[key].CH1QVASSY['AS'][jfo.dex(self.CH[key].CH1QVASSY,'B(T)',0)],'o',label=label)            
        ax.set_xlabel('B(T)')
        ax.set_ylabel(r'$\frac{X_AS}{X(B=0T)}$')
        if desrib==1:
            self.describ(ax,key)
    def plotXAS(self,key,ax=1,desrib=1,label=''):
        if ax==1:
            fig, ax = plt.subplots()
        ax.plot(self.CH[key].CH1QVASSY['B(T)'],self.CH[key].CH1QVASSY['AS'],'o',label=label)            
        ax.set_xlabel('B(T)')
        ax.set_ylabel(r'$\frac{X_AS}{X(B=0T)}$')
        if desrib==1:
            self.describ(ax,key)
            
    def plotNRdiff(self,key,ax=1,desrib=1,label='',**kwargs):
        if ax==1:
            fig, ax = plt.subplots()
        ax.plot(self.CH[key].CH1QVASSY['B(T)'][:-1],np.diff(self.CH[key].CH1VASSY['SY']/self.CH[key].CH1VASSY['SY'][jfo.dex(self.CH[key].CH1VASSY,'B(T)',0)]),'o',label=label,**kwargs) 
        ax.set_xlabel('B(T)')
        ax.set_ylabel(r'$\frac{d(R/R(B=0))}{dB}$')
        if desrib==1:
            self.describ(ax,key) 
            
    def plotNXdiff(self,key,ax=1,desrib=1,label='',**kwargs):
        if ax==1:
            fig, ax = plt.subplots()
        ax.plot(self.CH[key].CH1QVASSY['B(T)'][:-1],np.diff(self.CH[key].CH1QVASSY['SY']/self.CH[key].CH1QVASSY['SY'][jfo.dex(self.CH[key].CH1QVASSY,'B(T)',0)]),'o',label=label,**kwargs) 
        ax.set_xlabel('B(T)')
        ax.set_ylabel(r'$\frac{d(X/X(B=0))}{dB}$')
        if desrib==1:
            self.describ(ax,key)
            
    def plotRexp(self,key,ax=1,desrib=1,label=''):
        if ax==1:
            fig, ax = plt.subplots()
        ax.plot(self.CH[key].CH1QVASSY['B(T)'],self.CH[key].CH1VASSY['SY'],'o',label=label)            
        ax.set_xlabel('B(T)')
        ax.set_ylabel(r'R($\Omega$)')
        if desrib==1:
            self.describ(ax,key)  
            
    def paraperpRX(self,listdf,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):            
    
        fig = plt.figure(figsize=(36,32))
        gs0 = gridspec.GridSpec(2, 2, figure=fig,wspace=0.2,hspace=.20)
        
        ax0=fig.add_subplot(gs0[0,0])
        ax1=fig.add_subplot(gs0[0,1])
        ax2=fig.add_subplot(gs0[1,0])
        ax3=fig.add_subplot(gs0[1,1])
        
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].angle==0 and CH[i].F==70.19043 and round(CH[i].T)==2:
                self.plotR(i,ax0,desrib=1,label=CH[i].name) 
                ax0.set_title('Transversal')
                       
        
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].angle==90 and CH[i].F==70.19043 and round(CH[i].T)==2:
                self.plotR(i,ax1,desrib=0,label=CH[i].name) 
                ax1.set_title('logitudinal')    
                
                
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].angle==0 and CH[i].F==70.19043 and round(CH[i].T)==2:
                self.plotX(i,ax2,desrib=1) 
                ax0.set_title('Transversal')
              
                
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].angle==90 and CH[i].F==70.19043 and round(CH[i].T)==2:
                self.plotX(i,ax3,desrib=0,label=CH[i].name) 
                ax1.set_title('logitudinal')   
        ax1.legend(markerscale=2,frameon=False)
            
    def paraperp(self,listdf,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):
    
        fig = plt.figure(figsize=(36,32))
        gs0 = gridspec.GridSpec(2, 2, figure=fig,wspace=0.2,hspace=.20)
        
        ax0=fig.add_subplot(gs0[0,0])
        ax1=fig.add_subplot(gs0[0,1])
        ax2=fig.add_subplot(gs0[1,0])
        ax3=fig.add_subplot(gs0[1,1])
        
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].angle==0 and CH[i].F==70.19043 and round(CH[i].T)==2:
                self.plotRN(i,ax0,desrib=1) 
                ax0.set_title('Transversal')
                ax0.legend()
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].angle==90 and CH[i].F==70.19043 and round(CH[i].T)==2:
                self.plotRN(i,ax1,desrib=0,label=CH[i].name) 
                ax1.set_title('logitudinal')
        ax1.legend()    
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].angle==0 and CH[i].F==70.19043 and round(CH[i].T)==2:
                self.plotRdiff(i,ax2,desrib=0) 
                
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].angle==90 and CH[i].F==70.19043 and round(CH[i].T)==2:
                self.plotRdiff(i,ax3,desrib=0)
                
                
    def paraperpX(self,listdf,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):
    
        fig = plt.figure(figsize=(36,32))
        gs0 = gridspec.GridSpec(2, 2, figure=fig,wspace=0.2,hspace=.20)
        
        ax0=fig.add_subplot(gs0[0,0])
        ax1=fig.add_subplot(gs0[0,1])
        ax2=fig.add_subplot(gs0[1,0])
        ax3=fig.add_subplot(gs0[1,1])
        
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].angle==0 and CH[i].F==70.19043 and round(CH[i].T)==2:
                self.plotXN(i,ax0,desrib=1) 
                ax0.set_title('Transversal')
                ax0.legend()
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].angle==90 and CH[i].F==70.19043 and round(CH[i].T)==2:
                self.plotXN(i,ax1,desrib=0,label=CH[i].name) 
                ax1.set_title('logitudinal')
        ax1.legend()    
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].angle==0 and CH[i].F==70.19043 and round(CH[i].T)==2:
                self.plotXdiff(i,ax2,desrib=0) 
                
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].angle==90 and CH[i].F==70.19043 and round(CH[i].T)==2:
                self.plotXdiff(i,ax3,desrib=0) 
        ax1.legend(markerscale=2,frameon=False)       
                
    def roda(self,listdf,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):
    
        fig = plt.figure(figsize=(36,32))
        gs0 = gridspec.GridSpec(2, 2, figure=fig,wspace=0.2,hspace=.20)
        
        ax0=fig.add_subplot(gs0[0,0])
        ax1=fig.add_subplot(gs0[0,1])
        ax2=fig.add_subplot(gs0[1,0])
        ax3=fig.add_subplot(gs0[1,1])
        
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2:
                self.plotR(i,ax0,desrib=1) 
                
                
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2:
                self.plotX(i,ax1,desrib=0,label='%0.2f deg'%CH[i].angle) 
              
        ax1.legend(markerscale=2)    
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2:
                self.plotRdiff(i,ax2,desrib=0) 
                
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2:
                self.plotXdiff(i,ax3,desrib=0)
    def subtrairR(self,listdf,ax=1):
        Names=list(set([CH[i].name for i in listdf]))
        if ax==1:
            fig, ax = plt.subplots()
        for jj in Names:
            for i in listdf:
                if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==0 and CH[i].name==jj:
                    df0=self.CH[i].CH1VASSY
                elif CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==90 and CH[i].name==jj:
                    df90=self.CH[i].CH1VASSY
            ax.plot(df0['B(T)'],df0['SY']-df90["SY"])
            ax.set_xlabel('B(T)')
            ax.set_ylabel('R$_\perp$ - R$_\parallel$ ')
            
    def subtrairX(self,listdf,ax=1):
        Names=list(set([CH[i].name for i in listdf]))
        if ax==1:
            fig, ax = plt.subplots()
        for jj in Names:
            for i in listdf:
                if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==0 and CH[i].name==jj:
                    df0=self.CH[i].CH1QVASSY
                elif CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==90 and CH[i].name==jj:
                    df90=self.CH[i].CH1QVASSY
            ax.plot(df0['B(T)'],df0['SY']-df90["SY"])
            ax.set_xlabel('B(T)')
            ax.set_ylabel('X$_\perp$ - X$_\parallel$ ')
         
            
    def subtrairdR(self,listdf,ax=1):
        Names=list(set([CH[i].name for i in listdf]))
        if ax==1:
            fig, ax = plt.subplots()        
        for jj in Names:
            for i in listdf:
                if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==0 and CH[i].name==jj:
                    df0=self.CH[i].CH1VASSY
                elif CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==90 and CH[i].name==jj:
                    df90=self.CH[i].CH1VASSY
                    
            ax.plot(df0['B(T)'][:-1],np.diff(df0['SY']-df90['SY']),label=jj)
            ax.set_xlabel('B(T)')
            ax.set_ylabel('dR$_\perp$ - dR$_\parallel$ ')
    
    def subtrairdR(self,listdf,ax=1):
        Names=list(set([CH[i].name for i in listdf]))
        if ax==1:
            fig, ax = plt.subplots()        
        for jj in Names:
            for i in listdf:
                if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==0 and CH[i].name==jj:
                    df0=self.CH[i].CH1VASSY
                elif CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==90 and CH[i].name==jj:
                    df90=self.CH[i].CH1VASSY
                    
            ax.plot(df0['B(T)'][:-1],np.diff(df0['SY']-df90['SY']),label=jj)
            ax.set_xlabel('B(T)')
            ax.set_ylabel('dR$_\perp$ - dR$_\parallel$ ')

    def plotRlinearlog(self,key,ax=1,desrib=1,label=''):
        if ax==1:
            fig, ax = plt.subplots()
        ax.plot(self.CH[key].CH1QVASSY['B(T)'],self.CH[key].CH1VASSY['SY'],'o',label=label)            
        ax.set_xlabel('B(T)')
        ax.set_ylabel(r'R($\Omega$)')
        if desrib==1:
            self.describ(ax,key)
          

    def subtrainrdX(self,listdf,ax=1):
        Names=list(set([CH[i].name for i in listdf]))
        if ax==1:
            fig, ax = plt.subplots()
        for jj in Names:
            for i in listdf:
                if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==0 and CH[i].name==jj:
                    df0=self.CH[i].CH1QVASSY
                    
                elif CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==90 and CH[i].name==jj:
                    df90=self.CH[i].CH1QVASSY
                    
            ax.plot(df0['B(T)'][:-1],np.diff(df0['SY']/df0["SY"][jfo.dex(df0,"SY",0)]-df90['SY']/df90["SY"][jfo.dex(df90,"SY",0)]),label=jj)
            ax.set_xlabel('B(T)')
            ax.set_ylabel('X$_\perp$ - X$_\parallel$ ')
            
    def subtrainrdR(self,listdf,ax=1):
        Names=list(set([CH[i].name for i in listdf]))
        if ax==1:
            fig, ax = plt.subplots()        
        for jj in Names:
            for i in listdf:
                if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==0 and CH[i].name==jj:
                    df0=self.CH[i].CH1VASSY
                elif CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==90 and CH[i].name==jj:
                    df90=self.CH[i].CH1VASSY
                    
            ax.plot(df0['B(T)'][:-1],np.diff(df0['SY']/df0["SY"][jfo.dex(df0,"SY",0)]-df90['SY']/df90["SY"][jfo.dex(df90,"SY",0)]),label=jj)
            ax.set_xlabel('B(T)')
            ax.set_ylabel('R$_\perp$ - R$_\parallel$ ')
            
            
    def subtrairRchart(self,listdf,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):
    
        fig = plt.figure(figsize=(36,12))
        gs0 = gridspec.GridSpec(1, 3, figure=fig,wspace=0.25,hspace=.20)
        
        ax0=fig.add_subplot(gs0[0])
        ax1=fig.add_subplot(gs0[1])
        ax2=fig.add_subplot(gs0[2])
        
        
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==0:
                self.plotRN(i,ax0,desrib=1) 
        ax0.set_title('Transversal')        
                
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==90:
                self.plotRN(i,ax1,desrib=0,label=CH[i].name) 
        ax1.set_title('logitudinal')      
        ax1.legend(markerscale=2,frameon=False)    

        self.subtrairR(listdf,ax2) 
        
    def subtrairdRchart(self,listdf,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):
    
        fig = plt.figure(figsize=(36,12))
        gs0 = gridspec.GridSpec(1, 3, figure=fig,wspace=0.25,hspace=.20)
        
        ax0=fig.add_subplot(gs0[0])
        ax1=fig.add_subplot(gs0[1])
        ax2=fig.add_subplot(gs0[2])
        
        
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==0:
                self.plotRdiff(i,ax0,desrib=1) 
        ax0.set_title('Transversal')                    
            
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==90:
                self.plotRdiff(i,ax1,desrib=0,label=CH[i].name) 
              
        ax1.set_title('logitudinal')      
        ax1.legend(markerscale=2,frameon=False)          

        self.subtrairdR(listdf,ax2) 

        
    def subtrairXchart(self,listdf,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):
    
        fig = plt.figure(figsize=(36,12))
        gs0 = gridspec.GridSpec(1, 3, figure=fig,wspace=0.25,hspace=.20)
        
        ax0=fig.add_subplot(gs0[0])
        ax1=fig.add_subplot(gs0[1])
        ax2=fig.add_subplot(gs0[2])
        
        
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==0:
                self.plotXN(i,ax0,desrib=1) 
        ax0.set_title('Transversal')        
                
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==90:
                self.plotXN(i,ax1,desrib=0,label=CH[i].name) 
        ax1.set_title('logitudinal')      
        ax1.legend(markerscale=2,frameon=False)    

        self.subtrairX(listdf,ax2) 

        
    def subtrairdXchart(self,listdf,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):
    
        fig = plt.figure(figsize=(36,12))
        gs0 = gridspec.GridSpec(1, 3, figure=fig,wspace=0.25,hspace=.20)
        
        ax0=fig.add_subplot(gs0[0])
        ax1=fig.add_subplot(gs0[1])
        ax2=fig.add_subplot(gs0[2])
        
        
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==0:
                self.plotXdiff(i,ax0,desrib=1) 
        ax0.set_title('Transversal')                    
            
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==90:
                self.plotXdiff(i,ax1,desrib=0,label=CH[i].name) 
              
        ax1.set_title('logitudinal')      
        ax1.legend(markerscale=2,frameon=False)          

        self.subtrairdX(listdf,ax2) 
        
        
    def subtrairdNRchart(self,listdf,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):
    
        fig = plt.figure(figsize=(36,12))
        gs0 = gridspec.GridSpec(1, 3, figure=fig,wspace=0.25,hspace=.20)
        
        ax0=fig.add_subplot(gs0[0])
        ax1=fig.add_subplot(gs0[1])
        ax2=fig.add_subplot(gs0[2])
        
        
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==0:
                self.plotNRdiff(i,ax0,desrib=1) 
        ax0.set_title('Transversal')                    
            
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==90:
                self.plotNRdiff(i,ax1,desrib=0,label=CH[i].name) 
              
        ax1.set_title('logitudinal')      
        ax1.legend(markerscale=2,frameon=False)          

        self.subtrainrdR(listdf,ax2) 
        
    def subtrairdNXchart(self,listdf,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):
    
        fig = plt.figure(figsize=(36,12))
        gs0 = gridspec.GridSpec(1, 3, figure=fig,wspace=0.25,hspace=.20)
        
        ax0=fig.add_subplot(gs0[0])
        ax1=fig.add_subplot(gs0[1])
        ax2=fig.add_subplot(gs0[2])
        
        
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==0:
                self.plotNXdiff(i,ax0,desrib=1) 
        ax0.set_title('Transversal')                    
            
        for i in listdf:
            if CH[i].I==0.0001 and CH[i].F==70.19043 and round(CH[i].T)==2 and CH[i].angle==90:
                self.plotNXdiff(i,ax1,desrib=0,label=CH[i].name) 
              
        ax1.set_title('logitudinal')      
        ax1.legend(markerscale=2,frameon=False)          

        self.subtrainrdX(listdf,ax2) 
        
        
    def paraperp23(self,listdf0,listdf90,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):
    
        fig = plt.figure(figsize=(36,32))
        gs0 = gridspec.GridSpec(2, 2, figure=fig,wspace=0.2,hspace=.20)
        
        ax0=fig.add_subplot(gs0[0,0])
        ax1=fig.add_subplot(gs0[0,1])
        ax2=fig.add_subplot(gs0[1,0])
        ax3=fig.add_subplot(gs0[1,1])
        
        for i in listdf0:
            
            self.plotRN(i,ax0,desrib=0) 
       
            
        for i in listdf90:
            
            self.plotRN(i,ax1,desrib=0,label="%0.2f Hz"%CH[i].F) 
        ax1.set_title('logitudinal')
        ax1.legend()    
        for i in listdf0:
            
            self.plotRdiff(i,ax2,desrib=0) 
                
        for i in listdf90:
            
            self.plotRdiff(i,ax3,desrib=0)
            
        ax0.set_title('Transversal T=%i K'%CH[i].T)    
        ax1.set_title('logitudinal T=%i K'%CH[i].T)    
        
        
    def paraperp23X(self,listdf0,listdf90,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):
    
        fig = plt.figure(figsize=(36,32))
        gs0 = gridspec.GridSpec(2, 2, figure=fig,wspace=0.2,hspace=.20)
        
        ax0=fig.add_subplot(gs0[0,0])
        ax1=fig.add_subplot(gs0[0,1])
        ax2=fig.add_subplot(gs0[1,0])
        ax3=fig.add_subplot(gs0[1,1])
        
        for i in listdf0:
            
            self.plotXN(i,ax0,desrib=0) 
       
            
        for i in listdf90:
            
            self.plotXN(i,ax1,desrib=0,label="%0.2f Hz"%CH[i].F) 
        ax1.set_title('logitudinal')
        ax1.legend()    
        for i in listdf0:
            
            self.plotXdiff(i,ax2,desrib=0) 
                
        for i in listdf90:
            
            self.plotXdiff(i,ax3,desrib=0)
            
        ax0.set_title('Transversal T=%i K'%CH[i].T)    
        ax1.set_title('logitudinal T=%i K'%CH[i].T)  
        
        
        
    def paraperp23SYAS(self,listdf0,path=r'C:\Users\jaime\OneDrive\Documentos\Pos-doc\Data-plots\quadrature voltage\2ndHarmic'):
    
        fig = plt.figure(figsize=(36,32))
        gs0 = gridspec.GridSpec(2, 2, figure=fig,wspace=0.2,hspace=.20)
        
        ax0=fig.add_subplot(gs0[0,0])
        ax1=fig.add_subplot(gs0[0,1])
        ax2=fig.add_subplot(gs0[1,0])
        ax3=fig.add_subplot(gs0[1,1])
        
        for i in listdf0:
            
            
            self.plotRN(i,ax0,desrib=0) 
            self.plotRNAS(i,ax2,desrib=0) 
            
        for i in listdf0:
            
            self.plotXN(i,ax1,desrib=0,label="%0.2f Hz"%CH[i].F) 
            self.plotXNAS(i,ax3,desrib=0,label="%0.2f Hz"%CH[i].F) 
        ax1.legend()



         
#%%




def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The Ain idea behind this
    approach is to Ake for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.norAl(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import Atplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from Ath import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
        
    except (ValueError, msg):
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too sAll for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.At([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

#%%
def forceAspect(ax,aspect=1):
    im = ax.get_iAges()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

#%%
    
def saveQDPD(data,lista):
    '''
    Função para criar novos df from QD files 
    '''
    listanova=lista
    
    
    listanova=['T (K)' if x=='Temperature (K)' else x for x in listanova]
    listanova=['H (T)' if x=='Agnetic Field (Oe)' else x for x in listanova]
    listanova=['R (ohms)' if x=='Bridge 1 Resistance (Ohms)' else x for x in listanova]
    
    df=pd.DataFrame(columns=lista)
    
    for i in lista:
        df[i]=data[i]
    
    if lista.count('Agnetic Field (Oe)')==1:
        df['Agnetic Field (Oe)']=data['Agnetic Field (Oe)']/10000
    
    df.columns=listanova
    
    return df


#%%

def desin(A,mm):
    '''
    A dado em ohns/tesla e mm dados em milimetros
    resposta em cm^-3
    '''
    from scipy.constants import e
    
    return (1/(A*e*mm*10**-3))*10**-6 

#%%
    
def saveADR(data,T,H,X,Y,Vos,geo):
    df=pd.DataFrame(columns=['T','H','R','rho','sigA','x','y'])
    df['T']=data[T]
    df['R']=data[X]/(Vos*10E-3)
    df['rho']=(geo)*data[X]
    df['sigA']=1/df['rho']
    df['H']=(0.00291 + 0.03537*data[H])
    df['x']=data[X]
    df['y']=data[Y]
    return df
#%%
   

#%%


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
            função faz plots norAlizados dos dados em relação a um ponto em x
            '''
            yni=(np.abs(self.x-xn)).idxmin()
            plt.plot(self.x,self.y/self.y[yni],'*',**dic)
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
        '''
        Seleciona as colunas legais 
        '''
        dfn=pd.DataFrame(self.df[lista[0]])
        for i in np.arange(1,len(lista),1):
            dfn[lista[i]]=self.df[lista[i]]
        return dfn
    def inter(data,dadox,dadoy,xmin,xmax,npoints=1000,smoothing=1,K=3):
        ''' 
        função que faz a interpolção em de colunas de um dataframe
        '''
        import numpy as np
        from scipy import interpolate 
        data.sort_values(dadox,ascending=True,inplace=True)
        data.dropna(axis=1 ,how='all', inplace=True)
        data.drop_duplicates(dadox,inplace=True)
        df=data[(data[dadox]<xmax) & (data[dadox]>xmin)] 
    
    
        xs=np.linspace(xmin,xmax,npoints)
    
        spl = interpolate.UnivariateSpline(df[dadox], df[dadoy], k=K , s=smoothing,ext=3)
        return spl , xs
    def fit(f,x,y,limitx,plot=0):
        import numpy as np
        popt,pcov=curve_fit(f,x[(x>limitx[0])&(x<limitx[1])],y[(x>limitx[0])&(x<limitx[1])],bounds=(0, np.inf),maxfev = 10000)
        if plot==1:
            plt.plot(x[(x>limitx[0])&(x<limitx[1])],y[(x>limitx[0])&(x<limitx[1])],label='data')
            plt.plot(x[(x>limitx[0])&(x<limitx[1])],f(x[(x>limitx[0])&(x<limitx[1])],*popt),'+',label='fit')
            plt.legend()
        return popt,pcov
    def plot(*args):
        '''
        func para fazer os plots de dados e coloca os limites 10%aciA do colocado como utlimos parametros
        ex:
        dados limites em tuplas plot(*args,(xl,xu))
        '''
        plt.plot(*args[:-1])
        plt.xlim(*args[-1])
        col=args[0]
        xmin=args[-1][0]
        xmax=args[-1][1]
        x1=(np.abs(args[0]-xmin)).idxmin()
        x2=(np.abs(args[0]-xmax)).idxmin()
        if args[1][x1]<args[1][x2]:
            y1=(args[1][x1],args[1][x2])
        else:
            y1=(args[1][x2],args[1][x1])
        plt.ylim(y1)     
    def plotVRH(*args):
        '''
        func para fazer os plots de dados e coloca os limites 10%aciA do colocado como utlimos parametros
        ex:
        dados limites em tuplas plot(*args,(xl,xu))
        '''
        x=args[0]**(-1/4)
        y=np.log(args[1])
        plt.plot(x,y,*args[2:-1])
        plt.xlim(args[-1][1]**(-1/4),args[-1][0]**(-1/4))
        col=args[0]
        xmin=args[-1][0]
        xmax=args[-1][1]
        x1=(np.abs(args[0]-xmin)).idxmin()
        x2=(np.abs(args[0]-xmax)).idxmin()
        if (args[1][x1]**(-1/4))>(args[1][x2]**(-1/4)):
            y1=(args[1][x1],args[1][x2])
        else:
            y1=(args[1][x2],args[1][x1])
        plt.ylim(np.log(y1))
    
        
