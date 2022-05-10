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
os.chdir (work)


def pdinter(data1,dadox,dadoy,xmin,xmax,npoints=1000,smoothing=1,K=3):
    ''' 
    função que faz a interpolção em uma faixa 
    '''
    import numpy as np
    from scipy import interpolate
    import matplotlib.pyplot as plt
    data=data1
    data.sort_values(dadox,ascending=True,inplace=True)
    data.dropna(axis=1 ,how='all', inplace=True)
    data.drop_duplicates(dadox,inplace=True)
    
    
    
    xs=np.linspace(xmin,xmax,npoints)
    
    spl = interpolate.UnivariateSpline(data[dadox], data[dadoy], k=K , s=smoothing,ext=3)
    
    #plt.plot(data[dadox],data[dadoy])
    #plt.plot(xs,spl(xs))
    return spl , xs


def foldershow(path,ext):
    '''
    function to show all file to import files from an folder
    '''
    os.chdir (path)
    a= [name for name in os.listdir(".") if name.endswith("."+ ext)]
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
    df=pd.DataFrame(columns=['T','R','rho','sigma'])
    df['T']=data[T]
    df['R']=data[R]
    df['rho']=(0.05/1.221)*data[R]
    df['sigma']=1/df['rho']
    return df

def saveMR(data,T,R,F,rate):
    df=pd.DataFrame(columns=['T','R','rho','sigma','Field'])
    df['T']=data[T]
    df['R']=data[R]
    df['rho']=(0.05/1.221)*data[R]
    df['sigma']=1/df['rho']
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
        df1=pd.DataFrame(data={'T(K)':df['Temperature (K)'],'B(T)':df['Field (Oe)']/10000,'Angle (deg)':df['Sample Position (deg)'],'V(V)':df['In Phase Voltage Ampl Ch1 (V)'],'4V(V)':df['Quadrature Voltage Ch1 (V)'],'R(ohms)':df['Resistance Ch1 (Ohms)'],"Phase(deg)":df["Phase Angle Ch1 (deg)"]})
        df1.dropna(inplace=True)
        df1.reset_index(drop=True,inplace=True)
    
    except:
        df1=pd.DataFrame()
    try:
        df2=pd.DataFrame(data={'T(K)':df['Temperature (K)'],'B(T)':df['Field (Oe)']/10000,'Angle (deg)':df['Sample Position (deg)'],'V(V)':df['In Phase Voltage Ampl Ch2 (V)'],'4V(V)':df['Quadrature Voltage Ch2 (V)'],'R(ohms)':df['Resistance Ch2 (Ohms)'],"Phase(deg)":df["Phase Angle Ch2 (deg)"]})
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


def ASSYextractinter(dfinit,X='B(T)',Y='R(ohms)',Hmax=8):
    import jfodata
    
    inter=jfodata.pdinter(dfinit,X,Y,-Hmax,Hmax,smoothing=0.01,npoints=100,K=3)
    df=pd.DataFrame(columns=['T(K)','Angle (deg)',X,Y])
    for jj in range(len(inter[1])):
        df.loc[jj]=[dfinit['T(K)'][0],dfinit['Angle (deg)'][0],inter[1][jj],inter[0](inter[1][jj])]
    
    
    
    dfp=df[df[X]>0.1]
    dfp=dfp.sort_values(X)
    dfp.reset_index(drop=True,inplace=True)
    
    dfn=df[df[X]<-0.1]
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
    
    dff=pd.concat([dfn,dfp]).sort_values(X)
    return dff.reset_index(drop=True)







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
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
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
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
        
    except (ValueError, msg):
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

#%%
def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

#%%
    
def saveQDPD(data,lista):
    '''
    Função para criar novos df from QD files 
    '''
    listanova=lista
    
    
    listanova=['T (K)' if x=='Temperature (K)' else x for x in listanova]
    listanova=['H (T)' if x=='Magnetic Field (Oe)' else x for x in listanova]
    listanova=['R (ohms)' if x=='Bridge 1 Resistance (Ohms)' else x for x in listanova]
    
    df=pd.DataFrame(columns=lista)
    
    for i in lista:
        df[i]=data[i]
    
    if lista.count('Magnetic Field (Oe)')==1:
        df['Magnetic Field (Oe)']=data['Magnetic Field (Oe)']/10000
    
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
    df=pd.DataFrame(columns=['T','H','R','rho','sigma','x','y'])
    df['T']=data[T]
    df['R']=data[X]/(Vos*10E-3)
    df['rho']=(geo)*data[X]
    df['sigma']=1/df['rho']
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
            função faz plots normalizados dos dados em relação a um ponto em x
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
        func para fazer os plots de dados e coloca os limites 10%acima do colocado como utlimos parametros
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
        func para fazer os plots de dados e coloca os limites 10%acima do colocado como utlimos parametros
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
    
        
