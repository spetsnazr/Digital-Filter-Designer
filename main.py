# -*- coding: utf-8 -*-
import random
import math as m
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy import signal
import copy
import os
#print(os.getcwd())


#Vrsta filtera: 0-Lowpass
#               1-Highpass
#               2-Bandpass
#               3-Bandstop
# fw for 0 and 1 is the cutoff frequency
# for 2 and 3 it's a 2 element list with the border frequencies
# u Hz

def Filter(vrstaFilt,fHz,duz_filt,N,eps,T0,korak = 0.01,Fs = 2000, tipFilt = 0,K = 3,faktor_vel_okoline = 1,numer = 0,denumer = 0,):
   #lista = [[m.inf,m.inf]]*L;
   w = 0;
   if(tipFilt == 1): duz_filt = int(duz_filt/2);
   #num = [0.5,0.1,0.166,0.22];
   #num = [1]*duz_filt;
   num = np.random.uniform(-1,1,duz_filt);
   #den = [1]*duz_filt;
   den = 1;
   if(tipFilt != 0):
      # den = [0]*duz_filt;
      # den[0] = 1;
      den = np.random.uniform(-1,1,duz_filt);
   #tst = signal.dlti(num,den,dt=1./Fs);
   BestFilt = copy.deepcopy(num);
   bestDen = copy.deepcopy(den);
   #print(fHz[1]);
   err = 0;
   w = 0;
   if(numer != 0): 
       num = numer;       
       duz_filt = len(numer);
   if(denumer != 0):
       den = denumer;
   a,hn = signal.freqz(num,den,worN=Fs)
   #print(a);
   #Definisemo idealne filtere na osnovu kojih ćemo kreirati nas filter
   if(vrstaFilt == 0):
       w = max(fHz,0)*2*m.pi/Fs;
       def filt(w,wd):
           if(wd < w): return 1
           else: return 0;
   elif(vrstaFilt==1):
       w = max(fHz,0)*2*m.pi/Fs;
       def filt(w,wd):
           if(wd > w): return 1
           else: return 0;
   elif(vrstaFilt==2):
       w = [(fHz[0])*2*m.pi/Fs,(fHz[1])*2*m.pi/Fs];
       def filt(w,wd):
           if(wd > w[0] and wd < w[1]): return 1
           else: return 0;
   elif(vrstaFilt==3):
       w = [(fHz[0])*2*m.pi/Fs,(fHz[1])*2*m.pi/Fs];
       def filt(w,wd):
           if(wd < w[0] or wd > w[1]): return 1
           else: return 0;
   
    #grešku računamo kao sumu kvadrata pojednih grešaka za Fs različitih diskretnih frekvencija
    #od 0 do pi
   if(numer == 0):
        for i in range(Fs):
            err += (filt(w,a[i])-abs(hn[i]))*(filt(w,a[i])-abs(hn[i]));
   else:
        for i in range(Fs):
            err += abs((filt(w,a[i])-abs(hn[i])));
   print(err);       
   plt.figure();
   plt.plot(a,20*np.log10(abs(hn)));
   plt.show();   
   okol = [-korak,0,korak];
   i = 0;

   #print(err);
   while(i < N):
       #print(i);
       temperr = 0;
       a,hn = signal.freqz(num,den,worN=Fs)
       if(numer == 0):
            for p in range(Fs):
                temperr += (filt(w,a[p])-abs(hn[p]))*(filt(w,a[p])-abs(hn[p]));
       else:
            for p in range(Fs):
                temperr += abs((filt(w,a[p])-abs(hn[p])));
       l = 0;
       errrj = [];
       potrj = [];
       denpotrj = [];
       while(l < K):
           temp2err = 0;
           #generisanje neke tacke iz okoline
           for j in range(duz_filt*faktor_vel_okoline):
               temp2err = 0;
               temp = copy.deepcopy(num);
               dentemp = copy.deepcopy(den);
               for p in range(len(num)):
                    temp[p] += random.choice(okol);
                    if(tipFilt == 1):
                        dentemp[p] += random.choice(okol);
               #print(temp);
               x,temphn = signal.freqz(temp,dentemp,worN=Fs);
               if(numer == 0):
                    for p in range(Fs):
                        temp2err += (filt(w,a[p])-abs(temphn[p]))*(filt(w,a[p])-abs(temphn[p]));
               else:
                    for p in range(Fs):
                        temp2err += abs((filt(w,a[p])-abs(temphn[p])));
               if(temp2err < temperr):
                   potrj.append(temp);
                   denpotrj.append(dentemp);
                   errrj.append(temp2err);
               elif(m.exp(-abs(temperr-temp2err)/T0) > random.random()):
                   potrj.append(temp);
                   denpotrj.append(dentemp);
                   errrj.append(temp2err);
           #print(len(potrj));
           if(len(potrj) > 0):
                rannum = random.randint(0,len(potrj)-1);
                num = copy.deepcopy(potrj[rannum]);
                den = copy.deepcopy(denpotrj[rannum]);
                temperr = errrj[rannum];
           if(temperr < err):
               BestFilt  = copy.deepcopy(num);
               bestDen = copy.deepcopy(den);
               err = temperr;
           #print(i);  
           l+=1;
       T0 = 0.95*T0;
       if(T0 < eps): break;
       i+=1;
   a,hn = signal.freqz(BestFilt,bestDen,worN=Fs);
   print(err);
   plt.figure();
   plt.plot(a*Fs/(2*m.pi),20*np.log10(abs(hn)));
   plt.show();
   if(tipFilt == 1): rtn = bestDen;
   else: rtn = [1];
   return [BestFilt,rtn] ;
x = 0;
o = [300,500];

x = Filter(2,o,40,3000,1e-10,1000,0.001,2000,0,3,3,[-0.004406814030000001, -0.041759714499999996, -0.0317684973, 0.06974305010000001, 0.105837499, -0.02341463290000001, -0.154896872, -0.062230881800000006, 0.13779516399999997, 0.16919452799999996, -0.039947499444300004, -0.18959514100000002, -0.08336014970000001, 0.12301056509999998, 0.115137329, -0.02278915091, -0.0680948243, -0.0231891627, 0.01578269994, -0.006066464300000006, 0.016803293810000015, 0.044215653699999996, 0.025849964700000005, -0.05869460090000001, -0.0657177467, 0.02286155633999999, 0.052241209, 0.017698589599999995, -0.0196893256, -0.011278581799999987, -0.011547634300000014, -0.008674572400000012, 0.0054202144, 0.0286757709, 0.02706677669999999, -0.013781096910000007, -0.03300359749999999, 0.0014176363000000004, 0.011390863699999991, 0.0014949517000000002])
print(x);
