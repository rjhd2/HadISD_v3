#!/usr/local/sci/bin/python

import pylab
import numpy
import datetime
import matplotlib.pyplot as plt
import math


plt.clf()
plt.axes([0.1,0.1,0.87,0.85])
times=numpy.arange(27)

slps=(numpy.sin((numpy.arange(27)/24.)*2.*math.pi)*2.)+1010.

slps[14:17]=slps[14:17]+7

plt.plot(times[11:22],slps[11:22],'0.8')
plt.plot(times[11:22],slps[11:22],'bo',ms=10)

plt.annotate("",xy=(13.5,slps[13]),xycoords='data',xytext=(13.5,slps[14]),textcoords='data',arrowprops=dict(arrowstyle="<|-"))
plt.annotate("",xy=(13.5,slps[13]),xycoords='data',xytext=(13.5,slps[14]),textcoords='data',arrowprops=dict(arrowstyle="|-|"))
plt.text(times[13],(slps[13]+slps[14])/2.,'>t')
plt.annotate("",xy=(16.5,slps[16]),xycoords='data',xytext=(16.5,slps[17]),textcoords='data',arrowprops=dict(arrowstyle="<|-"))
plt.annotate("",xy=(16.5,slps[16]),xycoords='data',xytext=(16.5,slps[17]),textcoords='data',arrowprops=dict(arrowstyle="|-|"))
plt.text(times[17]+0.2,(slps[16]+slps[17])/2.,'>t',horizontalalignment='right')


plt.annotate("",xy=(12.5,slps[12]),xycoords='data',xytext=(12.5,slps[13]),textcoords='data',arrowprops=dict(arrowstyle="<|-"))
plt.annotate("",xy=(12.5,slps[12]),xycoords='data',xytext=(12.5,slps[13]),textcoords='data',arrowprops=dict(arrowstyle="|-|"))
#plt.annotate("",xy=(17.5,slps[17]),xycoords='data',xytext=(17.5,slps[18]),textcoords='data',arrowprops=dict(arrowstyle="<-")) # too small to show head of arrow
plt.annotate("",xy=(17.5,slps[17]),xycoords='data',xytext=(17.5,slps[18]),textcoords='data',arrowprops=dict(arrowstyle="|-|"))


plt.annotate("",xy=(14.5,slps[14]),xycoords='data',xytext=(14.5,slps[15]),textcoords='data',arrowprops=dict(arrowstyle="<|-"))
plt.annotate("",xy=(14.5,slps[14]),xycoords='data',xytext=(14.5,slps[15]),textcoords='data',arrowprops=dict(arrowstyle="|-|"))
plt.annotate("",xy=(15.5,slps[15]),xycoords='data',xytext=(15.5,slps[16]),textcoords='data',arrowprops=dict(arrowstyle="<|-"))
plt.annotate("",xy=(15.5,slps[15]),xycoords='data',xytext=(15.5,slps[16]),textcoords='data',arrowprops=dict(arrowstyle="|-|"))

plt.annotate("",xy=(15,1007),xycoords='data',xytext=(12.5,slps[13]),textcoords='data',arrowprops=dict(arrowstyle="<-",connectionstyle="arc3,rad=0.3",color='r'))
plt.annotate("",xy=(15,1007),xycoords='data',xytext=(14.5,slps[15]),textcoords='data',arrowprops=dict(arrowstyle="<-",connectionstyle="arc3,rad=0.1",color='r'))
plt.annotate("",xy=(15,1007),xycoords='data',xytext=(15.5,slps[16]),textcoords='data',arrowprops=dict(arrowstyle="<-",connectionstyle="arc3,rad=-0.1",color='r'))
plt.annotate("",xy=(15,1007),xycoords='data',xytext=(17.5,slps[17]),textcoords='data',arrowprops=dict(arrowstyle="<-",connectionstyle="arc3,rad=-0.3",color='r'))
props = dict(facecolor='white',boxstyle='round')
plt.text(15,1007,'<(t/2)',horizontalalignment='center', bbox=props)




plt.xlabel('hours')
plt.ylabel('slp (hPa)')
plt.ylim(1006,1018)
plt.xlim(11.5,20.5)
plt.title('Spike Check Example')

a = plt.axes([.67, .75, .3, .2])
plt.plot(times,slps,'k-')
plt.plot(times,slps,'bo')
plt.setp(a, ylim=[1000.,1020.], xlabel=('hours'),ylabel=('slp (hPa)'),yticks=[1005,1015],xticks=[0,10,20],xlim=(0,25))

plt.savefig('images/spike_example.eps')
plt.savefig('images/spike_example.png')




