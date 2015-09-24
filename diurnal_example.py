#!/usr/local/sci/bin/python

import pylab
import numpy
import datetime
import matplotlib.pyplot as plt
import math

# create schematic diagram for diurnal cycle test

plot_delta=0.07
plot_size=0.4

plt.clf()
plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')

# top left - schematic of temperatures and requirements

topleft=plt.axes([0.0+plot_delta,0.5+plot_delta,plot_size,plot_size])
plt.text(-5,10,'(a)')

hours=[0.,6.,12.,18.,24.]

# winter time temperatures - range >5

temperatures=[3,1,7,5,4]

plt.plot(hours,temperatures,'bo',ms=10)
plt.setp(topleft, ylim=[-1,10], xlabel=('hours'),ylabel=('Temperature (C)'),yticks=[0,5,10],xticks=[0,6,12,18,24],xlim=(-1,25),)
plt.plot([-1,25],[0,0],'k--')

plt.annotate("",xy=(9,temperatures[1]),xycoords='data',xytext=(9,temperatures[2]),textcoords='data',arrowprops=dict(arrowstyle="<-"))
plt.text(10,4,'>5C')

# top right - schematic of fitting sine curve and cost function

topright=plt.axes([0.5+plot_delta,0.5+plot_delta,plot_size,plot_size])
plt.text(-5,10,'(b)')
plt.plot(hours,temperatures,'bo',ms=10)
plt.setp(topright, ylim=[-1,10], xlabel=('hours'),ylabel=('Temperature (C)'),yticks=[0,5,10],xticks=[0,6,12,18,24],xlim=(-1,25),)
plt.plot([-1,25],[0,0],'k--')

fulltimes=numpy.arange(27)-1

for shift in [0,1,2,3]:
    sincurve=(numpy.sin((fulltimes-shift)*2.*math.pi/24.))*((max(temperatures)-min(temperatures))/2.)+((max(temperatures)+min(temperatures))/2.)
    exec("plt.plot(fulltimes,sincurve,'"+str(shift*0.27)+"')")

    exec("plt.text("+str(6+shift)+",max(sincurve)*1.1,'"+str(shift)+"',color='"+str(shift*0.27)+"')")
    

    for hr in range(len(hours)):
        exec("plt.plot([hours[hr]-0.3+(shift*0.2),hours[hr]-0.3+(shift*0.2)],[temperatures[hr],sincurve1[(6*hr)+1]],'"+str(shift*0.27)+"', ls=':')")


sincurve=(numpy.sin((fulltimes-9)*2.*math.pi/24.))*((max(temperatures)-min(temperatures))/2.)+((max(temperatures)+min(temperatures))/2.)
plt.plot(fulltimes,sincurve,'r-')

for hr in range(len(hours)):
    plt.plot([hours[hr]-0.1,hours[hr]-0.1],[temperatures[hr],sincurve[(6*hr)+1]],'r-')

plt.text(6+9,max(sincurve)*1.1,str(9),color='r')


# schematic of cost function for uncertainty


bottomleft=plt.axes([0.20+plot_delta,0.0+plot_delta,plot_size,plot_size])
plt.text(-4,25,'(c)')

costs=[18,17,16,14.5,12.5,10,9,8,6.5,6,6.5,7,7.5,8,9,10,11,12,13,14,15.5,16,16.5,17,18]

plt.plot(fulltimes[1:-1],costs,'go',ms=10)
plt.setp(bottomleft, ylim=[0,25], xlabel=('hours'),ylabel=('Cost Function'),yticks=[],xticks=[0,6,12,18,24],xlim=(-1,25),)

plt.text(9,1.8,str(9),color='r')

plt.plot([9,9],[0,25],'r-')
plt.plot([-1,25],[6,6],'k:')
plt.plot([-1,25],[10,10],'k:')
plt.plot([-1,25],[14,14],'k:')
plt.plot([-1,25],[18,18],'k:')

plt.plot([5,5],[0,25],'b-')
plt.plot([15,15],[0,25],'b-')
plt.text(5,1.8,str(5),color='b')
plt.text(15,1.8,str(15),color='b')

plt.annotate("",xy=(9,20),xycoords='data',xytext=(15,20),textcoords='data',arrowprops=dict(arrowstyle="<-"))
plt.annotate("",xy=(9,20),xycoords='data',xytext=(5,20),textcoords='data',arrowprops=dict(arrowstyle="<-"))

plt.plot([3,3],[0,25],'b--')
plt.text(3,1.8,str(3),color='b')
plt.annotate("",xy=(9,20),xycoords='data',xytext=(3,20),textcoords='data',arrowprops=dict(arrowstyle="<-"))


plt.savefig('diurnal_example.eps')
plt.savefig('diurnal_example.png')
