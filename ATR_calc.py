from __future__ import division
import csv, numpy, glob
import matplotlib.pyplot as plt


# This Funcion takes in the time from the CSV which is in the format
# HH:mm:ss and coverts it to seconds
def timeConv(tCSV):
    (h,m,s)     = tCSV.split(':')
    timD = float(h)*3600+float(m)*60+float(s)
    return timD


# function to find the slope of a data
def findslope(t,T):
    mc,res,_,_,_ = numpy.polyfit(t,T,1,full = True)
    return mc, res

# calculates the rsquared of a set of values given y and mean
def r2(y,yhat):
    y2 = []
    yhat2 = []
    yavg = numpy.mean(y)
    y2[:] = [(y[x]-yavg)**2 for x in range(len(y))]
    yhat2[:] = [(yhat[x]-yavg)**2 for x in range(len(yhat))]
    sst = sum(y2)
    ssr = sum(yhat2)
    rsq = ssr/sst
    return rsq

#Glob finds all csvs and opens them
for file in glob.glob("*.csv"):

#set these to zero as they will be appended later    
    timDec = []
    Traw =[]

#grabs data from the current csv, time converted to sec    
    with open(file,'rt') as csvfile:
        csvRead = csv.reader(csvfile,delimiter=',')
        csvfile.readline()
        for row in csvRead:
            tD = timeConv(row[0])
            Traw.append(float(row[1]))
            timDec.append(tD)
            
#finds max temp, as well as the integer position of Tmax
    Tmax = max(Traw)
    posTmax = [i for i,x in enumerate(Traw) if x == Tmax]
    posTmax = posTmax[0]
    
#these are inputs, grams of water, grams of polyol, polyol ew, cp, r,mt
#####future, get these from another csv
    Tamb = numpy.mean(Traw[0:20])
    Twork = 0
    Tcorr = []
    lnTfit = []
    pNCO = []
    r = 1
    gh2o = 1.718
    gpol = 42.96
    ewpol = 1169
    H  = 93900*(gpol/ewpol) + 125500*(gh2o/18) #j
    cp = 1.4 #j/g.K
    mt = 75
    
#fits the data to the area that dx/dt = 0 (after Tmax)
    fitt = timDec[posTmax+800:]
    Tt = Traw[posTmax+800:]
    Tt[:] = [x-Tamb for x in Tt]
    fitT = numpy.log(Tt)
    eq, R = findslope(fitt,fitT)
    lnTfit[:] = [fitt[x]*eq[0]+eq[1] for x in range(len(fitt))]
    rsq = r2(fitT,lnTfit)
    
#Calc theoretical max T from inputs, used to get conversion
    TadCalc = H/(cp*mt)
    
#Fill arrays for Tcorrected and pNO (conversion)    
    for k in range(len(Traw)):
        Tcorr.append(Traw[k] - eq[0]*(Traw[k]-Tamb)*timDec[k]- Twork)
        pNCO.append(100*r*(Tcorr[k]-Tamb)/(TadCalc-Tamb))

#Writing the SCSV output file with all new converted t and p
###########future, search for any csv's with CONVERTED in title and skip them    
    CSVout = file.split('.')[0] + "_Corrected.csv"    
    with open(CSVout,'w', newline='') as f: 
        head = ['t','Traw' ,'Tcorr' ,'P']
        csvwrite = csv.DictWriter(f,fieldnames = head)
        csvwrite.writeheader()
        for j in range(len(Traw)):
            csvwrite.writerow({'t': timDec[j], 'Traw': Traw[j],'Tcorr': Tcorr[j],'P': pNCO[j]})

#If you want to see plot change to Y    
    pltT = 'N'
    if pltT == 'Y' :
        
        fig = plt.figure(figsize=(10,12))
        
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)
        
        ax1.plot(fitt,fitT,'.', markersize=2)
        #ax1.text(10,12,r'$R^2$= '  +str(rsq),fontsize = 20)
        ax1.set_xlabel('time (s)',fontsize =20)
        ax1.set_ylabel('ln(Tt-Tamb)',fontsize =20)
        ax1.plot(fitt,lnTfit,'k')
        
        ax2.plot (timDec,Tcorr)
        ax2.plot (timDec,Traw)
        ax2.set_xlabel('time (s)',fontsize =20)
        ax2.set_ylabel('Temperature ($^\circ$C)',fontsize =20)
        ax2.set_xlim(0,500)
        ax2.set_ylim(bottom = 0)
        
        ax3.plot (timDec,pNCO)
        ax3.set_xlabel('time (s)',fontsize =20)
        ax3.set_ylabel('Conversion (%)',fontsize =20)
        ax3.set_xlim(0,500)
        ax3.set_ylim(0,100)
        
        plt.show()
        plt.clf()




