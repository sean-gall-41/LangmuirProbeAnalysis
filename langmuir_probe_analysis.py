# created: 4/19/2018
# most recent update: 12/14/2019
# authors: Jim Reardon and Sean Gallogly

import numpy as np
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
#from matplotlib.pyplot import cm
from scipy.optimize import curve_fit
import os
import pandas as pd
#A flag used for any and all while loops 
flag = False 
#A linear fitting function, used (potentially) for fitting ion saturation 
#region of IV Characteristic
def h(m,x,b):
    return m*x+b
#The exponential function fit for the exponential region of the IV 
#Characteristic. It includes a constant offset so it does not matter whether 
#values are negative.
def ExpFunc(x,A,B,C,D):
    return A *np.exp(B *(x - C))- D

#Begin the Analysis Algorithm 

#Import necessary files
path = os.getcwd()
files = os.listdir(path)
files_xls = [f for f in files if f[-3:] == 'xls']
#print files_xls
#Use below if files_xls is not in order. Current algorithm designed to deal 
#with imported file ordering of (numbered by run #): 4, 5, 6, 7, 8, 9, 10, 11, 
# 12, 13, 0, 1, 2, 3  
i = 0
while not flag:
    if i < 10: #May have to change this value depending on the files order
        files_xls.append(files_xls[0])
        files_xls.remove(files_xls[0])
        i+=1
    else:
        flag = True
flag = False
#Create electron energy list before running through all datasets s.t. it is not 
#reset when each new data set is analyzed
Te_List = []
#Reads each file from files_xls and creates pandas dataframe containing all data
for f in files_xls:
    runNum = files_xls.index(f) + 15
    print "\nThis is run",runNum
    df0 = pd.DataFrame()
    df1 = pd.DataFrame()
    data0 = pd.read_excel(f,'Voltage - Dev1_ai0')
    data1 = pd.read_excel(f,'Voltage - Dev1_ai1')
    df0 = df0.append(data0)
    df1 = df1.append(data1) 
    #initialize lists which are easier to work with
    Idata = []
    Vdata = []
    #putting all data in current run into lists
    for i in np.arange(len(df0)):
        #Avoids unecessary info imported from LabView Express at the time of
        #taking raw data
        if i < 7:
            continue
        else:
            Idata.append(df0.iat[i,1])
            Vdata.append(df1.iat[i,1])
    #calculation of acutal current and voltage values across the resistor
    resistance = 215.0
    amplification = 30.0
    for l in np.arange(len(Vdata)):
        Vdata[l] = amplification*Vdata[l] - Idata[l]
        Idata[l] = Idata[l]/resistance
    #Plot for the user, prompt for input
    plt.plot(Vdata,Idata,'mo')
    plt.show()
    userinput = raw_input("Use Data?(y/n): ")
    print " "
    while userinput not in ['y','n']:
        print "Invalid Input. Try Again."
        userinput = raw_input("\nUse Data?(y/n): ")
    #Below user may opt out of analysis w/o having to view every 
    # characteristic curve from every run
    if userinput in ['n']:
        userinput = raw_input("Are you finished analyzing data?(y/n): ")
        while userinput not in ['y','n']:
            print "Invalid Input. Try Again."
            userinput = raw_input("Use Data?(y/n): ")
            print " " 
        if userinput in['n']:
            continue
        else:
            break
    #Moving onto the data analysis
    else:
        if Idata and Vdata:
            #Initialize the initial guess list for parameters A-D in exp 
            # function above
            cutoffs = []
            voltCutOff = True
            inputparam = True
            while voltCutOff:
                cutoffs = []
                userinput = float(raw_input("Enter V_max: "))
                cutoffs.append(userinput)
                userinput = float(raw_input("Enter V_min: "))
                cutoffs.append(userinput)
                i = 1
                #Below allows user to input guesses for initial values. Could 
                # not figure out how to enumerate using Latin letters, so I 
                # have used the associations 1 -> A , 2 -> B , 3 -> C , 4 -> D 
                # for parameters A-D in def above If poor initial guess entered
                # the code will catch the RunTime Error arising from the 
                # leastsq method in scipy/optimize/minpack.py. If the method is 
                # not able to calculate the covariance matrix an optimize 
                # warning should pop up, But I have only seen this when the 
                # incorrect sign for parameter D is entered. Have not written 
                # code to catch this warning however as when I realized the 
                # correct sign for parameter D the leastsq method converged, 
                # though visually not to a very accurate result. I suspect this 
                # issue could be managed by allowing the user to define their 
                # own Vmax and Vmin and allow this to occur until the user is 
                # satisfied visually with the fit. 
                if inputparam:
                    initialGuess = []
                    while len(initialGuess) < 4:
                        print "Enter initial guess for parameter ",i,": "
                        userinput = float(raw_input())
                        initialGuess.append(userinput)
                        i+=1
                #print initialGuess
                Vforfit = []
                Iforfit = []
                #Vmax = 7.0
                #Vmin = -40.0
                Vmax = cutoffs[0]
                Vmin = cutoffs[1]
                for i in np.arange(len(Vdata)):
                    if Vdata[i] - Vmax < 0.0:
                        if Vdata[i] - Vmin > 0.0:
                            Vforfit.append(Vdata[i])
                            Iforfit.append(Idata[i])
                #Below plots the raw data and initial guess for the exponential 
                # model so user can adjust initial guess accordngly if least 
                # squares fit does not converge
                while not flag:
                    guessedFactors = [ExpFunc(x,*initialGuess) for x in Vforfit]
                    fig1 = plt.figure(1)
                    ax1 = fig1.add_subplot(1,1,1)
                    ax1.plot(Vforfit,Iforfit,linewidth=0.1,linestyle='',
                             marker='o',color='g',label="Raw Data")
                    ax1.plot(Vforfit,guessedFactors,linestyle='-',color='b',
                             label="initial guess")
                    ax1.legend(loc=0,title="graphs",fontsize=11)
                    ax1.set_ylabel("I(A)")
                    ax1.set_xlabel("V(V)")
                    ax1.grid()
                    plt.show()
                    #Here we try to fit the function to the data using the 
                    # initial guess and if a RuntimeError is raised, let user 
                    # input new initial guess. This process will repeat until 
                    # the fitting algorithmn converges or an optimize warning 
                    # is raised, which I have not taken account of as of 
                    # 5292018. Still looking into the theory of estimation of 
                    #the Jacobian and why the covariance matrix would be 
                    # indeterminate for certain initial guesses.
                    try: 
                        popt,pcov = curve_fit(ExpFunc,Vforfit,Iforfit,
                                              initialGuess)
                    except RuntimeError: #FIXME: except optimize warning 
                        print '''Could not fit function given input parameters. 
                              Try some others'''
                        for i in np.arange(4):
                            print "Enter initial guess for parameter ",i+1,": "
                            userinput = float(raw_input())
                            initialGuess[i] = userinput
                        #print initialGuess
                    else:
                        #We show the full covariance matrix and plot the 
                        # results as well as show the current electron temp 
                        # and a list f all calculated electron temperatures so 
                        # far
                        print "Here is the covariance matrix:\n",pcov
                        VContent = np.linspace(cutoffs[1],cutoffs[0],20000)
                        fittedData = [ExpFunc(x,*popt) for x in VContent]
                        fig2 = plt.figure(2)
                        ax2 = fig2.add_subplot(1,1,1)
                        ax2.plot(Vforfit,Iforfit,linewidth=0.1,linestyle='',
                                 marker='o',color='g',label="Raw Data")
                        #ax2.plot(Vdata[n-2000:n+1500],guessedFactors,
                        #         linestyle='',marker='^',color='b',
                        #         label="initial guess")
                        ax2.plot(VContent,fittedData,linestyle='-',
                                 color='#900000',
                                 label='''fit with ({0:0.2g},{1:0.2g},{2:0.2g},
                                                    {3:0.2g})'''.format(*popt))
                        ax2.legend(loc=0,title="graphs",fontsize=11)
                        ax2.set_ylabel("I(A)")
                        ax2.set_xlabel("V(V)")
                        ax2.grid()
                        plt.show()
                        flag = True
                flag = False
                userinput = raw_input('''Would you like to analyze this run
                                      again, with different voltage bounds?
                                      (y/n): ''')
                while userinput not in ['y','n']:
                    print "Invalid Input. Try Again."
                    userinput = raw_input('''Would you like to analyze this run 
                                          again, with different voltage bounds?
                                          (y/n): ''')
                if userinput in ['y']:
                    userinput = raw_input('''Would you also like to use initial 
                                          parameters just previously input?
                                          (y/n): ''')
                    if userinput in ['y']:
                        inputparam = False
                    else:
                        inputparam = True
                        continue
                else:
                    Te = 1.0/popt[1]
                    #only append if last time analyzing this dataset
                    Te_List.append(round(Te,2))
                    print "Electron temperature Te = ",Te,"eV"
                    print "Electron Temps calc so far: ",Te_List
                    voltCutOff = False
print "\n","Program Terminates." 

