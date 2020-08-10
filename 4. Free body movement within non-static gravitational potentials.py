# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 19:59:58 2019

@author: Pete
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Mar 03 11:38:13 2019

@author: Pete
"""
    
#import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
import math
style.use('_classic_test') # Changes the style of the graphs

def plot2(x, y,A,B, title, xaxislabel, yaxislabel):
    """This function plots a graph, labels axis and title.
    """
    #fig, ax = plt.subplots()
    #plt.figure(figsize=(10,20))
    fig = plt.figure(figsize=(15,20))
    ax1 = fig.add_subplot(222)
    ax1.axis('equal')
    #ax1.ticklabel_format(axis='both', style='sci')
    plt.xlabel(xaxislabel) # Labels x axis.
    plt.ylabel(yaxislabel) # Labels y axis.
    circle1 = plt.Circle((0, 0), radius = Mass1radius, fc='b')
    circle2 = plt.Circle((MoonPositionXMoving, MoonPositionYMoving), radius = Mass2radius, fc='y')
    plt.gca().add_patch(circle1)
    plt.gca().add_patch(circle2)
    plt.grid(True) 
    plt.title( title + "\n")
    #plt.legend()
    
    plt.plot(x,y,'r')
    plt.plot(A, B, 'g')
    plt.show()


def plot(x, y, title, xaxislabel, yaxislabel,section):
    """This function plots a graph, labels axis and title.
    """
    #fig, ax = plt.subplots()
    plt.figure(figsize=(10,20))
    fig = plt.figure(figsize=(10,20))
    ax1 = fig.add_subplot(212)
    ax1.axis('equal')
    ax1.ticklabel_format(axis='both', style='sci')
    plt.xlabel(str(xaxislabel)) # Labels x axis.
    plt.ylabel(str(yaxislabel)) # Labels y axis.
    circle1 = plt.Circle((0, 0), radius = Mass1radius, fc='b')
    plt.gca().add_patch(circle1)
    plt.grid(True)
    plt.title(str(title))
    
    if section != 'a':
        MoonPositionStaticX = 0 #Postion of M2 in x direction
        MoonPositionStaticY = 384400 *1000 #Postion of M2 in y direction
        circle2 = plt.Circle((MoonPositionStaticX, MoonPositionStaticY), radius = Mass2radius, fc='y')
        plt.gca().add_patch(circle2)
        
    plt.plot(x,y,'r')   
    plt.show()

def CreateWrite(x, y, t, Filename):
    """Creates a file, names it, transposes arrays and writes them to a file. 
    The order in the file is the same as the order input into the function."""
    np.savetxt(Filename, np.transpose([x,y,t]) , newline='\r\n', )
    
def CreateWrite2(x, t, Filename):
    """Creates a file, names it, transposes arrays and writes them to a file. 
    The order in the file is the same as the order input into the function."""
    np.savetxt(Filename, np.transpose([x,t]) , newline='\r\n', )
    


G = 6.67 * 10 **(-11)
M1 = 5.972 * 10**24
Mass1radius = 6371*1000
M2 = 7.3* 10 **22
Mass2radius = 1735.5 * 1000


def isfloat(x):
    """Checks if object *can* be converted to a float, if it can- it returns true."""
    try: 
        float(x)
    except ValueError: 
        return False
    return True

counter = []

def AcceqX2Mmoving(X,Y):
    #global MoonPositionX, MoonPositionY
    return -(G*M1*X/(X**2 + Y**2)**(3/2)) - (G*M2*(X-MoonPositionXMoving) / ((X - MoonPositionXMoving)**(2) + (Y - MoonPositionYMoving)**(2))**(3/2))

def AcceqY2Mmoving(X,Y):
    #global MoonPositionX, MoonPositionY
    return -(G*M1*Y/(X**2 + Y**2)**(3/2)) - (G*M2*(Y-MoonPositionYMoving) / ((X - MoonPositionXMoving)**(2) + (Y - MoonPositionYMoving)**(2))**(3/2))


def AcceqX2M(X,Y):

    #global MoonPositionStaticX, MoonPositionStaticY
    return -(G*M1*X/(X**2 + Y**2)**(3/2)) - (G*M2*(X-MoonPositionStaticX) / ((X - MoonPositionStaticX)**(2) + (Y - MoonPositionStaticY)**(2))**(3/2))

def AcceqY2M(X,Y):
    #global MoonPositionStaticX, MoonPositionStaticY
    return -(G*M1*Y/(X**2 + Y**2)**(3/2)) - (G*M2*(Y-MoonPositionStaticY) / ((X - MoonPositionStaticX)**(2) + (Y - MoonPositionStaticY)**(2))**(3/2))

def AcceqX(X,Y):
    return -G*M1*X/(X**2 + Y**2)**(3/2)

def AcceqY(X,Y):
    return -G*M1*Y/(X**2 + Y**2)**(3/2)

  

def RK(x0,y0,AccelerationEquationX,AccelerationEquationY,Veloc,Ang,section,MaxTime,h):
    
    #global OrbitCounter, G, M1, M2, Mass2radius, Mass1radius, MoonPositionX, MoonPositionY
    #OrbitCounter = -1
    #MaxTime = 27*24*3600/2
    Nsteps = int(MaxTime/h)

    global x,y,t,MoonPositionX,MoonPositionY,TotalEnergy
    x = np.zeros(Nsteps)
    y = np.zeros(Nsteps)
    vx = np.zeros(Nsteps)
    vy = np.zeros(Nsteps)
    t = np.zeros(Nsteps)
    TotalEnergy = np.zeros(Nsteps)

    x[0] = x0
    y[0] = y0
    vx[0] = Veloc * math.cos(Ang*math.pi/180)
    vy[0] = Veloc * math.sin(Ang*math.pi/180)
    t[0]=0
    TotalEnergy[0]= - G*M1*28881/(2*(x[0]**2+y[0]**2)**(1/2))
    
    Mphase= math.pi * 233/1000
    MoonPositionX = np.zeros(Nsteps)
    MoonPositionY = np.zeros(Nsteps)
    MV = 1035.34655 #2 pi r/T moons orbital velocity
    MoonPositionX[0] = 384400 *1000 * math.cos(((t[0])*MV/(384400 *1000))+Mphase)
    MoonPositionY[0] = 384400 *1000 * math.sin(((t[0])*MV/(384400 *1000))+Mphase)
    
    
    for i in range (1,Nsteps):
        
        t[i] = t[i-1] + h   
        
        #if section == 'c':
        global MoonPositionXMoving, MoonPositionYMoving
        if section == 'c':
            MoonPositionX[i] = 384400 *1000 * math.cos(((t[i])*MV/(384400 *1000))+Mphase)
            MoonPositionY[i] = 384400 *1000 * math.sin(((t[i])*MV/(384400 *1000))+Mphase)
            MoonPositionXMoving = MoonPositionX[i]
            MoonPositionYMoving = MoonPositionY[i]
            #global MoonPositionXMoving,MoonPositionYMoving
        elif section == 'b':
            global MoonPositionStaticX, MoonPositionStaticY
            MoonPositionStaticX = 0 #Postion of M2 in x direction
            MoonPositionStaticY = 384400 *1000 #Postion of M2 in y direction
            
        k1x = vx[i-1]
        k1y = vy[i-1]
        k1vx = AccelerationEquationX(x[i-1] , y[i-1])
        k1vy = AccelerationEquationY(x[i-1] , y[i-1])
        
        k2x = vx[i-1] + (h * (k1vx/2))
        k2y = vy[i-1] + (h * (k1vy/2))
        k2vx = AccelerationEquationX((x[i-1] + (h * (k1x/2))) , (y[i-1] + (h * (k1y/2))))
        k2vy = AccelerationEquationY((x[i-1] + (h * (k1x/2))) , (y[i-1] + (h * (k1y/2))))
        
        k3x = vx[i-1] + (h * (k2vx/2))
        k3y = vy[i-1] + (h * (k2vy/2))
        k3vx = AccelerationEquationX((x[i-1] + (h * (k2x/2))) , (y[i-1] + (h * (k2y/2))))
        k3vy = AccelerationEquationY((x[i-1] + (h * (k2x/2))) , (y[i-1] + (h * (k2y/2))))
        
        k4x = vx[i-1] + (h * k3vx)
        k4y = vy[i-1] + (h * k3vy)
        k4vx = AccelerationEquationX((x[i-1] + (h * k3x)) , (y[i-1] + (h * k3y)))
        k4vy = AccelerationEquationY((x[i-1] + (h * k3x)) , (y[i-1] + (h * k3y)))
        
        
        
        x[i] = x[i-1] + (h/6) * (k1x + 2 * k2x + 2 * k3x + k4x)
        y[i] = y[i-1] + (h/6) * (k1y + 2 * k2y + 2 * k3y + k4y)
        vx[i] = vx[i-1] + (h/6) * (k1vx + (2 * k2vx) + (2 * k3vx) + k4vx)
        vy[i] = vy[i-1] + (h/6) * (k1vy + (2 * k2vy) + (2 * k3vy) + k4vy)
        
        TotalEnergy[i]= - G*M1*28881/(2*(x[i]**2+y[i]**2)**(1/2))
        
        if ((x[i]**2 + y[i]**2)**(1/2)) <  Mass1radius or ((x[i-1]**2 + y[i-1]**2)**(1/2)) <  Mass1radius:
            print("The Rocket has crashed into Earth at: " + str(t[i]) + 'seconds.')
            break
        


        elif section == 'b':
            if (((x[i]-MoonPositionStaticX)**2 + (y[i]-MoonPositionStaticY)**2)**(1/2)) <  Mass2radius or (((x[0]-MoonPositionStaticX)**2 + (y[0]-MoonPositionStaticY)**2)**(1/2)) <  Mass2radius:
                print("The Rocket has crashed into Moon at: " + str(t[i]) + 'seconds.')
                break
        elif section =='c':
            if (((x[i]-MoonPositionX[i])**2 + (y[i]-MoonPositionY[i])**2)**(1/2)) <  Mass2radius or (((x[0]-MoonPositionX[i])**2 + (y[0]-MoonPositionY[i])**2)**(1/2)) <  Mass2radius:
                print("The Rocket has crashed into Moon at: " + str(t[i]) + 'seconds.')
                break
            
    x = np.trim_zeros(np.array(x), 'b')
    y = np.trim_zeros(np.array(y), 'b')
    vx = np.trim_zeros(np.array(vx), 'b')
    vy = np.trim_zeros(np.array(vy), 'b')
    t = np.trim_zeros(np.array(t), 'b')
    TotalEnergy =np.trim_zeros(np.array(TotalEnergy),'b')
    if section =='c':
        MoonPositionX = np.trim_zeros(np.array(MoonPositionX), 'b')
        MoonPositionY = np.trim_zeros(np.array(MoonPositionY), 'b')
        plot2(x,y,MoonPositionX,MoonPositionY,"title", "x position [m]", "y position [m]")
    elif section =='a':
        plot(x,y,"Plot of object orbiting Earth", "X Position in m", "Y Position in m",section)
        #CreateWrite(vx,vy,t, "Section " + str(MyInput) +" Initial velocity "+str(InitialVelocity) + "InitialAngle " + str(InitialAngle) + " TimeStep " + str(TimeStep) +"RocketVelocity"+ ".txt")
        
    elif section == 'b':
        plot(x,y,"Plot of object orbiting Earth and moon", "X Position in m", "Y Position in m",section)






MyInput = 0
while (MyInput != "q"):
    
    MyInput = input('Main menu: Enter a section, "a", "b", "c" Or "q" to quit: ')
    
    if MyInput == "a":
        
        XInitial = "P"
        YInitial = "P"
        InitialVelocity = "P"
        InitialAngle = "P"
        TimeStep = "P"
        MaxTime = "P"
        
        print("""\n This is part 'a' of the exercise.\n
              Here you will be able to simulate a rocket under the effects of earth as a fixed gravitational potential.""")
        
        while XInitial != "q" and YInitial != "q" and InitialVelocity != "q" and InitialAngle != "q" and TimeStep != "q": # Keeps user in the section, until they enter "q" into any/all the inputs below, to allow them to quit.
            
            XInitial = input(""" Please enter the starting X coordinate for your rocket in kilometers, e.g 0: """)
            YInitial = input(""" Please enter the starting Y coordinate for your rocket in kilometers, e.g -7000: """)
            InitialVelocity = input(""" Please enter the starting Velocity for your rocket in meters per second, e.g 7543.515664: """)
            InitialAngle = input(""" Please enter the starting angle for your rocket in degrees, e.g 0: """)
            TimeStep = input(""" Please enter a value for the time increment you wish to use for the calculations in seconds, e.g 10: """)
            MaxTime = input(""" Please enter a value for the maximum time you wish to simulate for in seconds, e.g 300000: """)            
            
            
            if isfloat(XInitial) and isfloat(YInitial)and isfloat(InitialVelocity)  and isfloat(InitialAngle) and isfloat(TimeStep) and isfloat(MaxTime): # Checks to see if inputs can be converted to a float, to avoid incorrect object types creating errors.
                
                XInitial = float(XInitial)
                YInitial = float(YInitial) 
                InitialVelocity = float(InitialVelocity)
                InitialAngle = float(InitialAngle)
                TimeStep = float(TimeStep)
                MaxTime = float(MaxTime)
                
                
                RK(XInitial*1000,YInitial*1000,AcceqX,AcceqY,InitialVelocity,InitialAngle,str(MyInput),MaxTime,TimeStep)
                
            elif  XInitial == "q" or YInitial == "q" or InitialVelocity == "q" or InitialAngle == "q" or TimeStep == "q":
                print("\nYou have chosen to exit this section.")
                
            else:
                print('\nThis is not a valid value.') 
                
            SaveQ = input("Would you like to save this data to a txt file? y/n ")
            if SaveQ  ==  'y':

                #CreateWrite(vx,vy,t, "Section" + str(MyInput) +"Initial velocity"+str(InitialVelocity) + "InitialAngle" + str(InitialAngle) + "TimeStep" + str(TimeStep) +"RocketVelocity"+ ".txt")
                CreateWrite(x,y,t, "Section" + str(MyInput) +"InitialXposition"+str(XInitial)+"InitialYposition" + str(YInitial)+ "InitialVelocity" + str(InitialVelocity) + "InitialAngle" + str(InitialAngle) + "TimeStep" + str(TimeStep) +"RocketMovement"+ ".txt")
                CreateWrite2(t,TotalEnergy, "Section" + str(MyInput) +"InitialXposition"+str(XInitial)+"InitialYposition" + str(YInitial)+ "InitialVelocity" + str(InitialVelocity) + "InitialAngle" + str(InitialAngle) + "TimeStep" + str(TimeStep) + "TE"+".txt")
                break
            
            else:
                break
        
        
        
    elif MyInput == "b":
        
        XInitial = "P"
        YInitial = "P"
        InitialVelocity = "P"
        InitialAngle = "P"
        TimeStep = "P"
        MaxTime = "P"
        
        print("""\n This is part 'b' of the exercise.\n
              Here you will be able to simulate a rocket under the effects of earth and the moon as two fixed gravitational potentials.""")
        
        while XInitial != "q" and YInitial != "q" and InitialVelocity != "q" and InitialAngle != "q" and TimeStep != "q": # Keeps user in the section, until they enter "q" into any/all the inputs below, to allow them to quit.
            
            XInitial = input(""" Please enter the starting X coordinate for your rocket in kilometers, e.g 0: """)
            YInitial = input(""" Please enter the starting Y coordinate for your rocket in kilometers, e.g -7000: """)
            InitialVelocity = input(""" Please enter the starting Velocity for your rocket in meters per second, e.g 10557.5: """)
            InitialAngle = input(""" Please enter the starting angle for your rocket in degrees, e.g 0: """)
            TimeStep = input(""" Please enter a value for the time increment you wish to use for the calculations in seconds, e.g 10: """)
            MaxTime = input(""" Please enter a value for the maximum time you wish to simulate for in seconds, e.g 900000: """)            
            
            
            if isfloat(XInitial) and isfloat(YInitial)and isfloat(InitialVelocity)  and isfloat(InitialAngle) and isfloat(TimeStep) and isfloat(MaxTime): # Checks to see if inputs can be converted to a float, to avoid incorrect object types creating errors.
                
                XInitial = float(XInitial)
                YInitial = float(YInitial) 
                InitialVelocity = float(InitialVelocity)
                InitialAngle = float(InitialAngle)
                TimeStep = float(TimeStep)
                MaxTime = float(MaxTime)
                
                
                RK(XInitial*1000,YInitial*1000,AcceqX2M,AcceqY2M,InitialVelocity,InitialAngle,str(MyInput),MaxTime,TimeStep)
                
            elif  XInitial == "q" or YInitial == "q" or InitialVelocity == "q" or InitialAngle == "q" or TimeStep == "q":
                print("\nYou have chosen to exit this section.")
                
            else:
                print('\nThis is not a valid value.') 

            SaveQ = input("Would you like to save this data to a txt file? y/n ")
            if SaveQ  ==  'y':

                CreateWrite(x,y,t, "Section" + str(MyInput) +"InitialXposition"+str(XInitial)+"InitialYposition" + str(YInitial)+ "InitialVelocity" + str(InitialVelocity) + "InitialAngle" + str(InitialAngle) + "TimeStep" + str(TimeStep) + "RocketMovement" ".txt")
                break
            
            else:
                break
        
        

    elif MyInput == "c":
        
        XInitial = "P"
        YInitial = "P"
        InitialVelocity = "P"
        InitialAngle = "P"
        TimeStep = "P"
        MaxTime = "P"
        
        print("""\n This is part 'c' of the exercise.\n
              Here you will be able to simulate a rocket under the effects of earth as a fixed gravitational potential and the moon as a moving gravitational potential.""")
        
        while XInitial != "q" and YInitial != "q" and InitialVelocity != "q" and InitialAngle != "q" and TimeStep != "q": # Keeps user in the section, until they enter "q" into any/all the inputs below, to allow them to quit.
            
            XInitial = input(""" Please enter the starting X coordinate for your rocket in kilometers, e.g 0: """)
            YInitial = input(""" Please enter the starting Y coordinate for your rocket in kilometers, e.g -7000: """)
            InitialVelocity = input(""" Please enter the starting Velocity for your rocket in meters per second, e.g 10600: """)
            InitialAngle = input(""" Please enter the starting angle for your rocket in degrees, e.g 0: """)
            TimeStep = input(""" Please enter a value for the time increment you wish to use for the calculations in seconds, e.g 10: """)
            MaxTime = input(""" Please enter a value for the maximum time you wish to simulate for in seconds, e.g 900000: """)            
            
            
            if isfloat(XInitial) and isfloat(YInitial)and isfloat(InitialVelocity)  and isfloat(InitialAngle) and isfloat(TimeStep) and isfloat(MaxTime): # Checks to see if inputs can be converted to a float, to avoid incorrect object types creating errors.
                
                XInitial = float(XInitial)
                YInitial = float(YInitial) 
                InitialVelocity = float(InitialVelocity)
                InitialAngle = float(InitialAngle)
                TimeStep = float(TimeStep)
                MaxTime = float(MaxTime)
                
                
                RK(XInitial*1000,YInitial*1000,AcceqX2Mmoving,AcceqY2Mmoving,InitialVelocity,InitialAngle,str(MyInput),MaxTime,TimeStep)
                
            elif  XInitial == "q" or YInitial == "q" or InitialVelocity == "q" or InitialAngle == "q" or TimeStep == "q":
                print("\nYou have chosen to exit this section.")
                
            else:
                print('\nThis is not a valid value.') 

            SaveQ = input("Would you like to save this data to a txt file? y/n ")
            if SaveQ  ==  'y':

                CreateWrite(x,y,t, "Section" + str(MyInput) +"InitialXposition"+str(XInitial)+"InitialYposition" + str(YInitial)+ "InitialVelocity" + str(InitialVelocity) + "InitialAngle" + str(InitialAngle) + "TimeStep" + str(TimeStep) +" RocketMovement" +".txt")
                CreateWrite(MoonPositionX,MoonPositionY,t, "Section" + str(MyInput)+"InitialXposition"+str(XInitial)+"InitialYposition" + str(YInitial)+ "InitialVelocity" + str(InitialVelocity) + "InitialAngle" + str(InitialAngle) + "TimeStep" + str(TimeStep) +" MoonMovement" +".txt")
                break
            
            else:
                break
        
        
    elif MyInput == 'q':
        print("\nYou have chosen to exit this program.")
        
    else:
        print('\nThis is not a valid value.') 

        




























    