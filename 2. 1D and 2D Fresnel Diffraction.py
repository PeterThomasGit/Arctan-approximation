# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 00:28:28 2017

@author: Pete
"""
import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
style.use('fivethirtyeight')

def plot(x,y,title,xaxislabel, yaxislabel):
    """This function plots a graph.
    """
    fig, ax = plt.subplots()
    plt.xlabel(xaxislabel)#Labels x axis.
    plt.ylabel(yaxislabel)#labels y axis.
    plt.grid(True) #adds crid to graph
    plt.title(title)#Labels title
    plt.plot(x,y)
    plt.show()
    
def plotsave(x,y,title,xaxislabel, yaxislabel, Filename):
    """This function plots a graph.
    """
    fig, ax = plt.subplots()
    plt.xlabel(xaxislabel)
    plt.ylabel(yaxislabel)    
    plt.grid(True)
    plt.title(title)
    plt.plot(x,y)
    fig.savefig(Filename)#saves the figure locally.
    
def isfloat(x):
    """Checks object for its type value, this case it checks if object is a float."""
    try: 
        float(x)
    except ValueError: 
        return False
    return True

def isint(x):
    """Checks object for its type value, this case it checks if object is an integer."""
    try: 
        int(x)
    except ValueError: 
        return False
    return True

def iseven(x):
    "Checks the if the number is even. If it is even, it passes."
    if x % 2 == 0:
        return True # Even 
    else:
        return False # Odd

def sinfunc(x):
    """test function used for section a."""
    y=math.sin(x)
    return y

def eq5(xprime):
    """This is subject of integrand in equation 5 from the exercise sheet, it describes the intensity of waves
    along the x axis. Depending on its wavelength, size of screen and distance
    from aperture."""
    k=2*math.pi/wavelength
    eq = cmath.exp(((1j*k)/(2*z))*(x-xprime)**2)
    return eq

def eq4(yprime):
    """This is equation 4 from the exercise sheet, it describes the intensity of waves
    along the x-y plane. Depending on its wavelength, size of aperture, size of screen and distance
    from aperture. It is integrated using the composite simpson's rule."""
    k=2*math.pi/wavelength
    eq = cmath.exp(((1j*k)/(2*z))*(y-yprime)**2)
    return eq*(CompSimp(yApertureMin,yApertureMax,Nvalue,eq5))

def eq4tri(yprime):
    """This is the equation 4 form the exercise sheet, however this time it has
     y aperature limits as a function of the x axis, describing a triangular shape.
    """
    k=2*math.pi/wavelength
    eq = cmath.exp(((1j*k)/(2*z))*(y-yprime)**2)
    return eq*(CompSimp((yApertureMin+yprime)/2,(yApertureMax-yprime)/2,Nvalue,eq5))

def eq4circ(yprime):
    """This is the equation 4 form the exercise sheet, however this time it has
     y aperature limits as a function of the x axis, describing a circular shape.
    """
    k=2*math.pi/wavelength
    eq = cmath.exp(((1j*k)/(2*z))*(y-yprime)**2)
    return eq*(CompSimp(-(yApertureMin**2 - yprime**2)**(1/2),(yApertureMax**2 - yprime**2)**(1/2),Nvalue,eq5))

def CompSimp(a,b,n,func):
    """Integrates given function using the composite simpson's rule.It takes
    upper and lower limits, the function to be integrated and the number of
    points to be evaluated.
    """
    a=float(a)#ensures lower limit is a float
    b=float(b)#ensures lower limit is a float
    n=int(n)#ensures number of evaluated points is an integer.
    h = (b-a)/n # creates step size for composite simpson's integration
    evenlist, oddlist, firstlast = [], [], [] #creates three empty lists for respective types.
    for i in range(0,n+1,2):#creates even values for the range of n points
        Sum=a+(h*i)
        evenlist.append(func(Sum))#adds function evaluated about the generated even value to evenlist.

    firstlast.append(evenlist[0])# adds first even number from list to a firstlast list
    firstlast.append(evenlist[-1])# adds last even number from list to a firstlast list
    evenlist.pop(0)
    evenlist.pop(-1)
    for i in range(1,n,2):#creates even values for the range of n points.
        Sum = a+(h*i)
        oddlist.append(func(Sum))#adds function evaluated about the generated odd value to oddlist.
        
    finalsum = (h/3*(np.sum(firstlast)+(4*np.sum(oddlist))+2*np.sum(evenlist)))#compiles all necessary values to evaluate the composite simpson's rule
    return(finalsum)

def lb():
    """Creates a line break segment"""
    print('---------------------------------------------------------------------------------------------------------------------------------------')
    
    
def TrapInt(a,b,n,func):
    """Integrates given function using the Trapezium rule.
    """
    a=float(a)#ensures lower limit is a float
    b=float(b)#ensures lower limit is a float
    n=int(n)#ensures number of evaluated points is an integer.
    h=(b-a)/n
    Midlist=[]
    firstlastsum=[]    
    for i in range(0,n+1):
        Sum=a+(h*i)
        Midlist.append(func(Sum))    
    firstlastsum.append(Midlist[0] + Midlist[-1])#Adds first and last values to seperate list for easier manipulating.
    Midlist.pop(0)#removes first and last values for easier manipulating.
    Midlist.pop(-1)#removes first and last values for easier manipulating.
    finalsum=(h/2*(np.sum(firstlastsum)+2*np.sum(Midlist)))
    return(finalsum)





print("""\nWelcome to exercise 2 of the Level 5 Laboratory: Computational Physics tasks set for second year physics students at the University of Bristol in 2018.\n""")
print(""" This python programme is designed to allow the user to integrate functions (sin and Fresnel equation) via the composite Simpon's rule. \n""")
print("""Graphs can be plotted for the 1 dimensional Fresnel equation along the x axis. Users will also be able to image Fresnel diffraction in a 2D plane.""")
print(""" Using 3 diffrent apertures; square, triangular and circular.""")
lb()
print("""Part (a1) of the exercise allows the user to integrate the sin function between two limits and n intervals in the composite Simpson's rule. \n""")
print("""Part (a2) of the exercise allows the user to integrate the sin function between two limits and n intervals in the Trapezium rule. \n""")
print("""Part (b) of the exercise allows the user to evaluate a 1D Fresnel equation and plot the intensity of the x axis as a function of x.\n""")
print("""Part (c) of the exercise allows the user to evaluate a 2D Fresnel equation and plot a 2D image of Fresnel diffraction, through a choice of 3 apertures.\n""")
lb()  
  
MyInput = 0
while (MyInput != "q"):

    MyInput = input('Main menu: Enter a section, "a1", "a2", "b", "c", or "q" to quit: ')
#-------------------------------------------------------------------------------------------------------------------------------         
    if MyInput == "a1": #Part a.  
        lb()
        Minlim = "P"#sets object to an aberitrary value to allow error checking in the following while loops.
        Maxlim = "P"
        Nvalue = "P"
        print("This is part a of the exercise. 1d integration using the composite simpson's rule. Enter 'q' in any input to quit after the run of inputs has completed.") 
        while Minlim != "q" and Maxlim != "q" and Nvalue != "q":#keeps user in the desired place until they decide to exit.
            
            Minlim = input("Enter the lower bound of the integral limit: ")
            Maxlim = input("Enter the upper bound of the integral limit: ")
            Nvalue = input("Enter the a number of even points the composite simpson's rule will evaluate: ")
            if  isint(Minlim) and isfloat(Minlim) and isint(Maxlim) and isfloat(Maxlim) and isint(Nvalue):#checks values to ensure they are the correct value type.
                Nvalue=int(Nvalue)
                if iseven(Nvalue):
                    print(CompSimp(Minlim,Maxlim,Nvalue,sinfunc))
                    
                else:
                    print('\nThis is not a valid value.')
            elif Nvalue =="q" or Minlim =="q" or Maxlim =="q":
                print("\nYou have chosen to exit this section.")
                Minlim = "q"#stops the function from looping and exits.
                Maxlim = "q"
                Nvalue = "q"
                break
    
            else:
                print('\nThis is not a valid value.')
                
#-------------------------------------------------------------------------------------------------------------------------------         
    elif MyInput == "a2": #Part a.  
        lb()
        Minlim = "P"#sets object to an aberitrary value to allow error checking in the following while loops.
        Maxlim = "P"
        Nvalue = "P"
        print("This is part a of the exercise. 1d integration using the trapezium rule. Enter 'q' in any input to quit after the run of inputs has completed.") 
        while Minlim != "q" and Maxlim != "q" and Nvalue != "q":#keeps user in the desired place until they decide to exit.
            
            Minlim = input("Enter the lower bound of the integral limit: ")
            Maxlim = input("Enter the upper bound of the integral limit: ")
            Nvalue = input("Enter the a number of points the trapezium rule will evaluate: ")
            if  isint(Minlim) and isfloat(Minlim) and isint(Maxlim) and isfloat(Maxlim) and isint(Nvalue):#checks values to ensure they are the correct value type.
                
                print(TrapInt(Minlim,Maxlim,Nvalue,sinfunc))
                    
            elif Nvalue =="q" or Minlim =="q" or Maxlim =="q":
                print("\nYou have chosen to exit this section.")
                Minlim = "q"#stops the function from looping and exits.
                Maxlim = "q"
                Nvalue = "q"
                break
    
            else:
                print('\nThis is not a valid value.')

  #-------------------------------------------------------------------------------------------------------------------------------       
    elif MyInput == "b": #Part b. 
        lb()   
        Minlim = "P"
        print("This is part b of the exercise. Enter 'q' in any input to quit after the run of inputs has completed.")
        screensize= "P"
        aperturesize = "P"
        Nvalue = "P"
        z="P"
        wavelength="P"
        while screensize != "q" and aperturesize != "q" and z != "q" and Nvalue != "q" and wavelength != "q":
            wavelength = input("Enter the wavelength of light in nm: ")
            screensize = input("Enter the size of the screen in mm: ")
            aperturesize = input("Enter the size of the aperture in mm: ")
            z = input("Enter the distance 'z' between the aperture and screen in m: ")
            Nvalue = input("Enter the a number of even points the composite simpson's rule will evaluate: ")
            if isfloat(screensize) and isfloat(aperturesize) and isfloat(z) and isint(Nvalue) and isfloat(wavelength):
                wavelength = float(wavelength)/1000000000 #unit conversion
                screenleft =  -(float(screensize)/2000)#unit conversion and halves the size.
                screenright = (float(screensize)/2000)
                a = -(float(aperturesize)/2000)
                b = (float(aperturesize)/2000)
                z=float(z)
                Nvalue = int(Nvalue)
                if iseven(Nvalue):
                    xList = np.linspace(screenleft, screenright, 100)
                    XmodsqList = []
                    for x in xList:
                        XmodsqList.append(abs(CompSimp(a,b,Nvalue,eq5))**2)

                    plot(xList,XmodsqList, 'Intensity', 'x', '|X(x)|^2')
                    screensize = "P"
                    Q=input("Would you like to save the graph and data locally? y/n ")
                    while Q!="q":
                        if (Q=="y"):
                            Filename = input('Name the graph: ' )
                            Filename2 = input('Name the data text file: ')
                            plotsave(xList,XmodsqList, 'Intensity', 'x', '|X(x)|^2', Filename+".png")
                            np.savetxt(Filename2 + '.txt', (xList,XmodsqList) , newline='\r\n'+'\r\n', )
                            print("Saved")
                            Q=="q"
                            break
                        elif(Q=="n"):
                            print("\nTaking you back to the start of section b.")
                            Q=="q"
                            break
                            
                        else:
                            print('\nThis is not within the set boundary conditions')
                    
                else:    
                    print('\nThis is not a valid value.')

            elif screensize == "q" or aperturesize == "q" or z == "q" or Nvalue == "q" or wavelength == "q":
                print("\nYou have chosen to exit this section.")
                break

            else:
                print('\nThis is not a valid value.')

    #-------------------------------------------------------------------------------------------------------------------------------       
    elif MyInput == "c": #Part c. 
        lb()
        Minlim = "P"
        print("This is part c of the exercise. Enter 'q' in any input to quit after the run of inputs has completed.")
        screensize= "P"

        aperturesize = "P"
        Nvalue = "P"
        z = "P"
        shape = "P"
        wavelength = "P"
        while shape != "q":
            shape="P"
            shape=input("What shape aperture would you like to simulate; square, triangle, or circle? ")
            if shape=="square" or shape=="triangle" or shape=="circle":
                while screensize != "q" and aperturesize != "q" and z != "q" and Nvalue != "q" and wavelength != "q":
                    wavelength = input("Enter the wavelength of light in nm: ")
                    screensize = input("Enter the size of the screen in mm: ")
                    aperturesize = input("Enter the size of the aperture in mm: ")
                    z = input("Enter the distance 'z' between the aperture and screen in mm: ")
                    Nvalue = input("Enter the a number of even points the composite simpson's rule will evaluate: ")
                    if isfloat(screensize) and isfloat(aperturesize) and isfloat(z) and isint(Nvalue) and isint(wavelength):
                        wavelength = float(wavelength)/1000000000
                        xScreenMin =  -(float(screensize)/2000)    
                        xScreenMax = (float(screensize)/2000)
                        yScreenMin =  -(float(screensize)/2000)    
                        yScreenMax = (float(screensize)/2000)
                        xApertureMin = -(float(aperturesize)/2000)
                        xApertureMax = (float(aperturesize)/2000)
                        yApertureMin = -(float(aperturesize)/2000)
                        yApertureMax = (float(aperturesize)/2000)
                        Nvalue = int(Nvalue)
                        z=float(z)/1000
                        if iseven(Nvalue):
                            
                            NumPoints = 100
                            dx=(xScreenMax-xScreenMin)/(NumPoints-1)
                            dy=(yScreenMax-yScreenMin)/(NumPoints-1)
                            intensity = np.zeros( (NumPoints,NumPoints) )         
                            k=2*math.pi/wavelength
                            E=1
                            for i in range(NumPoints):
                                x = xScreenMin + (i*dx)
                                for j in range(NumPoints):
                                    y= yScreenMin + (j*dy)
                                    if shape == "square":
                                        intensity[i,j]=abs(((k*E)/(2*math.pi*z))*CompSimp(xApertureMin,xScreenMax,Nvalue,eq4))**2
                                    elif shape == "triangle":
                                        intensity[i,j]=abs(((k*E)/(2*math.pi*z))*CompSimp(xApertureMin,xScreenMax,Nvalue,eq4tri))**2
                                    elif shape == "circle":
                                        intensity[i,j]=abs(((k*E)/(2*math.pi*z))*CompSimp(xApertureMin,xScreenMax,Nvalue,eq4circ))**2
                                        
                            plt.imshow(intensity)
                            plt.title('Fresnel Diffraction for a ' + shape + ' aperature.')
                            plt.show()   
                            Q=input("Would you like to save the image locally? y/n ")
                            
                            while Q!="q":
                                
                                if (Q=="y"):
                                    Filename = input('Name the file: ' )
                                    plt.imsave(Filename+".png", intensity)
                                    print("Saved")    
                                    Q=="q"
                                    break
                                
                                elif(Q=="n"):
                                    print("\nTaking you back to the start of section c.")
                                    Q=="q"
                                    break
                            
                                else:
                                    print('\nThis is not within the set boundary conditions')
                                    
                            break
                        else:
                            print('\nThis is not a valid value.')
                            
                    elif screensize == "q" or aperturesize == "q" or z == "q" or Nvalue == "q" or wavelength == "q":
                        print("\nYou have chosen to exit this section.") 
                        
                        shape = "q"
                        
                        break
                    else:
                        print('\nThis is not a valid value.')     
                        
            elif shape == "q" :
                        print("\nYou have chosen to exit this section.")
                        break
                    
            else:
                print('\nThis is not a valid value.') 
            

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        