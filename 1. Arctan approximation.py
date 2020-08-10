# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 15:09:32 2017
@author: Peter Thomas
"""
import matplotlib.pyplot as plt
import numpy as np
import math
import pylab
     
def MyArctanAll(Xvalue, Nvalue):
    """
    This function returns the sum of the arctan Taylor series expansion with N orders of terms about the argument x, for all x.
    Starts at N=0 and goes to Nvalue sums(the input).
    """
    N=0 #ensures the sum starts at N = 0
    a=0 # a is the the total running sum of the taylor series
    Xvalue=float(Xvalue)#the following sorts the arguments x to their required conditional equation based on its value.
    if (Xvalue>1) :
        Xvalue = Xvalue**(-1)# Converts x into a float and its inverse.
        while (N<=Nvalue):#runs the sum for each value of order from 0 to Nvalue so the result is the total sum of these orders.
            b = (((-1) ** N)/((2 * N) + 1))* (Xvalue **((2 * N) + 1)) #The taylor series for arctan where x is the argument x and Nvalue is the argument N.
            a+=b #adds each term of the taylor series onto the sum of the previous series.
            N+=1 #Increases N until it reaches the limit of Nvalue.
            c = ((math.pi/2) - a)            
        return c# Returns the total sum of the series.
        
    elif (Xvalue<-1):
        Xvalue = float(Xvalue)**(-1) # Converts x into a float and its inverse.
        while (N<=Nvalue):#runs the sum for each value of order from 0 to Nvalue so the result is the total sum of these orders.
            b = (((-1) ** N)/((2 * N) + 1))* (Xvalue **((2 * N) + 1)) #The taylor series for arctan where x is the argument x and Nvalue is the argument N.
            a+=b #adds each term of the taylor series onto the sum of the previous series.
            N+=1 #Increases N until it reaches the limit of Nvalue.
            c = -(math.pi/2) - a
        return c# Returns the total sum of the series.
    
    elif (Xvalue==0) :
        return 0
    
    elif (abs(Xvalue)<=1):
        Xvalue=float(Xvalue)
        while (N<=Nvalue):#runs the sum for each value of order from 0 to Nvalue so the result is the total sum of these orders.
            b = (((-1) ** N)/((2 * N) + 1))* (Xvalue **((2 * N) + 1)) #The taylor series for arctan where x is the argument x and Nvalue is the argument N.
            a+=b #adds each term of the taylor series onto the sum of the previous series.
            N+=1 #Increases N until it reachs the limit of Nvalue.
        return a # Returns the total sum of the series.

def plot(x, y, ytrue,Filename='results'):
    """This function plots and saves a graph, with one array of xvalues and 2 arrays of y values.
    """
    fig, ax = plt.subplots()
    plt.xlabel('x')
    plt.ylabel('arctan(x)')    
    plt.grid(True)
    plt.title('Taylor series approximation')
    plt.plot(x,y, label = 'Order '+ str(Nvalue))
    plt.plot(x,ytrue, label = 'arctan(x)')
    plt.legend()
    plt.show()
    fig.savefig(Filename)
    
def CreateWrite(x,y,z,Filename):
    """Creates a file, names it, and writes three numpy arrays to file. The order in the file is the same as the order input into the function."""
    np.savetxt(Filename, (x,y,z) , newline='\r\n'+'\r\n', )
    
def lb():
    """Creates a line break segment""" #this allows the console to be read easier and to distinguish between sections.
    print('---------------------------------------------------------------------------------------------------------------------------------------')

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

lb()
print("""\nWelcome to exercise 1 of the Level 5 Laboratory: Computational Physics tasks set for second year physics students at the University of Bristol in 2017.\n
This python programme is designed to allow the user to approximate the arctan(x) trigonometric function via a Taylor series expansion.
The function can be further used to plot graphs and save data for comparing the function against the built in arctan(x) function between the range -2<=x<=2,
and to find the lowest order of the expansion that can approximate Pi to 7 significant figures.\n""")
lb()
print("""Part (a) allows the user to produce approximations of arctan(x) using a Taylor series expansion, about the argument 'x' of order N. Where x belongs to any real number and N is user defined. \n""")
print("""Part (d) allows the user to produce approximations of arctan(x) using a Taylor series expansion, about the argument 'x' of order N . Where x is between the range -2<=x<=2 and N is user defined. This data is plotted into a graph and saved locally, the raw data can then be printed to the screen or saved to a file locally.\n""")
print("""Part (e) allows the user to use the Taylor series expansion of arctan(1), to find the lowest order of the expansion that can calculate Pi to 7 significant figures.
Given that arctan(1)=Pi/4. This will search through a range of N orders defined by the user (Part 1). There is also an option to save all the data locally so one can import said data into the graphing software of choice if further analyses is required(Part 2).\n""")

MyInput = 0
while (MyInput != "q"):#creates a loop so you stay in the menu regardless of incorrect inputs and allows the exit of the progromme via the input "q"
    lb()
    MyInput = input('Main menu: Enter a section, "a", "d", "e" or "q" to quit: ')
#-------------------------------------------------------------------------------------------------------------------------------         
    if MyInput == "a": #Part a.   
        print("""\nWelcome to part (a) of this excercise. To exit, enter "q".""")     
        #Takes user input values for arguments x and order N for extended taylor series of arctan.  
        Xvalue='a'
        Nvalue='a'
        while Xvalue != "q":             
            Xvalue = input('Please enter a value for the argument "x" where x is any real number (or "q" to return to the main menu): ')
            if isint(Xvalue) or isfloat(Xvalue): #checks to make sure input is number
                Xvalue=float(Xvalue)
                while Nvalue != "q":
                    Nvalue = input('Please enter an integer value "N" for the order of the expansion (or "q" to return to the main menu): ')#The user input for order N in the taylor series.
                    if isint(Nvalue):#checks to make sure input is integer
                        Nvalue = int(Nvalue)#Force N to be integer as the order can only be a whole number.
                        print('The sum of my arctan taylor series is = ', MyArctanAll(Xvalue, Nvalue))
                        print('The built in function of arctan evaluates to give =',pylab.arctan(Xvalue))
                        break #finishes section in order to repeat if required
                    elif Nvalue=='q':
                        print("\nYou have chosen to exit this section.")
                        Xvalue='q'#this stops the while loop from activating again when you want to exit
                        break
                    else:
                        print('\nThis is not within the set boundary conditions')
            elif Xvalue=='q':
                print("\nYou have chosen to exit this section.")
                break
            else:
                print('\nThis is not within the set boundary conditions')
#-------------------------------------------------------------------------------------------------------------------------------          
    elif MyInput == "d": #Part d.
        print("""\nWelcome to part (d) of this exerise. \n
Here you will be asked to give the order for which my defined taylor series runs for.\n
This is then plotted for a number of points, overlayed with the plot of the built in arctan(x) function over the same arguments for x. \n
The values for x and y in my defined function are also printed to screen or text if chosen and the graph saved.""")

        Nvalue="a"
        Npoints="a"#points on graph
        while Nvalue!="q":  
            y=[] #Empty array that is filled when evaluating MyArctanAll and for loop is active.
            Nvalue = input('Please enter a value for the order of the expansion, "N" (or "q" to return to the main menu): ')#The user input for order N in the taylor series.
            if isint(Nvalue):#checks to make sure input is number
               Nvalue = int(Nvalue)#Force N to be integer as the order can only be a whole number.
               while Npoints != "q":  
                   Npoints=input('How many data points would you like to evaluate?(or enter "q" to return to the main menu): ' )#for choosing number of points in the graph
                   if isint(Npoints):
                       Npoints=int(Npoints)
                       Lrange, Hrange = (-2), 2 #xlimits
                       x=np.linspace(Lrange, Hrange, Npoints)#creates array for Npoints' worth of values of x between the two ranges.
                       ytrue=pylab.arctan(x) #y values for built in function of arctan(x)
                       for Xvalue in x:   
                           y.append(MyArctanAll(Xvalue, Nvalue))#Appends y array with function outputs
                       Filename = input('Name the file: ' )
                       plot(x,y,ytrue,Filename)
                       #print('\n These are my values for x:\n',x,'\n','\n These are my values for y:\n',y)
                       Q = input('Would you like to save the data sets for x and y for my function and the built in y value to a text file? Please input y/n: ')
                       while Q!="q":
                           if (Q=="y"):
                               Filename = input('Name the file: ' )
                               CreateWrite(x,y,ytrue,Filename)
                               print("\nThe order of the data is: First, the x axis data. Second, the y axis data for my defined function for arctan(x). Lastly, the y axis values for the built-in arctan(x) function.")
                               Q=="q"#this stops the while loop from activating again when you want to exit
                               Nvalue="q"#this stops the while loop from activating again when you want to exit
                               Npoints="q"#this stops the while loop from activating again when you want to exit
                               break
                           elif(Q=="n"):
                               print("\nYou have chosen to exit this section.")
                               Q=="q"#this stops the while loop from activating again when you want to exit
                               Nvalue="q"#this stops the while loop from activating again when you want to exit
                               Npoints="q"#this stops the while loop from activating again when you want to exit
                               break
                           else:
                               print('\nThis is not within the set boundary conditions')
                   elif Npoints=='q':
                       print("\nYou have chosen to exit to the start of this section.")
                       Nvalue='q'#this stops the while loop from activating again when you want to exit
                       break
                   else:
                       print('\nThis is not within the set boundary conditions')

#-------------------------------------------------------------------------------------------------------------------------------            
    elif MyInput == "e": #Part e.
        print("""\nWelcome to part (e) of the exercise.\n""")
        lb()
        print("""Here you can do two things:\n
'Part1', use an algorithm to find the lowest order of the expansion that can calculate Pi to 7 significant figures. Given that arctan(1)=Pi/4. This will search through a range of N orders defined by the user.\n
'Part2', use another algorithm to write to a file, the N, output and difference values between a range of N orders defined by the user, about the argument x=1 for my defined function of arctan(1).""")
        lb()
        epartq='a'
        n='a'
        Nvalue='a'
        while epartq!="q":
            epartq=input('What section would you like to run? 1 or 2(or enter "q" to return to the main menu): ')
            if(epartq=='1'):#Part e, Part 1.
                while n!="q" or Nvalue!="q":
                    print("""Enter "q" to return to section (e)""")
                    n=input('Enter lowest limit: ')
                    Nvalue=input('Enter highest limit: ')
                    if isint(n) and isint(Nvalue):
                        Nvalue = int(Nvalue)
                        n=int(n)
                        tol=0.0000005#Tolerance allowed between MyArctanAll function multiplied by 4 and Pi to 7 sig fig.
                        if (abs(round(MyArctanAll(1,n)*4,6) - round((math.pi),6)))<tol: #checks to see if the limits applied contain a value of N that satisfies our search.
                            print('The first value for N satisfies our search, it is likely your lowest value is too high. Lower your lowest limit.')
                        elif(abs(round(MyArctanAll(1,Nvalue)*4,6) - round((math.pi),6)))>tol: #checks to see if the limits applied contain a value of N that satisfies our search.
                            print('There are no values for N in this range, raise your highest limit.')  
                        else:
                            Sum = MyArctanAll(1, n)
                            while(n<=Nvalue):  
                                rsum=round(Sum*4,6)
                                Diff =(rsum - round((math.pi),6))#difference between the rounded sum of the function multiplied by 4 to 7 significant figures and Pi to 7 significant figures.
                                if ((abs(Diff) < tol)):
                                    print('\nN-2 value is of order:',n-2)
                                    print('My function for N-2 multiplied by 4, gives Pi to be:', round(MyArctanAll(1,n-2)*4,6))
                                    print('\nN-1 value is of order:',n-1)
                                    print('My function for N-1 multiplied by 4, gives Pi to be:', round(MyArctanAll(1,n-1)*4,6))
                                    print('\nN value is of order: ',n)
                                    print('My function multiplied by 4, gives Pi to be:', round(Sum*4,6))
                                    print('Pi =', round((math.pi),6))
                                    n='q'#this stops the while loop from activating again when you want to exit
                                    Nvalue='q'#this stops the while loop from activating again when you want to exit
                                    break
                                else:
                                    n+=1 #required to increase the order of the expansion by 1 after each evaluation.
                                    Sum += (((-1) ** n)/((2 * n) + 1))* (1 **((2 * n) + 1))#this allows the process to be more efficient, rather than running ``MyArcanAll" for every N value.
                    elif n=='q' or Nvalue=='q':
                        print("\nYou have chosen to exit this section.")
                        epartq="q"#this stops the while loop from activating again when you want to exit
                        break
                    else:
                        print('\nThis is not within the set boundary conditions')

                        
            elif(epartq=='2'):#Part e, Part 2.
                 n='a'
                 Nvalue='a'
                 while n!="q" or Nvalue!="q":
                    print("""Enter "q" to return to section (e)""")
                    n=input('Enter lowest limit: ')
                    Nvalue=input('Enter highest limit: ')
                    if isint(n) and isint(Nvalue): #checks to make sure input is number
                        Filename = input('Name the file: ' )
                        n=int(n)
                        Nvalue=int(Nvalue)
                        tol=0.00000005
                        Y=[]#Empty array that fills with values as while loop is active
                        N=[]#Empty array that fills with values as while loop is active
                        D=[]#Empty array that fills with values as while loop is active
                        Sum = MyArctanAll(1, n)
                        while(n<=Nvalue): 
                            rsum=round(Sum*4,6)
                            Diff =(rsum - round((math.pi),6))
                            Y.append(rsum)#this fills the Y array with each output of the MyArctanAll function
                            N.append(n)#this fills the N array with each order of MyArctanAll function
                            D.append(Diff)#this fills the D array with the difference between pi and MyArctanAll function
                            n+=1 #required to increase the order of the expansion by 1 after each evaluation.
                            Sum += (((-1) ** n)/((2 * n) + 1))* (1 **((2 * n) + 1))  #this allows the process to be more efficient, rather than running ``MyArcanAll" for every N value.
                        CreateWrite(N,Y,D, Filename)
                        print("\nThe order of the data is: First, the N value data. Second, the y axis data for my defined function for arctan(x). Lastly, the difference between my function output multiplied by 4 and Pi to 7 significant figures")
                        print('\nComplete.')
                        epartq='q'#this stops the while loop from activating again when task is complete.
                        break
                    elif(epartq=="q"):
                        print("\nYou have chosen to exit this section.")
                        n='q'#this stops the while loop from activating again when you want to exit
                        Nvalue='q'#this stops the while loop from activating again when you want to exit
                        break
                    else:
                        print('\nThis is not within the set boundary conditions')
 #-------------------------------------------------------------------------------------------------------------------------------         
    elif MyInput != "q":
        print('This is not a valid choice')
        
print('\nYou have chosen to finish - goodbye.')
