#-* -  coding: utf - 8-* - 
"""
Created on Mon Feb 12 15:22:41 2018

@author: peter
"""
# ------------------------------------------------------------------------------ 
#imports and style
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style

style.use('fivethirtyeight') # Changes the style of the graphs

# ------------------------------------------------------------------------------ 
#Defined functions

def plot(x, y, title, xaxislabel, yaxislabel):
    """This function plots a graph, labels axis and title.
    """
    fig, ax = plt.subplots()
    plt.xlabel(xaxislabel) # Labels x axis.
    plt.ylabel(yaxislabel) # Labels y axis.
    plt.grid(True)
    plt.title(title)
    plt.plot(x,y)
    plt.show()

def isfloat(x):
    """Checks if object *can* be converted to a float, if it can- it returns true."""
    try: 
        float(x)
    except ValueError: 
        return False
    return True

def CreateWrite(x, y, Filename):
    """Creates a file, names it, transposes arrays and writes them to a file. 
    The order in the file is the same as the order input into the function."""
    np.savetxt(Filename, np.transpose([x,y]) , newline='\r\n', )



def FreeFallFunc(section,Y0,Yfunction, Vfunction,Ratio,PD=0):
    """This is a function that calculates position and veocity of an object in
    freefall, given its starting conditions, incremental moments in time and k/m ratio.
    This can be plotted to graphs or saved to file.
    
    'a' uses Eulers Method, 'b' uses an analytical prediciton, 'c' uses the modified Eulers Method,
    'd' uses Eulers Method and a function to evaulate airdensity as a function of height, 'd2' 
    is the same as 'd' with the addition of a parachute and 'e' calculates terminal velocity.
    """
    
    global V, Time, Yposition, YpositionMod, Vmod, Cd, A, Mass
    
    # Array creation.
    V = np.zeros(Nsteps) # Velocity array for use in Euler method.
    Time = np.zeros(Nsteps) # Time array for use in all methods.
    Yposition = np.zeros(Nsteps) # Position array for use in all methods.
    YpositionMod = np.zeros(Nsteps) # Position array for use in the modified euler method.
    Vmod = np.zeros(Nsteps) # Velocity array for use in the modified euler method.
    
    # Initial conditions.
    V[0] = 0 # Velocity array for Euler Method.
    Time[0] = 0 # Time array for all methods.
    Yposition[0] = Y0 # Position array for Euler Method.
    Vmod[0] = 0 # Velocity array for Modified Euler Method.
    YpositionMod[0] = Y0 # Position array for Modified Euler Method.


    if section  ==  "a" or section  ==  "b":
        for i in range (1,Nsteps):   
            
            # Iterative functions
            Time[i] = Time[i - 1] + dt # Increases time by dt and saves to an array for each iteration.
            Yposition[i] = Yfunction(i, Ratio) # Evaluates the method for position using the Euler method. Where "Ratio" is a scalar. Ratio used in 'b' but not 'a'.
            V[i] = Vfunction(i, Ratio, dt) # Evaluates the method for velocity using the Euler method. Where "Ratio" is a scalar.
            
            if Yposition[i]<0: # Checks to see if yposition negative? if yes go back 1 step and find time at Yposition = 0, replace values with corrected values.
                
                 # Remove all extra zeroes from end of array as the loop will end on the next iteration.
                V = np.trim_zeros(np.array(V), 'b')
                Time = np.trim_zeros(np.array(Time), 'b')
                Yposition = np.trim_zeros(np.array(Yposition), 'b')
                
                # Create the final values given change in step size that occurs due to function terminating on a time value thats not a multiple of dt.
                Time[i] = (Yposition[i - 1]/V[i - 1]) + Time[i - 1] # Calculates the time to get to Yposition  ==  zero from the (i - 1)th step, to the ith(final) step.
                Lastdt = Time[i] - Time[i - 1] # Finds the final time step size, given the termination statement.
                
                V[i] = Vfunction(i,Ratio,Lastdt) # CaLculates the final velocity, given the value of time stepsize has changed.
                Yposition[i] = 0 # Sets the final Yposition.
                break
            
            else:
                continue


        
    elif section == "c":  
        
        for i in range (1,Nsteps):     
           
            # Iterative functions
            Time[i] = Time[i - 1] + dt # Increases time by dt and saves to an array for each iteration.
            V[i] = Vfunction(i,Ratio,dt) # Evaluates the method for velocity using the Euler method. Where "Ratio" is a scalar.
            Vmod[i] = Vmod[i - 1] - dt/2*(( - g + Ratio*Vmod[i - 1]**2) + ( - g + Ratio*V[i]**2)) # Evaluates the velocity for the modified Euler method. Where "Ratio" is a scalar.
            YpositionMod[i] = YpositionMod[i - 1] - dt/2*(Vmod[i - 1] + V[i]) # Evaluates the method for position using the modified Euler Method.
     
            if YpositionMod[i]<0: # Checks to see if yposition negative? if yes go back 1 step and find time at Yposition = 0, replace values with corrected values.
                
                # Remove all extra zeroes from end of array as the loop will end on the next iteration.
                V = np.trim_zeros(np.array(V), 'b')
                Time=np.trim_zeros(np.array(Time), 'b')
                Vmod=np.trim_zeros(np.array(Vmod), 'b')
                YpositionMod=np.trim_zeros(np.array(YpositionMod), 'b')
                
                # Create the final values given change in step size that occurs due to function terminating on a time value thats not a multiple of dt.
                Time[i]=(YpositionMod[i-1]/Vmod[i-1])  +  Time[i - 1] # Calculates the time to get to Yposition  ==  zero from the (i - 1)th step, to the ith(final) step.
                Lastdt=Time[i] - Time[i - 1] # Finds the final time step size, given the termination statement.
                
                V[i] = Vfunction(i,Ratio,Lastdt) # CaLculates the final velocity using the EulerMethod, given the value of time stepsize has changed.
                Vmod[i] = Vmod[i - 1] - Lastdt/2*(( - g + Ratio*Vmod[i - 1]**2) + ( - g + Ratio*V[i]**2)) # CaLculates the final velocity using the ModifiedEulerMethod, given the value of time stepsize has changed.
                YpositionMod[i] = 0 # Sets the final Yposition.
                break
            
            else:
                continue
            
    elif section  ==  "d2" or section  ==  "d":
        
        for i in range (1,Nsteps):    
            
            # Iterative functions
            Time[i] = Time[i - 1] + dt # Increases time by dt and saves to an array for each iteration.
            Yposition[i] = Yfunction(i,Ratio(i,A,Cd,Mass)) # Evaluates the method for position using the Modified Euler method. Where "Ratio" is a function.
            V[i] = Vfunction(i,Ratio(i,A,Cd,Mass),dt) # Evaluates the method for velocity using the Modified Euler method. Where "Ratio" is a function.
            
            if Yposition[i]<0: # Checks if yposition is negative. If yes go back 1 step and find time at Yposition =0, replace values.
                
                # Remove all extra zeroes from end of array as the loop will end on the next iteration.
                V = np.trim_zeros(np.array(V), 'b')
                Time = np.trim_zeros(np.array(Time), 'b') 
                Yposition = np.trim_zeros(np.array(Yposition), 'b')
                
                # Create the final values given change in step size that occurs due to function terminating on a time value thats not a multiple of dt.
                Time[i] = (Yposition[i - 1]/V[i - 1]) + Time[i - 1] # Calculates the time to get to Yposition  ==  zero from the (i - 1)th step, to the ith(final) step.
                Lastdt = Time[i] - Time[i - 1] # Finds the final time step size, given the termination statement.
                
                V[i] = Vfunction(i,Ratio(i,A,Cd,Mass), Lastdt) # CaLculates the final velocity using the EulerMethod, given the value of time stepsize has changed.
                Yposition[i] = 0 # Sets the final Yposition.
                break
            
            elif Yposition[i]<=PD and section == "d2":# Checks to see if yposition is below minimum parachute height and deploys parachute, only if its in the correct section.
                Cd=1.75 # Change the drag coefficient to that of a parachute.
                A=25 # Change the cross sectional area to that of a parachute.
                
            else:
                continue           

#part A,C,D functions
def YfunctionACD(i, Ratio):
    """This is  equation 6 from the exercise, it is the Euler Method
    for determining the position of an object in free fall as function
    of velocity and time. This is an iterative method.
    """
    return Yposition[i - 1] - dt * V[i - 1]

def VfunctionACD(i,Ratio,dt):
    """This is equation 5 from the exerice, is the Euler Method for 
    determining velocity of an object in freefall as a function of time. 
    This is an iterative method.
    """
    return V[i - 1] - ( - g + Ratio*V[i - 1]**2)*dt

def KMFunction(i,A,Cd,Mass):
    """This is the equation for air density which is a function of height.
    """
    return (Cd*(FluidDensity*math.exp( - (Yposition[i])/h))*A)/(2*Mass)



#-----------------------------------------------------------------------------
# Part B functions.
def YfunctionB(i,Ratio):
    """This is equation 8 from the exercise, it is an analytical prediction
    of position for an object in freefall as a function of time.
    This is an iterative method.
    """
    return (Yposition[0]- ((Ratio**(-1))/(2))*np.log(np.cosh((g*Ratio)**(1/2)*Time[i])**2))
 
def VfunctionB(i,Ratio,dt):
    """This is equation 9 from the exercise, it is an analytical prediction
    of velocity for an object in freefall as a function of time.
    This is an iterative method.
    """
    return (g*Ratio**(-1))**(1/2)*np.tanh((Ratio*g)**(1/2)*Time[i])

    
#-----------------------------------------------------------------------------
# Part E function, this is an addition and does not represent section 'e' in the exercise.

    
def TerminalV(Maxheight,dy):
    """Calculates the terminal velocity as function of height, given a fluid density
    that is also a function of height. This is an iterative method.
    """
    global Height,Vterm
    
    Hsteps=int(Maxheight/dy) # Number of steps.
    # Array creation.
    Height=np.zeros(Hsteps)
    Vterm=np.zeros(Hsteps)
    
    Height[0]=0

    for i in range(1,Hsteps):
        Height[i]=Height[i - 1] + dy
        Vterm[i]= ((2*Mass*(g))/(A*Cd*(FluidDensity*math.exp( - (Height[i])/h))))**(1/2)
    
    return


def lb():
    """Creates a dashed line segment for clearer section seperation and easier reading"""
    print('------------------------------------------------------------------------------------------------------------------------------------')

#-----------------------------------------------------------------------------
#Constants and initial conditions

g=9.81 # Acceleration due to gravity in m/s^2.
h=7640 # Scale height of the atmosphere in m.

#-----------------------------------------------------------------------------
lb()

print("""Welcome to exercise 3 of the Level 5 Laboratory: Computational Physics tasks set for second year physics students at the University of Bristol in 2018.\n""")

print("""This python programme is designed to allow the user to simulate and predict the velocity and position of a mass in freefall through a fluid using Eulers Method and its modified counterpart.\n""")

print("""Graphs can be plotted for all simulations and prediction with the option to save a text file in order to further analyse the data. \n""")

print("""Users will be able to change the initial height, mass, density of the fluid and object, aswell as the the incremental time of evaluation. \n""")


lb()


MyInput = 0
while (MyInput != "q"):
    
    # Brief description of each section.
    print("""\nPart (a) of the exercise allows the user to simulate the freefall of a mass using Eulers method, through a fluid of constant density. \n""")
    
    print("""Part (b) of the exercise allows the user to use analytical prediction for the freefall of a mass, through a fluid of constant density. \n""")
    
    print("""Part (c) of the exercise allows the user to simulate the freefall of a mass using the Modified Eulers method, through a fluid of constant density. \n""")
    
    print("""Part (d) of the exercise allows the user to simulate the freefall of a mass using Eulers method, through air, of density which is a function of height. \n""")
    
    print("""Part (d2) of the exercise allows the user to simulate the freefall of a mass using Eulers method, through air, of density which is a function of height. Where a parachute will open at a user defined height. \n""")
    
    print("""Part (e) of the exercise allows the user to simulate the terminal velocity of a sphere as a function of height. \n""")

    lb()
    
    # Allow user to choose section to navigate to.
    
    MyInput = input('Main menu: Enter a section, "a", "b", "c", "d", "d2", "e". Or "q" to quit: ')

    if MyInput == "a":
        
        lb()
        
        # Here each value that the user can define is set to "P" to allow the while function to keep looping unless they wish to quit.
        Y0 = "P" # Initial height of freefall.
        Cd = "P" # Drag coeficient of object in freefall.
        Mass = "P" # Mass of the object in freefall.
        FluidDensity = "P" # Density of the fluid for which the object falls through.
        BodyDensity = "P" # Density of the object.
        dt = "P" # Increment for time step size.
        
        print("""\n This is part 'a' of the exercise.\n
              This simulates an object in freefall through a non - turbulent fluid, defined by its drag coefficient. The simulation is computed using the Euler's method. Enter 'q' into any user input to exit this section.""")
        
        while Cd !="q" and Y0 != "q" and Mass != "q" and FluidDensity != "q" and BodyDensity != "q" and dt != "q": # Keeps user in the section, until they enter "q" into any/all the inputs below, to allow them to quit.
            
            #Input values for the function FreeFallFunc.
            Y0 = input("Please enter a starting height for the object to begin its freefall, in metres: ")
            Cd = input("Please enter a drag coefficient for the mass. 0.47 can be used for a sphere, or if you wish to simulate the drag of a person, a value of '1.3' is appropriate: ")
            Mass = input("Please enter a mass in kg, for the object in freefall: ")
            FluidDensity = input("Please enter a value for the density of the fluid for which the object is in freefall, a value of '1.2' is suitable for 'air': ")
            BodyDensity = input("Please enter a value for the density of the object in free fall, a value of '1000' kg/m^3 is suitable for a human: ")
            dt = input("Please enter a time increment for frequency of data evaluation. For example a time of 0.01s is adequate for most purposes: ")
            #-----------------------------------------------------------------------------
            
            if isfloat(Y0) and isfloat(Cd)and isfloat(Mass)  and isfloat(FluidDensity) and isfloat(dt) and isfloat(BodyDensity): # Checks to see if inputs can be converted to a float, to avoid incorrect object types creating errors.
                
                # Convert object to float.
                Y0 = float(Y0)
                Cd = float(Cd) 
                Mass = float(Mass)
                FluidDensity = float(FluidDensity)
                dt = float(dt)
                BodyDensity = float(BodyDensity)
                #-----------------------------------------------------------------------------
                MaxTime= 10000 # Given an arbitrarily large value, such that the maxtime aloted is unlikely to be reached before the function selfterminates. But not so large it takes a long time to create.
                Nsteps =int(MaxTime/dt) # The number of data points per axis.
                #-----------------------------------------------------------------------------
                # Properties of object and constants.
                Volume = Mass/BodyDensity # Volume of object.
                A = math.pi*((3/(4*math.pi)*Volume)**(2/3)) # Cross sectional area of object.
                k = (Cd*FluidDensity*A)/2 # k constant.
                Ratio = k/Mass
                #-----------------------------------------------------------------------------
                # Runs the function and plots.
                FreeFallFunc("a", Y0, YfunctionACD, VfunctionACD, Ratio) 
                plot(Time, V, "Velocity vs Time Part: " + "a" + ", " +  str(dt)  +  "s" +  " "  +  ", km^-1=" + str(round(Ratio,6)),"Time(s)", "Velocity(m/s)")
                plot(Time, Yposition, "Y Position vs Time Part: "  + "a"  + ", "  +  str(dt)  +  "s" + ", km^-1=" + str(round(Ratio,6)), "Time(s)", "YPosition(m)")
                #-----------------------------------------------------------------------------
                print("The maximum velocity during the jump is: ", np.amax(V), "m/s")
                print(" The flight time of the object is: ", str(Time[ - 1])  + "s")
                #-----------------------------------------------------------------------------
                # Asks user if they wish to save the arrays to a text file. 
                SaveQ = input("Would you like to save this data to a txt file? y/n ")
                if SaveQ  ==  'y':
                    CreateWrite(Time,Yposition,"Y vs T Part " + "a" + " " +  str(dt)  +  "s" + " km^-1=" + str(round(Ratio,6)) + ".txt")
                    CreateWrite(Time,V,"V vs T Part " + "a"  + " " +  str(dt)  +  "s" + " km^-1=" + str(round(Ratio,6))  +  ".txt")
                    break
                
                else:
                    break
                
            elif Y0 == "q" or Cd == "q" or Mass == "q" or FluidDensity == "q" or BodyDensity == "q" or dt == "q":
                print("\nYou have chosen to exit this section.")
                
            else:
                print('\nThis is not a valid value.') 

 #-----------------------------------------------------------------------------
    elif MyInput == "b":
        
        lb()
        # Here each value that the user can define is set to "P" to allow the while function to keep looping unless they wish to quit.
        Y0="P" #Initial height of freefall
        Cd="P"# Drag coeficient of object in freefall
        Mass="P" #Mass of the object in freefall
        FluidDensity="P"#Density of the fluid for which the object falls through.
        BodyDensity ="P"#Density of the object
        dt = "P" # increment for frequency of evaluation
        
        print("""\n This is part 'b' of the exercise.\n
              This runs an analytical prediction for an object in freefall in a uniform gravitational field, in a non - turbulent fluid, defined by its drag coefficient. Enter 'q' into any user input to exit this section.""")
        
        while Cd !="q" and Y0 != "q" and Mass != "q" and FluidDensity != "q" and dt != "q" and BodyDensity != "q": #Keeps user in the program, until they enter "q" into any/all the inputs below.
            
            #Input values for the function FreeFallFunc.
            Y0 = input("Please enter a starting height for the sphere to begin its freefall, in metres: ")
            Cd = input("Please enter a drag coefficient for the mass. 0.47 can be used for a sphere, or if you wish to simulate the drag of a person a value of '1.3' is appropriate: ")
            Mass = input("Please enter a mass in kg, for the object in freefall: ")
            FluidDensity = input("Please enter a value for the density of the fluid for which the object is in freefall, a value of '1.2' is suitable for 'air': ")
            BodyDensity = input("Please enter a value for the density of the object in free fall, a value of '1000' kg/m^3 is suitable for a human: ")
            dt = input("Please enter a time increment for frequency of data evaluation. For example a time of 0.01s is adequate for most purposes: ")
            #-----------------------------------------------------------------------------
            
            if isfloat(Y0) and isfloat(Cd)and isfloat(Mass)  and isfloat(FluidDensity) and isfloat(dt) and isfloat(BodyDensity): #checks to see if inputs can be converted to a float.
                # Convert object to float.
                Y0 = float(Y0)
                Cd = float(Cd) 
                Mass = float(Mass)
                FluidDensity = float(FluidDensity)
                dt = float(dt)
                BodyDensity=float(BodyDensity)
                #-----------------------------------------------------------------------------
                MaxTime= 10000 # Given an arbitrarily large value, such that the maxtime aloted is unlikely to be reached before the function selfterminates.
                Nsteps =int(MaxTime/dt)
                #-----------------------------------------------------------------------------
                # Properties of object and constants.
                Volume=Mass/BodyDensity # Volume of object.
                A=math.pi*((3/(4*math.pi)*Volume)**(2/3)) # Cross sectional area of Volume/object.
                k=(Cd*FluidDensity*A)/2 # k constant.
                Ratio=k/Mass
                #-----------------------------------------------------------------------------
                # Runs the function and plots.
                FreeFallFunc("b",Y0, YfunctionB, VfunctionB,Ratio) 
                plot(Time, V, "Velocity vs Time Part: " + "b" + ", " +  str(dt)  +  "s" +  " "  +  ", km^-1=" + str(round(Ratio,6)),"Time(s)", "Velocity(m/s)")
                plot(Time, Yposition, "Y Position vs Time Part: "  + "b"  + ", "  +  str(dt)  +  "s" + ", km^-1=" + str(round(Ratio,6)), "Time(s)", "YPosition(m)")
                #-----------------------------------------------------------------------------
                print("The maximum velocity during the analytical model is: ", np.amax(V), "m/s")
                print(" The flight time of the object is: ", str(Time[ - 1])  + "s")
                #-----------------------------------------------------------------------------
                # Asks user if they wish to save the arrays to a text file. 
                SaveQ = input("Would you like to save this data to a txt file? y/n ")
                if SaveQ  ==  'y':
                    CreateWrite(Time,Yposition,"Y vs T Part " + "b" + " " +  str(dt)  +  "s" + " km^-1=" + str(round(Ratio,6)) + ".txt")
                    CreateWrite(Time,V,"V vs T Part " + "b"  + " " +  str(dt)  +  "s" + " km^-1=" + str(round(Ratio,6))  +  ".txt")
                    break
                else:
                    break
            elif Y0 == "q" or Cd == "q" or Ratio == "q" or FluidDensity == "q" or BodyDensity == "q" or dt == "q":
                print("\nYou have chosen to exit this section.")
            else:
                print('\nThis is not a valid value.') 
                
 #-----------------------------------------------------------------------------
    elif MyInput == "c":
        
        lb()
        # Here each value that the user can define is set to "P" to allow the while function to keep looping unless they wish to quit.
        Y0="P" #Initial height of freefall
        Cd="P"# Drag coeficient of object in freefall
        Mass="P" #Mass of the object in freefall
        FluidDensity="P"#Density of the fluid for which the object falls through.
        BodyDensity ="P"#Density of the object
        dt = "P" # increment for frequency of evaluation
        
        print("""\n This is part 'c' of the exercise.\n
              This simulates an object in freefall through a non - turbulent fluid, defined by its drag coefficient. The simulation is computed using the Modified Euler's method. Enter 'q' into any user input to exit this section.""")
        
        while Cd !="q" and Y0 != "q" and Mass != "q" and FluidDensity != "q" and dt != "q" and BodyDensity != "q": #Keeps user in the program, until they enter "q" into any/all the inputs below.
            
            #Input values for the function FreeFallFunc.
            Y0 = input("Please enter a starting height for the object to begin its freefall, in metres: ")
            Cd = input("Please enter a drag coefficient for the mass. 0.47 can be used for a sphere, or if you wish to simulate the drag of a person a value of '1.3' is appropriate: ")
            Mass = input("Please enter a mass in kg, for the object in freefall: ")
            FluidDensity = input("Please enter a value for the density of the fluid for which the object is in freefall, a value of '1.2' is suitable for 'air': ")
            BodyDensity = input("Please enter a value for the density of the object in free fall, a value of '1000' kg/m^3 is suitable for a human: ")
            dt = input("Please enter a time increment for frequency of data evaluation. For example a time of 0.01s is adequate for most purposes: ")
            #-----------------------------------------------------------------------------
            
            if isfloat(Y0) and isfloat(Cd)and isfloat(Mass)  and isfloat(FluidDensity) and isfloat(dt) and isfloat(BodyDensity): #checks to see if inputs can be converted to a float.
                
                # Convert object to float.
                Y0 = float(Y0)
                Cd = float(Cd) 
                Mass = float(Mass)
                FluidDensity = float(FluidDensity)
                dt = float(dt)
                BodyDensity=float(BodyDensity)
                #-----------------------------------------------------------------------------
                #Time and step sizes
                MaxTime= 10000 # Given an arbitrarily large value, such that the maxtime aloted is unlikely to be reached before the function selfterminates.
                Nsteps =int(MaxTime/dt) # The number of data points per axis.
                #-----------------------------------------------------------------------------
                #Properties of object and constants
                Volume=Mass/BodyDensity # Volume of object
                A=math.pi*((3/(4*math.pi)*Volume)**(2/3)) #Cross sectional area of Volume/object
                k=(Cd*FluidDensity*A)/2 # k constant
                Ratio=k/Mass
                #-----------------------------------------------------------------------------
                # Runs the function and plots.
                FreeFallFunc("c",Y0, YfunctionACD, VfunctionACD,Ratio) 
                plot(Time, Vmod, "Velocity vs Time Part: " + "c" + ", " +  str(dt)  +  "s" +  " "  +  ", km^-1=" + str(round(Ratio,6)),"Time(s)", "Velocity(m/s)")
                plot(Time, YpositionMod, "Y Position vs Time Part: "  + "c"  + ", "  +  str(dt)  +  "s" + ", km^-1=" + str(round(Ratio,6)), "Time(s)", "YPosition(m)")
                #-----------------------------------------------------------------------------
                print("The maximum velocity during the jump is: ", np.amax(Vmod), "m/s")
                print(" The flight time of the object is: ", str(Time[ - 1])  + "s")
                #-----------------------------------------------------------------------------
                # Asks user if they wish to save the arrays to a text file. 
                SaveQ = input("Would you like to save this data to a txt file? y/n ")
                if SaveQ  ==  'y':
                    CreateWrite(Time,YpositionMod,"Y vs T Part " + "c" + " " +  str(dt)  +  "s" + " km^-1=" + str(round(Ratio,6)) + ".txt")
                    CreateWrite(Time,Vmod,"V vs T Part " + "c"  + " " +  str(dt)  +  "s" + " km^-1=" + str(round(Ratio,6))  +  ".txt")
                    break
                
                else:
                    break
                
            elif Y0 == "q" or Cd == "q" or Mass == "q" or FluidDensity == "q" or BodyDensity == "q" or dt == "q":
                print("\nYou have chosen to exit this section.")
                
            else:
                print('\nThis is not a valid value.')  
                
 #-----------------------------------------------------------------------------
    elif MyInput == "d":
        
        lb()
        # Here each value that the user can define is set to "P" to allow the while function to keep looping unless they wish to quit.
        Y0="P" #Initial height of freefall
        Cd="P"# Drag coeficient of object in freefall
        Mass="P" #Mass of the object in freefall
        BodyDensity ="P"#Density of the object
        dt = "P" # Frequency of evaulation.
        
        print("""\n This is part 'd' of the exercise.\n
              This simulates an object in freefall through non - turbulent air. Where the air density is a function of height. The simulation is computed using the Euler's method. Enter 'q' into any user input to exit this section.""")
        
        while Cd !="q" and Y0 != "q" and Mass != "q" and dt != "q" and BodyDensity != "q": #Keeps user in the program, until they enter "q" into any/all the inputs below.
            
            #Input values for the function FreeFallFunc.
            Y0 = input("Please enter a starting height for the sphere to begin its freefall, in metres. Felix's height was 39045m.: ")
            Cd = input("Please enter a drag coefficient for the mass. 0.47 can be used for a sphere, or if you wish to simulate the drag of a person a value of '1.3' is appropriate: ")
            Mass = input("Please enter a mass in kg, for the object in freefall. Felix's Mass was approximately 73kg: ")
            BodyDensity = input("Please enter a value for the density of the object in free fall, a value of '1000' kg/m^3 is suitable for a human: ")
            dt = input("Please enter a time increment for frequency of data evaluation. For example a time of 0.01s is adequate for most purposes: ")
            #-----------------------------------------------------------------------------
            
            if isfloat(Y0) and isfloat(Cd) and isfloat(dt) and isfloat(Mass) and isfloat(BodyDensity): #checks to see if inputs can be converted to a float.
                
                #Convert object to float. 
                Y0 = float(Y0)
                Cd = float(Cd) 
                Mass = float(Mass)
                dt = float(dt)
                BodyDensity=float(BodyDensity)
                #-----------------------------------------------------------------------------
                # Time and step sizes.
                MaxTime= 10000 # Given an arbitrarily large value, such that the maxtime aloted is unlikely to be reached before the function selfterminates.
                Nsteps =int(MaxTime/dt) # The number of data points per axis.
                #-----------------------------------------------------------------------------
                # Properties of object and constants.
                FluidDensity = 1.2
                FluidDensity = float(FluidDensity)
                Volume=Mass/BodyDensity # Volume of object.
                A=math.pi*((3/(4*math.pi)*Volume)**(2/3)) # Cross sectional area of Volume/object.
                #-----------------------------------------------------------------------------
                # Runs the function and plots.
                FreeFallFunc("d",Y0, YfunctionACD, VfunctionACD,KMFunction) 
                plot(Time, V, "Velocity vs Time Part: " + "d" + ", " +  str(dt)  +  "s" +  " " + " Alpha " + str(round(Cd*A/Mass,6)) ,"Time(s)", "Velocity(m/s)")
                plot(Time, Yposition, "Y Position vs Time Part: "  + "d"  + ", "  +  str(dt)  +  "s" +  " " + " Alpha " + str(round(Cd*A/Mass,6)), "Time(s)", "YPosition(m)")
                #-----------------------------------------------------------------------------
                print("The maximum velocity during the jump is: ", np.amax(V), "m/s")
                print(" The flight time of the object is: ", str(Time[ - 1])  + "s")
                #-----------------------------------------------------------------------------
                # Asks user if they wish to save the arrays to a text file. 
                SaveQ = input("Would you like to save this data to a txt file? y/n ")
                if SaveQ  ==  'y':
                    CreateWrite(Time,Yposition,"Y vs T Part " + "d" + " " +  str(dt)  +  "s, Height " + str(Y0)  + " Alpha " + str(round(Cd*A/Mass,6))  +  ".txt")
                    CreateWrite(Time,V,"V vs T Part " + "d"  + " " +  str(dt)  +  "s, Height " + str(Y0)  + " Alpha " + str(round(Cd*A/Mass,6))  +  ".txt")
                    break
                
                else:
                    break
                
            elif Y0 == "q" or Cd == "q" or Mass == "q" or BodyDensity == "q" or dt == "q":
                print("\nYou have chosen to exit this section.")
                
            else:
                print('\nThis is not a valid value.')                 
 #-----------------------------------------------------------------------------
    elif MyInput == "d2":
        
        lb()
        # Following values are set to "P" to keep the user in a loop until they wish to quit.
        Y0="P" # Initial height of freefall.
        Cd="P" # Drag coeficient of object in freefall.
        Mass="P" # Mass of the object in freefall.
        BodyDensity ="P" # Density of the object.
        PD = "P" # Height above ground where the parachute deploys.
        FluidDensity = 1.2 # Density of air.
        dt = "P" # Size of time step for eache evaulation.
        #-----------------------------------------------------------------------------
        
        print("""\n This is part 'd2' of the exercise.\n
              This simulates an object in freefall through non - turbulent air. Where the air density is a function of height. The simulation is computed using the Euler's method. Enter 'q' into any user input to exit this section.""")
        print("This simulation also has a parachute function, where by at a set height parachute will open at 1500m (the height at which Felix opened his 'chute).")
        while Cd !="q" and Y0 != "q" and Mass != "q" and dt != "q" and PD != "q" and BodyDensity != "q": #Keeps user in the program, until they enter "q" into any/all the inputs below.
            #Input values for the function FreeFallFunc.
            Y0 = input("Please enter a starting height for the sphere to begin its freefall, in metres. Felix's height was 39045m.: ")
            Cd = input("Please enter a drag coefficient for the mass. 0.47 can be used for a sphere, or if you wish to simulate the drag of a person a value of '1.3' is appropriate: ")
            Mass = input("Please enter a mass in kg, for the object in freefall. Felix's Mass was approximately 73kg: ")
            BodyDensity = input("Please enter a value for the density of the object in free fall, a value of '1000' kg/m^3 is suitable for a human: ")
            PD = input("Please enter a height for when Felix's parachute deploys, in his jump it deploys at 1524m: ")
            dt = input("Please enter a time increment for frequency of data evaluation. For example a time of 0.01s is adequate for most purposes: ")
            #-----------------------------------------------------------------------------
            
            if isfloat(Y0) and isfloat(Cd) and isfloat(dt) and isfloat(Mass) and isfloat(PD) and isfloat(BodyDensity): #checks to see if inputs can be converted to a float.
                     
                #Convert object to float.   
                Y0 = float(Y0)
                Cd = float(Cd) 
                Mass = float(Mass)
                dt = float(dt)
                PD =float(PD)
                BodyDensity = float(BodyDensity)
                #-----------------------------------------------------------------------------
                #Time and step sizes
                MaxTime= 10000 # Given an arbitrarily large value, such that the maxtime aloted is unlikely to be reached before the function selfterminates.
                Nsteps =int(MaxTime/dt) # The number of data points per axis.
                #-----------------------------------------------------------------------------
                #Properties of object and constants
                FluidDensity = 1.2 # Density of air.
                FluidDensity = float(FluidDensity)
                Volume=Mass/BodyDensity # Volume of object.
                A=math.pi*((3/(4*math.pi)*Volume)**(2/3)) # Cross sectional area of Volume/object.        
                #-----------------------------------------------------------------------------
                # Runs the function and plots.
                FreeFallFunc("d2",Y0, YfunctionACD, VfunctionACD,KMFunction,PD) 
                plot(Time, V, "Velocity vs Time Part: " + "d2" + ", " +  str(dt)  +  "s" +  " " + " Alpha " + str(round(Cd*A/Mass,6)) ,"Time(s)", "Velocity(m/s)")
                plot(Time, Yposition, "Y Position vs Time Part: "  + "d2"  + ", "  +  str(dt)  +  "s" +  " " + " Alpha " + str(round(Cd*A/Mass,6)), "Time(s)", "YPosition(m)")
                
                print("The maximum velocity during the jump is: ", np.amax(V), "m/s")
                print(" The flight time of the object is: ", str(Time[ - 1])  + "s")
                # Asks user if they wish to save the arrays to a text file.
                SaveQ = input("Would you like to save this data to a txt file? y/n ")
                if SaveQ  ==  'y':
                    CreateWrite(Time,Yposition, "Y vs T Part " + "d2" + " " +  str(dt)  +  "s, Height " + str(Y0)  + " Alpha " + str(round(Cd*A/Mass,6))  +  ".txt")
                    CreateWrite(Time, V, "V vs T Part " + "d2"  + " " +  str(dt)  +  "s, Height " + str(Y0)  + " Alpha " + str(round(Cd*A/Mass,6))  +  ".txt")
                    break
                else:
                    break
            elif Y0 == "q" or Cd == "q" or Mass == "q" or BodyDensity == "q" or dt == "q" or PD == "q":
                print("\nYou have chosen to exit this section.")
            else:
                print('\nThis is not a valid value.') 
                
 #-----------------------------------------------------------------------------
    elif MyInput == "e":
        lb()
        # Following values are set to "P" to keep the user in a loop until they wish to quit.
        Maxheight = "P" # Maxheight to test terminal velocity from 0 to Maxheight.
        Cd = "P" # Drag coeficient of object in freefall.
        Mass = "P" # Mass of the object in freefall.
        BodyDensity = "P" # Density of the object.
        dy = "P" # Size of distance step for eache evaulation.
        #-----------------------------------------------------------------------------
        print("\n This section allows you to create a graph for Terminal velocity of an object as a function of height and other characteristics for air.")
        
        while Cd !="q" and Maxheight != "q" and Mass != "q" and BodyDensity != "q" and dy!= "q":
            Maxheight = input("Please enter a max height- in metres- To calculate terminal velocity for: ")
            Cd = input("Please enter a drag coefficient for the mass. 0.47 can be used for a sphere, or if you wish to simulate the drag of a person a value of '1.3' is appropriate: ")
            Mass = input("Please enter a mass in kg, for the object in freefall: ")
            BodyDensity = input("Please enter a value for the density of the object in free fall, a value of '1000' kg/m^3 is suitable for a human: ")
            dy = input("Please enter a distance increment for frequency of data evaluation. For example a time of 5m is adequate for most purposes: ")
            #-----------------------------------------------------------------------------
            if isfloat(Maxheight) and isfloat(Mass) and isfloat(Cd) and isfloat(BodyDensity) and isfloat(dy):
                Maxheight =float(Maxheight)
                Mass = float(Mass)
                Cd = float(Cd)
                BodyDensity = float(BodyDensity)
                FluidDensity = 1.2
                FluidDensity = float(FluidDensity)
                dy=float(dy)
                #-----------------------------------------------------------------------------
                #Properties of object and constants
                Volume=Mass/BodyDensity # Volume of object
                A=math.pi*((3/(4*math.pi)*Volume)**(2/3)) #Cross sectional area of Volume/object
                
                #-----------------------------------------------------------------------------
                TerminalV(Maxheight,dy)
                plot(Height, Vterm, "Terminal Velocity vs Position Part: " + "d" + ", " +  str(dy)  +  "m" +  " " ,"Position(m)", "Velocity(m/s)")
                
                SaveQ = input("Would you like to save this data to a txt file? y/n ")
                if SaveQ  ==  'y':
                    CreateWrite(Height,Vterm,"TermV vs position Part " + "e"  + " " +  str(dy)  +  "m, Height " + str(Maxheight)  + " Alpha " + str(round(Cd*A/Mass,6))  +  ".txt")
                    break
                
                else:
                    break
                
            elif Maxheight == "q" or Cd == "q" or Mass == "q" or BodyDensity == "q" or dy == "q":
                print("\nYou have chosen to exit this section.")
                
            else:
                print('\nThis is not a valid value.') 

    elif MyInput == 'q':
        print("\nYou have chosen to exit this program.")
        
    else:
        print('\nThis is not a valid value.') 
































