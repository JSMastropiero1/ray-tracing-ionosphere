import numpy as np

def elec_dens_linear_day(r):
    '''Linear model for electron density distribution in Ionosphere'''
    if (r < 301):
        a = (4e11 - 6e8)/(300 - 50)
        b = 6e8 - a*50
        return a*r + b
    
    else:
        a = (9e10-4e11)/(1000-300)
        b = 9e10 - a*1000
        return a*r + b

####################################

def elec_dens_linear_night(r):
    '''Linear model for electron density distribution in Ionosphere at night'''
    if (r < 501):
        a = (8.5e10 - 6e8)/(500 - 70)
        b = 6e8 - a*70
        return a*r + b
    
    else:
        a = (5e10-8.5e10)/(1000-500)
        b = 8.5e10 - a*500
        return a*r + b

####################################
    
def elec_dens_quadratic_day(r):
    '''Second order polynomial function approximation for Ionesphere's electro density'''
    if (r < 301):
        x = 50.
        y = 6.e8
        xmax = 300 # Maximum electron density at 300 km
        ymax = 4.e11 # Maximum density during day
        # having a and c in func of b, using an extra point we can get b. The expression is
        b = (y - ymax)/(-1*x*x/(2.*xmax) + x - xmax/2)  
        c = ymax - b*xmax/2 # for condition given to a b has be this
        a = -1.*b/(2.*xmax) # condition to get the max in this point
        return a*r*r + b*r + c
    else :
        x = 1000.
        y = 9.e10
        xmax = 300 # Maximum electron density at 300 km
        ymax = 4.e11 # Maximum density during day
        # having a and c in func of b, using an extra point we can get b. The expression is
        b = (y - ymax)/(-1*x*x/(2.*xmax) + x - xmax/2)  
        c = ymax - b*xmax/2 # for condition given to a b has be this
        a = -1.*b/(2.*xmax) # condition to get the max in this point
        return a*r*r + b*r + c

####################################
    
def elec_dens_quadratic_night(r):
    '''Second order polynomial function approximation for Ionesphere's 
    electro density at night'''
    if (r < 501 ):
        x = 70.
        y = 6.e8
        xmax = 500 # Maximum electron density at 500 km
        ymax = 8.5e10 # Maximum density during night
        b = (y - ymax)/(-1*x*x/(2.*xmax) + x - xmax/2)
        c = ymax - b*xmax/2
        a = -1.*b/(2.*xmax)
        return a*r*r + b*r + c
    else :
        x = 1000.
        y = 5.e10
        xmax = 500 # Maximum electron density at 300 km
        ymax = 8.5e10 # Maximum density during day
        # having a and c in func of b, using an extra point we can get b. The expression is
        b = (y - ymax)/(-1*x*x/(2.*xmax) + x - xmax/2)  
        c = ymax - b*xmax/2 # for condition given to a b has be this
        a = -1.*b/(2.*xmax) # condition to get the max in this point
        return a*r*r + b*r + c

###############################

def ref_index(freq,elec_dens,parameter = [-1,-1],delta = 0):
    total_y = 1200
    total = 5000
    N = np.zeros(total) # Refraction index

    # First 80 km the refraction index will be 1 (air ref. index approx 1)
    
    for i in range(0,50,1):
        N[i] = 1

    vacuum_permittivity = 8.854e-12 # Permittivity of free space in Faradays over meters
    e_charge = -1.602e-19 # Electron charge in Coulombs
    e_mass = 9.10938356e-31 # Electron mass in Kilograms
    
# The electron density has to be given in units of meters^(-3)

    for i in range(50,total_y,1):
        if (i>=parameter[0] and i<=parameter[1]):
            pl_freq = np.sqrt(((elec_dens(i)+delta)*e_charge*e_charge)/(e_mass*vacuum_permittivity))
        else:
            pl_freq = np.sqrt((elec_dens(i)*e_charge*e_charge)/(e_mass*vacuum_permittivity))
        
        if (pl_freq >= freq):
            N[i] = 0
        else:
            N[i] = np.sqrt(1 - (pl_freq/freq)**2)

    return N


###############################


def GetRay(angle,n):
    c = 1
    #We choose the initial position and angle of the ray:
    ray=[[0,6370]]
    #We just do iterations:
    step_t=0.5 #constant, time passing between iterations
    current_index=0 #To see in which pixel we are
    sum_h=0 #The sum of the height, to see if it passes to the next pixel
    delta_y = 1
    #We track te iterations number too:
    i=0

    #A boolean to stop the iterations:
    boolean=True

    while boolean: 
        #The n index depends on the height:
        n0=n[current_index]
        x=step_t*c/n0*np.cos(angle)
        y=step_t*c/n0*np.sin(angle)
        # print("El pinche y:", y)
        #We need to see if we have reached the next height with this last step
        #To do so, we check if the sum is bigger than the grid size:
        '''
        print("sum of height: ",sum_h)
        print("refraction index: ",n0)
        print("angle: ", angle)
        print("current index: ", current_index)
        '''
        # print("PREVIOUS sum of height: ",sum_h0)
        if (abs(sum_h+y)>=delta_y):
            k=1
            while((abs(sum_h+y)-k*delta_y)>delta_y):
                k+=1
            n1=np.mean(n[current_index+1:current_index+k+1])
            # print("n1: ",n1)
            current_index+=k #We go to the next grid
            #See if changes refraction index:
            if (n1==n0):
                # print("n1 is equal to n0")
                sum_h+= y - delta_y #Goes to next cell normal
                ray.append([ray[i-1][0]+x,ray[i-1][1]+y])
            else:
                #Snell law to calculate new angle.
                #CAUTION: Snell angle is taken from the normal and ours is from the horizontal axis:
                #We need to do sin(90-theta)=cos(theta)
                angle0=angle
                #CAREFUL: arcsin input has to be between -1 and 1. If it exceedes, add/substract 1
                if (((n0*np.cos(angle0))/n1)>=1 or n1 == 0):
                    #There has been an internal reflection, the wave can't pass. We stop
                    boolean=False
                else:
                    # print("x expected: ",x)
                    # print("y expected: ",y)
                    # print("sum of y: ",sum_h)
                    # print("n0: ",n0)
                    # print("angle: ",angle0)
                    angle=(np.pi/2.)-np.arcsin((n0*np.cos(angle0))/n1)
                    #Now we need to calculate the height.
                    #It divides into two heights:
                    #1)height in previous cell: delta_y-sum_h
                    #2)height in new cell with new angle and n:
                    y1=(delta_y-sum_h)
                    #vprint("y1: ",y1)
                   # print("THE FUCKING Y1: ",y1)
                    y2=(c/n1)*np.sin(angle)*(step_t-((delta_y-sum_h)/(c*np.sin(angle0)/n0)))
                    y=y1+y2
                    #vprint("y2: ",y2)
                    #The y can be negative, but the x is always positive (it moves forward)
                    x=abs(y1/np.tan(angle0))+abs(y2/np.tan(angle))
                    #print("x real: ",x)
                    #print("y real: ",y)
                    #print("n1: ",n1)
                    #print("angle: ",angle)
                    #We put the sum back to current height in next cell, to count again:
                    sum_h=y2 
                    #We add the new total position of our ray:
                    ray.append([ray[i-1][0]+x,ray[i-1][1]+y])

        else:
            #If it doesn't jump to next cell the ray will travel normally:
            ray.append([ray[i-1][0]+x,ray[i-1][1]+y])
            # print("El pinche y:", y)
            # print("The previous sum_h: ",sum_h)
            # print("la suma: ",sum_h+y)
            #sum_h=sum_h+y
            sum_h = sum_h + y
            # print("The new sum_h: ",sum_h)

        #Next iteration:
        i+=1
        # print(boolean)
    return ray
