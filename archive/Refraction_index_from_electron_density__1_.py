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
    
    for i in range(0,80,1): #Estaba 50 pero así los night daban negativo
        N[i] = 1

    vacuum_permittivity = 8.854e-12 # Permittivity of free space in Faradays over meters
    e_charge = -1.602e-19 # Electron charge in Coulombs
    e_mass = 9.10938356e-31 # Electron mass in Kilograms
    
# The electron density has to be given in units of meters^(-3)

    for i in range(80,total_y,1):
        if ((((elec_dens(i)+delta)*e_charge*e_charge)/(e_mass*vacuum_permittivity))<=0):
            N[i]=1
        else:
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
    i=1

    #A boolean to stop the iterations:
    boolean=True

    while boolean: 
        #The n index depends on the height:
        n0=n[current_index]
        x=step_t*c/n0*np.cos(angle)
        y=step_t*c/n0*np.sin(angle)
        #We need to see if we have reached the next height with this last step
        #To do so, we check if the sum is bigger than the grid size:
        
        
        if (abs(sum_h+y)>=delta_y):
            k=1
            while((abs(sum_h+y)-k*delta_y)>delta_y):
                k+=1
            n1=np.mean(n[current_index+1:current_index+k+1])
            current_index+=k #We go to the next grid
            #See if changes refraction index:
            if (n1==n0):
                
                sum_h+= y - delta_y #Goes to next cell normal
                ray.append([ray[i-1][0]+x,ray[i-1][1]+y])
            else:
                #Snell law to calculate new angle.
                #CAUTION: Snell angle is taken from the normal and ours is from the horizontal axis:
                #We need to do sin(90-theta)=cos(theta)
                angle0=angle
                #CAREFUL: arcsin input has to be between -1 and 1. If it exceedes, add/substract 1
                if (n1==0):
                    #If the ray has passed the plasma without bouncing back
                    #We can't accept this ray, it's null:
                    if (current_index>=1200):
                        ray.append([-1,-1])
                    boolean=False
                    break
                elif (abs(((n0*np.cos(angle0))/n1)>=1)):
                    #There has been an internal reflection, the wave can't pass. We stop
                    boolean=False
                    break
                else:
                    
                    angle=(np.pi/2.)-np.arcsin((n0*np.cos(angle0))/n1)
                    #Now we need to calculate the height.
                    #It divides into two heights:
                    #1)height in previous cell: delta_y-sum_h
                    #2)height in new cell with new angle and n:
                    y1=(delta_y-sum_h)
                    
                    y2=(c/n1)*np.sin(angle)*(step_t-((delta_y-sum_h)/(c*np.sin(angle0)/n0)))
                    y=y1+y2
                    
                    #The y can be negative, but the x is always positive (it moves forward)
                    x=abs(y1/np.tan(angle0))+abs(y2/np.tan(angle))
                    
                    #We put the sum back to current height in next cell, to count again:
                    sum_h=y2 
                    #We add the new total position of our ray:
                    ray.append([ray[i-1][0]+x,ray[i-1][1]+y])

        else:
            #If it doesn't jump to next cell the ray will travel normally:
            ray.append([ray[i-1][0]+x,ray[i-1][1]+y])
            sum_h = sum_h + y

        #Next iteration:
        i=i+1
        # print(boolean)
    return ray

def GetRayPolar(theta0,n):
    c = 1
    delta_r = 1
    
    #The radius of the earth:
    R_Earth=6370
    # initial angle, our initial point is in the "top" of the earth
    phi = np.pi/2 
    #We choose the initial position of the ray:
    
    # Getting the coordinates in (r,phi) polar coordinates
    # Should be easier
    ray=[[R_Earth , phi]]
    
    #We just do iterations:
    step_t=0.5 #constant, time passing between iterations
    current_index=0 #To see in which pixel we are
    sum_r=0 #The sum of the distance from the centre, to see if it passes to the next cell

    #We track te iterations number too:
    i=1

    #A boolean to stop the iterations:
    boolean=True
    
    # Angle of rotation of each step of the Snell law
    ubre = np.pi/2 - phi #Can be negative(don't know)

    while boolean: 
        #The n index depends on the height:
        n0=n[current_index]
        x0=step_t*(c/n0)*np.cos(theta0)
        y0=step_t*(c/n0)*np.sin(theta0)
        
        #As we are in a sphere, we need to rotate these values clockwise
        x = x0 * np.cos(ubre) - y0*np.sin(ubre)
        y = x0 * np.sin(ubre) + y0*np.cos(ubre)
        

        #We need to see if we have reached the next height with this last step
        #To do so, we check if the difference in modulus is bigger than the grid size:
        mod_r0 =ray[i-1][0]
        angle_0=ray[i-1][1]
        
        # Vector coordinates
        r1 = [mod_r0*np.cos(angle_0) + x, mod_r0*np.sin(angle_0)  + y] 
        #The modulus:
        mod_r1 = np.sqrt(r1[0]**2 + r1[1]**2)
    
    
        if (abs(sum_r + (mod_r1 - mod_r0)>=delta_r)):
            k=1
            while(abs(sum_r+(mod_r1 - mod_r0)-k*delta_r)>delta_r):
                k+=1
            n1=np.mean(n[current_index+1:current_index+k+1])
           
            current_index+=k #We go to the next grid
            #See if changes refraction index:
            if (n1==n0):
        
                sum_r=sum_r+ (mod_r1-mod_r0) - k*delta_r #Goes to next cell normal
                phi=np.arctan2(r1[1],r1[0])
               
                #There is no need for rotation because theta is still expressed in previous axis
                
                ray.append([mod_r1,phi])
    
            else:
                #Snell law to calculate new angle.
                #CAUTION: Snell angle is taken from the normal and ours is from the horizontal axis:
                
                alpha0=np.pi/2-theta0
                #CAREFUL: arcsin input has to be between -1 and 1. If it exceedes, add/substract 1
                if (n1==0):
                    #If the ray has passed the plasma without bouncing back
                    #We can't accept this ray, it's null:
                    if (current_index>=1200):
                        ray.append([-1,-1])
                    boolean=False
                    break
                elif (abs((n0*np.sin(alpha0))/n1)>=1):
                    #There has been an internal reflection, the wave can't pass. We stop
                    boolean=False
                    break
                else:
                    alpha=np.arcsin((n0*np.sin(alpha0))/n1)
                    theta=np.pi/2-alpha
                    # We divide the total step into 2:
                    
                    y0=delta_r-sum_r
                    x0=y0/np.tan(theta0)
                    
                    y1=c/n1*np.sin(theta)*(step_t-(c/n0)*y0/np.sin(theta0))
                    x1=y1/np.tan(theta0)
                    
                    sum_r=y1
                    
                    #Now rotate each of them manualy:
                    x_prev0=ray[i-1][0]*np.cos(ray[i-1][1])
                    y_prev0=ray[i-1][0]*np.sin(ray[i-1][1])
                    phi_0=np.arctan2(y_prev0,x_prev0)
                    ubre_0=np.pi/2-phi_0
                    
                    x0_rot = x0 * np.cos(ubre_0) - y0*np.sin(ubre_0)
                    y0_rot = x0 * np.sin(ubre_0) + y0*np.cos(ubre_0)
                    
                    # To do the next rotation we need to find the position after first step:
                    x_step0=x_prev0+x0_rot
                    y_step0=y_prev0+y0_rot
                    
                    phi_1=np.arctan2(y_step0,x_step0)
                    ubre_1=np.pi/2-phi_1
                    
                    #We add this little step to our ray, for resolution purposes:
                    r_step0=np.sqrt(x_step0**2+y_step0**2)
                    
                    
                    ray.append([r_step0,phi_1])
                    
                    #We need to add 1 more iteration to compensate this extra point:
                    i=i+1
                    
                    x1_rot = x1 * np.cos(ubre_1) - y1*np.sin(ubre_1)
                    y1_rot = x1 * np.sin(ubre_1) + y1*np.cos(ubre_1)
                    
                    #finally, we need the final angle of the ray
                    #To know the rotation of the next points:
                    
                    x_step1=x_step0+x1_rot
                    y_step1=y_step0+y1_rot
                    
                    phi_fin=np.arctan2(y_step1,x_step1)
                    r_fin=np.sqrt(x_step1**2+y_step1**2)
                    #We need to see if after Snell the ray is higher or lower than before:
                    if (r_fin<ray[i-2][0]):
                        boolean=False
                        break
                    ray.append([r_fin,phi_fin])
                    
                    #For next iterations, the rotation will be given by the angle of the last point:
                    ubre=np.pi/2-phi_fin
                    
                    #The new theta will be the last one:
                    theta0=theta
                    
            
        else:
            #If it doesn't jump to next cell the ray will travel normally:
            # The quantity that advance is the difference of r1 and r0
            # The angle is the one of the advanced ray
            
            phi=np.arctan2(r1[1],r1[0])
            ray.append([mod_r1, phi])
            sum_r+=(mod_r1-mod_r0)

            #There is no need for rotation because theta is still expressed in previous axis
            
            
        #Next iteration:
        i=i+1
        
    return ray