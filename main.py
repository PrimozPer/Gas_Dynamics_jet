import numpy as np
import matplotlib.pyplot as plt
def get_user_in():
    N_chars=input("Insert deisced number of characteristics: ")
    #test if N_chars is an integer
    try:
        N_chars=int(N_chars)
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return None
    return N_chars

def compute_PM_angle(M, gamma):
    if M < 1:
        raise ValueError("Mach number must be greater than 1 for Prandtl-Meyer expansion.")
    nu = np.sqrt((gamma+1)/(gamma-1)) * np.arctan(np.sqrt((gamma-1)*(M**2 - 1)/(gamma+1))) - np.arctan(np.sqrt(M**2 - 1))
    
    return nu

def compute_mach(PM,gamma):
    for M in np.arange(1.0, 3.0, 0.01):
        nu=compute_PM_angle(M,gamma)
        if np.isclose(nu, PM, atol=1e-2):
            Mach=M
            break
    if Mach is not None:
        return Mach
    else:
        raise ValueError("No Mach number found for given Prandtl-Meyer angle.")
    

def do_MOC_plus(phi1=None, nu1=None, phi2=None, nu2=None):
    '''Method of Characteristics
    Inputs: three of phi1, nu1, phi2, nu2
    Outputs: missing input
    calculates the missing input using the MOC equations
    '''
    if phi1 is not None and nu1 is not None and phi2 is not None:
        nu2 = nu1 - phi1 + phi2
        return nu2
    elif phi1 is not None and nu1 is not None and nu2 is not None:
        phi2 = phi1 - nu1 + nu2
        return phi2
    elif phi1 is not None and phi2 is not None and nu2 is not None:
        nu1 = nu2 - phi2 + phi1
        return nu1
    elif nu1 is not None and phi2 is not None and nu2 is not None:
        phi1 = phi2 - nu2 + nu1
        return phi1
    pass

#############################


#############################
#####INITIALIZATION##########
#############################

pa=101325 #Pa
pe=2*pa
Me=2
gamma=1.4 #air
philist=[]
nulist=[]
###inside the fan###
philist_fan = []
nulist_fan = []

#flow total pressure
p_tot = pe*(1+(gamma-1)/2*Me**2)**(gamma/(gamma-1))

#boundary mach number
M_boundary = np.sqrt((2/(gamma-1))*((p_tot/pa)**((gamma-1)/gamma)-1))
PM_outlet = compute_PM_angle(Me, gamma)
phi_outlet = 0
PM_boundary = compute_PM_angle(M_boundary, gamma)
philist.append(phi_outlet)
nulist.append(PM_outlet)
nulist.append(PM_boundary)
shockwave=False



######get the user input#######

N_chars = get_user_in()
print("Input ", N_chars, "is of type", type(N_chars))
if N_chars is not None:
    print("Valid input received.")
else:
    print("No valid input received.")
    
#############################
#####MAIN LOOP###############
#############################

#get states in uniform region
while shockwave==False:
    #fill in missing values
    if len(philist)<len(nulist):
        #get phi
        phi_new=do_MOC_plus(nu1=nulist[-2], phi1=philist[-1], nu2=nulist[-1])
        philist.append(phi_new)

    if len(nulist)<len(philist):
        #get nu
        nu_new=do_MOC_plus(nu1=nulist[-2], phi1=philist[-1], nu2=nulist[-1])
        nulist.append(nu_new)
    if len(philist)==len(nulist):
        shockwave=True
        
for i in range(len(nulist)-1): ##for each fan

    dphi=(philist[i]-philist[i+1])/(N_chars-1)
    #create a temporary instance of in-fan values
    philist_fan.append(philist[i])
    nulist_fan.append(nulist[i])
    for j in range(1,N_chars): ##for each char in the fan
        phi_new=philist[i]-j*dphi
        nulist_fan.append(do_MOC_plus(phi1=philist_fan[-1], nu1=nulist_fan[-1], phi2=phi_new))
        philist_fan.append(phi_new)
    plot_list=np.arange(len(nulist_fan))



# print("philist: ", np.degrees(philist))
# print("nulist: ", np.degrees(nulist))
# print("philist_fan: ", np.degrees(philist_fan))
# print("nulist_fan: ", np.degrees(nulist_fan))

a = 5  # starting y coordinate (point at x=0, y=a)

# Plot setup
plt.figure(figsize=(6,6))
plt.axhline(0, color='black', linewidth=1)  # x-axis
plt.axvline(0, color='gray', linestyle='--')  # y-axis reference
plt.scatter(0, a, color='red', label=f"Start (0,{a})")

# Loop over angles
for i in range(len(philist_fan)):
    # Adjust angle so that we substract mach angle
    mach_angle = np.arcsin(1 / compute_mach(nulist_fan[i], gamma))
    print("Mach angle (deg): ", np.degrees(mach_angle))
    print("Phi (deg): ", np.degrees(philist_fan[i]))
    
    theta = mach_angle -  philist_fan[i]  # angle in radians
    print("Theta (deg): ", np.degrees(theta))
    # Parametric line: (x,y) = (0,a) + t*(cosθ, sinθ), with t >= 0
    # We want intersection with y=0:
    #   0 = a + t*sinθ  =>  t = -a/sinθ
    if np.sin(theta) != 0:
        t = a / np.sin(theta)
        if t > 0:  # only forward rays
            x_end = t * np.cos(theta)
            plt.plot([0, x_end], [a, 0], label=f"{np.degrees(philist_fan[i])   }°")
        else:
            print("Ray does not intersect y=0 in the positive direction for Phi =", np.degrees(philist_fan[i]))
    else:
        # Vertical line case (phi=90 or 270 after shift)
        # It will never hit y=0 unless a=0
        print("Vertical line, no intersection with y=0 at Phi =", np.degrees(philist_fan[i]))
        pass

plt.legend()

plt.xlabel("x")
plt.ylabel("y")
plt.title("Lines from (0,a) at given angles")
plt.grid(True)
plt.show()