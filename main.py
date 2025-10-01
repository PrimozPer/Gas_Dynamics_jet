import numpy as np
import matplotlib.pyplot as plt
import time
def get_user_in():
    N_chars=input("Insert deisced number of characteristics: ")
    #test if N_chars is an integer
    try:
        N_chars=int(N_chars)
    except ValueError:
        print("Invalid input. Please enter an integer.")
        #try again
        return get_user_in()
    return N_chars

def compute_PM_angle(M, gamma):
    if M < 1:
        raise ValueError("Mach number must be greater than 1 for Prandtl-Meyer expansion.")
    nu = np.sqrt((gamma+1)/(gamma-1)) * np.arctan(np.sqrt((gamma-1)*(M**2 - 1)/(gamma+1))) - np.arctan(np.sqrt(M**2 - 1))
    
    return nu

def compute_mach(PM, gamma, tol=1e-6, max_iter=100):
    def f(M):
        return compute_PM_angle(M, gamma) - PM
    
    def df(M):
        # derivative of Prandtl-Meyer angle wrt M
        return np.sqrt(M**2 - 1) / (1 + 0.5*(gamma-1)*M**2) / M
    
    M = 2.0  # initial guess
    for _ in range(max_iter):
        M_new = M - f(M)/df(M)
        if abs(M_new - M) < tol:
            return M_new
        M = M_new
    
    raise ValueError("Did not converge")

    

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

def do_MOC_minus(phi1=None, nu1=None, phi2=None, nu2=None):
    '''Method of Characteristics
    Inputs: three of phi1, nu1, phi2, nu2
    Outputs: missing input
    calculates the missing input using the MOC equations
    '''
    if phi1 is not None and nu1 is not None and phi2 is not None:
        nu2 = nu1 + phi1 - phi2
        return nu2
    elif phi1 is not None and nu1 is not None and nu2 is not None:
        phi2 = phi1 + nu1 - nu2
        return phi2
    elif phi1 is not None and phi2 is not None and nu2 is not None:
        nu1 = nu2 + phi2 - phi1
        return nu1
    elif nu1 is not None and phi2 is not None and nu2 is not None:
        phi1 = phi2 + nu2 - nu1
        return phi1
    pass

def compute_mach_angle(nulist_fan, j, gamma, down):
    mach_angle = np.arcsin(1 / compute_mach(nulist_fan[j], gamma))
    if down == True: #computing angle of gamma+ chars
        theta = mach_angle + philist_fan[j]
    else:
        theta = mach_angle - philist_fan[j]
    return theta

def compute_fan_gamma_minus(theta, x_start, y_start, N_chars, philist_fan, nulist_fan, start_points, reflected,plot_list):
    if np.sin(theta) != 0:
                t = y_start / np.sin(theta)
                if t > 0:
                    x_end = x_start + t * np.cos(theta)
                    y_end = y_start - t * np.sin(theta)

                    plot_list.append((theta,x_start, y_start, x_end, y_end, 1, i, j)) #type 1 is characteristic line, type 2 is shear line

                    new_start_points.append((x_end, y_end))
                    new_reflected.append(True)  # mark as bounced
                else:
                    new_start_points.append((x_start, y_start))
                    new_reflected.append(reflected[j])
    else:
        new_start_points.append((x_start, y_start))
        new_reflected.append(reflected[j])
    return new_reflected

def compute_fan_gamma_plus(theta, x_start, y_start, N_chars, philist_fan, nulist_fan, start_points, reflected, shear_anchor,plot_list):
    # -------- UPWARD case: intersect shear line --------
    length = a * 5  # finite length
    x_edge = shear_anchor[0] + length * np.cos(philist_fan[j])
    y_edge = shear_anchor[1] + length * np.sin(philist_fan[j])

    plot_list.append((theta,shear_anchor[0], shear_anchor[1], x_edge, y_edge, 2, i, j))  # type 2 is shear line, i and j are indices of fan and char

    #get the flow angle for the current region in the fan

    phi_flow = philist_fan[j]
    
    
    if debug:
        print("Fan number: ", i, "Char number: ", j, "Phi flow (deg): ", np.degrees(phi_flow), 'Theta: ', np.degrees(theta))
        print("Previous end points:", shear_anchor)
    #make new line equation for shear line
    A = np.array([[np.cos(theta), -np.cos(phi_flow)],
                    [np.sin(theta), -np.sin(phi_flow)]])
    #start the new shar equation at the ennd of the previous line
    b = np.array([shear_anchor[0]-x_start,shear_anchor[1]-y_start])


    try:
        t, s = np.linalg.solve(A, b)
        if t > 0 and s > 0:
            x_end = x_start + t * np.cos(theta)
            y_end = y_start + t * np.sin(theta)


            plot_list.append((theta,x_start, y_start, x_end, y_end, 1,i,j)) #type 1 is characteristic line, type 2 is shear line

            new_start_points.append((x_end, y_end))
            shear_anchor=(x_end,y_end)
            new_reflected.append(True)
        else:
            new_start_points.append((x_start, y_start))
            new_reflected.append(reflected[j])
    except np.linalg.LinAlgError:
        new_start_points.append((x_start, y_start))
        new_reflected.append(reflected[j])
    return shear_anchor

def plotting_routine(plot_list):
    '''Plotting routine for the characteristics and shear line with the style depending on the type of line'''
    for line in plot_list:
        theta, x_start, y_start, x_end, y_end, type, fan, char = line

        #type 1 is characteristic line, type 2 is shear line
        if type == 1:
            plt.plot([x_start, x_end], [y_start, y_end], '-', color='blue')
        elif type == 2:
            plt.plot([x_start, x_end], [y_start, y_end], '--', color='gray')
            

    plt.legend()
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.title("Characteristics and Shear Lines")
    plt.grid()
    plt.show()

#############################


#############################
#####INITIALIZATION##########
#############################


debug=0

pa=101325 #Pa
pe=2*pa
Me=2
gamma=1.4 #air
philist=[]
nulist=[]
###inside the fan###
philist_fan = []
nulist_fan = []
a = 1  # starting y coordinate (point at x=0, y=a)
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

plot_list=[]

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
start_time = time.time()
print("Computing...")
#get states in uniform region
while shockwave==False:
    #fill in missing values
    

    if len(philist)<len(nulist): #nu is known, phi is unknown (outer edge)
        #get phi
        phi_new=do_MOC_plus(nu1=nulist[-2], phi1=philist[-1], nu2=nulist[-1]) #following gamma +
        philist.append(phi_new)

    if len(nulist)<len(philist):#phi is known, nu is unknown (inner edge)
        #get nu
        nu_new=do_MOC_minus(nu1=nulist[-1], phi1=philist[-2], phi2=philist[-1]) #following gamma -
        nulist.append(nu_new)

    if len(philist)>2:
        shockwave=True
    elif len(philist)==len(nulist) and philist[-1]!=0:
        
        philist.append(0) #we are on the outer edge, moving inwards next
    elif len(philist)==len(nulist) and philist[-1]==0:
        nulist.append(PM_boundary) #we are on the inner edge, moving outwards next
    
    




# Initialize first fan start point(s)
# anchor point for shear line (starts at nozzle edge)
shear_anchor = (0, a)
start_points = [shear_anchor] * N_chars   # nozzle lip
reflected = [False] * N_chars       # track if ray has bounced

for i in range(len(nulist) - 1):  # for each fan
    
    # ## Shear Line ##
    # length = a * 5  # finite length
    # x_edge = length * np.cos(philist[i])
    # y_edge = a + length * np.sin(philist[i])
    # plt.plot([0, x_edge], [a, y_edge], '--', color='gray')

    
    # Divide the fan into N_chars characteristics
    dphi = (philist[i] - philist[i+1]) / (N_chars - 1)


    # Create temporary instance of in-fan values
    philist_fan = [philist[i]]
    nulist_fan  = [nulist[i]]

    if philist[i+1]>philist[i]:
        down=False
        for j in range(1, N_chars):  # build characteristic angles
            phi_new = philist[i] - j * dphi
            nulist_fan.append(
                do_MOC_plus(phi1=philist_fan[-1], nu1=nulist_fan[-1], phi2=phi_new))
            philist_fan.append(phi_new)
    elif philist[i+1]<philist[i]:
        if debug:
            print("Downwards expansion")
        down=True
        for j in range(1, N_chars):  # build characteristic angles
            phi_new = philist[i] - j * dphi
            nulist_fan.append(
                do_MOC_minus(phi1=philist_fan[-1], nu1=nulist_fan[-1], phi2=phi_new))
            philist_fan.append(phi_new)
        
    new_start_points = []
    new_reflected = []
    
    if debug:
        print('Fan_phi', philist_fan)
        print("Nulist fan", nulist_fan)
    for j in range(N_chars):
        x_start, y_start = start_points[j]

        theta = compute_mach_angle(nulist_fan, j, gamma, down)  # angle of characteristic line

        # refl = -1 → before bounce, refl = +1 → after bounce
        refl = 1 if reflected[j] else -1  

        if refl == -1:
            # -------- DOWNWARD case: intersect y=0 --------
            compute_fan_gamma_minus(theta, x_start, y_start, N_chars, philist_fan, nulist_fan, start_points, reflected,plot_list)

        else:
            # -------- UPWARD case: intersect shear line --------
            shear_anchor=compute_fan_gamma_plus(theta, x_start, y_start, N_chars, philist_fan, nulist_fan, start_points, reflected, shear_anchor, plot_list)
        #for reflection from the shear line, the new start point also depends on the next, downwards facing characteristic and essentially, the flow expands downwards again, making a new phi
        

    # update
    start_points = new_start_points
    reflected = new_reflected#reverse to maintain bottom-to-top order


plot_list=np.array(plot_list)

if debug:
    print("_______________________")
    print("Plotting list storing all lines to be plotted (theta, x_start, y_start, x_end, y_end, type (1 for char 2 for shear), fan, char):")
    print("Plot list shape: ", plot_list.shape)
    print(plot_list)

print('Computed, see graph.')
end_time = time.time()
print(f"Execution time: {end_time - start_time} seconds")
plotting_routine(plot_list)
#end timing

print(f"Execution time: {end_time - start_time:.4f} seconds")

