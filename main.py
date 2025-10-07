import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
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

def compute_mach(PM, gamma, tol=1e-8, max_iter=100):
    """Find M such that compute_PM_angle(M,gamma) == PM.
    PM must be in radians.
    Robust Newton with numerical derivative.
    """
    def f(M):
        return compute_PM_angle(M, gamma) - PM

    # initial guess: if PM small use 1.2, otherwise 2.0
    M = 2.0 if PM > 0.2 else 1.2
    for _ in range(max_iter):
        # ensure M > 1
        if M <= 1.0:
            M = 1.0001
        # numerical derivative (centered)
        eps = 1e-6 * max(1.0, M)
        f_plus = f(M + eps)
        f_minus = f(M - eps)
        df_num = (f_plus - f_minus) / (2 * eps)
        if df_num == 0:
            raise ValueError("Derivative zero in Newton iteration")
        M_new = M - f(M) / df_num
        if abs(M_new - M) < tol:
            return float(M_new)
        M = M_new
    raise ValueError("compute_mach: did not converge")

    

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

def compute_mach_angle(nulist_fan, philist_fan, j, gamma, down):
    mach_angle = np.arcsin(1 / compute_mach(nulist_fan[j], gamma))
    if down == True: #computing angle of gamma+ chars
        theta = mach_angle + philist_fan[j]
    else:
        theta = mach_angle - philist_fan[j]
    return theta


def get_entry(plot_list, fan, char):
    for entry in plot_list:
        if entry["fan"] == fan and entry["char"] == char:
            return entry
    raise KeyError(f"No entry found with fan={fan}, char={char}")


def get_max_entry(plot_list, fan, char):
    hold_list=[]
    for entry in plot_list:
        if entry["fan"] == fan and entry["char"] == char:
            hold_list.append(entry)
            
    if hold_list==[]:
        raise KeyError(f"No entry found with fan={fan}, char={char}")
    print(max(hold_list, key=lambda x: x["pair_key"]))
    return max(hold_list, key=lambda x: x["pair_key"])
def get_min_entry(plot_list, fan, char):
    hold_list=[]
    for entry in plot_list:
        if entry["fan"] == fan and entry["char"] == char:
            hold_list.append(entry)
            
    if hold_list==[]:
        raise KeyError(f"No entry found with fan={fan}, char={char}")
    return min(hold_list, key=lambda x: x["pair_key"])
    
    

def compute_fan_gamma_minus(theta, x_start, y_start,  philist_fan, nulist_fan, reflected,plot_list):

    if np.sin(theta) != 0:
                t = y_start / np.sin(theta)
                if t > 0:
                    x_end = x_start + t * np.cos(theta)
                    y_end = y_start - t * np.sin(theta)

                    entry = {
                        "theta": theta,
                        "x0": x_start, "y0": y_start,
                        "x1": x_end, "y1": y_end,
                        "type": 1,
                        "fan": i,
                        "char": j,
                        "phi": philist_fan[j],
                        "nu": nulist_fan[j],
                        "pair_key": 0,
                        "merged": False
                    }
                    plot_list.append(entry)
                    new_start_points.append((x_end, y_end))
                    new_reflected.append(True)  # mark as bounced
                else:
                    new_start_points.append((x_start, y_start))
                    new_reflected.append(reflected[j])
    else:
        new_start_points.append((x_start, y_start))
        new_reflected.append(reflected[j])
    return plot_list

def store_upwards_simple_region(i,j,philist_fan,nulist_fan,theta,plot_list,xstart=None,ystart=None):
    entry = {
        "theta": theta,  # to be filled later
        "x0": xstart, "y0": ystart,
        "x1": None, "y1": None,
        "type": 1,
        "fan": i,
        "char": j,
        "phi": philist_fan[j],
        "nu": nulist_fan[j],
        "pair_key": 0,
        "merged": False
    }
    plot_list.append(entry)
    return plot_list

def store_downwards_simple_region(i,j,philist_fan,nulist_fan,theta,plot_list,xstart=None,ystart=None):
    entry = {
        "theta": theta,  # to be filled later
        "x0": xstart, "y0": ystart,
        "x1": None, "y1": None,
        "type": 1,
        "fan": i,
        "char": j,
        "phi": philist_fan[j],
        "nu": nulist_fan[j],
        "pair_key": 0,
        "merged": False
    }
    plot_list.append(entry)
    return plot_list
def overwrite_entry(plot_list,fan,char,xend,yend):
    entry_kp1_0 = get_max_entry(plot_list, fan, char)
    entry_kp1_0["x1"] = xend
    entry_kp1_0["y1"] = yend
    # Overwrite if entry exists, otherwise append
    found = False
    for idx, entry in enumerate(plot_list):
        if entry is entry_kp1_0:
            plot_list[idx] = entry_kp1_0
            found = True
            break
    if not found:
        plot_list.append(entry_kp1_0)
    return plot_list

def compute_reflection_from_middle(plot_list,k):
    #fan is the inbound fan number (starting with 0 )
    for j in range(N_chars):
        print("new J")
        previous_char=get_max_entry(plot_list,k,j)
        theta = previous_char["theta"]
        x_start = previous_char["x0"]
        y_start = previous_char["y0"]
        if np.sin(theta) != 0:
            t = y_start / np.sin(theta)
            if t > 0:
                x_end = x_start + t * np.cos(theta)
                y_end = y_start - t * np.sin(theta)
            else:
                x_end = x_start
                y_end = y_start
        else:
            x_end = x_start
            y_end = y_start
        print(x_end,y_end)    
        plot_list=overwrite_entry(plot_list,k,j,x_end,y_end)
        previous_char=get_max_entry(plot_list,k,j)
        #get the flow characteristics with moc minus ( phi = 0 )
        #use the endpoints of 'previous_char' as start points, leave edpoints blank,
        #save this as fan (k+1) char 0 pair 1
        # get the flow characteristics with moc minus (phi = 0)
        # use the endpoints of 'previous_char' as start points, leave endpoints blank,
        # save this as fan (k+1) char 0 pair 1
        phi_new = 0
        latest_endpoints=[previous_char["x1"], previous_char["y1"]]
        nu_new = do_MOC_minus(phi1=previous_char["phi"], nu1=previous_char["nu"], phi2=phi_new)
        theta=compute_mach_angle([nu_new],[phi_new],0,1.4,True)
        if j!=N_chars-1:
            entry = {
                "theta": theta,  
                "x0": latest_endpoints[0], "y0": latest_endpoints[0],
                "x1": None, "y1": None,
                "type": 1,
                "fan": k+1,
                "char": j,
                "phi": phi_new,
                "nu": nu_new,
                "pair_key": 1,
                "merged": False
            }
            plot_list.append(entry)
        
        
        
            
            for i in range(j+1,N_chars):
                print(i-1,j)
                if i!=N_chars-1:
                    # get the k fan max pair entry (from previous loop)
                    prev_entry = get_max_entry(plot_list, k, i)
                    # get the k+1 fan max pair entry (from just before)
                    curr_entry = get_max_entry(plot_list, k+1, j)
                    # compute intersection with linalg - fill in endpoints (don't overwrite)
                    theta_prev = prev_entry["theta"]
                    theta_curr = curr_entry["theta"]
                    x0_prev, y0_prev = prev_entry["x0"], prev_entry["y0"]
                    x0_curr, y0_curr = curr_entry["x0"], curr_entry["y0"]
                    A = np.array([[np.cos(theta_prev), -np.cos(theta_curr)],
                                    [np.sin(theta_prev), -np.sin(theta_curr)]])
                    b = np.array([x0_curr - x0_prev, y0_curr - y0_prev])
                    try:
                        t, s = np.linalg.solve(A, b)
                        x_int = x0_prev + t * np.cos(theta_prev)
                        y_int = y0_prev + t * np.sin(theta_prev)
                    except np.linalg.LinAlgError:
                        x_int, y_int = x0_prev, y0_prev  # fallback if no intersection

                    overwrite_entry(plot_list,k+1,j,x_int,y_int)
                    overwrite_entry(plot_list,k,i,x_int,y_int)
                    latest_endpoints=[x_int,y_int]
                    nu_new=0.5*(prev_entry['nu']+curr_entry['nu']+prev_entry['phi']-curr_entry['phi'])
                    phi_new=0.5*(prev_entry['nu']-curr_entry['nu']+prev_entry['phi']+curr_entry['phi'])
                    theta=compute_mach_angle([nu_new],[phi_new],0,1.4,True)
                    entry = {
                        "theta": theta,  # to be filled later
                        "x0": latest_endpoints[0], "y0": latest_endpoints[0],
                        "x1": None, "y1": None,
                        "type": 1,
                        "fan": k,
                        "char": i,
                        "phi": phi_new,
                        "nu": nu_new,
                        "pair_key": j+1,
                        "merged": False
                    }
                    plot_list.append(entry)
                    theta=compute_mach_angle([nu_new],[phi_new],0,1.4,False)
                    entry = {
                        "theta": theta,  # to be filled later
                        "x0": latest_endpoints[0], "y0": latest_endpoints[0],
                        "x1": None, "y1": None,
                        "type": 1,
                        "fan": k+1,
                        "char": j,
                        "phi": phi_new,
                        "nu": nu_new,
                        "pair_key": i+1,
                        "merged": False
                    }
                    plot_list.append(entry)
                    
                    #get flow characteristics (moc minus and get them from intersection of charactersitics) for fan k pair 1 & save with start points as the previous endpoints and no endpoints
                    #get flow characteristics (moc plus and get them from intersection of charactersitics) for fan k+1 pair 2 & save with start points as the previous endpoints and no endpoints
                else:
                    

                    # Upwards char j in fan k+1
                    prev_char_entry = get_max_entry(plot_list, k+1, j)
                    entry= get_min_entry(plot_list,k+1,j)
                    x_end,y_end=prev_char_entry['x1'],prev_char_entry['y1']
                    entry['x0']=x_end
                    entry['y0']=y_end
                    for idx, en in enumerate(plot_list):
                        if en["fan"] == k and en["char"] == j and en["pair_key"] == 0:
                            plot_list[idx] = entry
                            break
                        else:
                            # If not found, append
                            plot_list.append(entry)
                        
                    
                    
        else:
            entry_kp1_0 = get_entry(plot_list, k+1, j)
            print('theentry',entry_kp1_0)

            entry_kp1_0['x0']=x_end
            entry_kp1_0['y0']=y_end
            # Overwrite based on fan and char identity
            for idx, entry in enumerate(plot_list):
                if entry["fan"] == k and entry["char"] == j and entry["pair_key"] == 0:
                    plot_list[idx] = entry_kp1_0
                    break
            else:
                # If not found, append
                plot_list.append(entry_kp1_0)
    
    return plot_list

def compute_reflection_from_shear(plot_list,k):
    #fan is the inbound fan number (starting with 0 )
    return plot_list

def compute_reflection_from_shear_line(i,j,plot_list,shear_anchor,theta_now):
    #get first shear line from phi and last start point
    #the shear anchor is the start of the last characteristic treated 
    
    #since the upwards fan is just computed and stored, we need to ste the index down by 1
    i=i-1
    # i is the index of the last upward fan
    # i+1 is the index of the current downward fan
    # i-1 is the previous downward fan
    
    if j==0: #special case for 1st char in fan
        #get start of last char in previous fan
        region = get_entry(plot_list, i-1, N_chars-1)
        shear_anchor = (region['x0'], region['y0'])
        
        #get the flow angle for the current region in the fan
        shear_line_angle = region['phi']
        nu_down_line=[region['nu']]
        if debug:
            print("shear angle for 1st char in fan:", shear_line_angle)
            print('previous line: ' , region)
    else:
        #shear anchor is the end of the previous char in the current fan
        shear_anchor = get_entry(plot_list, i, j-1)
        shear_anchor = (shear_anchor['x1'], shear_anchor['y1'])
        up_line = get_entry(plot_list, i, j-1)
        #the shear is now the phi after the previous char (downwards) in the current fan
        down_line= get_entry(plot_list, i+1, j-1)
        shear_line_angle = 0.5*(up_line['phi']+down_line['phi']-up_line['nu']+down_line['nu'])

        if debug:
            print("shear anchor for char ", j, " in fan:", shear_anchor)
    #get the endpoint of the downward char in the previous fan
    prev_char_start = get_entry(plot_list, i-1, j)
    
    x_start = prev_char_start['x1']
    y_start = prev_char_start['y1'] #y should be 0

    theta = get_entry(plot_list, i, j)['theta']
    if debug:
        print("Fan number: ", i, "Char number: ", j, "Shear angle (deg): ", np.degrees(shear_line_angle), 'Theta: ', np.degrees(theta))
        print("Previous end points:", shear_anchor)
    #solve for intersection of shear and characteristic line
    A = np.array([[np.cos(theta), -np.cos(shear_line_angle)],
                    [np.sin(theta), -np.sin(shear_line_angle)]])
    #start the new shear equation at the end of the previous line
    b = np.array([shear_anchor[0]-x_start,shear_anchor[1]-y_start])

    try:
        t, s = np.linalg.solve(A, b)
        if t > 0 and s > 0:
            x_end = x_start + t * np.cos(theta)
            y_end = y_start + t * np.sin(theta)
            #update the previously emptu upwards entry in plot list
            entry = get_entry(plot_list, i, j)
            entry["theta"] = theta
            entry["x0"] = x_start
            entry["y0"] = y_start
            entry["x1"] = x_end
            entry["y1"] = y_end
            entry["type"] = 1
            entry["fan"] = i
            entry["char"] = j
            entry["phi"] = entry["phi"] #dont change
            entry["nu"] = entry["nu"] #dont change
            entry["pair_key"] = j
            entry["merged"] = False
            plot_list.append(entry)
        else:
            if debug:
                print("No valid intersection found for fan ", i, " char ", j, "t or s < 0")
            x_end = x_start + t * np.cos(theta)
            y_end = y_start + t * np.sin(theta)
            #update the previously emptu upwards entry in plot list
            entry = get_entry(plot_list, i, j)
            entry["theta"] = theta
            entry["x0"] = x_start
            entry["y0"] = y_start
            entry["x1"] = x_end
            entry["y1"] = y_end
            entry["type"] = 1
            entry["fan"] = i
            entry["char"] = j
            entry["phi"] = entry["phi"] #dont change
            entry["nu"] = entry["nu"] #dont change
            entry["pair_key"] = j
            entry["merged"] = False
            plot_list.append(entry)
    except np.linalg.LinAlgError:
        if debug:
            print("No valid intersection found for fan ", i, " char ", j)
    #now compute the downward char reflection
    #shear anchor is at the end of the just computed upward char
    #the new start point is the end of the previous upward char
    x_start = x_end
    y_start = y_end
    theta = theta_now
    if debug:
        print("Fan number: ", i+1, "Char number: ", j, "Shear angle (deg): ", np.degrees(shear_line_angle), 'Theta: ', np.degrees(theta))
        print("Previous end points:", shear_anchor)
    #solve for intersection with y=0
    if np.sin(theta) != 0:
        t = y_start / np.sin(theta)
        if t > 0:
            x_end = x_start + t * np.cos(theta)
            y_end = y_start - t * np.sin(theta)

            entry = {
                "theta": theta,
                "x0": x_start, "y0": y_start,
                "x1": x_end, "y1": y_end,
                "type": 1,
                "fan": i+1,
                "char": j,
                "phi": philist_fan[j],
                "nu": nulist_fan[j],
                "pair_key": j,
                "merged": False
            }
            plot_list.append(entry)
            new_start_points.append((x_end, y_end))
            new_reflected.append(True)  # mark as bounced
        else:
            new_start_points.append((x_start, y_start))
            new_reflected.append(reflected[j])
    #store the shear line (starts at at shear anchor, ends at start of downward char)
    shear_entry = {
        "theta": shear_line_angle,  # shear follows phi, not theta
        "x0": shear_anchor[0], "y0": shear_anchor[1],
        "x1": x_start, "y1": y_start,
        "type": 2,
        "fan": i+1,
        "char": j,
        "phi": shear_line_angle,
        "nu": nulist_fan[j],
        "pair_key": j,
        "merged": False
    }
    plot_list.append(shear_entry)
    
    


def plotting_routine(plot_list):
    '''Plotting routine for the characteristics and shear line with the style depending on the type of line'''
    for line in plot_list:
        theta = line["theta"]
        x_start = line["x0"]    
        y_start = line["y0"]
        x_end = line["x1"]
        y_end = line["y1"]
        type = line["type"]
        fan = line["fan"]
        char = line["char"]
        phi = line["phi"]
        nu = line["nu"]
        pair_key = line["pair_key"]
        merged = line["merged"]

        #type 1 is characteristic line, type 2 is shear line
        if type == 1:
            plt.plot([x_start, x_end], [y_start, y_end], '-', color='blue')
        elif type == 2:
            plt.plot([x_start, x_end], [y_start, y_end], '--', color='gray')
            
    plt.axis('equal')
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


debug=1

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

    if len(philist)>5:
        shockwave=True
    elif len(philist)==len(nulist) and philist[-1]!=0:
        
        philist.append(0) #we are on the outer edge, moving inwards next
    elif len(philist)==len(nulist) and philist[-1]==0:
        nulist.append(PM_boundary) #we are on the inner edge, moving outwards next
    
inter_mat=np.zeros([N_chars,N_chars,len(philist)])




# Initialize first fan start point(s)
# anchor point for shear line (starts at nozzle edge)
shear_anchor = (0, a)
start_points = [shear_anchor] * N_chars   # nozzle lip
reflected = [False] * N_chars       # track if ray has bounced

for i in range(len(nulist) - 1):  # for each fan
    
    # Divide the fan into N_chars characteristics
    dphi = (philist[i] - philist[i+1]) / (N_chars - 1)


    # Create temporary instance of in-fan values
    philist_fan = [philist[i]]
    nulist_fan  = [nulist[i]]

    if i%2==0:
        down=False
        for j in range(1, N_chars):  # build characteristic angles
            phi_new = philist[i] - j * dphi
            nulist_fan.append(
                do_MOC_plus(phi1=philist_fan[-1], nu1=nulist_fan[-1], phi2=phi_new))
            philist_fan.append(phi_new)
    elif i%2==1:
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
        x_start, y_start = [0, a]  # reset start point for each char
        theta = compute_mach_angle(nulist_fan, philist_fan, j, gamma, down)  # angle of characteristic line
        if debug:
            print('CURRENT FAN:', i, 'CHAR:', j, 'THETA (deg):', np.degrees(theta), 'PHI (deg):', np.degrees(philist_fan[j]), 'NU (deg):', np.degrees(nulist_fan[j]))
        #consider the 1st fan seperately since it originates from one point
        if i == 0:
            if debug:
                print("First fan, j=", j)
            plot_list=store_downwards_simple_region(i,j,philist_fan,nulist_fan,theta,plot_list,0,a)
            #instead of computing each up and down case, compute an up and down case and do each characterstic across both cases
        #for upward case, store the flow phi and pm angle for each char
        elif i%2==1: #odd fan, upward case
            #store current characterisitc flow properties, dont do any geometry with start and end points yet
            plot_list=store_upwards_simple_region(i,j,philist_fan,nulist_fan,theta,plot_list)
            if debug:
                print("Upward case, j=", j)
        elif i%2==0: #even fan, downward case
            plot_list=store_downwards_simple_region(i,j,philist_fan,nulist_fan,theta,plot_list)  
            if debug:
                print("Downward case, j=", j)
        else:
            print("Error in fan counting")
            exit()
        
plot_list_df = pd.DataFrame(plot_list)
print(plot_list_df)
#compute non_simple_regions

for i in range(len(nulist) - 1):  # for each fan
    if i==0: #fan at the middle
        plot_list=compute_reflection_from_middle(plot_list=plot_list,k=i)
        

if debug:
    plot_list_df = pd.DataFrame(plot_list)
    print("_______________________")
    print("Plotting list storing all lines to be plotted (theta, x_start, y_start, x_end, y_end, type (1 for char 2 for shear), fan, char):")
    print(plot_list_df)
    print(inter_mat)
    print("Do you wish to save it as a csv? (y/n)")
    save_csv = input().lower()
    if save_csv == 'y':
        plot_list_df.to_csv('plot_list.csv', index=False)
        print("Saved as plot_list.csv")
    else:
        print("Not saved.")
    
