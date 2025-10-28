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

def get_entry_by_key(plot_list, fan, char, pair_key):
    for entry in plot_list:
        if (entry["fan"] == fan and 
            entry["char"] == char and 
            entry["pair_key"] == pair_key):
            return entry
    raise KeyError(f"No entry found with fan={fan}, char={char}, pair_key={pair_key}")


def get_max_entry(plot_list, fan, char):
    hold_list=[]
    for entry in plot_list:
        if entry["fan"] == fan and entry["char"] == char:
            hold_list.append(entry)
            
    if hold_list==[]:
        raise KeyError(f"No entry found with fan={fan}, char={char}")
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
def replace_entry(plot_list, fan, char, pair_key, new_entry):
    """
    Replaces an entry in plot_list that matches fan, char, and pair_key.
    If no such entry exists, appends the new one.
    """

    for idx, entry in enumerate(plot_list):
        if (entry["fan"] == fan and 
            entry["char"] == char and 
            entry["pair_key"] == pair_key):
            plot_list[idx] = new_entry
            break
    else:
        # Only runs if no break occurs → entry not found
        plot_list.append(new_entry)

    return plot_list


def compute_reflection_from_middle(plot_list,k):
    #fan is the inbound fan number (starting with 0 )
    for j in range(N_chars):
        entry=get_max_entry(plot_list,k,j) #get i=j just before the ground C-
        # Compute intersection of characteristic line with y=0
        x_start = entry['x0']
        y_start = entry['y0']
        theta = entry['theta']
        if np.sin(theta) != 0:
            t = y_start / np.sin(theta)
            if t > 0:
                x_end = x_start + t * np.cos(theta)
                y_end = y_start - t * np.sin(theta)
                entry["x1"] = x_end
                entry["y1"] = y_end
                plot_list=replace_entry(plot_list,entry['fan'],entry['char'],entry['pair_key'],entry)
        
        
        #treat last char seperately
        if j!=N_chars-1:
            print("Running inner loop")
            phi_loc=0
            region_above=get_max_entry(plot_list,k,j)
            nu_loc=region_above['phi']+region_above['nu'] #since phi=0
            theta=compute_mach_angle([nu_loc],[phi_loc],0,1.4,False)
            entry = {
                "theta": theta,  # to be filled later
                "x0": x_end, "y0": y_end,
                "x1": None, "y1": None,
                "type": 1,
                "fan": k+1,
                "char": j,
                "phi": phi_loc,
                "nu": nu_loc,
                "pair_key": 1,
                "merged": False
            }
            plot_list.append(entry)
            
            for i in range (j+1,N_chars):
                if i!=N_chars-1:
                    print(i-1,j)
                    entry_above=get_max_entry(plot_list,k,i) #gamma -
                    entry_below=get_max_entry(plot_list,k+1,j) #gamma +
                    print("entry above",entry_above)
                    print("entry Below",entry_below)
                    # Find intersection of entry_above and entry_below characteristic lines
                    theta_above = entry_above['theta']
                    theta_below = entry_below['theta']
                    x0_above, y0_above = entry_above['x0'], entry_above['y0']
                    x0_below, y0_below = entry_below['x0'], entry_below['y0']

                    A = np.array([
                        [np.cos(theta_above), -np.cos(theta_below)],
                        [np.sin(theta_above), -np.sin(theta_below)]
                    ])
                    b = np.array([x0_below - x0_above, y0_below - y0_above])

                    try:
                        t_above, t_below = np.linalg.solve(A, b)
                        x_inter = x0_above + t_above * np.cos(theta_above)
                        y_inter = y0_above + t_above * np.sin(theta_above)
                        # Update entry_above and entry_below with intersection point
                        entry_above['x1'] = x_inter
                        entry_above['y1'] = y_inter
                        entry_below['x1'] = x_inter
                        entry_below['y1'] = y_inter
                        plot_list = replace_entry(plot_list, entry_above['fan'], entry_above['char'], entry_above['pair_key'], entry_above)
                        plot_list = replace_entry(plot_list, entry_below['fan'], entry_below['char'], entry_below['pair_key'], entry_below)
                    except np.linalg.LinAlgError:
                        print("Error in linalg")
                        pass
                    
                    
                    phi_loc=0.5*(entry_above['phi']+entry_above['nu']+entry_below['phi']-entry_below['nu'])
                    nu_loc=0.5*(entry_above['phi']+entry_above['nu']-entry_below['phi']+entry_below['nu'])
                    theta_up=compute_mach_angle([nu_loc],[phi_loc],0,1.4,True)
                    theta_down=compute_mach_angle([nu_loc],[phi_loc],0,1.4,False)
                    
                    #add upwards_line
                    entry = {
                        "theta": theta_up,  # to be filled later
                        "x0": x_inter, "y0": y_inter,
                        "x1": None, "y1": None,
                        "type": 1,
                        "fan": k+1,
                        "char": j,
                        "phi": phi_loc,
                        "nu": nu_loc,
                        "pair_key": entry_below['pair_key']+1,
                        "merged": False
                    }
                    plot_list.append(entry)
                    #add downwards_line
                    entry = {
                        "theta": theta_down,  # to be filled later
                        "x0": x_inter, "y0": y_inter,
                        "x1": None, "y1": None,
                        "type": 1,
                        "fan": k,
                        "char": i,
                        "phi": phi_loc,
                        "nu": nu_loc,
                        "pair_key": entry_above['pair_key']+1,
                        "merged": False
                    }
                    plot_list.append(entry)
                    
                else:
                    entry_above=get_max_entry(plot_list,k,i) #gamma -
                    entry_below=get_max_entry(plot_list,k+1,j) #gamma +
                    print("entry abouve",entry_above)
                    print("entry Below",entry_below)
                    # Find intersection of entry_above and entry_below characteristic lines
                    theta_above = entry_above['theta']
                    theta_below = entry_below['theta']
                    x0_above, y0_above = entry_above['x0'], entry_above['y0']
                    x0_below, y0_below = entry_below['x0'], entry_below['y0']

                    A = np.array([
                        [np.cos(theta_above), -np.cos(theta_below)],
                        [np.sin(theta_above), -np.sin(theta_below)]
                    ])
                    b = np.array([x0_below - x0_above, y0_below - y0_above])

                    try:
                        t_above, t_below = np.linalg.solve(A, b)
                        x_inter = x0_above + t_above * np.cos(theta_above)
                        y_inter = y0_above + t_above * np.sin(theta_above)
                        # Update entry_above and entry_below with intersection point
                        entry_above['x1'] = x_inter
                        entry_above['y1'] = y_inter
                        entry_below['x1'] = x_inter
                        entry_below['y1'] = y_inter
                        plot_list = replace_entry(plot_list, entry_above['fan'], entry_above['char'], entry_above['pair_key'], entry_above)
                        plot_list = replace_entry(plot_list, entry_below['fan'], entry_below['char'], entry_below['pair_key'], entry_below)
                    except np.linalg.LinAlgError:
                        print("Error in linalg")
                        pass
                    
                    
                    phi_loc=0.5*(entry_above['phi']+entry_above['nu']+entry_below['phi']-entry_below['nu'])
                    nu_loc=0.5*(entry_above['phi']+entry_above['nu']-entry_below['phi']+entry_below['nu'])
                    theta_up=compute_mach_angle([nu_loc],[phi_loc],0,1.4,True)
                    theta_down=compute_mach_angle([nu_loc],[phi_loc],0,1.4,False)
                    #add downwards_line
                    entry = {
                        "theta": theta_down,  # to be filled later
                        "x0": x_inter, "y0": y_inter,
                        "x1": None, "y1": None,
                        "type": 1,
                        "fan": k,
                        "char": i,
                        "phi": phi_loc,
                        "nu": nu_loc,
                        "pair_key": entry_above['pair_key']+1,
                        "merged": False
                    }
                    plot_list.append(entry)
                    check_entry=get_entry_by_key(plot_list,k+1,j,0)
                    print("CHECK HERE: ",check_entry['theta'],theta_up)
                    entry = {
                    "theta": theta_up,  # to be filled later
                    "x0": x_inter, "y0": y_inter,
                    "x1": None, "y1": None,
                    "type": 1,
                    "fan": k+1,
                    "char": j,
                    "phi": phi_loc,
                    "nu": nu_loc,
                    "pair_key": 0, #set to 0 for simple region
                    "merged": False
                    }
                    plot_list=replace_entry(plot_list,check_entry['fan'],check_entry['char'],0,entry)
                    
            
        else:
            entry=get_entry_by_key(plot_list,k+1,j,0)
            entry['x0']=x_end
            entry['y0']=y_end
            plot_list=replace_entry(plot_list,k+1,j,0,entry)
            print("Running else")
            
        plot_list_df = pd.DataFrame(plot_list)
        print(plot_list_df)
    
    
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
        
plotting_routine(plot_list)
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
    
