############################################
##
## Ray Tracing and Optical Collection Simulations
##
## Group A - Underwater Optical Communication
## Sponsored by Northrup Grumman
##
## Version 5.0.1
##
############################################

##all distances are in mm
##all angles are in radians


##Constants
pi = 3.1415926535897
a_small = 1.5
R = 3000
THETA_STEP = .001
MAX_ANGLE = (int)(pi/2/THETA_STEP)
PATHS_PER_ANGLE = 1000
N = 2


##Cone Geometry

##Wall 1 : m = -tan(theta_cone)
##         g = -a_small*tan(theta_cone/2)
##Wall 2 : m =  tan(theta_cone)
##         g = -a_small*tan(theta_cone/2)


############################################
## Bounce Calculations
############################################


def find_intersection(m, g, m1, g1, m2, g2, *intersection):
    x1 = (g1 - g) / (m - m1)
    x2 = (g2 - g) / (m - m2)
    y1 = m1*x1 + g1
    y2 = m2*x2 + g2

    if -g/m < a_small and -g/m > -a_small: ##If the ray intersects the collector (y=0) within the -a_small < x < a_small range
        intersection[0] = -g/m
        intersection[1] = 0
    elif y1 < 0 :
        intersection[0] = x1
        intersection[1] = y1

    elif y2 < 0:
        intersection[0] = x2
        intersection[1] = y2

    elif y1 > y2:
        intersection[0] = x2
        intersection[1] = y2

    else:
        intersection[0] = x1
        intersection[1] = y1

    return



def bounce(m,g):
    theta = pi/2 - atan(m)
    m1 = -tan(theta_cone)
    m2 = -m1
    g1 = -a_small*tan(theta_cone/2)
    g2 = g1
    find_intersection(m, g, m1, g1, m2, g2, intersection)

    ## 0 -> x    1 -> y

    angle_to_surf = theta + theta_cone/2
    theta_new = theta + theta_cone/2
    m_new = tan(pi/2 - theta_new);
    g_new = intersection[1] - (a_big*m_new)

    if intersection[1] <= 0:
        return true

    if intersection[1] >= L:
        return false

    return bounce (m_new, g_new)

############################################
##Utility Functions
############################################

def radtodeg(rad):
    return rad*180/pi

def degtorad(deg):
    return deg/180*pi

############################################
##Main Program
############################################

def main():

    f = open('Output.txt','w')
    f.write('Starting Cone Simulation...\n\n')
    theta_max = degtorad(25)

    ############################################
    ##Calculate the Geometry and Parameters
    ############################################

    a_big = a_small/(sin(theta_max))
    L = (a_small+a_big)/(tan(theta_max))
    theta_cone = 2*atan((a_big - a_small)/L)
    b = L/(a_big*a_big - a_small*a_small)
    n = a_big/a_small
    G = n*n
    theta_cutoff = atan(a_big/R)

    print ('Theta Max: %d\n'  % theta_max)
    print ('Small Aperture:  %d\n'  % a_small)
    print ('Big Aperture: %d\n'  % a_big)
    print ('Distance from Source:  %d\n'  % R)
    print ('Collector Length:  %d\n'  % L)
    print ('Theta Cone:  %d\n'  % radtodeg(theta_cone))
    print ('Parabola Coefficient (b): %d\n'  % b)
    print ('"Aperture Ratio: %d\n'  % n)
    print ('Area Ratio: %d\n'  % G)
    print ('Theta Cutoff: %d\n'  % radtodeg(theta_cutoff))

    f.write('All Parameters Calculated\n\n')

    ############################################
    ##Normalize The Power
    ############################################

    ##i*THETA_STEP is the actual angle

    f.write('Normalizing Power Model...\n\n')

    sum = 0;
    power[MAX_ANGLE];
    weight[MAX_ANGLE];

    for i in range(0, MAX_ANGLE):
        power[i] = pow(cos(i*THETA_STEP),N)
        sum +=pow(cos(i*THETA_STEP),N)


    for i in range(0, MAX_ANGLE):
        power[i] /= sum;


    for i in range(0, MAX_ANGLE):
        f.write('Angle : %d    Power: %d\n\n' % (i*THETA_STEP , power[i]))

    ############################################
    ##Simulate (Ray-Tracing)
    ############################################

    ## y = mx + g

    ## m = tan(pi/2 - theta) * a_big
    ## g = L - (a_big*m) + path_height

    f.write('Starting Ray Tracing...')

    for i in range(MAX_ANGLE): ##For Each Theta
        theta = i*THETA_STEP; ## theta incident
        total=0;
        for j in range(PATHS_PER_ANGLE):  ##For Each (of the 1000) paths per Theta i

            if(i*THETA_STEP!=0):
                m = tan(pi/2 - theta)
                g = L - (a_big*m) + j/(2*a_big*m)
                    if(bounce(m,g)):
                        total+=1

            else:
                total+=1;

        weight[i] = ((double)total)/PATHS_PER_ANGLE

    ############################################
    ##Integrate Results (Sum of Multiplication)
    ############################################

    f.write('Summing the Multiplication...')

    total_power = 0

    for i in range(MAX_ANGLE):
        power_captured[i] = power[i]*weight[i]
        total_power += power_captured[i]


    ############################################
    ##Interrupt Results
    ############################################


    ##power over all transmitted
    ##power over power without collector
    ##power over power incident upon collector

    print "Power Captured: %d" % total_power

    ##These are linear powers for comparison only
    ##For accurate 3D Power Simulations we need to consider an integral

    f.close()
    return
