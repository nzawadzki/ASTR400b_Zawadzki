# Homework 4 
# Center of Mass Position and Velocity

# import modules
import numpy as np
import astropy.units as u
from ReadFile import Read


class CenterOfMass:
   
    def __init__(self, filename, ptype):
        # read in the file and particle type
        self.time, self.total, self.data = Read(filename)
            
        #create an array to store indexes of particles of desired Ptype
        self.index = np.where(self.data['type'] == ptype)
    
        # store the mass, positions, velocities of only the particles of the given type
        self.m = self.data['m'][self.index]
    
        ##### PLACE other particle properties here: x,y,z,vx,vy,vz #####
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]

        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]


    def total_mass(self, m):
        #Note: you can add other keyword arguments into the function, but 'self' must be first
        return np.sum(m)
        
    ##### PLACE OTHER FUNCTIONS BELOW #####
    #function that returns the coordinates of the center of mass
    def COMdefine(self, i, j, k, mass):
        xCOM = np.sum(i*mass) / self.total_mass(mass)
        yCOM = np.sum(j*mass) / self.total_mass(mass)
        zCOM = np.sum(k*mass) / self.total_mass(mass)

        return xCOM, yCOM, zCOM

    #function to refine COM
    def COM_P(self, d):
        #initializing tolerance and position variables
        delta = 100
        xi = self.x
        yi = self.y
        zi = self.z
        mi = self.m

        #first guess of Center of Mass
        COMxi, COMyi, COMzi = self.COMdefine(xi,yi,zi,mi)
        rCOM = np.sqrt(COMxi**2 + COMyi**2 + COMzi**2)
        
        #change into center of mass frame
        #and find magnitude of particles from center of mass
        rnew = np.sqrt((xi-COMxi)**2 + (yi-COMyi)**2 + (zi-COMzi)**2)
        #find max distance--needed when using smaller volume to recalculate COM
        rmax = np.amax(rnew)

        #loop to converge center of mass calculation
        while delta > d:
            #index to find particles within smaller volume
            index = np.where(rnew < (rmax/2))
            #store particles that fall within smaller volume
            xnew = xi[index]
            ynew = yi[index]
            znew = zi[index]
            mnew = mi[index]

            #calculate center of mass again with smaller volume
            xCOM2, yCOM2, zCOM2 = self.COMdefine(xnew, ynew, znew, mnew)

            #print xCOM2, yCOM2, zCOM2

            #storing magnitude of vector  
            rCOM2 = np.sqrt(xCOM2**2 + yCOM2**2 + zCOM2**2)

            #find difference to see if COM calculation has converged
            delta = abs(rCOM-rCOM2)
            #print delta

            #redefine variables and go back to beginning of loop
            xnew = xi - xCOM2
            ynew = yi - yCOM2
            znew = zi - zCOM2
            rnew = np.sqrt(xnew**2 + ynew**2 + znew**2)
            rmax = rmax/2

            rCOM = rCOM2
            xCOM = xCOM2
            yCOM = yCOM2
            zCOM = zCOM2

        return xCOM, yCOM, zCOM

    def COM_V(self):
        #temporary variables to store center of mass position
        tempx, tempy, tempz = self.COM_P(1)
        xcom = self.x - tempx
        ycom = self.y - tempy
        zcom = self.z - tempz
        r = np.sqrt(xcom**2 + ycom**2 + zcom**2)
        
        #store velocities of particles within 15kpc of the COM
        #index to find particles
        index = np.where(r < 15)

        #velocity arrays
        VX = self.vx[index]
        VY = self.vy[index]
        VZ = self.vz[index]
        M = self.m[index]

        #calculate COM velocity using COMdefine
        vxCOM, vyCOM, vzCOM = self.COMdefine(VX, VY, VZ, M)

        return vxCOM, vyCOM, vzCOM
   
# EXAMPLE OF USING A CLASS
##########################
#Answers to Question 1
# Create a Center of mass object for the MW
print ""
MWCOM = CenterOfMass("MW_000.txt", 2)
# Calculate quantities for MW data
COM_MWx, COM_MWy, COM_MWz = MWCOM.COM_P(1) #calculated with tolerance of 3 pc
print "1. MW Disk CoM:", (COM_MWx, COM_MWy, COM_MWz)*u.kpc 

vCOM_MWx, vCOM_MWy, vCOM_MWz = MWCOM.COM_V()
print "MW Disk CoM Velocity:", (vCOM_MWx, vCOM_MWy, vCOM_MWz)*u.km/u.s

#Create a Center of mass object for M31
M31COM = CenterOfMass("M31_000.txt", 2)
print ""
#Calculate quantities for M31 data
COM_M31x, COM_M31y, COM_M31z = M31COM.COM_P(1)
print "M31 Disk CoM:", (COM_M31x, COM_M31y, COM_M31z)*u.kpc

vCOM_M31x, vCOM_M31y, vCOM_M31z = M31COM.COM_V()
print "M31 Disk CoM Velocity:", (vCOM_M31x, vCOM_M31y, vCOM_M31z)*u.km/u.s

#Create a Center of mass object for M33
M33COM = CenterOfMass("M33_000.txt", 2)
print ""
#Calculate quantities for M33 data
COM_M33x, COM_M33y, COM_M33z = M33COM.COM_P(1)
print "M33 Disk CoM:", (COM_M33x, COM_M33y, COM_M33z)*u.kpc

vCOM_M33x, vCOM_M33y, vCOM_M33z = M33COM.COM_V()
print "M33 Disk CoM Velocity:", (vCOM_M33x, vCOM_M33y, vCOM_M33z)*u.km/u.s
print ""

#Answers to Question 2
rCOM_MW = np.sqrt(COM_MWx**2 + COM_MWy**2 + COM_MWz**2)
rCOM_M31 = np.sqrt(COM_M31x**2 + COM_M31y**2 + COM_M31z**2)
P_diff = abs(rCOM_MW - rCOM_M31)
print "2. Current MW/M31 position separation:", P_diff*u.kpc
vCOM_MW = np.sqrt(vCOM_MWx**2 + vCOM_MWy**2 + vCOM_MWz**2)
vCOM_M31 = np.sqrt(vCOM_M31x**2 + vCOM_M31y**2 + vCOM_M31z**2)
V_diff = abs(vCOM_MW - vCOM_M31)
print "Current MW/M31 velocity separation:", V_diff*u.km/u.s
print ""

#Answer to Question 3
#Take the diffference of each component, then find the magnitude of that 
psepx = abs(COM_M31x - COM_M33x)
psepy = abs(COM_M31y - COM_M33y)
psepz = abs(COM_M31z - COM_M33z)
rsep = np.sqrt(psepx**2 + psepy**2 + psepz**2)
print "3. Current M31/M33 position separation:", rsep*u.kpc
vsepx = abs(vCOM_M31x - vCOM_M33x)
vsepy = abs(vCOM_M31y - vCOM_M33y)
vsepz = abs(vCOM_M31z - vCOM_M33z)
vsep = np.sqrt(vsepx**2 + vsepy**2 + vsepz**2)
print "Current M31/M33 velocity separation:", vsep*u.km/u.s
print ""
#Answer to Question 4
print '4. The iterative process is important because when the galaxies begin to colllide they will change shape due to tidal forces. This will change where the center of mass is located.'
print ""
