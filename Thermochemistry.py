import sys
import math
from subprocess import Popen, STDOUT, PIPE

#
# Thermochemistry 1.0
#
# Author: Dr. Pere Miro
# Contact: pere.miro@gmail.com
#
# Execute: python Thermochemistry.py input.dat checkpoint.chk temp(k) press(atm) solid replace transitionstate
#
# The name of all these files can be different but the order matters
#
# input.dat = Contains system variables
# checkpoint.chk = Gaussian checkpoint file
# output.out = Gaussian output
#

#
# Function that calculate the thermochemistry
#

def thermo(T,P,mass,frequencies,Spin,inx,iny,inz):

  global TotUvib
  global TotSvib
  global TotZPE
  global Srot
  global Urot
  global Strans
  global Utrans
  global Selec
  global Uelec

# pi constant
  pi=3.14159265359

# Boltzman constant kJ/K
  boltzman=1.380648813*(10**(-23))

# h in kJ.s
  planck=6.6260695729*(10**(-34))

# Combination of constants 
#  - Avogadro's number in molecules/mol
#  - h in kJ.s
#  - speed of light in cm/2
  variable1=0.011962656594658400

# Combination of constants (theta)
#  - h in kJ.s
#  - speed of light in cm/s
#  - Boltzman constant kJ/K     
  variable2=1.438776947066190
      
  for a in range (0,int(len(frequencies))):

# Calculation of the vibrational contribution to the internal energy 
    if(frequencies[a]!=0.0):
      Uvib=0.008314462175*((variable2*(float(frequencies[a])))*(0.5+(1/(math.exp((variable2*(float(frequencies[a]))/T))-1))))
    if(frequencies[a]==0.0):
      Uvib=0.0
          
# Calculation of the vibrational contribution to the enthropy 
    if(frequencies[a]!=0.0):
      Svib=0.008314462175*(((variable2*(float(frequencies[a]))/T)/(math.exp((variable2*(float(frequencies[a]))/T))-1))-(math.log(1-math.exp(-(variable2*(float(frequencies[a]))/T)))))
    if(frequencies[a]==0.0):
      Svib=0.0 
      
# Calculation of the zero point energy contribution of the given frequency
    ZPE=0.5*((float(frequencies[a]))*variable1)

# Sums of the total vibrational contribution to the internal energy and the entropy 
    TotUvib=TotUvib+Uvib
    TotSvib=TotSvib+Svib
    
# Sums of the total zero point energy  
    TotZPE=TotZPE+ZPE

# Mass from a.m.u. to Kg
  mass=mass*1.66053892*(10**(-27))
         
# Calculation of the rotational partition function and the rotational contribution to the entropy
# The eigenvalues of the moment of inertia matrix are transformed from a.m.u. and Bohrs to Kg and meters
  inx=inx*1.66053892*(10**(-27))*5.29177249*(10**(-11))*5.29177249*(10**(-11))
  iny=iny*1.66053892*(10**(-27))*5.29177249*(10**(-11))*5.29177249*(10**(-11))
  inz=inz*1.66053892*(10**(-27))*5.29177249*(10**(-11))*5.29177249*(10**(-11))
  intertia=math.sqrt(pi*inx*iny*inz)
  qr=intertia*(((8*pi*pi*boltzman*T)/(planck*planck))**(1.5))  
  
# Calculation of the rotational contribution to entropy   
  Srot=0.008314462175*(math.log(qr)+1.5)

# Calculation of the rotational contribution to the internal energy 
  Urot=0.0124716932625*T

# Calculation of the translational partition function and the translational contribution to the entropy
  qt=(((2*pi*mass*boltzman*T)/(planck*planck))**(1.5))*((boltzman*T)/P)
  Strans=0.008314462175*((math.log(qt))+1+1.5)
  
# Calculation of the translational contribution to the internal energy 
  Utrans=0.0124716932625*T

# Calculation of the electronic contribution to the entropy  
  Selec=0.008314462175*(math.log(Spin))

  return

#
# Function to calculate the moment of inertia and the molecular weight
#

def intertia(natoms):

  global mass
  global inx
  global iny
  global inz

# Library with the atomic masses
  atomicmass= {'H':'1.00783','O':'15.99491','Ni':'57.93535','C':'12.0000','Co':'58.93320','Pd':'105.90320','Rh':'102.9048','Ti':'47.94795','Cu':'62.92960','N':'14.00307'}

# Open the file with the geometry
  file = open('geometry.xyz', 'r+')

#Initialize the variables that will contain the coordinates and atom label
  x=[]
  y=[]
  z=[]
  name=[]

#Initialize the variables that will contain the center of mass
  centerofmassx=0
  centerofmassy=0
  centerofmassz=0

#Loop to read the geometry and store the x, y and z components
  for a in range (0,int(natoms)):
    line = file.readline()
    columns = line.split()
    name.append(columns[0])
    x.append(float(columns[1]))
    y.append(float(columns[2]))
    z.append(float(columns[3]))

# Calculate the center of mass
    centerofmassx=centerofmassx+(float(atomicmass[columns[0]])*float(columns[1]))
    centerofmassy=centerofmassy+(float(atomicmass[columns[0]])*float(columns[2]))
    centerofmassz=centerofmassz+(float(atomicmass[columns[0]])*float(columns[3]))

  centerofmassx=centerofmassx/mass
  centerofmassy=centerofmassy/mass
  centerofmassz=centerofmassz/mass

# Close the file with the geometry
  file.close()

# Translate the center of mass to (0,0,0)
  for a in range (0,int(natoms)):
    x[a]=x[a]-centerofmassx
    y[a]=y[a]-centerofmassy
    z[a]=z[a]-centerofmassz
   
# Initialize the moment of inertia 
  Ixx=0
  Ixy=0
  Ixz=0
  Iyx=0
  Iyy=0
  Iyz=0
  Izx=0
  Izy=0
  Izz=0

# Compute the values of each component in the moment of inertia matrix
# The 1.8897 is the Angstroms to Bohr conversion
  for a in range (0,int(natoms)):
    Ixx=Ixx+(float(atomicmass[name[a]])*((y[a]*y[a]*1.8897*1.8897)+(z[a]*z[a]*1.8897*1.8897)))
    Ixy=Ixy+(float(atomicmass[name[a]])*(x[a]*y[a]*1.8897*1.8897))
    Ixz=Ixz+(float(atomicmass[name[a]])*(x[a]*z[a]*1.8897*1.8897))
    Iyx=Iyx+(float(atomicmass[name[a]])*(y[a]*x[a]*1.889725989*1.889725989))
    Iyy=Iyy+(float(atomicmass[name[a]])*((x[a]*x[a]*1.889725989*1.889725989)+(z[a]*z[a]*1.8897*1.8897)))
    Iyz=Iyz+(float(atomicmass[name[a]])*(y[a]*z[a]*1.889725989*1.889725989))
    Izx=Izx+(float(atomicmass[name[a]])*(z[a]*x[a]*1.889725989*1.889725989))
    Izy=Izy+(float(atomicmass[name[a]])*(z[a]*y[a]*1.889725989*1.889725989))
    Izz=Izz+(float(atomicmass[name[a]])*((x[a]*x[a]*1.889725989*1.889725989)+(y[a]*y[a]*1.8897*1.8897)))

# Put the correct signs into the moment of inertia matrix
  Ixx=Ixx
  Ixy=-Ixy
  Ixz=-Ixz
  Iyx=-Iyx
  Iyy=Iyy
  Iyz=-Iyz
  Izx=-Izx
  Izy=-Izy
  Izz=Izz

# Open the file with the moment of the inertia matrix
  file = open("momentinertia.dat","wb")

# Write the inertia matrix into the file
  file.write(str(Ixx)+'\n')
  file.write(str(Ixy)+'\n')
  file.write(str(Ixz)+'\n')
  file.write(str(Iyx)+'\n')
  file.write(str(Iyy)+'\n')
  file.write(str(Iyz)+'\n')
  file.write(str(Izx)+'\n')
  file.write(str(Izy)+'\n')
  file.write(str(Izz)+'\n')

# Close the file with the moment of the inertia matrix
  file.close()

# Call the fortran subroutine that diagonalize the moment of inertia matrix
  cmd = './diagonalization.exe'
  ret = Popen([cmd], shell=True, stdout=PIPE, stderr=STDOUT)
  (stdout, stderr) = ret.communicate()

# Open the file with the eigenvalues
  file = open('eigenvalues.dat', 'r+')

# Read the eigenvalues
  inx = float(file.readline())
  iny = float(file.readline())
  inz = float(file.readline())

# Close the file with the eigenvalues
  file.close()

  return

# Get the input variables
T=float(sys.argv[3])
P=float(sys.argv[4])*101325
solid=int(sys.argv[5])

# Load Gaussian in Quest
cmd = "module load Gaussian"
ret = Popen([cmd], shell=True, stdout=PIPE, stderr=STDOUT)
(stdout,stderr) = ret.communicate()

# Use newzmat to generate the xyz
cmd = "newzmat -ichk -step 9999 -oxyz "+ str(sys.argv[2]) + " geometry.xyz"
ret = Popen([cmd], shell=True, stdout=PIPE, stderr=STDOUT)
(stdout,stderr) = ret.communicate()

# Get the total electronic energy
Energy=0.0
cmd = "grep -i 'SCF Done:  ' "+ str(sys.argv[1])+ "| tail -1 | awk '{ print $5 }'"
ret = Popen([cmd], shell=True, stdout=PIPE, stderr=STDOUT)
(stdout,stderr) = ret.communicate()
Energy=float(stdout)

# Get the number of atoms
natoms=0
cmd = "wc -l geometry.xyz | awk '{ print $1 }' "
ret = Popen([cmd], shell=True, stdout=PIPE, stderr=STDOUT)
(stdout,stderr) = ret.communicate()
natoms=int(stdout)

# Get the multiplicity
Spin=0
cmd = "grep -i 'Multiplicity' "+ str(sys.argv[1]) + "| tail -1 | awk '{ print $6}' "
ret = Popen([cmd], shell=True, stdout=PIPE, stderr=STDOUT)
(stdout,stderr) = ret.communicate()
Spin=float(stdout)

# Get the frequencies
temp1=[]
frequencies=[]
cmd = "grep -i 'Frequencies --  ' "+ str(sys.argv[1]) + "| awk '{ print $3}' "
ret = Popen([cmd], shell=True, stdout=PIPE, stderr=STDOUT)
(stdout,stderr) = ret.communicate()
temp1=stdout.split()
for a in range(0,len(temp1)):
  frequencies.append(float(temp1[a]))
temp1=[]
cmd = "grep -i 'Frequencies --  ' "+ str(sys.argv[1]) + "| awk '{ print $4}' "
ret = Popen([cmd], shell=True, stdout=PIPE, stderr=STDOUT)
(stdout,stderr) = ret.communicate()
temp1=stdout.split()
for a in range(0,len(temp1)):
  frequencies.append(float(temp1[a]))
temp1=[]
cmd = "grep -i 'Frequencies --  ' "+ str(sys.argv[1]) + "| awk '{ print $5}' "
ret = Popen([cmd], shell=True, stdout=PIPE, stderr=STDOUT)
(stdout,stderr) = ret.communicate()
temp1=stdout.split()
for a in range(0,len(temp1)):
  frequencies.append(float(temp1[a]))

# Order the frequencies from small to large

frequencies.sort()

Replace=float(sys.argv[6])
TransitionState=int(sys.argv[7])

# Loop that replaces low frequencies. It can be deactivated above.  

count=0  
for a in range (0,int(len(frequencies))):  
 if(Replace != 0):
   if(frequencies[a] < Replace):
     if(TransitionState == 0):
        frequencies[a]=Replace 
	count=count+1
     if((TransitionState == 1)and(a==0)):
        frequencies[a]=0.0 
     if((TransitionState == 1)and(a!=0)):
        frequencies[a]=Replace 
	count=count+1   
 if(Replace == 0):
   if((TransitionState == 1)and(a==0)):
      frequencies[a]=0.0

for a in range (0,int(len(frequencies))):
 if(Replace == 0):
   if(frequencies[a] < 0.0):
    print ''
    print 'There are negative frequencies and the replacements are not on.'
    print 'You  MUST turn the replacements on to continue.'
    print ''
    sys.exit()

# Initialization of variables
mass=0.0
inx=0.0
iny=0.0
inz=0.0
TotUvib=0
TotSvib=0
TotZPE=0
Srot=0
Urot=0
Strans=0
Utrans=0
Selec=0
Uelec=0

# Get the mass
mass=0
cmd = "grep -i 'Molecular mass' "+ str(sys.argv[1]) + "| awk '{ print $3}' "
ret = Popen([cmd], shell=True, stdout=PIPE, stderr=STDOUT)
(stdout,stderr) = ret.communicate()
mass=float(stdout)

# Execute the function inertia
intertia(natoms)

# Execute the function thermo
thermo(T,P,mass,frequencies,Spin,inx,iny,inz)

# Calculate the correections, remember that RT is required to get the enthalpy
Enthaply=0.0
Entropy=0.0

if solid==0:
  Enthalpy=((TotUvib+Urot+Utrans+Uelec+(0.0083144621*T)))/(4.184*27.211*23.06)
  Entropy=-((TotSvib+Srot+Strans+Selec)*T)/(4.184*27.211*23.06)
if solid==1:
  Enthalpy=((TotUvib+Uelec+(0.0083144621*T)))/(4.184*27.211*23.06)
  Entropy=-((TotSvib+Selec)*T)/(4.184*27.211*23.06)
  
# Print section
print ''
print T,'Kelvin'
print P,'Pascals'
print ''

if solid==0:
  print 'Gas phase species'
if solid!=0:
  print 'Solid species. qrot=0 and qtrans=0'

if Replace==0:
  print 'No low frequencies were replaced'
  for a in range (0,int(len(frequencies))):  
    if(frequencies[a] < 25.0):
     if(TransitionState == 0):
       print 'Warining frequency number ',a+1,'is ', frequencies[a],' cm-1'  
     if(TransitionState != 0):
       if(a != 0):
         print 'Warining frequency number ',a+1,'is ', frequencies[a],' cm-1'     
if (Replace != 0):
  print count,' low frequencies were replaced to ',Replace,'  cm-1 ones'
print ''
if (TransitionState == 0):
  print 'This is a minima'
if (TransitionState == 1):
  print 'This is a transition state'
print ''
 
print 'Enthalpy (hartree)'
print Energy+Enthalpy

if solid==0:
  print 'Gibbs Free Energy with Full Rot/Trans Entropy (hartree)'
  print Energy+Enthalpy+Entropy
  Entropy=-((TotSvib+((((Srot+Strans)*2)/3))+Selec)*T)/(4.184*27.211*23.06)

if solid==1:
  print 'Gibbs Free Energy (hartree)'
  print Energy+Enthalpy+Entropy
print ''




