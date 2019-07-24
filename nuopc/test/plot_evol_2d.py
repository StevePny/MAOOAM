import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

ndim = 36     # total system dimension (number of modes)
na = 10       # Atmospheric modes
no = 8        # Oceanic modes
n = 1.5       # 2*Ly/Lx :: aspect ratio
Ly = 5.0*10.0**3 # pi * L (km) :: y-dim length scale
Lx = 3*Ly
nx = 300
ny = 100

# Define the basis functions:
def FA(x,y,P):
  return np.sqrt(2.0) * np.cos(P*y)

def FK(x,y,M,P):
  return 2.0 * np.cos(M*n*x) * np.sin(P*y)

def FL(x,y,H,P):
  return 2.0 * np.sin(H*n*x) * np.sin(P*y)

def phi(x,y,H0,P0):
  return 2.0 * np.sin(0.5*H0*n*x) * np.sin(P0*y)

# Define the ordering of the modes
def F(i,x,y):
  i = i+1
  # psi variables:
  if (i==1):
    return FA(x,y,P=1)
  elif (i==2):
    return FK(x,y,M=1,P=1)
  elif (i==3):
    return FL(x,y,H=1,P=1)
  elif (i==4):
    return FA(x,y,P=2)
  elif (i==5):
    return FK(x,y,M=1,P=2)
  elif (i==6):
    return FL(x,y,H=1,P=2)
  elif (i==7):
    return FK(x,y,M=2,P=1)
  elif (i==8):
    return FL(x,y,H=2,P=1)
  elif (i==9):
    return FK(x,y,M=2,P=2)
  elif (i==10):
    return FL(x,y,H=2,P=2)

  # theta variables:
  elif (i==11):
    return FA(x,y,P=1)
  elif (i==12):
    return FK(x,y,M=1,P=1)
  elif (i==13):
    return FL(x,y,H=1,P=1)
  elif (i==14):
    return FA(x,y,P=2)
  elif (i==15):
    return FK(x,y,M=1,P=2)
  elif (i==16):
    return FL(x,y,H=1,P=2)
  elif (i==17):
    return FK(x,y,M=2,P=1)
  elif (i==18):
    return FL(x,y,H=2,P=1)
  elif (i==19):
    return FK(x,y,M=2,P=2)
  elif (i==20):
    return FL(x,y,H=2,P=2)

  # A variables (ocean streamfunction)
  elif (i==21):
    return phi(x,y,H0=0.5,P0=1)
  elif (i==22):
    return phi(x,y,H0=0.5,P0=2)
  elif (i==23):
    return phi(x,y,H0=0.5,P0=3)
  elif (i==24):
    return phi(x,y,H0=0.5,P0=4)
  elif (i==25):
    return phi(x,y,H0=1.0,P0=1)
  elif (i==26):
    return phi(x,y,H0=1.0,P0=2)
  elif (i==27):
    return phi(x,y,H0=1.0,P0=3)
  elif (i==28):
    return phi(x,y,H0=1.0,P0=4)

  # T variables (ocean temperature)
  elif (i==29):
    return phi(x,y,H0=0.5,P0=1)
  elif (i==30):
    return phi(x,y,H0=0.5,P0=2)
  elif (i==31):
    return phi(x,y,H0=0.5,P0=3)
  elif (i==32):
    return phi(x,y,H0=0.5,P0=4)
  elif (i==33):
    return phi(x,y,H0=1.0,P0=1)
  elif (i==34):
    return phi(x,y,H0=1.0,P0=2)
  elif (i==35):
    return phi(x,y,H0=1.0,P0=3)
  elif (i==36):
    return phi(x,y,H0=1.0,P0=4)

  else:
    exit('F:: error: invalid index %d'%i)

def F2xy(s,xs,ys):

  # Define xs and ys, arrays of x and y coordinates
  nx = len(xs)
  ny = len(ys)

  # Initialize coordinates to 0
  psi = np.zeros((nx,ny))
  theta = np.zeros((nx,ny))
  A = np.zeros((nx,ny))
  T = np.zeros((nx,ny))

  for k in range(na):
    for i,x in enumerate(xs):
      for j,y in enumerate(ys):
#       print('k,i,j = ',k,i,j)
#       print('s[k] = ',s[k])
#       print('x = ',x)
#       print('y = ',y)
#       print('F(k,x,y) = ', F(k,x,y))
        psi[i,j]   = psi[i,j]   + s[k]    * F(k,x,y)
        theta[i,j] = theta[i,j] + s[k+na] * F(k+na,x,y)

  for k in range(no):
    for i,x in enumerate(xs):
      for j,y in enumerate(ys):
        A[i,j] = A[i,j] + s[k+2*na]    * F(k+2*na,x,y)
        T[i,j] = T[i,j] + s[k+2*na+no] * F(k+2*na+no,x,y)
   
  return psi,theta,A,T

# Test
def xy2F(psi,theta,A,T,xs,ys):

  # Apply (multiply) a basis function to the value of the function at each grid point and add them up
  s = np.zeros(ndim+1) #(first element is the time)
  dotest = False       # Should equal 1 in all fields (probably won't be exact due to numerical roundoff)
  for k in range(na):
    for i,x in enumerate(xs): 
      for j,y in enumerate(ys):
        if (dotest):
          a = F(k,x,y)
          b = F(na+k,x,y)
        else:
          a = psi[i,j]
          b = theta[i,j]
        s[1+k]    =  s[1+k]    + a * F(k,x,y)
        s[1+na+k] =  s[1+na+k] + b * F(na+k,x,y)

  for k in range(no):
    for i,x in enumerate(xs): 
      for j,y in enumerate(ys):
        if (dotest):
          a = F(2*na+k,x,y)
          b = F(2*na+no+k,x,y)
        else:
          a = A[i,j]
          b = T[i,j]
        s[1+2*na+k]    = s[1+2*na+k]    + a * F(2*na+k,x,y)
        s[1+2*na+no+k] = s[1+2*na+no+k] + b * F(2*na+no+k,x,y)

  # Average fields to get final value, since the spatial average Fi*Fi = 1
  s = s/(nx*ny)

  low_thresh = 0.1e-16
  idx = np.abs(s) < low_thresh
  s[idx] = 0

  return s

#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------

# Read in data
#infile='../SingleModelProto/evol_field.dat'
infile='../../fortran/evol_field.dat'
with open(infile, 'r') as f:
  x0 = f.read().splitlines()

S0 = np.zeros((len(x0),ndim))
for i,s0 in enumerate(x0):
  s0_ = s0.split()
  s = list(map(float,s0_[1:]))
  print('s = ')
  print (s)
  S0[i,:] = s

# Convert from spectral basis to x,y coordinates:
xs = np.arange(nx) / (nx-1) * (2.0*np.pi/n)
ys = np.arange(ny) / (ny-1) * (np.pi)
xg, yg = np.meshgrid(xs, ys)
print ('xs = ')
print (xs)
print ('ys = ')
print (ys)
psi,theta,A,T = F2xy(s,xs,ys)

# Compute the inverse field:
s2 = xy2F(psi,theta,A,T,xs,ys)
print ('s2 = ')
print (s2)

# Check inversion procedure:
#psi,theta,A,T = F2xy(s2,xs,ys)

# Compute the vector fields
ua = -np.gradient(psi.T,axis=0) # Charney and Straus (1980) set u = - d(psi)/dy
va = np.gradient(psi.T,axis=1)
uo = -np.gradient(A.T,axis=0)
vo = np.gradient(A.T,axis=1)

stride=int(np.round(ny/10))
print ('stride = ', stride)
skip=(slice(None,None,stride),slice(None,None,stride))

# Plot the data
fig = plt.figure()
ax = fig.gca()
plt.xlabel('x-dimension')
plt.ylabel('y-dimension')
#plt.imshow(psi)
c = ax.contourf(xg, yg, psi.T, cmap=cm.BrBG)
cbar = fig.colorbar(c)
q = ax.quiver(xg[skip], yg[skip], ua[skip], va[skip], units='width', pivot='mid') #, width=0.022, scale=1 / 0.15)
#q = ax.quiver(xg,yg,ua,va)
plt.title('psi (atmosphere streamfunction)')
plt.show()

fig = plt.figure()
ax = fig.gca()
plt.xlabel('x-dimension')
plt.ylabel('y-dimension')
#plt.imshow(theta)
c = ax.contourf(xg, yg, theta.T, cmap=cm.coolwarm)
cbar = fig.colorbar(c)
plt.title('theta (atmosphere temperature delta)')
plt.show()

fig = plt.figure()
ax = fig.gca()
plt.xlabel('x-dimension')
plt.ylabel('y-dimension')
#plt.imshow(A)
c = ax.contourf(xg, yg, A.T, cmap=cm.BrBG)
cbar = fig.colorbar(c)
q = ax.quiver(xg[skip], yg[skip], uo[skip], vo[skip], units='width', pivot='mid') #, width=0.022, scale=1 / 0.15)
plt.title('A (ocean streamfunction)')
plt.show()

fig = plt.figure()
ax = fig.gca()
plt.xlabel('x-dimension')
plt.ylabel('y-dimension')
#plt.imshow(T)
c = ax.contourf(xg, yg, T.T, cmap=cm.coolwarm)
cbar = fig.colorbar(c)
plt.title('T (ocean temperature delta)')
plt.show()
