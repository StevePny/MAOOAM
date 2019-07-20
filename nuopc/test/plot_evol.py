import matplotlib.pyplot as plt
import numpy as np

ndim = 36

# Read in data

#(a) add single coupled model (black)
infile='../SingleModelProto/evol_field.dat'
with open(infile, 'r') as f:
  x0 = f.read().splitlines()

S0 = np.zeros((len(x0),ndim))
for i,s0 in enumerate(x0):
  s0_ = s0.split()
  s = list(map(float,s0_[1:]))
  print('s = ')
  print (s)

  S0[i,:] = s
# plt.plot(S0[i,:],'ko')

#plt.axis([0, len(s), min(s), max(s)])
#plt.show()
  

#(b) add AtmOcnProto
infile='../AtmOcnProto/evol_atmos.dat'
with open(infile, 'r') as f:
  x1a = f.read().splitlines()

S1a = np.zeros((len(x1a),ndim))
for i,s1a in enumerate(x1a):
  s1a_ = s1a.split()
  s = list(map(float,s1a_[1:]))
  print('s = ')
  print (s)

  S1a[i,:] = s
  s = S1a[i,:] - S0[i,:]
  print('s = ')
  print (s)

  plt.plot(s,'bo')

plt.axis([0, len(s), min(s), max(s)])
plt.show()
  

infile='../AtmOcnProto/evol_ocean.dat'
with open(infile, 'r') as f:
  x1o = f.read().splitlines()

S1o = np.zeros((len(x1o),ndim))
for i,s1o in enumerate(x1o):
  s1o_ = s1o.split()
  s = list(map(float,s1o_[1:]))
  print('s = ')
  print (s)

  S1o[i,:] = s
  s = S1o[i,:] - S0[i,:]
  print('s = ')
  print (s)

  plt.plot(s,'go')

plt.axis([0, len(s), min(s), max(s)])
plt.show()

#(c) add AtmOcnMedProto
infile='../AtmOcnMedProto/evol_atmos.dat'
with open(infile, 'r') as f:
  x2a = f.read().splitlines()

S2a = np.zeros((len(x2a),ndim))
for i,s2a in enumerate(x2a):
  s2a_ = s2a.split()
  s = list(map(float,s2a_[1:]))
  print('s = ')
  print (s)

  S2a[i,:] = s
  s = S2a[i,:] - S0[i,:]
  print('s = ')
  print (s)

  plt.plot(s,'bo')

plt.axis([0, len(s), min(s), max(s)])
plt.show()
  
infile='../AtmOcnMedProto/evol_ocean.dat'
with open(infile, 'r') as f:
  x2o = f.read().splitlines()

S2o = np.zeros((len(x2o),ndim))
for i,s2o in enumerate(x2o):
  s2o_ = s2o.split()
  s = list(map(float,s2o_[1:]))
  print('s = ')
  print (s)

  S2o[i,:] = s
  s = S2o[i,:] - S0[i,:]
  print('s = ')
  print (s)

  plt.plot(s,'go')

plt.axis([0, len(s), min(s), max(s)])
plt.show()

#(d) add ingest proto
infile='../AtmOcnMedIngestFromConfigProto/evol_atmos.dat'
with open(infile, 'r') as f:
  x3a = f.read().splitlines()

S3a = np.zeros((len(x3a),ndim))
for i,s3a in enumerate(x3a):
  s3a_ = s3a.split()
  s = list(map(float,s3a_[1:]))
  print('s = ')
  print (s)

  S3a[i,:] = s
  s = S3a[i,:] - S0[i,:]
  print('s = ')
  print (s)

  plt.plot(s,'bo')

plt.axis([0, len(s), min(s), max(s)])
plt.show()
  
infile='../AtmOcnMedIngestFromConfigProto/evol_ocean.dat'
with open(infile, 'r') as f:
  x3o = f.read().splitlines()

S3o = np.zeros((len(x3o),ndim))
for i,s3o in enumerate(x3o):
  s3o_ = s3o.split()
  s = list(map(float,s3o_[1:]))
  print('s = ')
  print (s)

  S3o[i,:] = s
  s = S3o[i,:] - S0[i,:]
  print('s = ')
  print (s)

  plt.plot(s,'go')

plt.axis([0, len(s), min(s), max(s)])
plt.show()


