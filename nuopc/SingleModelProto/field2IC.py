import re

# Read in either the evol_field.dat or mean_field.dat files and read in the state vectors
# Output in the format of the IC.dat file
# https://stackabuse.com/read-a-file-line-by-line-in-python/

# Set index for output:
oidx = 0

datfile='evol_field.dat'
with open(datfile, 'r') as fp:
  for cnt, line in enumerate(fp):
    print("Line {}: {}".format(cnt, line))
    
    # If this is the index we want, then make it the new Initial Condition (IC) in IC.nml
    if (cnt == oidx):
      IC = [float(s) for s in line.split()]
      print('using index oidx = ', oidx)
      print('IC = ')
      print(IC)
      print('len(IC) = ', len(IC))
      break

# Read in IC.nml file and replace IC definitions with new IC values
ICfile='IC.nml.orig'
newfile='IC.nml.new'
#  IC(1) =   -0.0000000000000000         ! typ= A, Nx= 0.0, Ny= 1.0
# https://regexone.com/problem/matching_decimal_numbers
p = re.compile(r'-?\d+\.\d+')
idx=0
fo = open(newfile,'w')
with open(ICfile, 'r') as fp:
  for cnt, line in enumerate(fp):
    print("Line {}: {}".format(cnt, line))
    
    # If this is the index we want, then make it the new Initial Condition (IC) in IC.nml
    ICstr = 'IC('+str(idx+1)+')'
    print('check ICstr = ', ICstr)
    if (ICstr in line):
      new = "{:30.29e}".format(IC[idx])
      print('new = ', new)
      line = p.sub(new,line,1)
      print('new line = ')
      print(line)      
      idx = idx + 1

    # Write new line to file
    fo.write(line)

fo.close()
