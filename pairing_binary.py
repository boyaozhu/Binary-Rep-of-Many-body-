
import math
import numpy as np
from sys import argv
import time as tm
import itertools

# bitCount() counts the number of bits set (not an optimal function)

def bitCount(int_type):
    """ Count bits set in integer """
    count = 0
    while(int_type):
        int_type &= int_type - 1
        count += 1
    return(count)


# testBit() returns a nonzero result, 2**offset, if the bit at 'offset' is one.

def testBit(int_type, offset):
    mask = 1 << offset
    return(int_type & mask) >> offset

# setBit() returns an integer with the bit at 'offset' set to 1.

def setBit(int_type, offset):
    mask = 1 << offset
    return(int_type | mask)

# clearBit() returns an integer with the bit at 'offset' cleared.

def clearBit(int_type, offset):
    mask = ~(1 << offset)
    return(int_type & mask)

# toggleBit() returns an integer with the bit at 'offset' inverted, 0 -> 1 and 1 -> 0.

def toggleBit(int_type, offset):
    mask = 1 << offset
    return(int_type ^ mask)

# binary string made from number

def bin0(s):
    return str(s) if s<=1 else bin0(s>>1) + str(s&1)

def bin(s, L = 0):
    ss = bin0(s)
    if L > 0:
        return '0'*(L-len(ss)) + ss
    else:
        return ss
    
    

class Slater:
    """ Class for Slater determinants """
    def __init__(self, x = int(0)):
        self.word = int(x)

    def create(self, j):
        #print "c^+_" + str(j) + " |" + bin(self.word) + ">  = ",
        # Assume bit j is set, then we return zero.
        s = 0
        # Check if bit j is set.
        isset = testBit(self.word, j)
        if isset == 0:
            bits = bitCount(self.word & ((1<<j)-1))
            s = pow(-1, bits)
            self.word = setBit(self.word, j)
        else:
          self.word = 0
          return self.word
            
        #print str(s) + " x |" + bin(self.word) + ">"
        return self.word
        
    def annihilate(self, j):
        #print "c_" + str(j) + " |" + bin(self.word) + ">  = ",
        # Assume bit j is not set, then we return zero.
        s = 0
        # Check if bit j is set.
        isset = testBit(self.word, j)
        if isset == 1:
            bits = bitCount(self.word & ((1<<j)-1))
            s = pow(-1, bits)
            self.word = clearBit(self.word, j)
        else:
          self.word = 0
          return self.word

        #print str(s) + " x |" + bin(self.word) + ">"
        return self.word

def basis(particles, Nstates):
  N_sp   = []
  states = []
  temp   = []
  M      = []
  for i in range(Nstates):
    N_sp.append(i)
  
  for basis in itertools.combinations(N_sp, particles):
    temp.append(basis)
  for i in temp:
    phi = Slater()
    m   = 0
    for j in i:
      if j%2 == 0:
        m += 1
      else:
        m -= 1
      a = phi.create(j)
    states.append(a)
    M.append(m)

  return states, M
'''
def basis3(particles, Nstates, Slater):
  states = []
  for i in range(Nstates):
    phi = Slater()
    a = phi.create(i)
    states.append(a)
  k = 0
  print (states)
  temp = states
  print (temp)
  while 1:
    k += 1
    for i in temp:
      for j in range(len(states)):
        states.append(states[j]+i)
    temp = states
    print (states)
    if k == particles:
      break
  f = []
  for i in range(len(states)):
    if bitCount(states[i] == particles):
      f.append(states[i])
  return f
# create basis
def basis(particles, Nstates, Slater):
  states = []
  M      = []
  spin = 1
  for i in range(Nstates-particles+1):
      
    for j in range(i+1,Nstates-particles+2):
        
      for k in range(j+1,Nstates-particles+3):
        for l in range(k+1,Nstates-particles+4):
          phi = Slater()
          phi.create(i)
          phi.create(j)
          phi.create(k)
          a = phi.create(l)
          states.append(a)
          if i%2 == 0:
            m1 = spin
          else:
            m1 = -spin
          if j%2 == 0:
            m2 = spin
          else:
            m2 = -spin
          if k%2 == 0:
            m3 = spin
          else:
            m3 = -spin
          if l%2 == 0:
            m4 = spin
          else:
            m4 = -spin
          
          temp = m1 + m2 + m3 + m4
          M.append(temp)


  return states, M

'''
# collect states based on M
def M_scheme(states, M):
  states_0 = []
  states_2 = []
  states_2_= []
  states_4 = []
  states_4_= []
  for i in range(len(M)):
    
    if M[i] == 0:
      states_0.append(states[i])
    if M[i] == 2:
      states_2.append(states[i])
    if M[i] == -2:
      states_2_.append(states[i])
    if M[i] == 4:
      states_4.append(states[i])
    if M[i] == -4:
      states_4_.append(states[i])

  return states_0, states_2, states_2_, states_4, states_4_

# For demonstration
def states_binary(states_0):
  states_binary = []
  for i in states_0:
    states_binary.append(bin(i))
  return states_binary

# indexing of basis
def construct_index(particles, Nstates, basis, M_scheme):
  
  states, M = basis(particles, Nstates)
  states_0, states_2, states_2_, states_4, states_4_ = M_scheme(states, M)
  index = { }
  
  for i, state in enumerate(states_0):
    index[state] = i
  return index


def Hamiltonian(delta, g, f, user_data):
  
  index     = user_data["index"]
  particles = user_data["particles"]
  Nstates   = user_data["Nstates"]
  basis     = user_data["basis"]
  Slater    = user_data["Slater"]
  M_scheme  = user_data["M_scheme"]
  
  states, M = basis(particles, Nstates)
  states_0, states_2, states_2_, states_4, states_4_ = M_scheme(states, M)

  # H = delta * /sum(a'a) -  g * /sum(a'a'aa) - f * /sum(a'a'ab + H.c.)
  H = np.zeros((len(states_0),len(states_0)))
  
  # one - body interaction
  for i in range(len(states_0)):
    N = 0
    j = 0
    while 1:
      a = testBit(states_0[i],j)
      if a == 1:
        H[i,i] += delta*np.floor_divide(j, 2)
        N += 1
      if N == particles:
        break
      j += 1

  # pick up two pair states
  stat_0 = []
  for i in range(len(states_0)):
    for l in range(Nstates):
      if l%2 == 1:
        continue
    
      a = testBit(states_0[i],l)
      b = testBit(states_0[i],l+1)
      if a == 1 and b == 1:
        if l == 6:
          stat_0.append(states_0[i])
        continue
      if a == 0 and b == 0:
        if l == 6:
          stat_0.append(states_0[i])
        continue
      else:
        break

  # pairing interaction
  for i in range(len(stat_0)):
    for j in range(Nstates):
      if j%2 == 1:
        continue
      phi = Slater(stat_0[i])
      a = phi.annihilate(j)
      a = phi.annihilate(j+1)
      if a == 0:
        continue
      for k in range(Nstates):
        if k%2 == 1:
          continue
        ph = Slater(a)
        b = ph.create(k)
        b = ph.create(k+1)
        if b == 0:
          continue
        #print (b)
        for m in range(len(stat_0)):
          if b == stat_0[m]:
            H[index[stat_0[m]],index[stat_0[i]]] -= 0.5 * g
            break
          else:
            continue

  ## particle - hole interaction
  for i in range(len(states_0)):
    
    for q in range(Nstates):
      if q%2 == 1:
        continue
      phi = Slater(states_0[i])
      a = phi.annihilate(q)
      if a == 0:
        continue
      
      for r in range(Nstates):
        if r%2 == 0 or r == q+1:
          continue
        phh = Slater(a)
        c = phh.annihilate(r)
        if c == 0:
          continue
      
        for p in range(Nstates):
          if p%2 == 1:
            continue
          ph = Slater(c)
          b = ph.create(p)
          b = ph.create(p+1)
          if b == 0:
            continue
       
          
          for m in range(len(states_0)):
            if b == states_0[m]:
              H[index[states_0[m]],index[states_0[i]]] -= 0.5 * f
              H[index[states_0[i]],index[states_0[m]]] -= 0.5 * f
              break
            else:
              continue
  return H
'''
  for i in range(len(states_0)):
    for p in range(Nstates):
      if p%2 == 1:
        continue
      phi = Slater(states_0[i])
      d = phi.annihilate(p)
      d = phi.annihilate(p+1)
      if d == 0:
        continue
      for q in range(Nstates):
        if q%2 == 1:
          continue
        ph = Slater(d)
        e = ph.create(q)
        if e == 0:
          continue
        for r in range(Nstates):
          if r%2 == 0 or r == q+1:
            continue
          phh = Slater(e)
          n = phh.create(r)
          if n == 0:
            continue
          for m in range(len(states_0)):
            if n == states_0[m]:
              H[index[states_0[m]],index[states_0[i]]] -= 0.5 * f
              break
            else:
              continue

  return H

'''


def main():
  # grab delta and g from the command line
  delta     = 1.0    #float(argv[1])
  g         = 0.5    #float(argv[2])
  f         = 0.3    #float(argv[3])
    
  particles = 6
    
  # number of single particle state
  Nstates   = 16
  
  
  
  index = construct_index(particles, Nstates, basis, M_scheme)
  
  
  
  user_data = {
      
      "particles":  particles,
      "Nstates":    Nstates,
      "index":      index,
      "basis":      basis,
      "Slater":     Slater,
      "M_scheme":   M_scheme
  
  }
  #index = construct_index(user_data)
  states, M = basis(particles, Nstates)
  #print (states)
  #print (M)
  states_0, states_2, states_2_, states_4, states_4_ = M_scheme(states, M)
  #print (states_0)
  print (len(states_0))

  

  
  #a = states_binary(states_0)
  #print (a)

  #a = testBit(states_0[0],2)
  #print (a)

  #a = bitCount(states_0[0])
  #print (a)


  H = Hamiltonian(delta, g, f, user_data)
  print (H)
  np.savetxt("matrix.txt", H, fmt = "%.2f")
  
  #H1, H2, H3, H4= Hamiltonian(delta, g, Slater, basis, M_scheme, user_data)
  #print (H1)
  #print ("*"*50)
  #print (H2)
  #print ("*"*50)
  #print (H3)
  #print ("*"*50)
#print (H4)
'''
  tbg = tm.time()
  H = Hamiltonian(delta, g, f, user_data)

  eigva = np.linalg.eigvals(H)
  ted = tm.time()
#print (eigva)

  print (ted - tbg)

  np.savetxt("matrix.txt", H, fmt = "%.2f")
'''

if __name__ == "__main__":
  main()
    

