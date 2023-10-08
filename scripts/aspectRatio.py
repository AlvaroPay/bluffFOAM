import numpy as np

def calculateAR(bluff):
  x_ref_length = 1
  y_ref_length = int(bluff[1]+bluff[2])
  AoA = np.radians(int(bluff[3])*1.5)
  
  if y_ref_length <= 10:
      return 0.1*y_ref_length*x_ref_length, (0.1*y_ref_length*x_ref_length)*np.cos(AoA)+np.sin(AoA)
  else:
      return (y_ref_length/5-1)*x_ref_length, ((y_ref_length/5-1)*x_ref_length)*np.cos(AoA)+np.sin(AoA)
