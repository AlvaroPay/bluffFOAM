def calculateAR(bluff):
  x_ref_length = 1
  y_ref_length = int(bluff[1]+bluff[2])
  if y_ref_length <= 5:
      return 0.2*y_ref_length*x_ref_length
  else:
      return (y_ref_length/2 - 1.5)*x_ref_length