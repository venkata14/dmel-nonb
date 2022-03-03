import numpy as np
import pandas as pd

def removeOutliers(x, outlierConstant=1.5):
  a = np.array(x)
  upper_quartile = np.percentile(a, 75)
  lower_quartile = np.percentile(a, 25)
  IQR = (upper_quartile - lower_quartile) * outlierConstant
  quartileSet = (lower_quartile - IQR, upper_quartile + IQR)
  resultList = []
  for y in a.tolist():
    if y >= quartileSet[0] and y <= quartileSet[1]:
      resultList.append(y)
  return resultList