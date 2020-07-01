import sklearn as sk
from sklearn.metrics import confusion_matrix
import numpy as np
import matplotlib.pyplot as plt

#Overall
TP = float(1138)
FP = float(13175)
FN = float(62)

sens = TP/(TP+FN)
ppv = TP/(TP+FP)

print(sens)
print(ppv)
f1 = 2*sens*ppv/(sens+ppv)
print(f1)

colors = (0,0)

plt.scatter(ppv, sens,c=colors, alpha=0.5)
plt.title('Toy Simulated Data')
plt.xlabel('PPV')
plt.ylabel('Sens')
plt.show()
