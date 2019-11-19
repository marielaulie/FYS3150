

import numpy as np
import matplotlib.pyplot as plt
import sys

c = []
E = []

text = sys.argv[1]

f = open(text)
text = sys.argv[1]
F = np.genfromtxt(text, delimiter="")
c = F[:,0]
E = F[:,1]

f.close()

print (E)
print (c)

n, bins, patches = plt.hist(E, bins=100)
plt.title("Probability distribution ", fontweight='bold')
plt.suptitle("MCcycles = 100.000, L=20, ordered initial matrix, T=1.0 [kT/J]")
plt.xlabel("E [J]")
plt.ylabel("P(E)")
plt.savefig('Histo_T1_O_100bins')

plt.show()
