
import random 
import numpy as np

def listToStringWithoutBrackets(list1):
    return str(list1).replace('[','').replace(']','')

num_elements = 10**4

#elemsx = np.linspace(0.0, 129.99999999, num_elements)
#elemsy = np.array([0.1 for i in range(0, num_elements)])
#elemsz = np.array([0.1 for i in range(0, num_elements)])
#elements = (elemsx, elemsy, elemsz)
elements = np.random.uniform(low=0.0, high=129.9999999999, size=(num_elements, 3))


print elements

file_name = '/scratch/chasonnscr/known/DM_knownResearch.dat'
target_file = open(file_name, 'w')

np.savetxt(target_file, elements, fmt='%lf')


#new_elements = listToStringWithoutBrackets(elements)
#new_elements = str(elements).replace('[','').replace(']','')

#print new_elements
#print >> target_file, new_elements
#for item in elements:
#  print >> target_file, str(item).replace('[','').replace(']','')

target_file.close()

#print(elements)


