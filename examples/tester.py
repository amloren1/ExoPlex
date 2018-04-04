import sys
import os

path = os.path.dirname(os.path.realpath(__file__))
path = path[0:-9]
print path
sys.path.append(path)
print sys.path


import exoplex.run
