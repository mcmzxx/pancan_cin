import os
from itertools import combinations 
import re
os.chdir('/data/zhang/pancan_cin/step5_pathway/code')
import os
import pandas as pd
import numpy as np
i=1
for bt in list(range(1,5,1)):
  r_cmd="Rscript ssgsea{}.R".format(str(bt))
  shell_cmd=r'screen -dmS "session${}" sh -c "{}; exec bash"'.format(i,r_cmd)
  print(shell_cmd)
  os.system(shell_cmd)
  i=i+1
