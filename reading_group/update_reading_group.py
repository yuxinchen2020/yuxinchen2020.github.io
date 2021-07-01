import sys
import os
from subprocess import Popen


child = Popen(['./submit.sh'],shell=True);
child.wait()
if not child.returncode == 0: 
  print 'some error!! not updating'
  sys.exit(-1);
