import os
from glob import glob

# Get current path of this file
curpath = os.path.dirname(os.path.abspath(__file__))
# Go up three directories
curpath = os.path.abspath(os.path.join(curpath,os.pardir))
curpath = os.path.abspath(os.path.join(curpath,os.pardir))
curpath = os.path.abspath(os.path.join(curpath,os.pardir))
# Go to jobs directory
PATH_TO_JOBS = os.path.abspath(os.path.join(curpath,'jobs'))

assert os.path.isdir(PATH_TO_JOBS)

jobfiles = glob('%s/*.sh' % PATH_TO_JOBS)
JOBS = dict( ('.'.join(os.path.basename(f).split('.')[:-1]), f) for f in jobfiles )
