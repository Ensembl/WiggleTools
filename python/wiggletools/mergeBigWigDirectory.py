#!/usr/bin/env python
# Copyright 2014 EMBL-EBI                                                                             
# 
# Licensed under the Apache License, Version 2.0 (the "License");                                     
# you may not use this file except in compliance with the License.                                    
# You may obtain a copy of the License at                                                             
# 
# http://www.apache.org/licenses/LICENSE-2.0
#       
# Unless required by applicable law or agreed to in writing, software                                 
# distributed under the License is distributed on an "AS IS" BASIS,                                   
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.                            
# See the License for the specific language governing permissions and                                 
# limitations under the License.

import sys
import subprocess
import glob
import shutil
import os

directory = sys.argv[1]

command = ['catBigWigs', directory] + glob.glob(os.path.join(directory + "x", '*.bw'))

p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

return_code = p.wait()

if return_code != 0:
	out, err = p.communicate()
	print "Failed: %s" % " ".join(command)
	print "Stdout:"
	print out
	print "Stderr:"
	print err
	sys.exit(1)

shutil.rmtree(directory + "x")
