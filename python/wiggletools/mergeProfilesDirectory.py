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
import wigglePlots
import os
import shutil

file=sys.argv[1]
if (subprocess.call("sort -k1,1 -k2,2n -k3,3n -m %sx/* > %s" % (file, file), shell=True)):
	print 'Error processing directory %sx' % file
	sys.exit(1)
shutil.rmtree(file + "x")
wigglePlots.make_profiles_matrix(file, file + ".pdf")
