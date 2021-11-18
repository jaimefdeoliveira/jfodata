# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 16:06:59 2020

@author: jaime
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May 27 07:52:07 2020

@author: jaime
"""

import pyvisa
import time

rm = pyvisa.ResourceManager()
#it.close()
it=rm.open_resource(rm.list_resources()[0])

it.write_termination='\n'
time.sleep(1)
it.read_termination = '\n'
time.sleep(1)
it.timeout=2
time.sleep(1)
print(it.query('*IDN?',delay=1))

it.close()