from lantz import Feat
from lantz.messagebased import MessageBasedDriver
from time import sleep

class arduino(MessageBasedDriver):
    
    DEFAULTS = {'COMMON': {'write_termination': '\n',
                           'read_termination': '\n'},
                'ASRL':   {'write_termination': '\n',
                         'read_termination': '\n',
                         'baud_rate': 9600}
                }

     
    @Feat()
    def idn(self):
        

        return self.query('*IDN?')

    

    

            
            
if __name__ == '__main__':
    
    from time import sleep
    from datetime import datetime
    import pyvisa
    import lantz.messagebased
    
    lantz.messagebased._resource_manager = pyvisa.ResourceManager('@py')
    
    inst = arduino('ASRLCOM3::INSTR') #Digite seu endereço no lugar do COM10 
    sleep(1)
    inst.initialize()
    inst.timeouttimeout = 1
    sleep(1)
    inst.query_delay=2
    sleep(1)

    print(inst.idn)
    inst.finalize()