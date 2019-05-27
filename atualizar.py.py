from shutil import copyfile
import os 

works=os.getcwd().split('\\')[0]+'\\' +os.getcwd().split('\\')[1]+'\\'+os.getcwd().split('\\')[2]+'\\'


copyfile(works+r'OneDrive\Documentos\python proj\jfodata\jfodata.py',works+r'Anaconda3\Lib\jfodata.py')