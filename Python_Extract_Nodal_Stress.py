# This simple python script extracts the Von Mises Stress # at selected Node set
#(current Nodes Set in current view in Abaqus ODB viewer) for multiple ODB files
# and multiple frames (each frame was saved as a single file)
# ------------------------------------------------------------------------------
# This file was created with the help of Abaqus Macro ('File- Macro Manager')
# 1. Create a new macro and start recording 
# 2. Perform operations you want in GUI (e.g., save field output, create model, etc.
# 3. Finish recording and go to your current working directory, open the newly created
#    python macro file
# 4. This macro actually creates a user-defined function, you can delete this declaration
#    and just keep the codes in the defined function (pay attention to the identations!)
# 5. Modify the codes as you wish, for example, add for loops, make your own file names..
# 6. Save the file, and go to Abaqus ->File->Run Script, run your code
#
# PS: I am still a rookie Python user, so this code is not optimized. However, I do hope it gives
#     you a quick start of how to combine the Abaqus Macro and basic python coding to achieve
#     your specific goal!
#   
# by Yun Peng, PhD, Massachusetts General Hospital/ Harvard Medical School
# 09/07/2017 

from abaqus import *
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import shutil
import os
for ODBname in os.listdir("E:/FEM_THA/ODB"): # this is where all your ODB files are located
    SUBname = ODBname[3:6]  # My subject ID is saved in the ODB name - this helps me create the file
    print('Current File: '+ODBname)
    ODBnamefull='E:/FEM_THA/ODB/'+ODBname   # Full path to the ODB file, otherwise abaqus tries to find it in default work directory
    o1 = session.openOdb(name=ODBnamefull)  # open the ODB file
    odb = session.odbs[ODBnamefull]         
    file_pre = 'G:/GoogleDrive/1 - Research/1 - Journal/04_FEM_THA_Wear/Data/VonMises/'+SUBname+'/' # prefix of my file to write
    if os.path.exists(file_pre):
        shutil.rmtree(file_pre)
    os.makedirs(file_pre)
    for x in range(0,101):  # frame number = 101
        session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=x)
        session.writeFieldReport(
            fileName=file_pre+'SV_'+str(x).zfill(3)+'.txt',
            append=OFF, sortItem='Node Label', odb=odb, step=1, frame=x,
            outputPosition=NODAL, variable=(('S', INTEGRATION_POINT, ((INVARIANT,'Mises'), )), ))
