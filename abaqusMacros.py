# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__


#def Area():
#    import section
#    import regionToolset
#    import displayGroupMdbToolset as dgm
#    import part
#    import material
#    import assembly
#    import step
#    import interaction
#    import load
#    import mesh
#    import optimization
#    import job
#    import sketch
#    import visualization
#    import xyPlot
#    import displayGroupOdbToolset as dgo
#    import math
#    import numpy as np
##    import connectorBehavior
#    
#    #initialising variables
#    nodeLabels = []
#    elementLabels = []
#    # selecting elements and nodes from 3 point planes
#    p = mdb.models['Model-1'].parts['Part-1']
#    # creating plane
#    e2 = p.edges
#    p.DatumPlaneByThreePoints(point1=p.InterestingPoint(edge=e2[1], rule=MIDDLE), 
#        point2=p.InterestingPoint(edge=e2[3], rule=MIDDLE), 
#        point3=p.InterestingPoint(edge=e2[8], rule=MIDDLE)) 
#    c = p.cells
#    d = p.datums
#    print(p.datums.keys())
#    #partitioning part into two sections
#    p.PartitionCellByDatumPlane(datumPlane=d[p.datums.keys()[0]], cells=p.cells[0])
#    
#    size =2
#    #generating mesh
#    p.seedPart(size=size, deviationFactor=0.1, minSizeFactor=0.1)
#    pickedRegions = c.getSequenceFromMask(mask=('[#2 ]', ), )
#    p.generateMesh(regions=pickedRegions)
#    
#    #creating set from face at the plane
#    p.Set('facename',faces=p.faces.findAt(((0,0,10),),)) 
#    a = mdb.models['Model-1'].rootAssembly
#    
#    #creating instance
#    a1 = mdb.models['Model-1'].rootAssembly
#    
#    a1.Instance(name='Part-1-1', part=p, dependent=ON)
#    e21 = a.instances['Part-1-1'].edges
#    a1.DatumCsysByThreePoints(name='Datum csys-1', coordSysType=CARTESIAN, origin=(
#        0.0, 0.0, 0.0), point1=a.instances['Part-1-1'].InterestingPoint(
#        edge=e21[1], rule=MIDDLE), 
#        point2=a.instances['Part-1-1'].InterestingPoint(edge=e21[10], 
#        rule=MIDDLE))
#    p = mdb.models['Model-1'].parts['Part-1']
#    #creating element and nodeset using face set at the plane
#    p.Set('facename',faces=p.faces.findAt(((0,0,10),),)) 
#    
#    #creating instance
#    a1 = mdb.models['Model-1'].rootAssembly
##    a1.DatumCsysByDefault(CARTESIAN)
##    
##    a1.Instance(name='Part-1-1', part=p, dependent=ON)
##    
#    #creating element and nodeset using face set at the plane
#    elementSet = a1.instances['Part-1-1'].sets['facename'].elements
#    nodeSet = a1.instances['Part-1-1'].nodes
##    
#    print(len(elementSet))
#    # Storing node labels at the plane
#    for i in range(0,len(nodeSet)):
#        nodeLabels.insert(0,nodeSet[i].label)
#        
#    for i in range(0,len(elementSet)):
#        elementLabels.insert(0,elementSet[i].label)
##   # creating set from node labels to visualise plane
#    a1.SetFromNodeLabels(name='nodeSet', nodeLabels=(('Part-1-1',nodeLabels),))
#    a1.SetFromElementLabels(name='elementSet', elementLabels=(('Part-1-1',elementLabels),))
#    
#    for i in range(0,len(nodeSet)):
#        node = []
#        node.insert(0,nodeSet[i].label)
#        a1.SetFromNodeLabels(name='node{0}'.format(node[0]), nodeLabels=(('Part-1-1', node),))
#   # initialising to creat sets
#    j = 0
#    k = 0 
#    z = 0
#    for i in range(0,len(elementSet[1].connectivity)):
#        if elementSet[1].connectivity[i] in nodeLabels:
#            k += 1
#    # initialising arrary for storing connected nodes to each element
#    connectedNode = np.empty([len(elementLabels),k])
#    
#    # element set creation and creating connected node array
#    for i in range(0,len(elementSet)):
#        element = []
#        element.insert(0,elementSet[i].label)
#        a1.SetFromElementLabels(name='element{0}'.format(element[0]), elementLabels=(('Part-1-1', element),))
#        for b in elementSet[i].connectivity:
#            if b in nodeLabels:
#                
#                connectedNode[j,z] = b
#                z += 1
#        j +=1
#        z = 0
#    
#    area = []
#    # creating co-ordiante array and calculating area of each element
#    coord = []
#    for j in range(0,len(elementLabels)):
#        coord = []
#        for i,k in enumerate(nodeLabels):
#            if k in connectedNode[j,:]:
#                # appending coordinates for each element surface
#                coord.append(nodeSet[i].coordinates)
#
#    # area calculation 
#    
#        coord = np.array(coord)
#        pq = coord[1,:] - coord[2,:]
#        pr = coord[1,:] - coord[3,:]
#        cros = np.cross(pq,pr)
#        a = np.linalg.norm(cros)  
#            # for tet element
#        if a < 2*size**2 and len(coord) == 4:
#            area.append(a)
#        elif a < 2*size**2 and len(coord) == 3:
#            a = a/2
#            area.append(a)
##            
##    for i in range(0,len(area)):
##        if area[i] > 10:
##            area[i] = 0        
#    print(area)
#    print(connectedNode[0,:])
def fEDC():
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
    import numpy as np
    import math
    import time
    from collections import Counter
    import csv
    import os
    p = mdb.models['Model-1'].parts['Part-1']
    a1 = mdb.models['Model-1'].rootAssembly
    
    # timing script
    t1 = time.time()
    # creating plane
    area2 = {}
    
#    for z in range(-10,15,1):
    lambRecord = []
    s = np.array([0,50,0])
    q = np.array([0,50,6])
    r = np.array([50,50,0])
    
    
    p.DatumPlaneByThreePoints(point1=s, 
        point2=q, 
        point3=r) 
    
    pq = q - s
    pr = r - s
    
    perp = np.cross(pq,pr)
    d = -np.dot(perp,-s)
#    nodeSet = a1.instances['Part-1-1'].nodes
#    elementSet = a1.instances['Part-1-1'].elements
#    
    nodeSet = p.nodes
    elementSet = p.elements
    #select elements
    vec = np.zeros((12,3))
    vertexj = np.zeros((12,3))
    vertexi = np.zeros((12,3))
    
    
    #finding area of intercepts in the general case of any plane interception on the element
    for e in range(0,len(elementSet)):
        lamb = []
        connected = []
        
        connected = elementSet[e].connectivity
        label = elementSet[e].label
        # sort node index to keep track of lines and vertex
        sort = np.sort(connected)

        # finding any vertices and edges on a particular element
        #parallel vertices
        if abs(np.dot(perp,np.array(nodeSet[sort[0]].coordinates)) - d)/1 <= abs(d/10):
            for i in range(0,4):
                vertexj[i,:] = np.array(nodeSet[sort[2*i]].coordinates)
                vertexi[i,:] = np.array(nodeSet[sort[2*i+1]].coordinates)
    
            #parallel vertex
            for i in range(4,8):
                vertexj[i,:] = np.array(nodeSet[sort[i]].coordinates)
                vertexi[i,:] = np.array(nodeSet[sort[i-4]].coordinates)
            #parallel vertex
            vertexj[8,:] = np.array(nodeSet[sort[2]].coordinates)
            vertexj[9,:] = np.array(nodeSet[sort[3]].coordinates)
            vertexj[10,:] =np.array(nodeSet[sort[6]].coordinates)
            vertexj[11,:] =np.array(nodeSet[sort[7]].coordinates)
       
            vertexi[8,:] = np.array(nodeSet[sort[0]].coordinates)
            vertexi[9,:] = np.array(nodeSet[sort[1]].coordinates)
            vertexi[10,:] =np.array( nodeSet[sort[4]].coordinates)
            vertexi[11,:] =np.array( nodeSet[sort[5]].coordinates)
            
    
                # Algorithm for lines and finding lambda (gradient of intercept)
                # refer to paper
            for i in range(0,len(vec)):
                vec[i,:] = vertexi[i,:] - vertexj[i,:]
                if np.dot(perp, vec[i,:]) != 0:
                    lamb.append((d - np.dot(perp, vertexj[i,:]))/(np.dot(perp, vec[i,:]))) 
                    
                else: # if the dot product of gradient of line and normal of plane is 0 add -1 to lambda
                    lamb.append(-1)
            k = 0  
            intercept = np.zeros((6,3))
            count = Counter(lamb)
            # checking for every value of lambda if there is an intercept
            for i in range(0,len(lamb)):
                if lamb[i] <= 1 and lamb[i] > 0:
                    #if there is an intercept find the coordinates using parametric method
                    temp = (np.multiply(vec[i,:],lamb[i]) + vertexj[i,:])
                    intercept[k,:] = temp
                    #if lambda = 1 check if there are over lapping intersections
                    if lamb[i] == 1:
                        for j in range(0,len(intercept)):
                            if j < k and np.array_equal(intercept[j,:], intercept[k,:]):
                                k -= 1
                    k += 1
                # check if intersection is at the nodes and if the itersections lie
                # diagonally directly at the nodes so we have lambda = 1 and 0 
                # the same number of times
                if lamb[i] == 0 and count[0] == count[1]:
                    temp = (np.multiply(vec[i,:],lamb[i]) + vertexj[i,:])
                    intercept[k,:] = temp
                    if lamb[i] == 0:
                        for j in range(0,len(intercept)):
                            if j < k and np.array_equal(intercept[j,:], intercept[k,:]):
                                k -= 1
                    k += 1
                    
                if lamb[i] <=1 and lamb[i] >= 0:
                    lambRecord.append(lamb[i])
    
        
            #find the areas by sorting the coordinates of intersections
            norm = 0
            angles = {}
    
            #for just three points of intersection find area of triangle
            if k   == 3:
                pq = intercept[0,:] - intercept[1,:]
                pr = intercept[0,:] - intercept[2,:]
                cros = np.cross(pq,pr)
                norm += np.linalg.norm(cros)/2
                area2[label] = norm
            #things become harder with more than 3 points
            elif k  > 4 or k == 4:
                #define reference vector
                ref = intercept[0,:] - intercept[1,:]
                #set the reference angle of the reference vector to zero
                angles[1] = 0
                for i in range(2,k):
                    #find all other vectors with respect to the first point
                    #find their angles with respect to the reference angle
                    #then arrange from smallest to largest angle
                    ref2 = intercept[0,:] - intercept[i,:]
                    cros = np.cross(ref,ref2)
                    refNorm = np.linalg.norm(cros)
                    product = np.linalg.norm(ref)*np.linalg.norm(ref2)
                    dot = np.dot(ref,ref2)
                    #if statement added as to test which side the point lies
                    #on reference vector. If the z componenet of the cross
                    #product of the reference vector and the new vector
                    #is negative then the point lies on one side, if positive
                    #point lies on other side
                    if cros[2] < 0:
                        sin = -refNorm/product
                    else:
                        sin = refNorm/product
                    cos = dot/product
                    #between 0 and 90
                    if sin >= 0 and cos >= 0:
                        theta = (180*math.acos(cos))/math.pi
                        angles[i] = theta
                    #between -90 and 0
                    elif sin <= 0 and cos >= 0:
                        theta = (180*math.asin(sin))/math.pi
                        angles[i] = theta
                    #between 90 and 180
                    elif sin >= 0 and cos <= 0:
                        theta = 180-((180*math.asin(sin))/math.pi)
                        angles[i] = theta
                    #between -180 and -90 
                    elif sin <= 0 and cos <= 0:
                        theta = -180 + ((180*math.asin(sin))/math.pi)
                        angles[i] = theta
                sorted_angles = sorted(angles.items(), key=lambda x: x[1])
                sorted_index = np.array(sorted_angles)[:,0]
                for i in range(0,len(sorted_index)-1):
                    pq = intercept[0,:] - intercept[sorted_index[i],:]
                    pr = intercept[0,:] - intercept[sorted_index[i+1],:]
                    cros = np.cross(pq,pr)
                    norm += np.linalg.norm(cros)/2
                area2[label] = norm
           #5.7 seconds when added to only consider elements near the plane and 28.4 seconds without modifycation 
    if lambRecord[0:] == [1] * len(lambRecord) or lambRecord[0:] == [0] * len(lambRecord):
        area2 = {}
        
    print('Area calculated from intersecting planes: {0}'.format(area2))
#    print('Area calculated for special case where plane intersected nodes: {0}'.format(area))
    t2 = time.time()
    print('Run time: {0}'.format(t2-t1))
    
    
    # extracting material orientation
    Name = 'csys-plane'
    #creating coordinate system for plane
    p.DatumCsysByThreePoints( origin = s, name = Name,coordSysType = CARTESIAN,point1 = q, point2 = s + perp)
   #extracting coordinate axis
    planeAxis1 = p.datums[p.features[Name].id].axis1
    planeAxis2 = p.datums[p.features[Name].id].axis2
    planeAxis3 = p.datums[p.features[Name].id].axis3
    
    
    nLayups = len(p.compositeLayups.keys())
    Rx = np.zeros((3,3))
    Ry = np.zeros((3,3))
    Rz = np.zeros((3,3))
    materialAngle = {}
    #for each defined ply layup
    for i in range(0,nLayups):
        key = p.compositeLayups.keys()[i]
        sys = p.compositeLayups[key].orientation.localCsys
        refAxis = str(p.compositeLayups[key].orientation.axis)
        addAngle = p.compositeLayups[key].orientation.angle
        #extract layup axis
        plyaxis1 = p.datums[sys].axis1
        plyaxis2 = p.datums[sys].axis2
        plyaxis3 = p.datums[sys].axis3
        nPlies = len(p.compositeLayups[key].plies)
        #for every ply
        for j in range(0,nPlies):
            plyOrientation = p.compositeLayups[key].plies[j].orientation
            plyAngle = p.compositeLayups[key].plies[j].orientationValue
            plyOrient = str(p.compositeLayups[key].plies[j].axis)
            plyOrientType = str(p.compositeLayups[key].plies[j].orientationType)
            
            #if no material orientation is defined for each individual py then 
            #orientatio will be that of the layup
            if plyOrientation == None:
                direction1 = plyaxis1.direction
                direction2 = plyaxis2.direction
                direction3 = plyaxis3.direction
                totalAngle = addAngle + plyAngle
                
                if plyOrientType == 'ANGLE_45':
                    totalAngle = addAngle + 45
                elif plyOrientType == 'ANGLE_90':
                    totalAngle = addAngle + 90
                elif plyOrientType == 'ANGLLE_NEG45':
                    totalAngle = addAngle - 45
                    
                    
                if totalAngle != 0: 
                    if refAxis == 'AXIS_1':
                #rotation wrt x axis
                        Rx = np.array([[1,0,0],[0,math.cos(totalAngle*math.pi/180),math.sin(totalAngle*math.pi/180)],[0,-math.sin(totalAngle*math.pi/180),math.cos(totalAngle*math.pi/180)]])
                        direction2 = np.dot(Rx,plyaxis2.direction)
                        direction3 = np.dot(Rx,plyaxis3.direction)
                       
                    elif refAxis == 'AXIS_2':
                        Ry = np.array([[math.cos(totalAngle*math.pi/180),0,-math.sin(totalAngle*math.pi/180)],[0,1,0],[math.sin(totalAngle*math.pi/180),0,math.cos(45*math.pi/180)]])
                        direction1 = np.dot(Ry,plyaxis1.direction)
                        direction3 = np.dot(Ry,plyaxis3.direction)
                   
                    elif refAxis == 'AXIS_3':
                        Rz = np.array([[math.cos(totalAngle*math.pi/180),-math.sin(totalAngle*math.pi/180),0],[-math.sin(totalAngle*math.pi/180),math.cos(totalAngle*math.pi/180),0],[0,0,1]])
                        direction1 = np.dot(Rz,plyaxis1.direction)
                        direction2 = np.dot(Rz,plyaxis2.direction)
            else:
                #if there is a defines Csys for each ply then change the axis defined
                plyaxis1 = p.datums[plyOrientation].axis1
                plyaxis2 = p.datums[plyOrientation].axis2
                plyaxis3 = p.datums[plyOrientation].axis3
            
                direction1 = plyaxis1.direction
                direction2 = plyaxis2.direction
                direction3 = plyaxis3.direction
#               
                totalAngle = plyAngle
                if plyOrientType == 'ANGLE_45':
                    totalAngle = 45
                elif plyOrientType == 'ANGLE_90':
                    totalAngle = 90
                elif plyOrientType == 'ANGLLE_NEG45':
                    totalAngle = -45
                    
                if plyOrient == 'AXIS_1':
                #rotation wrt x axis
                    Rx = np.array([[1,0,0],[0,math.cos(totalAngle*math.pi/180),math.sin(totalAngle*math.pi/180)],[0,-math.sin(totalAngle*math.pi/180),math.cos(totalAngle*math.pi/180)]])
                    direction2 = np.dot(Rx,plyaxis2.direction)
                    direction3 = np.dot(Rx,plyaxis3.direction)
                   
                elif plyOrient == 'AXIS_2':
                    Ry = np.array([[math.cos(totalAngle*math.pi/180),0,-math.sin(totalAngle*math.pi/180)],[0,1,0],[math.sin(totalAngle*math.pi/180),0,math.cos(45*math.pi/180)]])
                    direction1 = np.dot(Ry,plyaxis1.direction)
                    direction3 = np.dot(Ry,plyaxis3.direction)
                    
                elif plyOrient == 'AXIS_3':
                    Rz = np.array([[math.cos(totalAngle*math.pi/180),-math.sin(totalAngle*math.pi/180),0],[-math.sin(totalAngle*math.pi/180),math.cos(totalAngle*math.pi/180),0],[0,0,1]])
                    direction1 = np.dot(Rz,plyaxis1.direction)
                    direction2 = np.dot(Rz,plyaxis2.direction)

#            
            #finding angle between direction of plane axis and fibre axis    
            
            dotx = np.dot(direction1,planeAxis1.direction)
            doty = np.dot(direction1,planeAxis2.direction)
            dotz = np.dot(direction1,planeAxis3.direction)
            productx = np.linalg.norm(direction1)*np.linalg.norm(planeAxis1.direction)
            producty = np.linalg.norm(direction1)*np.linalg.norm(planeAxis2.direction)
            productz = np.linalg.norm(direction1)*np.linalg.norm(planeAxis3.direction)
            
            cosx = dotx/productx
            cosy = doty/producty
            cosz = dotz/productz
            
            theta = [180*math.acos(cosx)/math.pi, 180*math.acos(cosy)/math.pi, 180*math.acos(cosz)/math.pi]
            
            for l in range(0,len(theta)):
                if theta[l] > 90:
                    theta[l] = 180 - theta[l]
                    
            #assigning angles to elements
            plyRegion = p.compositeLayups[key].plies[j].region[0]
            plyElements = p.sets[plyRegion].elements
            for k in range(0,len(plyElements)):
                label = plyElements[k].label
                if label in area2.keys():
                    materialAngle[label] = theta[1]
    
    print('Angle output for each element with respect to x, y and z axis of cut plane: {0}'.format(materialAngle))
    
    #calculating toughness of composite in relation to fibre orientation
    #Then storing the toughness of composite in a database if existing 
    #orientation already exist
    materialG = {}
    angles_record = []
    #intralaminar fracture toughness in kJ/m2
    G90 = 0.21
    #interlaminar fracture toughness in kJ/m2
    G0 = 159
    for element in area2.keys():
        theta = materialAngle[element]
        #writing  and reading to database
        if os.path.exists('./composite_toughness_database.csv') == True:
            with open('composite_toughness_database.csv',mode='r') as csv_file:
                csv_reader = csv.DictReader(csv_file)
                for row in csv_reader:
                    if float(row['angle']) not in angles_record:
                        angles_record.append(float(row['angle']))
                    if round(theta,1) == round(float(row['angle']),1):
                        materialG[element] = float(row['toughness'])
                        break
                if round(theta,1) not in np.around(angles_record,decimals = 1):
                    materialG[element] = G0 * math.cos(theta*math.pi/180) + G90 * math.sin(theta*math.pi/180)
                    with open('composite_toughness_database.csv','a') as csv_f:
                        fieldnames = ['angle','toughness']
                        csv_writer = csv.DictWriter(csv_f, fieldnames = fieldnames)
                        csv_writer.writerow({fieldnames[0]:theta,fieldnames[1]:materialG[element]})
        else:
            materialG[element] = G0 * math.cos(theta*math.pi/180) + G90 * math.sin(theta*math.pi/180)
            with open('composite_toughness_database.csv',mode='w') as csv_file:
                fieldnames = ['angle','toughness']
                csv_writer = csv.DictWriter(csv_file, fieldnames = fieldnames)
                csv_writer.writeheader()
                csv_writer.writerow({fieldnames[0]:theta,fieldnames[1]:materialG[element]})   
        #calculate energy dissipation
    totalU = np.dot(materialG.values(),area2.values())
    print('Total energy dissipated from failure: {0} kJ'.format(totalU))
        
    
    
    
    