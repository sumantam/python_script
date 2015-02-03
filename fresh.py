from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import *
from pprint import *
from operator import *

##################################################################
#some global variables and constants which could be used
##################################################################

sixsixseven = 0.6666667
threethreethree = 0.33333
G_dict_face2edges = dict()
############################## Functions are defined below ###############


def dump_obj(obj, level=0):
    for key, value in obj.__dict__.items():
        if not isinstance(value, types.InstanceType):
             print " " * level + "%s -> %s" % (key, value)
        else:
            dump_obj(value, level + 2)


def create_session(wid, ht):
    session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=wid,
    height=ht)
    session.viewports['Viewport: 1'].makeCurrent()
    session.viewports['Viewport: 1'].maximize()
    executeOnCaeStartup()
    session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)

    return session

def get_parts(model):
    p = model.parts['Part-1']
    return p

def get_model():
    m = mdb.models['Model-1']
    return m


def create_sketch():
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
        sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.sketchOptions.setValues(viewStyle=AXISYM)
    s.setPrimaryObject(option=STANDALONE)
    s.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    s.FixedConstraint(entity=g.findAt((0.0, 0.0)))

    return s

def create_part(s, point1, point2):

    # Creating a rectangle
    s.rectangle(point1=point1, point2=point2)
    p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=AXISYMMETRIC,
        type=DEFORMABLE_BODY)
    mm = get_model()
    p = get_parts(mm)

    p.BaseShell(sketch=s)

    ############################################################
    # NOTE NOTE NOTE :
    # This is not a good place to unset the primary object.
    # A better place has to be found where the unsetting should
    # be done
    ############################################################
    s.unsetPrimaryObject()

    mm = get_model()
    p = get_parts(mm)
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']
    return p

##############################################################################
# Could need to relook at this function later for a better way of doing things
##############################################################################

def create_face_coord(point1, point2):

    ret_coord_x = 0;
    ret_coord_y = 0;

    dist = tuple(map(lambda x, y: abs(x - y) , point1, point2))

    ret_coord_x = (point1[0] + sixsixseven*dist[0]) if (point1[0] < point2[0]) \
        else (point1[0] - threethreethree*dist[0]);
    ret_coord_y = (point1[1] + sixsixseven*dist[1]) if (point1[1] < point2[1]) \
        else (point1[1] - threethreethree*dist[1]);

    coord2 = (ret_coord_x, ret_coord_y, 0.0)

    return coord2

def create_origin_coord(point1, point2):
    coord = tuple(map(lambda x, y: x*0.5 + y*0.5 , point1, point2))
    return (coord[0], coord[1], 0.0)


def get_x(object):
    return object.pointOn[0][0]

def get_y(object):
    return object.pointOn[0][1]

def check_horizontal_membership(root_assembly, edge):

    vv = edge.getVertices()
    point1 = root_assembly.vertices[vv[0]]
    point2 = root_assembly.vertices[vv[1]]

    if ((point1.pointOn[0][1]) == (point2.pointOn[0][1])):
        ee = edge
    else :
        ee = None

    return ee



def bifurcate(point1, point2, ratio = 0.5):

    coord = (0.0, 0.0)
    pp_1 = point1
    pp_2 = point2

    if (point1[0] > point2[0]):
        pp_1 = point2;
        pp_2 = point1;
    elif ((point1[0] == point2[0]) and (point1[1] > point2[1])):
        pp_1 = point2;
        pp_2 = point1;


    if (point1[1] == point2[1]):
        coord = (ratio*pp_1[0]+ (1-ratio)*pp_2[0], pp_1[1]);

    elif (point1[0] == point2[0]):
        coord = (pp_1[0], ratio*pp_1[1]+ (1-ratio)*pp_2[1]);

    else:
        print "There is an error in the points, pls check"

    return coord

#def unset_and_delete_sketch(p, s):

def create_new_sketch_transform(p, face_coord, origin):

    f = p.faces
    t = p.MakeSketchTransform(sketchPlane=f.findAt(coordinates= face_coord,
    normal= normal), sketchPlaneSide=SIDE1, origin=origin)

    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
        sheetSize=100, gridSpacing=1.11, transform=t)

    s1.setPrimaryObject(option=SUPERIMPOSE)
    p = mdb.models['Model-1'].parts['Part-1']
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)

    return s1

def bifurcate_edges(p, s1, face_number = 0):

    edges = p.faces[face_number].getEdges()
    face = p.faces[face_number]
    face_edges = face.getEdges()
    points = []

    for e in edges :
        ##############################################
        #Check if the edge is in the list of vertices
        ###############################################

        if not(e in face_edges):
            continue

        vertices = p.edges[e].getVertices()
        vertex_1 = vertices[0]
        vertex_2 = vertices[1]


        v_1 = (p.vertices[vertex_1])
        v_2 = (p.vertices[vertex_2])

        if (v_1.pointOn[0][1] == v_2.pointOn[0][1]):
            coord = bifurcate(v_1.pointOn[0], v_2.pointOn[0], 0.3)
            print " The coord, vertex 1 and vertex 2  is "
            print (coord, v_1, v_2)
            points.append(coord)

    if (len(points)==2):
        print " U r in the correct place"
        print points
    else:
        print " U r in the wrong place"

    s1.Line(point1=points[0], point2=points[1])
    centre_point_of_line = bifurcate(points[0], points[1])

    g = s1.geometry

    s1.VerticalConstraint(entity=g.findAt(centre_point_of_line), addUndoState=False)
    s1.setPrimaryObject(option=SUPERIMPOSE)
    p = mdb.models['Model-1'].parts['Part-1']
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    p = mdb.models['Model-1'].parts['Part-1']
    f = p.faces


    print "The faces available are the following"
    print f[face_number].pointOn
    #pickedFaces = f.findAt(((14.0, 8.0, 0.0), ))
    print " The face coordinate is given as"
    print f[face_number].pointOn

    pickedFaces = f.findAt((f[face_number].pointOn[0], ))
    e1, d2 = p.edges, p.datums
    p.PartitionFaceBySketch(faces=pickedFaces, sketch=s1)

    s1.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']



def create_material_prop(prop):

    mat_list = mdb.models['Model-1'].materials
    mat_list_len = len(mat_list)

    mat_name = 'Mat-' + str(mat_list_len+1)

    #print "The name of the material list is ", mat_name
    mdb.models['Model-1'].Material(name=mat_name)
    mdb.models['Model-1'].materials[mat_name].Elastic(table=(prop, ))

    mm = mdb.models['Model-1'].materials[mat_name]
    return mm

def section_create(mat):

    section_list = mdb.models['Model-1'].sections.keys()
    section_list_len = len(section_list)

    mat_name = mat.name
    sect_name = 'Section-' + str(section_list_len+1)
    ss = mdb.models['Model-1'].HomogeneousSolidSection(name=sect_name,
        material=mat_name, thickness=None)
    return ss

def region_create(section, face):

    p = mdb.models['Model-1'].parts['Part-1']
    region_list = p.sets.keys()
    region_list_len = len(region_list)
    region_name = 'Set-' + str(region_list_len + 1)

    #This is an ugly way but needs to look up later
    ###########################################################
    f = p.faces
    faces = f.findAt((face.pointOn[0], ))


    section_name = section.name

    region = p.Set(faces=faces, name=region_name)
    p.SectionAssignment(region=region, sectionName=section_name, offset=0.0,
        offsetType=MIDDLE_SURFACE, offsetField='',thicknessAssignment=FROM_SECTION)


def create_assembly():
    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    a = mdb.models['Model-1'].rootAssembly
    a.DatumCsysByThreePoints(coordSysType=CYLINDRICAL, origin=(0.0, 0.0, 0.0),
        point1=(1.0, 0.0, 0.0), point2=(0.0, 0.0, -1.0))
    p = mdb.models['Model-1'].parts['Part-1']
    a.Instance(name='Part-1-1', part=p, dependent=ON)


def create_steps():

    list_step = ()
    list_step = session.viewports['Viewport: 1'].assemblyDisplay.step
    len_list_step = len(list_step)
    step_name = 'Set-' + str(len_list_step + 1)

    session.viewports['Viewport: 1'].assemblyDisplay.setValues(
        adaptiveMeshConstraints=ON)
    mdb.models['Model-1'].StaticStep(name=step_name, previous='Initial')
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step=step_name)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON,
        predefinedFields=ON, connectors=ON, adaptiveMeshConstraints=OFF)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step=step_name)

def get_face_boundary(facenumber, top=TRUE):

    hor_list = list()
    a = mdb.models['Model-1'].rootAssembly
    edges = a.instances['Part-1-1'].edges
    e_index_list = a.instances['Part-1-1'].faces[facenumber].getEdges()

    for i in e_index_list :
        e = edges[i]
        hh = check_horizontal_membership(a.instances['Part-1-1'], e)

        if hh is None :
            continue
        else:
            hor_list.append(hh)


    if (top == TRUE):
        side = hor_list[0];
        for m in hor_list:
            print "The value of m is ", m.pointOn[0]
            if (get_y(m) < get_y(side)):
                side = m

    else :
        side = hor_list[0];
        for m in hor_list:
            if (get_y(m) > get_y(side)):
                side = m

    #print "The vertices in the hor_list  ", hor_list

    return side



def print_faces():
    p = mdb.models['Model-1'].parts['Part-1']
    faces = p.faces

    for f in faces:
        print f


def clear_and_delete_db():
    print " The database needs to be cleared and deleted"


def create_and_run_jobs(
                        session,
                        point1,
                        point2,
                        normal,
                        boundary_face_number,
                        fixed_face_number,
                        mat_prop_1,
                        mat_prop_2
                        ):


    s = create_sketch()
    p = create_part(s, point1, point2)

    f, e, d1 = p.faces, p.edges, p.datums
    face_coord = create_face_coord(point1, point2)
    origin = (0.0, 0.0, 0.0)

    s1 = create_new_sketch_transform(p, face_coord, origin)
    bifurcate_edges(p, s1)

    s2 = create_new_sketch_transform(p, face_coord, origin)
    bifurcate_edges(p, s2)

    mm1 = create_material_prop(mat_prop_1)
    mm2 = create_material_prop(mat_prop_2)

    ####################################################################
    # create same number of sections as the material number
    ####################################################################

    ss1 = section_create(mm1)
    ss2 = section_create(mm2)

    region_create(ss1, p.faces[0])
    region_create(ss2, p.faces[1])
    region_create(ss1, p.faces[2])

    create_assembly()
    create_steps()


    face_bound_3 = get_face_boundary(boundary_face_number, TRUE)
    face_bound_1 = get_face_boundary(fixed_face_number, FALSE)

    a = mdb.models['Model-1'].rootAssembly

    edge_List = a.instances['Part-1-1'].edges

    region = a.Surface(side1Edges=edge_List.findAt((face_bound_3.pointOn[0],)), name='Surf-1')

    ll = len(mdb.models['Model-1'].steps.keys())
    name = mdb.models['Model-1'].steps.keys()[-1]

    ############## Create Boundary Conditions #########################################
    mdb.models['Model-1'].Pressure(name='Load-1', createStepName=name,
            region=region, distributionType=UNIFORM, field='', magnitude=-5000.0,
            amplitude=UNSET)

    region = a.Set(edges=edge_List.findAt((face_bound_1.pointOn[0],)), name='Set-1')

    mdb.models['Model-1'].DisplacementBC(name='BC-1', createStepName=name,
        region=region, u1=0.0, u2=0.0, ur3=0.0, amplitude=UNSET, fixed=OFF,
        distributionType=UNIFORM, fieldName='', localCsys=None)

    ###################################  Meshing ##########################

    p = mdb.models['Model-1'].parts['Part-1']
    p.seedPart(size=2.2, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    a1 = mdb.models['Model-1'].rootAssembly
    a1.regenerate()

    ########################### Job Submission ############################

    job_list = mdb.jobs.keys()
    job_list_len = len(job_list)

    job_name = 'Job-' + str(job_list_len+1)

    mdb.Job(name=job_name, model='Model-1', description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1,
        numGPUs=0)
    mdb.jobs[job_name].submit(consistencyChecking=OFF)


    return mdb.jobs[job_name]


################################################################
#   it is important to provide the order of the points in the same
#   way to the create_parts since a change in the order will change
#   the order of the faces.
#   This seems to be a very quirky part of the tool but I am not
#   sure how to overcome that just now 22/1/2015
################################################################

############## Beginning of the User Input ###########################

point1=(5.0, 20.0)
point2=(50.0, 10.0)
normal=(0.0, 0.0, 1.0)
boundary_face_number = 2 # The index of the faces are starting from 0, so essentially the 3rd face
fixed_face_number = 0
mat_prop_1 = (1200, 0.3)    # Youngs modulus , poisson ratio
mat_prop_2 = (1500, 0.15)   # Youngs modulus , poisson ratio

############# End of the User Input ###############################


##############################################################
#   This is the start of the main part of the code
#   Above all functions are called
##############################################################

cliCommand("""session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)""")

session = create_session(435, 184)

job = create_and_run_jobs(session,
                    point1,point2,
                    normal,boundary_face_number,
                    fixed_face_number,mat_prop_1,
                    mat_prop_2)

job.waitForCompletion()

output_path = 'C:/Temp/' + job.name + '.odb'

print "The output path of the job is given as ", output_path

o1 = session.openOdb(name=output_path)
session.viewports['Viewport: 1'].setValues(displayedObject=o1)

def tt():
    session.viewports['Viewport: 1'].setValues(displayedObject=odb)
    session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
        CONTOURS_ON_DEF, ))

    session.Path(name='Path-1', type=POINT_LIST, expression=((5.0, 15.0, 0.0), (
        50.0, 15.0, 0.0)))
    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
        variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT,
        'S11'))
    xyp = session.XYPlot('XYPlot-1')
    chartName = xyp.charts.keys()[0]
    chart = xyp.charts[chartName]
    pth = session.paths['Path-1']
    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=True,
        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10,
        projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE)
    c1 = session.Curve(xyData=xy1)
    chart.setValues(curvesToPlot=(c1, ), )
    session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
    pth = session.paths['Path-1']
    session.XYDataFromPath(name='XYData-1', path=pth, includeIntersections=True,
        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10,
        projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE)
    x0 = session.xyDataObjects['XYData-1']
    session.writeXYReport(fileName='abaqus.rpt', xyData=(x0, ))






