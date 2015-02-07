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
    m = get_model()
    s = m.ConstrainedSketch(name='__profile__',
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
    mm = get_model()
    p = mm.Part(name='Part-1', dimensionality=AXISYMMETRIC,
        type=DEFORMABLE_BODY)

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

def check_vertical_membership(root_assembly, edge):

    vv = edge.getVertices()
    point1 = root_assembly.vertices[vv[0]]
    point2 = root_assembly.vertices[vv[1]]

    if ((point1.pointOn[0][0]) == (point2.pointOn[0][0])):
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

    s1 = get_model().ConstrainedSketch(name='__profile__',
        sheetSize=100, gridSpacing=1.11, transform=t)

    s1.setPrimaryObject(option=SUPERIMPOSE)
    p = mdb.models['Model-1'].parts['Part-1']
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)

    return s1

def get_face_edges(face_number):

    a = get_model().rootAssembly
    p = get_model().parts['Part-1']
    face = a.faces[face_number]
    face_edges = face.getEdges()

    ver = a.vertices()

    for e in face_edges :
        ##############################################
        #Check if the edge is in the list of vertices
        ###############################################


        vertices = a.edges[e].getVertices()
        vertex_1 = vertices[0]
        vertex_2 = vertices[1]

        print "The verticesw are given as ", vertex_1, vertex_2

###########################################################################
# This is a bad algorithm and I am sure there are better ways to do it
# But I do not have a handle to the data structure and so it is painful
###########################################################################

#def get_face_number(point1, point2):
#
#    face_list = p.faces
#
#    for face in face_list:
#        edges = face.getEdges()
#
#        for e in edges :



#############################################################
# Have to know how to overload a function in python
#############################################################


def split_edges(p, s1, point1, point2, face_number=0):

    s1.Line(point1, point2)
    centre_point_of_line = bifurcate(point1, point2)

    g = s1.geometry

    s1.VerticalConstraint(entity=g.findAt(centre_point_of_line), addUndoState=False)
    s1.setPrimaryObject(option=SUPERIMPOSE)
    p = get_model().parts['Part-1']
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    p = get_model().parts['Part-1']
    f = p.faces

    pickedFaces = f.findAt((point1[0], point1[1], 0),)
    p.PartitionFaceBySketch(faces=pickedFaces, sketch=s1)



def bifurcate_edges(p, s1, ratio = 0.5, face_number = 0):

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
            coord = bifurcate(v_1.pointOn[0], v_2.pointOn[0], ratio)
            #print " The coord, vertex 1 and vertex 2  is "
            #print (coord, v_1, v_2)
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
    p = get_model().parts['Part-1']
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    p = get_model().parts['Part-1']
    f = p.faces


    print "The faces available are the following"
    print f[face_number].pointOn
    #pickedFaces = f.findAt(((14.0, 8.0, 0.0), ))
    print " The face coordinate is given as"
    print f[face_number].pointOn

    pickedFaces = f.findAt((f[face_number].pointOn[0], ))
    e1, d2 = p.edges, p.datums
    p.PartitionFaceBySketch(faces=pickedFaces, sketch=s1)



def delete_sketch(sketch):
    sketch.unsetPrimaryObject()
    del get_model().sketches['__profile__']


def create_material_prop(prop):

    mat_list = mdb.models['Model-1'].materials
    mat_list_len = len(mat_list)

    mat_name = 'Mat-' + str(mat_list_len+1)

    #print "The name of the material list is ", mat_name
    get_model().Material(name=mat_name)
    get_model().materials[mat_name].Elastic(table=(prop, ))

    mm = mdb.models['Model-1'].materials[mat_name]
    return mm

def section_create(mat):

    section_list = get_model().sections.keys()
    section_list_len = len(section_list)

    mat_name = mat.name
    sect_name = 'Section-' + str(section_list_len+1)
    ss = get_model().HomogeneousSolidSection(name=sect_name,
        material=mat_name, thickness=None)
    return ss

def region_create(section, face):

    p = get_model().parts['Part-1']
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
    a = get_model().rootAssembly
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

def get_key(item):
    return item.pointOn[0][0]

def get_v_face_boundary(facenumber, left = TRUE):

    ver_list = list()
    a = get_model().rootAssembly
    edges = a.instances['Part-1-1'].edges

    face_list = a.instances['Part-1-1'].faces

    sorted_list = sorted(face_list, key=get_key)

    e_index_list = sorted_list[facenumber].getEdges()

    for i in e_index_list :
        e = edges[i]
        hh = check_vertical_membership(a.instances['Part-1-1'], e)

        if hh is None :
            continue
        else:
            ver_list.append(hh)

    if (left == TRUE):
        side = ver_list[0];
        for m in ver_list:
           # print "The value of m is ", m.pointOn[0]
            if (get_x(m) < get_x(side)):
                side = m

    else :
        side = ver_list[0];
        for m in ver_list:
            if (get_x(m) > get_x(side)):
                side = m

    #print "The vertices in the hor_list  ", hor_list

    return side



def get_face_boundary(facenumber, top=TRUE):

    hor_list = list()
    a = get_model().rootAssembly
    edges = a.instances['Part-1-1'].edges

    face_list = a.instances['Part-1-1'].faces

    sorted_list = sorted(face_list, key=get_key)

    e_index_list = sorted_list[facenumber].getEdges()

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
            # print "The value of m is ", m.pointOn[0]
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
    p = get_model().parts['Part-1']
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
                        mat_prop_2,
                        ratio_one,
                        ratio_second,
                        Load,
                        pipe_pressure,
                        mesh_size = 0.02,
                        option = FALSE
                        ):


    s = create_sketch()
    p = create_part(s, point1, point2)

    f, e, d1 = p.faces, p.edges, p.datums
    face_coord = create_face_coord(point1, point2)
    origin = (0.0, 0.0, 0.0)

    ##################################################################################
    # Since most users would not want to think of the complications of ratios
    # we would want to put that independent and have the coordinates at which
    # the partition being done given at one shot sometimes.
    # In case not,  an option is being provided to do the old way
    ##################################################################################


    if (option == TRUE):
        s1 = create_new_sketch_transform(p, face_coord, origin)
        delete_sketch(s1)
        bifurcate_edges(p, s1, ratio_one)

        s2 = create_new_sketch_transform(p, face_coord, origin)
        delete_sketch(s2)
        bifurcate_edges(p, s2, ratio_second)

    else :

        s1 = create_new_sketch_transform(p, face_coord, origin)

        div_point_one = point1
        div_point_one_prime = (point2[0], point1[1], 0)

        div_point_two = point2
        div_point_two_prime = (point1[0], point2[1], 0)

        pp1 = bifurcate(div_point_one ,  div_point_one_prime, ratio_one)
        pp2 = bifurcate(div_point_two, div_point_two_prime, ratio_one)

        split_edges(p, s1, pp1, pp2)

        pp1 = bifurcate(div_point_one ,  div_point_one_prime, ratio_second)
        pp2 = bifurcate(div_point_two, div_point_two_prime, ratio_second)

        split_edges(p, s1, pp1, pp2)

        delete_sketch(s1)


    mm1 = create_material_prop(mat_prop_1)
    mm2 = create_material_prop(mat_prop_2)

    ####################################################################
    # create same number of sections as the material number
    ####################################################################

    ss1 = section_create(mm1)
    ss2 = section_create(mm2)

    face_list = p.faces
    sorted_list = sorted(face_list, key=get_key)

    region_create(ss1, sorted_list[0])
    region_create(ss2, sorted_list[1])
    region_create(ss1, sorted_list[2])

    create_assembly()
    create_steps()

    face_bound_0 = get_v_face_boundary(fixed_face_number, TRUE)
    face_bound_3 = get_face_boundary(boundary_face_number, TRUE)
    face_bound_1 = get_face_boundary(fixed_face_number, FALSE)

    mm = get_model()
    a = mm.rootAssembly

    edge_List = a.instances['Part-1-1'].edges

    region = a.Surface(side1Edges=edge_List.findAt((face_bound_3.pointOn[0],)), name='Surf-1')

    ll = len(mm.steps.keys())
    name = mm.steps.keys()[-1]

    ############## Create Boundary Conditions #########################################
    mm.Pressure(name='Load-1', createStepName=name,
            region=region, distributionType=UNIFORM, field='', magnitude=Load,
            amplitude=UNSET)

    #print " the value of the vertical face is ", face_bound_0.pointOn[0]
    region = a.Surface(side1Edges=edge_List.findAt((face_bound_0.pointOn[0],)), name='Surf-2')

    mm.Pressure(name='Load-2', createStepName=name,
            region=region, distributionType=UNIFORM, field='', magnitude=pipe_pressure,
            amplitude=UNSET)


    region = a.Set(edges=edge_List.findAt((face_bound_1.pointOn[0],)), name='Set-1')


    mm.DisplacementBC(name='BC-1', createStepName=name,
        region=region, u1=0.0, u2=0.0, ur3=0.0, amplitude=UNSET, fixed=OFF,
        distributionType=UNIFORM, fieldName='', localCsys=None)

    ###################################  Meshing ##########################

    p = mm.parts['Part-1']
    p.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    a1 = mm.rootAssembly
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

    mdb.jobs[job_name].waitForCompletion()

    return mdb.jobs[job_name]

def generate_path(job,
                  point1,
                  point2,
                  session
                ):

    path_list = session.viewports.keys()
    path_list_len = len(path_list)
    path_name = 'Path-' + str(path_list_len)

    div_point_one = point1
    div_point_one_prime = (point2[0], point1[1], 0)

    div_point_two = point2
    div_point_two_prime = (point1[0], point2[1], 0)

    pp1 = bifurcate(div_point_one ,  div_point_one_prime, ratio_one)
    pp2 = bifurcate(div_point_two, div_point_two_prime, ratio_one)

    #print "The generated points for the paths are", pp1, pp2

    pt1 = (pp1[0], pp1[1], 0)
    pt2 = (pp2[0], pp2[1], 0)
    path = session.Path(name=path_name, type=POINT_LIST, expression=(pt1, pt2))

    return path

def normalize_with_stress( list_ll, applied_stress ):

    #print "The list is ",list_ll
    F = list(map (lambda x: (float(x[0])/applied_stress, x[1]), list_ll))
    return F

def print_output(
                job,
                point1,
                point2,
                stress_plot_way,
                session,
                Load=1
                ):

    output_path = 'C:/Temp/' + job.name + '.odb'

    viewport_list = session.viewports.keys()
    viewport_list_len = len(viewport_list)
    job_name = 'Viewport: ' + str(viewport_list_len)


    o1 = session.openOdb(name=output_path)
    session.viewports[job_name].setValues(displayedObject=o1)

    pth = generate_path(job, point1, point2, session)

    xyplot_len = len(session.xyPlots)
    plot_name = 'XYPlot-'+ str(xyplot_len+1)

    path_name = pth.name
    xyp = session.XYPlot(plot_name)

    var = ('S', INTEGRATION_POINT, ((COMPONENT, 'S12'),),)

    xydata_list_s12 = session.XYDataFromPath(
                                                name='xyName',
                                                path=pth,
                                                includeIntersections=True,
                                                pathStyle=PATH_POINTS,
                                                numIntervals=10,
                                                labelType=TRUE_DISTANCE_Y,
                                                shape=UNDEFORMED,
                                                variable=var
                                            )

    xx1 = xydata_list_s12/Load

    session.writeXYReport(
            'C:/Temp/sumo.txt',
            xydata_list_s12
    )

    session.writeXYReport(
            'C:/Temp/sumo_1.txt',
            xx1
    )

    session.viewports[job_name].odbDisplay.setPrimaryVariable(
        variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT,
        'S12')
    )

    chartName = xyp.charts.keys()[0]
    chart = xyp.charts[chartName]
    pth = session.paths[path_name]

    c1 = session.Curve(xyData=xx1)
    chart.setValues(curvesToPlot=(c1, ), )

    session.viewports[job_name].setValues(displayedObject=xyp)

    return xydata_list_s12

################################################################
#   it is important to provide the order of the points in the same
#   way to the create_parts since a change in the order will change
#   the order of the faces.
#   This seems to be a very quirky part of the tool but I am not
#   sure how to overcome that just now 22/1/2015
#   problems in choosing material radius and inner dia of the pipe
# Make all of the material variables as a data structure and pass
################################################################

############## Beginning of the User Input ###########################

print " ###############################################################################################"
print "          WELCOME    TO  THE  SCRIPT                                                            "
print " Give the dimensions in metres for the geometries and use SI units for other physical quantities"
print "################################################################################################"


point1=(0.5, 20.0)
point2=(2.0, 10.0)
normal=(0.0, 0.0, 1.0)
boundary_face_number = 2 # The index of the faces are starting from 0, so essentially the 3rd face
fixed_face_number = 0
mat_prop_1 = (120000, 0.3)    # Youngs modulus , poisson ratio
mat_prop_2 = (1500, 0.15)   # Youngs modulus , poisson ratio

pipe_pressure = 3

ratio_one = 0.5            # These two need some more explanation or writing
ratio_two = 0.498            # ( This needs to be understood as something happening as (1-0.8)*(1 - 0.3)

stress_plot_way = ratio_one*0.5 + ratio_two*0.5
Load = -50




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
                    mat_prop_2,
                    ratio_one,
                    ratio_two,
                    Load,
                    pipe_pressure)


aa = print_output(job, point1, point2, stress_plot_way, session, Load)




