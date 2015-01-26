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


def get_x(object_pointOn):
    return pointOn[0][0]

def get_y(object_pointOn):
    return pointOn[0][1]

def check_horizontal_membership(face, edges):
    for e in edges  :
        face_coord = face.pointOn
        vv = e.getVertices()


        #if (get_x(face.pointOn) )



def populate_global_dict(p, points_0, points_1, pickedFaces, vertical = True):

    faces = p.faces
    edges = p.edges

    for f in faces :

        new_list = list()
        print "The faces are given as "
        print f

        if (vertical == True):
            print "In the correct point"

            face_x = get_x(f.pointOn)
            face_y = get_y(f.pointOn)

            print "The face x and y are ", (face_x, face_y)

            #if ()

        elif (vertical != true):
            print "Needs to be solved later"


    for e in edges :
        vv = e.getVertices()
        print "The edges are given as"
        print p.vertices[vv[0]], p.vertices[vv[1]]


    print "The points are given as "
    print points_0, points_1

    print " <<<<<<<<<<<<<< END >>>>>>>>>>>>>>>>"


    #A face is picked and hence all its edges
    #will either be in one hash of the face or another
    #hash of another face.
    #After that remove the pickedFaces


    #del G_dict_face2edges[pickedFaces]




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

    print "~~~~~~~~~~~~~~~~~~~~~~"
    print "The coordinates are given as :"
    print coord, ratio
    print "~~~~~~~~~~~~~~~~~~~~~~"
    return coord

#def unset_and_delete_sketch(p, s):

def create_new_sketch_transform(p, face_coord, origin):

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

    #populate_global_dict(p, points[0], points[1], pickedFaces)
    s1.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']


def print_faces():
    p = mdb.models['Model-1'].parts['Part-1']
    faces = p.faces

    for f in faces:
        print f








##############################################################
#   This is the start of the main part of the code
#   Above all functions are called
##############################################################

cliCommand("""session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)""")

session = create_session(435, 184)
s = create_sketch()

################################################################
#   it is important to provide the order of the points in the same
#   way to the create_parts since a change in the order will change
#   the order of the faces.
#   This seems to be a very quirky part of the tool but I am not
#   sure how to overcome that just now 22/1/2015
################################################################


point1=(5.0, 20.0)
point2=(50.0, 10.0)
normal=(0.0, 0.0, 1.0)

p = create_part(s, point1, point2)

f, e, d1 = p.faces, p.edges, p.datums
face_coord = create_face_coord(point1, point2)
origin = (0.0, 0.0, 0.0)

s1 = create_new_sketch_transform(p, face_coord, origin)
bifurcate_edges(p, s1)

s2 = create_new_sketch_transform(p, face_coord, origin)
bifurcate_edges(p, s2)

#bifurcate_edges(p, s1, )

print_faces()