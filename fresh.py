from abaqus import *
from abaqusConstants import *
from caeModules import *
#from driverUtils import executeOnCaeStartup
from driverUtils import *
from pprint import *
from operator import *

##################################################################
#some global variables and constants which could be used
##################################################################

sixsixseven = 0.6666667
threethreethree = 0.33333

############################## Functions are defined below ###############


def dump_obj(obj, level=0):
    for key, value in obj.__dict__.items():
        if not isinstance(value, types.InstanceType):
             print " " * level + "%s -> %s" % (key, value)
        else:
            dump_obj(value, level + 2)


def create_session(wid, ht):
#    session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=435.0,
#    height=184.0)
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

    #print " Before : The values of the geometry are %s" % (g)
    s.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    #print()
    #print " After : The values of the geometry are  %s " % (str(g))

    # Fixing up the origin
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

    print " The first point is :"
    print point1
    print "The second point is :"
    print point2

    dist = tuple(map(lambda x, y: abs(x - y) , point1, point2))
    print " The value of the coordinates are ------>"
    print (dist)

    ret_coord_x = (point1[0] + sixsixseven*dist[0]) if (point1[0] < point2[0]) \
        else (point1[0] - threethreethree*dist[0]);
    ret_coord_y = (point1[1] + sixsixseven*dist[1]) if (point1[1] < point2[1]) \
        else (point1[1] - threethreethree*dist[1]);

    #coord2 = tuple(map(lambda x, y : y + .6666667*x, dist, point1))
    print " The value of the coordinates are <-----"
    coord2 = (ret_coord_x, ret_coord_y, 0.0)
    print (coord2)

    return coord2

def create_origin_coord(point1, point2):
    coord = tuple(map(lambda x, y: x*0.5 + y*0.5 , point1, point2))
    print "The value of the origin coordinate is "
    print coord
    return (coord[0], coord[1], 0.0)

#def create_line_with_constraint(point1, point2):
#
#    g = s1.geometry
#    v = s1.vertices
#    d = s1.dimensions
#    c = s1.constraints
#
#    s1.Line(point1=point1, point2=point2)
#    s1.VerticalConstraint(entity=g.findAt((-1.801388, 0.0)), addUndoState=False)
#    s1.PerpendicularConstraint(entity1=g.findAt((0.0, 3.0)), entity2=g.findAt((
#        -1.801388, 0.0)), addUndoState=False)
#    s1.CoincidentConstraint(entity1=v.findAt((-1.801388, 3.0)), entity2=g.findAt((
#        0.0, 3.0)), addUndoState=False)
#    s1.CoincidentConstraint(entity1=v.findAt((-1.801388, -3.0)), entity2=g.findAt((
#        0.0, -3.0)), addUndoState=False)

def bifurcate(point1, point2):

    coord = (0.0, 0.0)
    if (point1[1] == point2[1]):
        coord = (0.5*point1[0]+ 0.5*point2[0], point1[1]);

    elif (point1[0] == point2[0]):
        coord = (point1[0], 0.5*point1[1]+ 0.5*point2[1]);

    else:
        print "There is an error in the points, pls check"

    return coord


def bifurcate_edges(p, s1):
    edges = p.edges

    points = []


    for e in edges :
        vertices = e.getVertices()
        vertex_1 = vertices[0]
        vertex_2 = vertices[1]


        v_1 = (p.vertices[vertex_1])
        v_2 = (p.vertices[vertex_2])

        if (v_1.pointOn[0][1] == v_2.pointOn[0][1]):
            coord = bifurcate(v_1.pointOn[0], v_2.pointOn[0])
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

    print "The value of the centre_point_line is"
    print centre_point_of_line

    g = s1.geometry
    print "The geometry is given as "
    print g[2], g[3], g[4], g[5], g[6], g[7]

#    s1.Line(point1=(-1.80138778686523, 3.0), point2=(-1.80138778686523, -3.0))
    s1.VerticalConstraint(entity=g.findAt(centre_point_of_line), addUndoState=False)
#    s1.PerpendicularConstraint(entity1=g.findAt((0.0, points[0][1])), entity2=g.findAt((
#        points[1][0], 0.0)), addUndoState=False)
#    s1.CoincidentConstraint(entity1=v.findAt(points[0]), entity2=g.findAt((
#        0.0, 3.0)), addUndoState=False)
#    s1.CoincidentConstraint(entity1=v.findAt((-1.801388, -3.0)), entity2=g.findAt((
#        0.0, -3.0)), addUndoState=False)
    s1.setPrimaryObject(option=SUPERIMPOSE)
    p = mdb.models['Model-1'].parts['Part-1']
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    p = mdb.models['Model-1'].parts['Part-1']
    f = p.faces
    pickedFaces = f.findAt(((14.0, 8.0, 0.0), ))
    e1, d2 = p.edges, p.datums
    p.PartitionFaceBySketch(faces=pickedFaces, sketch=s1)
    s1.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']

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
#origin = create_origin_coord(point1, point2)
origin = (0.0, 0.0, 0.0)

t = p.MakeSketchTransform(sketchPlane=f.findAt(coordinates= face_coord,
    normal= normal), sketchPlaneSide=SIDE1, origin=origin)

s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
        sheetSize=100, gridSpacing=1.11, transform=t)

#g = s1.geometry
#v = s1.vertices
#d = s1.dimensions
#c = s1.constraints

s1.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Part-1']
p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
print "The vertices are the following ~~~~~~~"


bifurcate_edges(p, s1)

#p = mdb.models['Model-1'].parts['Part-1']



#s1.Line(point1=point1, point2=point2)
#s1.Line(point1=(-1.80138778686523, 3.0), point2=(-1.80138778686523, -3.0))
#s1.VerticalConstraint(entity=g.findAt((-1.801388, 0.0)), addUndoState=False)
#s1.PerpendicularConstraint(entity1=g.findAt((0.0, 3.0)), entity2=g.findAt((
#    -1.801388, 0.0)), addUndoState=False)
#s1.CoincidentConstraint(entity1=v.findAt((-1.801388, 3.0)), entity2=g.findAt((
#    0.0, 3.0)), addUndoState=False)
#s1.CoincidentConstraint(entity1=v.findAt((-1.801388, -3.0)), entity2=g.findAt((
#    0.0, -3.0)), addUndoState=False)


#g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints

###############################################
#create_line_with_constraint()
################################################

#s1.setPrimaryObject(option=SUPERIMPOSE)
#p = mdb.models['Model-1'].parts['Part-1']
#p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)

#p = mdb.models['Model-1'].parts['Part-1']

#print " faces = " % p.faces
#dump_obj(p)
# create_part(s)


#session.viewports['Viewport: 1'].setValues(displayedObject=p)