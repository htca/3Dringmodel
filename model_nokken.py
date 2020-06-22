import math
import re
from os import getcwd

############################################################################################
## Class definitions
############################################################################################

class tunnel_element:
    def __init__(self, name):
        self.name = name
        self.ring = int(re.split(" |_", name)[-2])
        self.elnr = int(re.split(" |_", name)[-1])
        self.faces = faces(name)
        self.BBox = boundingBox(name)
        self.center = [(self.BBox[0]+self.BBox[1])/2, (self.BBox[2]+self.BBox[3])/2, (self.BBox[4]+self.BBox[5])/2]
        angle = math.atan2(self.center[2], self.center[1])
        if angle < 0: angle += 2 * math.pi
        self.angle = angle
        self.outer_face = []
        self.inner_face = []
        self.long_face = []
        self.trans_face = []

############################################################################################
## Function definitions
############################################################################################

def distance(coor1, coor2):
    return (math.sqrt((coor1[0]-coor2[0])**2 + (coor1[1]-coor2[1])**2 + (coor1[2]-coor2[2])**2))

def get_rings(elements, ringnumber):
    sel_elements=[]
    for element in elements:
        if ringnumber == element.ring:
            sel_elements.append(element)
    return(sel_elements)

############################################################################################
## Parameter input
############################################################################################
# Switches
create_ascending_bedding = True
variable_outside_loading = True
create_dummy_interface_long = True
create_dummy_interface_trans = True
create_analysis = True
create_mesh = True
run_analysis_linsta = False
run_analysis_nlsta = False
mc_int_trans = True
mc_int_long = mc_int_trans

n_rings = 5
d_inner = 11
t_ring = 0.4
d_outer = d_inner + 2 * t_ring
n_segment = 7
l_ring = 2
b_nok = 0.3

Es = 4000E6
H_water = 20  # waterlevel above centerline

# MC interface vars
mu_conc_to_conc = 0.4
MC_phi_fric = math.atan(mu_conc_to_conc)
MC_phi_dil = MC_phi_fric
MC_cohesion = 1
MC_tension_cut_off = 0.2*10**6


# k = 0.0123457 N/mm3
k_bedding =  0.0123457 *1E9
k_bedding = Es/(d_outer/2)
F_nok = 14E6
#bedding with ascend
k_bedding_low = 300
k_bedding_high = 30000
n_ring_k_bedding_low = 1
n_ring_k_bedding_var = 3
n_ring_90degree_k_bedding = 1
h_ring_90degree_k_bedding = d_outer * ((2**(1/2)) / 2)
x_start = 0
x_k_bedding_high = n_ring_90degree_k_bedding * l_ring
x_k_bedding_low = x_k_bedding_high + n_ring_k_bedding_var * l_ring
x_k_bedding_end = x_k_bedding_low + n_ring_k_bedding_low * l_ring

#soil/water
cover = 15.95 #m on top of tunnel
P_0   = 10E3    # surface load in N/m^2
gamma_sat = 18E3 # density of saturated soil kN/m3
gamma_water = 10E3  # density of saturated soil kN/m3
z_w = 5 # waterlevel under surface
# gamma_dry = 16E3 # density of dry soil kN/m3, not used
sigma_top = P_0 + gamma_sat*z_w + (gamma_sat-gamma_water)*(cover-z_w)
if variable_outside_loading:
    K_0 = 0.5
else:
    phi_soil = 37.5  # degrees
    K_0 = 1 - math.sin(math.radians(phi_soil))
H_water = cover-z_w
# variable outside loading
n_stressless_ring = 1
n_grout_ring = 3
n_water_soil_ring = 1
pressure_grout_center = 2*10**5  #N/m^2
gamma_grout = 16*10**3  #N/m^3
pressure_grout_bottom = pressure_grout_center - gamma_grout * d_outer / 2
pressure_grout_top = pressure_grout_center + gamma_grout * d_outer / 2



Meshsize = 0.2

colors = ["#ff0000", "#00ff00", "#ffff00", "#ff00ff", "#00ffff", "#0000ff", "#d5a6bd"]

############################################################################################
## Initialize project
############################################################################################

closeProject()
if "908644" in getcwd():
    newProject(r"C:\Users\908644\Documents\Projecten\2020_boortunnel\github_workfolder", 1000)
else:
    newProject( r"D:\projects_d\3D ringmodel\test", 1000 )
setViewerEnabled(False)
setModelAnalysisAspects(["STRUCT"])
setModelDimension("3D")
setDefaultMeshOrder("QUADRATIC")
setDefaultMesherType("HEXQUAD")
setDefaultMidSideNodeLocation("ONSHAP")


############################################################################################
## Creating data
############################################################################################

addElementData( "NL-interfaces" )
setParameter( DATA, "NL-interfaces", "./NOGEOM", [] )
setParameter( DATA, "NL-interfaces", "./INTEGR", [] )
setParameter( DATA, "NL-interfaces", "INTEGR", "HIGH" )

addElementData( "LIN-interfaces" )
setParameter( DATA, "LIN-interfaces", "./NOGEOM", [] )
############################################################################################
## Creating shapes
############################################################################################

createCylinder("Cylinder 1", [-l_ring, 0, 0], [1, 0, 0], d_inner/2, l_ring)
createCylinder("Cylinder 2", [-l_ring, 0, 0], [1, 0, 0], d_outer/2, l_ring)
subtract("Cylinder 2", ["Cylinder 1"], False, True)


createSheet("Sheet 1", [[-1-l_ring, 0,-(d_outer/2+1)],[1,0, -(d_outer/2+1)],[1,0, d_outer/2+1],[-1-l_ring,0, d_outer/2+1]])
arrayCopy(["Sheet 1"], [0, 0, 0], [0, 0, 0], [2*math.pi/n_segment, 0, 0], 1)

subtract("Cylinder 2", ["Sheet 2", "Sheet 1"], False, True)


renameShape("Cylinder 2_3", "element")
for shape in shapes():
    if "Cylinder" in shape:
        removeShape([shape])

#create ring
arrayCopy(["element"], [0, 0, 0], [0, 0, 0], [2*math.pi/n_segment, 0, 0], 6)
addSet(SHAPESET, "ring 0")
i = 1

for shape in namesIn(SHAPESET, "Shapes"):
    renameShape(shape, "element " + str(0) + "_" + str(i))
    i +=1
moveToShapeSet(namesIn(SHAPESET, "Shapes"),"ring 0")
ring_0_elements = namesIn(SHAPESET, "ring 0")

for i_ring in range(1, n_rings+1):
    if (i_ring % 2) == 0:
        alpha = math.pi/n_segment
    else:
        alpha = 0
    arrayCopy(namesIn(SHAPESET, "ring 0"),[l_ring * (i_ring), 0, 0], [0, 0, 0], [alpha, 0, 0], 1)
    addSet(SHAPESET, "ring " + str(i_ring))
    for shape in namesIn(SHAPESET, "ring 0"):
        if not shape in ring_0_elements:
            moveToShapeSet([shape],"ring " + str(i_ring))
    i = 1
    for shape in namesIn(SHAPESET, "ring " + str(i_ring)):
        renameShape(shape, "element " + str(i_ring) + "_" + str(i))
        setShapeColor(colors[i-1], ["element " + str(i_ring) + "_" + str(i)])
        i += 1

remove(SHAPESET, ["ring 0"])

for i_ring in range(1, n_rings+1):
    targets = namesIn(SHAPESET, "ring " + str(i_ring-1))
    tools = namesIn(SHAPESET, "ring " + str(i_ring))
    for target in targets:
        for tool in tools:
            imprintIntersection(target, tool, True)


#create nok
createSheet("nok 1", [[-1, -0.5*b_nok, (d_inner/2-0.1)],[-1, 0.5*b_nok, (d_inner/2-0.1)],[-1, 0.5*b_nok, (d_outer/2+0.1)],[-1, -0.5*b_nok, (d_outer/2+0.1)]])
rotate(["nok 1"], [0, 0, 0], [1, 0, 0], 2*math.pi/(n_segment)/3)
mirror(["nok 1"], [0, 0, 0], [False, True, False], True)
arrayCopy(["nok 1", "nok 2"], [0, 0, 0], [0, 0, 0], [2*math.pi/(n_segment),0, 0 ], n_segment-1)

nok=[]
for sh in shapes():
    if "nok" in sh:
        nok.append(sh)

nok_angles = []
for n in nok:
    angle= math.atan2((boundingBox(n)[5]+boundingBox(n)[4]),(boundingBox(n)[3]+boundingBox(n)[2]))
    if angle < 0 : angle += 2*math.pi
    nok_angles.append(angle)
#in case nokken between elements are needed
# for sh in shapes():
#     if "element" in sh:
#         projection(sh, nok, [1, 0, 0], True)
translate(nok, [l_ring*n_rings+1, 0, 0])
for sh in shapes():
    if "element "+str(n_rings) in sh:
        projection(sh, nok, [-1, 0, 0], True)
removeShape(nok)

tunnel_elements = []
for shape in shapes():
    tunnel_elements.append(tunnel_element(shape))

for i_ring in range(1,n_rings+1):
    sel_elements = get_rings(tunnel_elements,i_ring)
    for t_e in sel_elements:
        for face in t_e.faces:
            if distance([t_e.BBox[0],0,0],[face[0],0,0]) < 1E-3 or distance([t_e.BBox[1],0,0],[face[0],0,0]) < 1E-3:
                t_e.trans_face.append(face)
            else:
                if (abs(math.sqrt(face[1]**2+face[2]**2)-d_outer/2) < 1E-3):
                    t_e.outer_face.append(face)
                elif (abs(math.sqrt(face[1]**2+face[2]**2)-(d_inner/2)) < 1E-3):
                    t_e.inner_face.append(face)
                else:
                    t_e.long_face.append(face)


k_if = 100*(30E9 * 1E-3)/(1E-3*d_outer*math.pi/n_segment)
k_if_dummy = k_if * 0.0001
#longitudinal interfaces
if mc_int_long:
    addMaterial("interface long", "INTERF", "FRICTI", [])
    setParameter("MATERIAL", "interface long", "LINEAR/ELAS6/DSNZ", k_if)
    setParameter("MATERIAL", "interface long", "LINEAR/ELAS6/DSSX", k_if / 10)
    setParameter("MATERIAL", "interface long", "LINEAR/ELAS6/DSSY", k_if / 10)
    setParameter("MATERIAL", "interface long", "COULOM/COHESI", MC_cohesion)
    setParameter("MATERIAL", "interface long", "COULOM/PHI", MC_phi_fric)
    setParameter("MATERIAL", "interface long", "COULOM/PSI", MC_phi_dil)
    setParameter("MATERIAL", "interface long", "COULOM/OPNTYP", "TECTOF")
    setParameter("MATERIAL", "interface long", "COULOM/TENSTR", MC_tension_cut_off)
else:
    addMaterial("interface long", "INTERF", "NONLIF", [])
    setParameter(MATERIAL, "interface long", "LINEAR/ELAS6/DSNZ", k_if)
    setParameter(MATERIAL, "interface long", "LINEAR/ELAS6/DSSX", k_if/10)
    setParameter(MATERIAL, "interface long", "LINEAR/ELAS6/DSSY", k_if/10)





if create_dummy_interface_long:
    addMaterial("dummy interface long", "INTERF", "ELASTI", [])
    setParameter(MATERIAL, "dummy interface long", "LINEAR/ELAS6/DSNZ", k_if_dummy)
    setParameter(MATERIAL, "dummy interface long", "LINEAR/ELAS6/DSSX", k_if_dummy/10)
    setParameter(MATERIAL, "dummy interface long", "LINEAR/ELAS6/DSSY", k_if_dummy/10)

found_long =[]
ifound = 0
for element in tunnel_elements:
    for face1 in element.long_face:
        found_already = False
        for found in found_long:
            if distance(face1,found) < 1E-3:
                found_already = True
        if not found_already:
            for element2 in tunnel_elements:
                if element.name != element2.name and element.ring == element2.ring :
                    for face2 in element2.long_face:
                        if distance(face1, face2) < 1E-3:
                            found_long.append(face1)
                            ifound+=1
                            con_name = "long " + element.name+  " " + element2.name
                            createConnection(con_name , "INTER", SHAPEFACE, SHAPEFACE)
                            setParameter(GEOMETRYCONNECTION, con_name, "MODE", "CLOSED")
                            setElementClassType(GEOMETRYCONNECTION, con_name, "STPLIF")
                            assignMaterial("interface long", GEOMETRYCONNECTION, con_name)
                            setParameter(GEOMETRYCONNECTION, con_name, "FLIP", False)
                            attachTo(GEOMETRYCONNECTION, con_name, "SOURCE", element.name,
                                       [face1])
                            attachTo(GEOMETRYCONNECTION, con_name, "TARGET", element2.name,
                                       [face2])
                            assignElementData("NL-interfaces", GEOMETRYCONNECTION, con_name)
                            if create_dummy_interface_long:
                                con_name = "dummy long " + element.name + " " + element2.name
                                createConnection(con_name, "INTER", SHAPEFACE, SHAPEFACE)
                                setParameter(GEOMETRYCONNECTION, con_name, "MODE", "CLOSED")
                                setElementClassType(GEOMETRYCONNECTION, con_name, "STPLIF")
                                assignMaterial("dummy interface long", GEOMETRYCONNECTION, con_name)
                                setParameter(GEOMETRYCONNECTION, con_name, "FLIP", False)
                                attachTo(GEOMETRYCONNECTION, con_name, "SOURCE", element.name,
                                         [face1])
                                attachTo(GEOMETRYCONNECTION, con_name, "TARGET", element2.name,
                                         [face2])
                                assignElementData("LIN-interfaces", GEOMETRYCONNECTION, con_name)


#transverse interfaces
if mc_int_trans:
    addMaterial("interface trans", "INTERF", "FRICTI", [])
    setParameter("MATERIAL", "interface trans", "LINEAR/ELAS6/DSNZ", k_if)
    setParameter("MATERIAL", "interface trans", "LINEAR/ELAS6/DSSX", k_if / 10)
    setParameter("MATERIAL", "interface trans", "LINEAR/ELAS6/DSSY", k_if / 10)
    setParameter("MATERIAL", "interface trans", "COULOM/COHESI", MC_cohesion)
    setParameter("MATERIAL", "interface trans", "COULOM/PHI", MC_phi_fric)
    setParameter("MATERIAL", "interface trans", "COULOM/PSI", MC_phi_dil)
    setParameter("MATERIAL", "interface trans", "COULOM/OPNTYP", "TECTOF")
    setParameter("MATERIAL", "interface trans", "COULOM/TENSTR", MC_tension_cut_off)
else:
    addMaterial("interface trans", "INTERF", "NONLIF", [])
    setParameter(MATERIAL, "interface trans", "LINEAR/ELAS6/DSNZ", k_if)
    setParameter(MATERIAL, "interface trans", "LINEAR/ELAS6/DSSX", k_if/10)
    setParameter(MATERIAL, "interface trans", "LINEAR/ELAS6/DSSY", k_if/10)
if create_dummy_interface_trans:
    addMaterial("dummy interface trans", "INTERF", "ELASTI", [])
    setParameter(MATERIAL, "dummy interface trans", "LINEAR/ELAS6/DSNZ", k_if_dummy)
    setParameter(MATERIAL, "dummy interface trans", "LINEAR/ELAS6/DSSX", k_if_dummy/10)
    setParameter(MATERIAL, "dummy interface trans", "LINEAR/ELAS6/DSSY", k_if_dummy/10)

found_trans =[]
ifound = 0
for element in tunnel_elements:
    for face1 in element.trans_face:
        mindist = 1e6
        found_already = False
        for found in found_trans:
            if distance(face1,found) < 1E-3:
                found_already = True
        if not found_already:
            for element2 in tunnel_elements:

                if abs(element.ring - element2.ring) == 1:
                    for face2 in element2.trans_face:
                        if distance(face1, face2)< mindist and abs(face1[0] - face2[0]) < 1E-3:
                            mindist = distance(face1, face2)
                            face_target = face2
                            target_elm = element2
            if face1[0] > 0.1*l_ring and face1[0] < n_rings * l_ring - 0.1*l_ring:
                # print(element.name, element2.name)
                found_trans.append(face1)
                found_trans.append(face_target)

                ifound+=1
                con_name = "trans " + element.name+  " " + target_elm.name
                createConnection(con_name , "INTER", SHAPEFACE, SHAPEFACE)
                setParameter(GEOMETRYCONNECTION, con_name, "MODE", "CLOSED")
                setElementClassType(GEOMETRYCONNECTION, con_name, "STPLIF")
                assignMaterial("interface trans", GEOMETRYCONNECTION, con_name)
                setParameter(GEOMETRYCONNECTION, con_name, "FLIP", False)
                attachTo(GEOMETRYCONNECTION, con_name, "SOURCE", element.name,
                           [face1])
                attachTo(GEOMETRYCONNECTION, con_name, "TARGET", target_elm.name,
                           [face_target])
                assignElementData("NL-interfaces", GEOMETRYCONNECTION, con_name)
                if create_dummy_interface_trans:
                    con_name = "dummy trans " + element.name + " " + target_elm.name
                    createConnection(con_name, "INTER", SHAPEFACE, SHAPEFACE)
                    setParameter(GEOMETRYCONNECTION, con_name, "MODE", "CLOSED")
                    setElementClassType(GEOMETRYCONNECTION, con_name, "STPLIF")
                    assignMaterial("dummy interface trans", GEOMETRYCONNECTION, con_name)
                    setParameter(GEOMETRYCONNECTION, con_name, "FLIP", False)
                    attachTo(GEOMETRYCONNECTION, con_name, "SOURCE", element.name,
                             [face1])
                    attachTo(GEOMETRYCONNECTION, con_name, "TARGET", target_elm.name,
                             [face_target])
                    assignElementData("LIN-interfaces", GEOMETRYCONNECTION, con_name)


#bedding interfaces
if not create_ascending_bedding:
    addMaterial("outer interface", "INTERF", "ELASTI", [])
    setParameter(MATERIAL, "outer interface", "LINEAR/ELAS6/DSNZ", k_bedding)
    setParameter(MATERIAL, "outer interface", "LINEAR/ELAS6/DSSX", k_bedding/10)
    setParameter(MATERIAL, "outer interface", "LINEAR/ELAS6/DSSY", k_bedding/10)
else:
    # setFunctionValues("bedding",
    #                   [x_start, x_k_bedding_low, x_k_bedding_high, x_k_bedding_high+0.001, x_k_bedding_end],
    #                   [],
    #                   [0, h_ring_90degree_k_bedding, h_ring_90degree_k_bedding+0.001, d_outer],
    #                   [k_bedding_low, k_bedding_low, k_bedding_high, k_bedding_high, k_bedding_high,
    #                    k_bedding_low, k_bedding_low, k_bedding_high, k_bedding_high, k_bedding_high,
    #                    k_bedding_low, k_bedding_low, k_bedding_high, 0, 0,
    #                    k_bedding_low, k_bedding_low, k_bedding_high, 0, 0])
    setFunctionValues("bedding",
                      [x_start, x_k_bedding_high, x_k_bedding_high+0.001, x_k_bedding_low, x_k_bedding_end],
                      [],
                      [0, h_ring_90degree_k_bedding, h_ring_90degree_k_bedding+0.001, d_outer],
                      [k_bedding_high, k_bedding_high, k_bedding_high, k_bedding_low, k_bedding_low,
                       k_bedding_high, k_bedding_high, k_bedding_high, k_bedding_low, k_bedding_low,
                       0, 0, k_bedding_high, k_bedding_low, k_bedding_low,
                       0, 0, k_bedding_high, k_bedding_low, k_bedding_low])

    addMaterial("outer interface", "INTERF", "ELASTI", [])
    setParameter(MATERIAL, "outer interface", "LINEAR/ELAS6/DSNZ", 1)
    setParameter(MATERIAL, "outer interface", "LINEAR/ELAS6/DSSX", 1/10)
    setParameter(MATERIAL, "outer interface", "LINEAR/ELAS6/DSSY", 1/10)
    setMaterialFunction("outer interface", "LINEAR/ELAS6/DSNZ", "bedding")
    setMaterialFunction("outer interface", "LINEAR/ELAS6/DSSX", "bedding")
    setMaterialFunction("outer interface", "LINEAR/ELAS6/DSSY", "bedding")
addSet(GEOMETRYSUPPORTSET, "soilsprings")
createSurfaceSupport("total",  "soilsprings")
setParameter(GEOMETRYSUPPORT, "total", "AXES", [1, 2])
setParameter(GEOMETRYSUPPORT, "total", "TRANSL", [1, 1, 1])
setParameter(GEOMETRYSUPPORT, "total", "ROTATI", [0, 0, 0])

createConnection("bedding", "BOUNDA", SHAPEFACE)
setParameter(GEOMETRYCONNECTION, "bedding", "MODE", "CLOSED")
setElementClassType(GEOMETRYCONNECTION, "bedding", "STPLIF")
assignMaterial("outer interface", GEOMETRYCONNECTION, "bedding")
setParameter(GEOMETRYCONNECTION, "bedding", "FLIP", False)
assignElementData("LIN-interfaces", GEOMETRYCONNECTION, "bedding")
for element in tunnel_elements:
    for face in element.outer_face:
        attach(GEOMETRYSUPPORT, "total", element.name, [face])
        attachTo(GEOMETRYCONNECTION, "bedding", "SOURCE", element.name, [face])

addGeometry("axiaal coor", "SOLID", "STRSOL", [])
setParameter(GEOMET, "axiaal coor", "AXIAL", True)
setParameter(GEOMET, "axiaal coor", "AXIAL/CYLIN", [0, 0, 0, 1, 0, 0])

addMaterial("concrete", "CONCR", "LEI", [])
setParameter(MATERIAL, "concrete", "LINEAR/ELASTI/YOUNG", 3e+10)
setParameter(MATERIAL, "concrete", "LINEAR/ELASTI/YOUNG", 3e+10)
setParameter(MATERIAL, "concrete", "LINEAR/ELASTI/POISON", 0.15)
setParameter(MATERIAL, "concrete", "LINEAR/ELASTI/POISON", 0.15)
setParameter(MATERIAL, "concrete", "LINEAR/MASS/DENSIT", 2500)
setParameter(MATERIAL, "concrete", "LINEAR/ELASTI/YOUNG", 3e+10)
setElementClassType(SHAPE, shapes(), "STRSOL")
assignMaterial("concrete", SHAPE, shapes())
assignGeometry("axiaal coor", SHAPE, shapes())


setElementSize(shapes(), Meshsize, -1, True)
setMesherType(shapes(), "HEXQUAD")
clearMidSideNodeLocation(shapes())


############################################################################################
## Creating supports
############################################################################################

addSet(GEOMETRYSUPPORTSET, "axiaal")
createSurfaceSupport("Support 1", "axiaal")
setParameter(GEOMETRYSUPPORT, "Support 1", "AXES", [1, 2])
setParameter(GEOMETRYSUPPORT, "Support 1", "TRANSL", [1, 0, 0])
setParameter(GEOMETRYSUPPORT, "Support 1", "ROTATI", [0, 0, 0])
# attach(GEOMETRYSUPPORT, "Support 1", "element 1_7", [[0, 5.5996158, -0.89424583]])
sel_elements = get_rings(tunnel_elements,1)
for element in sel_elements:
    for face in element.trans_face:
        if face[0] < 0.5*l_ring:
            attach(GEOMETRYSUPPORT, "Support 1", element.name, [face])

############################################################################################
## Creating Loads
############################################################################################

addSet(GEOMETRYLOADSET, "nokken")
q_nok = F_nok/(t_ring * b_nok)

#find axial loading points
inok = 1
for te in tunnel_elements:
    if te.ring == n_rings:
        for face in te.trans_face:
            if face[0] > l_ring*n_rings - 0.5*l_ring:
                angle = math.atan2(face[2],face[1])
                if angle<0: angle +=2*math.pi
                mindist=2*math.pi
                for nok in nok_angles:
                    if abs(angle-nok) < 0.01:
                        createSurfaceLoad("Load nok "+str(inok), "nokken")
                        setParameter(GEOMETRYLOAD, "Load nok "+str(inok), "FORCE/VALUE", -q_nok)
                        attach(GEOMETRYLOAD, "Load nok "+str(inok), te.name, [face])
                        inok+=1

# create function
z_top = cover - 1
z_cent = cover + d_outer/2
z_bot = cover + d_outer + 1
load_ver_top = P_0 + gamma_sat * z_w + (gamma_sat - gamma_water) * (z_top - z_w)
load_ver_center = P_0 + gamma_sat * z_w + (gamma_sat - gamma_water) * (z_cent - z_w)
load_ver_bot = P_0 + gamma_sat * z_w + (gamma_sat - gamma_water) * (z_bot -z_w)
load_hor_top = K_0 * load_ver_top
load_hor_center = K_0 * load_ver_center
load_hor_bot = K_0 * load_ver_bot

# soil load
setFunctionValues("verticalsoilload", [], [-d_outer/2-1, -0.001, 0, d_outer/2+1],
                                             [-d_outer/2-1, -0.001, 0, d_outer/2+1],
                                             [load_ver_bot,     load_ver_bot,     load_ver_bot,     load_ver_bot,
                                              load_ver_center,  load_ver_center,  load_ver_center,  load_ver_center,
                                             -load_ver_center, -load_ver_center, -load_ver_center, -load_ver_center,
                                             -load_ver_top,    -load_ver_top,    -load_ver_top,    -load_ver_top]
                                            )

if variable_outside_loading:
    # Grout load
    setFunctionValues("hydro_grout", [], [], [-d_outer, d_outer], [pressure_grout_bottom, pressure_grout_top])
    grout_elements = [tunnel_element
                      for tunnel_element in tunnel_elements
                      if n_water_soil_ring < tunnel_element.ring <= n_water_soil_ring + n_grout_ring]
    addSet(GEOMETRYLOADSET, "grout")
    for te in grout_elements:
        createSurfaceLoad(te.name + " grout", "grout")
        setParameter(GEOMETRYLOAD, te.name + " grout", "FORCE/VALUE", 1)
        attach(GEOMETRYLOAD, te.name + " grout", te.name, te.outer_face)
        setValueFunction(GEOMETRYLOAD, te.name + " grout", "hydro_grout")

    # water and soil load
    water_soil_elements = [tunnel_element
                           for tunnel_element in tunnel_elements
                           if tunnel_element.ring <= n_water_soil_ring]

    addSet(GEOMETRYLOADSET, "buitenbelasting verticaal")
    for te in water_soil_elements:
        createSurfaceLoad(te.name + " buiten vert", "buitenbelasting verticaal")
        setParameter(GEOMETRYLOAD, te.name + " buiten vert", "FORCE/VALUE", 1)
        attach(GEOMETRYLOAD, te.name + " buiten vert", te.name, te.outer_face)
        setValueFunction(GEOMETRYLOAD, te.name + " buiten vert", "verticalsoilload")
        setParameter(GEOMETRYLOAD, te.name + " buiten vert", "FORCE/DIRECT", 3)

    setFunctionValues("horizontalsoilload", [], [-d_outer/2-1, -0.001, 0, d_outer/2+1],
                                                   [-d_outer/2-1, -0.001, 0, d_outer/2+1],
                                                   [load_hor_bot,    load_hor_bot,    -load_hor_bot,    -load_hor_bot,
                                                     load_hor_center, load_hor_center, -load_hor_center, -load_hor_center,
                                                     load_hor_center, load_hor_center, -load_hor_center, -load_hor_center,
                                                     load_hor_top,    load_hor_top,    -load_hor_top,    -load_hor_top])
    addSet(GEOMETRYLOADSET, "buitenbelasting horizontaal")
    for te in water_soil_elements:
        createSurfaceLoad(te.name + " buiten hor", "buitenbelasting horizontaal")
        setParameter(GEOMETRYLOAD, te.name + " buiten hor", "FORCE/VALUE", 1)
        attach(GEOMETRYLOAD, te.name + " buiten hor", te.name, te.outer_face)
        setValueFunction(GEOMETRYLOAD, te.name + " buiten hor", "horizontalsoilload")
        setParameter(GEOMETRYLOAD, te.name + " buiten hor", "FORCE/DIRECT", 2)
    # water load
    setFunctionValues("water", [], [], [-H_water, H_water], [2*H_water*10*1000 , 0])
    addSet(GEOMETRYLOADSET, "waterbelasting")
    for te in water_soil_elements:
        createSurfaceLoad(te.name + " water", "waterbelasting")
        setParameter(GEOMETRYLOAD, te.name + " water", "FORCE/VALUE", -1)
        attach(GEOMETRYLOAD, te.name + " water", te.name, te.outer_face)
        setValueFunction(GEOMETRYLOAD, te.name + " water", "water")
else:
    addSet(GEOMETRYLOADSET, "buitenbelasting verticaal")
    for te in tunnel_elements:
        createSurfaceLoad(te.name + " buiten vert", "buitenbelasting verticaal")
        setParameter(GEOMETRYLOAD, te.name + " buiten vert", "FORCE/VALUE", 1)
        attach(GEOMETRYLOAD, te.name + " buiten vert", te.name, te.outer_face)
        setValueFunction(GEOMETRYLOAD, te.name + " buiten vert", "verticalsoilload")
        setParameter(GEOMETRYLOAD, te.name + " buiten vert", "FORCE/DIRECT", 3)

    setFunctionValues("horizontalsoilload", [], [-d_outer/2-1, -0.001, 0, d_outer/2+1],
                                                   [-d_outer/2-1, -0.001, 0, d_outer/2+1],
                                                   [load_hor_bot,    load_hor_bot,    -load_hor_bot,    -load_hor_bot,
                                                     load_hor_center, load_hor_center, -load_hor_center, -load_hor_center,
                                                     load_hor_center, load_hor_center, -load_hor_center, -load_hor_center,
                                                     load_hor_top,    load_hor_top,    -load_hor_top,    -load_hor_top])
    addSet(GEOMETRYLOADSET, "buitenbelasting horizontaal")
    for te in tunnel_elements:
        createSurfaceLoad(te.name + " buiten hor", "buitenbelasting horizontaal")
        setParameter(GEOMETRYLOAD, te.name + " buiten hor", "FORCE/VALUE", 1)
        attach(GEOMETRYLOAD, te.name + " buiten hor", te.name, te.outer_face)
        setValueFunction(GEOMETRYLOAD, te.name + " buiten hor", "horizontalsoilload")
        setParameter(GEOMETRYLOAD, te.name + " buiten hor", "FORCE/DIRECT", 2)

    # water load
    setFunctionValues("water", [], [], [-H_water, H_water], [2*H_water*10*1000 , 0])
    addSet(GEOMETRYLOADSET, "waterbelasting")
    for te in tunnel_elements:
        createSurfaceLoad(te.name + " water", "waterbelasting")
        setParameter(GEOMETRYLOAD, te.name + " water", "FORCE/VALUE", -1)
        attach(GEOMETRYLOAD, te.name + " water", te.name, te.outer_face)
        setValueFunction(GEOMETRYLOAD, te.name + " water", "water")

#selfweight
addSet(GEOMETRYLOADSET, "Selfweight")
createModelLoad("Selfweight", "Selfweight")


############################################################################################
## Creating loadcombinations
############################################################################################
if variable_outside_loading:
    LC = 7
else:
    LC = 6
setDefaultGeometryLoadCombinations()
addGeometryLoadCombination("")
setGeometryLoadCombinationFactor(f"Geometry load combination {LC}", "buitenbelasting verticaal", 1)
setGeometryLoadCombinationFactor(f"Geometry load combination {LC}", "buitenbelasting horizontaal", 1)
setGeometryLoadCombinationFactor(f"Geometry load combination {LC}", "waterbelasting", 1)
setGeometryLoadCombinationFactor(f"Geometry load combination {LC}", "Selfweight", 1)
setGeometryLoadCombinationFactor(f"Geometry load combination {LC}", "nokken", 0)
if variable_outside_loading:
    setGeometryLoadCombinationFactor(f"Geometry load combination {LC}", "grout", 1)


rename(GEOMETRYLOADCOMBINATION, "Geometry load combination 1", "nokken")
if variable_outside_loading:
    rename(GEOMETRYLOADCOMBINATION, "Geometry load combination 2", "grout")
    i = 3
else:
    i = 2
rename(GEOMETRYLOADCOMBINATION, f"Geometry load combination {i}", "buitenbelasting verticaal")
i += 1
rename(GEOMETRYLOADCOMBINATION, f"Geometry load combination {i}", "buitenbelasting horizontaal")
i += 1
rename(GEOMETRYLOADCOMBINATION, f"Geometry load combination {i}", "waterbelasting")
i += 1
rename(GEOMETRYLOADCOMBINATION, f"Geometry load combination {i}", "Selfweight")
i += 1
rename(GEOMETRYLOADCOMBINATION, f"Geometry load combination {i}", "radial and tangitional")
############################################################################################
## Creating analysis
############################################################################################
if create_analysis:
    addAnalysis("Analysis1")
    addAnalysisCommand("Analysis1", "LINSTA", "Structural linear static")


    addAnalysis("Analysis2")
    addAnalysisCommand("Analysis2", "NONLIN", "Structural nonlinear")
    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "TYPE/GEOMET", True)
    addAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(1)/LOAD/LOADNR")
    if variable_outside_loading:
        setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(1)/LOAD/LOADNR", 7)
    else:
        setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(1)/LOAD/LOADNR", 6)
    if mc_int_trans:
        setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(1)/LOAD/STEPS/EXPLIC/SIZES", "0.100000(2) 0.0500000(6) 0.0100000(46) 0.00500000(8)")
    else:
        setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(1)/LOAD/STEPS/EXPLIC/SIZES",
                                 "0.1000(2) 0.0500(6) 0.01(50)")
    # setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(1)/ITERAT/METHOD/NEWTON/TYPNAM", "MODIFI")
    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(1)/ITERAT/METHOD/METNAM", "SECANT")
    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(1)/ITERAT/METHOD/SECANT/TYPNAM", "BFGS")
    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(1)/ITERAT/MAXITE", 20)
    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(1)/ITERAT/CONVER/ENERGY/NOCONV", "CONTIN")
    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(1)/ITERAT/CONVER/DISPLA", False)
    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(1)/ITERAT/CONVER/FORCE", False)
    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(1)/ITERAT/CONVER/ENERGY", True)
    if mc_int_trans:
        setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(1)/ITERAT/CONVER/ENERGY/TOLCON", 0.001)
    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(1)/ITERAT/LINESE", True)


    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT/EXETYP", "LOAD")
    addAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(2)/LOAD/LOADNR")
    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(2)/LOAD/LOADNR", 1)
    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(2)/LOAD/STEPS/EXPLIC/SIZES", "0.2(5)")

    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(2)/ITERAT/MAXITE", 10)
    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(2)/ITERAT/CONVER/ENERGY", True)
    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(2)/ITERAT/CONVER/DISPLA", False)
    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(2)/ITERAT/CONVER/FORCE", False)
    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(2)/ITERAT/CONVER/ENERGY/TOLCON", 0.001)
    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(2)/ITERAT/LINESE", True)
    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(2)/ITERAT/CONVER/ENERGY/NOCONV", "CONTIN")

    copyAnalysisCommandDetail("Analysis2", "Structural nonlinear", "OUTPUT(1)", "")
    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "OUTPUT(2)/SELTYP", "USER")
    addAnalysisCommandDetail("Analysis2", "Structural nonlinear", "OUTPUT(2)/USER")
    addAnalysisCommandDetail("Analysis2", "Structural nonlinear", "OUTPUT(2)/USER/STRESS(1)/TOTAL/CAUCHY/LOCAL")
    addAnalysisCommandDetail("Analysis2", "Structural nonlinear", "OUTPUT(2)/USER/STRESS(2)/TOTAL/CAUCHY/PRINCI")
    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "OUTPUT(2)/USER/STRESS(1)/TOTAL/CAUCHY/LOCAL/LOCATI", "INTPNT")
    setAnalysisCommandDetail("Analysis2", "Structural nonlinear", "OUTPUT(2)/USER/STRESS(2)/TOTAL/CAUCHY/PRINCI/LOCATI", "INTPNT")
    addAnalysisCommandDetail("Analysis2", "Structural nonlinear", "OUTPUT(2)/USER/DISPLA")


    renameAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(1)", "radial and tangentional")
    renameAnalysisCommandDetail("Analysis2", "Structural nonlinear", "EXECUT(2)", "nokken")
    renameAnalysis("Analysis1", "LINSTA")
    renameAnalysis("Analysis2", "NLSTA")
if create_mesh:
    if "Shapes" in names('SHAPESET'):
        remove(SHAPESET, ["Shapes"])
    generateMesh([])


setViewerEnabled(True)
draw()
hide( ELEMENTSET, [ "bedding" ] )

if run_analysis_linsta:
    runSolver(["LINSTA"])
if run_analysis_nlsta:
    runSolver(["NLSTA"])

print("Python script is finished.")
