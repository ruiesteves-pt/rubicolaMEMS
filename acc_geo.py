# MEMS Accelerometer 3D Geometry
# @ruiesteves

# Imports
import pygmsh as pg     # geometry & meshing definition
import meshio           # meshing export


#Functions
def build(suspension_beam_width,proof_mass_length):

    # Geometric parameters
    cl = 55
    proof_mass_cl = 150
    scale = 1e-6
    suspension_beam_length = 3300*scale
    beam_thickness = 69*scale
    small_beam_length = 500*scale
    proof_mass_thickness = 320*scale
    beam_dist = 500*scale
    beam_l = 122.5 * scale
    beam_h = 177.5 * scale
    beam_dist2 = beam_dist + suspension_beam_length
    beam_dist3 = proof_mass_length + suspension_beam_width + small_beam_length
    beam_dist_final = beam_dist2 - beam_dist3
    beam_lower = (proof_mass_thickness - beam_thickness)/2
    beam_to_mass = (beam_dist_final+suspension_beam_width+small_beam_length + proof_mass_length) - suspension_beam_length

    # Geometry build
    geom = pg.opencascade.Geometry()

    # Beam1_1
    p1 = [0,beam_dist_final,beam_lower]
    p2 = [suspension_beam_length,suspension_beam_width,beam_thickness]
    beam1_1 = geom.add_box(p1,p2,char_length=cl*scale)

    # Beam1_2
    p3 = [suspension_beam_length-suspension_beam_width,beam_dist_final + suspension_beam_width,beam_lower]
    p4 = [suspension_beam_width,small_beam_length,beam_thickness]
    beam1_2 = geom.add_box(p3,p4,char_length=cl*scale)

    # Beam1 complete
    beam1 = geom.boolean_union([beam1_1,beam1_2])

    # Proof mass
    p1 = [beam_dist_final+suspension_beam_width+small_beam_length,beam_dist_final+suspension_beam_width+small_beam_length,0]
    p2 = [proof_mass_length,proof_mass_length,proof_mass_thickness]
    proof_mass = geom.add_box(p1,p2,char_length=proof_mass_cl*scale)

    # Beam2_1
    p1 = [beam_dist_final,beam_dist_final+suspension_beam_width+small_beam_length+beam_to_mass,beam_lower]
    p2 = [suspension_beam_width,suspension_beam_length,beam_thickness]
    beam2_1 = geom.add_box(p1,p2,char_length=cl*scale)

    # Beam2_2
    p1 = [beam_dist_final+suspension_beam_width,beam_dist_final+suspension_beam_width+small_beam_length+beam_to_mass,beam_lower]
    p2 = [small_beam_length,suspension_beam_width,beam_thickness]
    beam2_2 = geom.add_box(p1,p2,char_length=cl*scale)

    # Beam 2
    beam2 = geom.boolean_union([beam2_1,beam2_2])


    # Beam3_1
    p1 = [beam_dist_final+suspension_beam_width+small_beam_length+proof_mass_length,beam_dist_final+suspension_beam_width+small_beam_length+proof_mass_length - beam_to_mass - suspension_beam_width,beam_lower]
    p2 = [small_beam_length,suspension_beam_width,beam_thickness]
    beam3_1 = geom.add_box(p1,p2,char_length=cl*scale)

    # Beam3_2
    p1 = [beam_dist_final+suspension_beam_width+small_beam_length+proof_mass_length+small_beam_length,0,beam_lower]
    p2 = [suspension_beam_width,suspension_beam_length,beam_thickness]
    beam3_2 = geom.add_box(p1,p2,char_length=cl*scale)

    # Beam3
    beam3 = geom.boolean_union([beam3_1,beam3_2])


    # Beam4_1
    p1 = [beam_dist_final+suspension_beam_width+small_beam_length+beam_to_mass,beam_dist_final+suspension_beam_width+small_beam_length+proof_mass_length,beam_lower]
    p2 = [suspension_beam_width,small_beam_length,beam_thickness]
    beam4_1 = geom.add_box(p1,p2,char_length=cl*scale)

    # Beam4_2
    p1 = [beam_dist_final+suspension_beam_width+small_beam_length+beam_to_mass,beam_dist_final+suspension_beam_width+small_beam_length+proof_mass_length+small_beam_length,beam_lower]
    p2 = [suspension_beam_length,suspension_beam_width,beam_thickness]
    beam4_2 = geom.add_box(p1,p2,char_length=cl*scale)

    # Beam 4
    beam4 = geom.boolean_union([beam4_1,beam4_2])


    #Complete Union
    final = geom.boolean_union([beam1,proof_mass,beam2,beam3,beam4])

    mesh = pg.generate_mesh(geom,gmsh_path="/home/ruiesteves/Documents/Tese/MechanicalModel/gmsh-4.5.2-Linux64/bin/gmsh")   # Be sure to change the gmsh_path to the installed folder
    meshio.write("accelerometer.xml",mesh)
