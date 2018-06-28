import numpy as np
from numpy import linalg as la
from modflow import get_fd_mesh, ALGORITHM, LENGTH_UNIT, TIME_UNIT
import unicodedata

def to_range(indices):
    indices.sort()
    rng = [(indices[0], indices[0])]
    for i in indices[1:]:
        if i == (rng[-1][1]+1):
            rng[-1] = (rng[-1][0], i)
        else:
            rng.append((i, i))
    return rng

def redu_range(values):
    cv = values[0][0]
    idxs = [values[0][1]]
    redu = []
    
    for v, k in values[1:]:
        if v == cv:
            idxs.append( k )
        else:
            redu.append( ( cv, idxs ) )
            cv = v
            idxs = [k]
    
    redu.append( ( cv, idxs ) )
    return redu

def write_node( fd, nump, shape, node_str ):
    # NODE
    fd.write("NODE\n")
    
    
    all_elems = []
    
    for lay in xrange(shape[0]):
        ini_lay = lay*(shape[1]+1)*(shape[2]+1)
        nxt_lay = (lay+1)*(shape[1]+1)*(shape[2]+1)
        for row in xrange(shape[1]):
            ini_row = row*(shape[2]+1)
            nxt_row = (row+1)*(shape[2]+1)
            for col in xrange(shape[2]):
                elem = [ini_lay+nxt_row+col+1,
                        ini_lay+ini_row+col+1,
                        ini_lay+ini_row+col,
                        ini_lay+nxt_row+col,
                        nxt_lay+nxt_row+col+1,
                        nxt_lay+ini_row+col+1,
                        nxt_lay+ini_row+col,
                        nxt_lay+nxt_row+col]
                
                all_elems.append(elem)
                
                fd.write( "".join( map(node_str, elem) ) )
                fd.write( "\n" )
    
    return all_elems

def write_coor( fd, shape, dX, dY, X_inf, Y_sup ):
    # COOR
    fd.write("COOR\n")
    xy_points = []
    for row in xrange(shape[1]+1):
        for col in xrange(shape[2]+1):
            point = (X_inf+col*dX,Y_sup-row*dY)
            # This layer is the top.
            xy_points.append(point)
    
    for i in xrange(len(xy_points)): 
        fd.write( "%21.14le," % xy_points[i][0] )
        if (i+1)%12 == 0:
            fd.write("\n")
    
    if (i+1)%12 != 0:
        fd.write("\n")
    
    for i in xrange(len(xy_points)): 
        fd.write( "%21.14le," %  xy_points[i][1] )
        if (i+1)%12 == 0:
            fd.write("\n")
    
    if (i+1)%12 != 0:
        fd.write("\n")
    
    return xy_points

def write_elev_i( fd, shape, dX, dY, X_inf, Y_sup, Z_top, Z_bottoms, node_str ):
    fd.write("ELEV_I\n")
    all_points = []
    node_lists = { }
    node_lists["top"] = []
    slice_values = []
    conv_range = lambda r: str(r[0]+1) if r[0]==r[1] else "%s-%s" % (r[0]+1, r[1]+1)
    
    def write_slice( ):
        for Z, idxs in slice_values:
            f = "%21.14le %s\n" if Z > 0 else " %21.14le %s\n"
            fd.write(f % ( Z, " ".join(map(conv_range,to_range(idxs) ) )))
    
    for row in xrange(shape[1]+1):
        for col in xrange(shape[2]+1):
            cnt = 0
            Z = 0
            
            for i, j in [( -1, -1 ), (-1, 1), (1, -1), (1, 1)]:
                if row+i < 0 or row+i >= shape[1]:
                    continue
                if col+j < 0 or col+j >= shape[2]:
                    continue
                cnt += 1
                Z += Z_top[row+i,col+j]
            
            Z /= cnt
            
            point = (X_inf+col*dX,Y_sup-row*dY, Z)
            
            # This layer is the top.
            node_lists["top"].append(len(all_points))
            all_points.append(point)
            idx = len(slice_values)
            slice_values.append( (Z, idx) )
    
    slice_values.sort()
    slice_values = redu_range( slice_values )
    write_slice()
    
    node_lists["top"] = to_range(node_lists["top"])
    node_lists["bottom"] = []
    node_lists["front"]  = []
    node_lists["back"]   = []
    node_lists["left"]   = []
    node_lists["right"]  = []
    
    for lay in xrange(shape[0]):
        slice_values = []
        for row in xrange(shape[1]+1):
            for col in xrange(shape[2]+1):
                cnt = 0
                Z = 0
                for i, j in [( -1, -1 ), (-1, 1), (1, -1), (1, 1)]:
                    if row+i < 0 or row+i >= shape[1]:
                        continue
                    if col+j < 0 or col+j >= shape[2]:
                        continue
                    Z += Z_bottoms[lay,row+i,col+j]
                    cnt += 1
                
                Z /= cnt
                point = (X_inf+col*dX,Y_sup-row*dY, Z)
                
                if lay+1 == shape[0]:
                    node_lists["bottom"].append(len(all_points))
                elif row == 0:
                    node_lists["back"].append(len(all_points))
                elif row == shape[1]:
                    node_lists["front"].append(len(all_points))
                elif col == 0:
                    node_lists["left"].append(len(all_points))
                elif col == shape[2]:
                    node_lists["right"].append(len(all_points))
                
                all_points.append(point)
                idx = len(slice_values)
                slice_values.append( (Z, idx) )

        fd.write( "%d\n" % (lay+2) ) 
        slice_values.sort()
        slice_values = redu_range( slice_values )
        write_slice()
    
    node_lists["bottom"] = to_range(node_lists["bottom"])
    node_lists["front"] =  to_range(node_lists["front"])
    node_lists["back"] =   to_range(node_lists["back"])
    node_lists["left"] =   to_range(node_lists["left"])
    node_lists["right"] =  to_range(node_lists["right"])

    return node_lists, all_points

def write_xyzcoor( fd, shape, dX, dY, X_inf, Y_sup, Z_top, Z_bottoms ):
    # XYZCOOR
    fd.write("XYZCOOR\n")
    all_points = []
    node_lists = { }
    node_lists["top"] = []
    
    for row in xrange(shape[1]+1):
        for col in xrange(shape[2]+1):
            cnt = 0
            Z = 0
            
            for i, j in [( -1, -1 ), (-1, 1), (1, -1), (1, 1)]:
                if row+i < 0 or row+i >= shape[1]:
                    continue
                if col+j < 0 or col+j >= shape[2]:
                    continue
                cnt += 1
                Z += Z_top[row+i,col+j]
            
            Z /= cnt
            
            point = (X_inf+col*dX,Y_sup-row*dY, Z)
            
            # This layer is the top.
            node_lists["top"].append(len(all_points))
            all_points.append(point)
            fd.write( "%21.14le,%21.14le,%21.14le\n" % point )
    
    node_lists["top"] = to_range(node_lists["top"])
    node_lists["bottom"] = []
    node_lists["front"]  = []
    node_lists["back"]   = []
    node_lists["left"]   = []
    node_lists["right"]  = []
    
    for lay in xrange(shape[0]):
        for row in xrange(shape[1]+1):
            for col in xrange(shape[2]+1):
                cnt = 0
                Z = 0
                for i, j in [( -1, -1 ), (-1, 1), (1, -1), (1, 1)]:
                    if row+i < 0 or row+i >= shape[1]:
                        continue
                    if col+j < 0 or col+j >= shape[2]:
                        continue
                    Z += Z_bottoms[lay,row+i,col+j]
                    cnt += 1
                
                Z /= cnt
                point = (X_inf+col*dX,Y_sup-row*dY, Z)
                
                if lay+1 == shape[0]:
                    node_lists["bottom"].append(len(all_points))
                elif row == 0:
                    node_lists["back"].append(len(all_points))
                elif row == shape[1]:
                    node_lists["front"].append(len(all_points))
                elif col == 0:
                    node_lists["left"].append(len(all_points))
                elif col == shape[2]:
                    node_lists["right"].append(len(all_points))
                
                all_points.append(point)
                fd.write( "%21.14le,%21.14le,%21.14le\n" % tuple(point) )
    
    node_lists["bottom"] = to_range(node_lists["bottom"])
    node_lists["front"] =  to_range(node_lists["front"])
    node_lists["back"] =   to_range(node_lists["back"])
    node_lists["left"] =   to_range(node_lists["left"])
    node_lists["right"] =  to_range(node_lists["right"])

    return node_lists, all_points

def create_feflow_input(name, model, units_data, 
                        length_units=LENGTH_UNIT.METERS, 
                        rows=100, cols=100, layers=100, 
                        bbox=None, angle=20, dz_min=1.0, 
                        time_units=TIME_UNIT.SECONDS, 
                        algorithm=ALGORITHM.REGULAR, 
                        faults_data={},
                        faults_method=u'regular'):
    
    layers, Z_top, Z_bottoms, dY, dX, X_inf, Y_sup, I_bound, chani_var, K_hor, K_ver, K_anisotropy_hor, bbox = get_fd_mesh(model, units_data, length_units, 
                                                                                                                           rows, cols, layers, bbox, angle, 
                                                                                                                           dz_min, time_units, algorithm, 
                                                                                                                           faults_data, faults_method)
    fd = open( name+".fem", "w" )
    
    # PROBLEM
    fd.write( "PROBLEM: " )
    fd.write( name )
    fd.write( "\n" )
    
    # CLASS
    fd.write( "CLASS: (v.7.007.14984)\n" )
    #    2    1    0    3    1    0    8    8    0    0
    fd.write( "%(ic1)4d %(ic2)4d %(proj)4d %(ndm)4d %(n_layers)4d %(ic0)4d %(save_fsize_rreal)4d %(save_fsize_creal)4d %(n_species)4d %(has_nonvertical_join_faces)4d\n" % 
              {"ic1": 0, "ic2": 0, "proj": 0, "ndm": 3, "n_layers": layers, "ic0": 0, "save_fsize_rreal": 8, "save_fsize_creal": 8, "n_species": 0, "has_nonvertical_join_faces": 1}  )
    
    # DIMENS
    fd.write('DIMENS\n')

    shape = Z_bottoms.shape
    nump = (shape[0]+1)*(shape[1]+1)*(shape[2]+1)
    ne = shape[0]*shape[1]*shape[2]
    
    first  = "%(nump)6d %(ne)6d %(nbn)6d %(numb_dt)6d " % {"nump": nump, "ne": ne, "nbn": 8, "numb_dt": 1}
    second = "%(icrank)6d %(upwind)6d %(obs)6d %(optim)6d %(phreatic)6d " % {"icrank": 0, "upwind": 0, "obs": 0, "optim": 0, "phreatic": 0}
    third  = "%(nwca)6d %(npcor)6d %(adaptive_mesh)6d %(special_fem_process_id)6d " % {"nwca": 2, "npcor": 0, "adaptive_mesh": 0, "special_fem_process_id": 0}
    fourth = "%(sorption_type)6d %(reaction_type)6d %(dispersion_type)6d " % {"sorption_type": 0, "reaction_type": 0, "dispersion_type": 0}
    fifth  = "%(reaction_popt)6d %(dual_nodes)6d " % {"reaction_popt": 0, "dual_nodes": 0}
    
    fd.write( first + second + third + fourth + fifth + '\n' )
    
    # SCALE
    fd.write("SCALE\n")
    PIX = 1000
    bdiff = [ bbox[5]-bbox[2], bbox[4]-bbox[1], bbox[3]-bbox[0] ]
    maxd = max( bdiff )
    fd.write("%(scale_fac)13.6e,%(scale_l)13.6e,%(scale_ratio)13.6e,%(eps)13.6e,%(x_o)13.6e,%(y_o)13.6e\n" % 
             {"scale_fac": PIX/maxd, "scale_l": PIX, "scale_ratio": 1, "eps": 1.0e-7, "x_o": 0, "y_o": 0})
    
    if nump < 10000:
        num_columns = 80/6
        node_str = lambda n: "%5d" % (n+1)
    elif nump < 1000000:
        num_columns = 80/8
        node_str = lambda n: "%7d" % (n+1)
    else:
        num_columns = 80/10
        node_str = lambda n: "%9d" % (n+1)
    
    #NODE
    all_elems = write_node( fd, nump, shape, node_str )
    
    #XYZCOOR
    # node_lists, all_points = write_xyzcoor( fd, shape, dX, dY, X_inf, Y_sup, Z_top, Z_bottoms )
    write_coor( fd, shape, dX, dY, X_inf, Y_sup )
    node_lists, all_points = write_elev_i( fd, shape, dX, dY, X_inf, Y_sup, Z_top, Z_bottoms, node_str )
    
    # FLOW_I_BC.
    fd.write("FLOW_I_BC\n")
    
    # INIT_I_FLOW.
    fd.write("INIT_I_FLOW\n")
    fd.write("%(head_ini_level)21.14le 1-%(nump)s\n" % { "head_ini_level": 0.0, "nump": nump})
    
    # MAT_I_FLOW
    first_elems = (shape[1]+1)*(shape[2]+1)
    elem_units = []
    for i, e in enumerate(all_elems):
        if model.params.get('map', 'basic') == 'soils' and i < first_elems:
            mid = [0, 0, 0]
            for n in e[:4]:
                for i in range(3):
                    mid[i] += all_points[n][i]
            for i in range(3):
                mid[i]=mid[i]/4.0
            h = model.height( (mid[0], mid[1]) )
            unit = model.closest( (mid[0], mid[1], h) )[0]
        else:
            mid = [0, 0, 0]
            for n in e:
                for i in range(3):
                    mid[i] += all_points[n][i]
            for i in range(3):
                mid[i]=mid[i]/8.0
            unit = model.closest( mid )[0]
        elem_units.append(unit)
    
    unit_indices = {}
    for i in xrange(len(all_elems)):
        if elem_units[i] in unit_indices:
            unit_indices[elem_units[i]].append(i)
        else:
            unit_indices[elem_units[i]] = [i]
    
    ranges = {}
    for u in unit_indices:
        indices = unit_indices[u]
        ranges[u] = to_range(indices)
    
    mtidxs = [101, 103, 107 ,137]
    comments = ["K_xx", "K_yy", "K_zz", "Inactive"]
    for u in units_data:
        ud = units_data[u]
        # units_data[u] = (ud[1], ud[1]*ud[2], ud[3], 1.0-ud[0])
        units_data[u] = (ud[0], ud[0]*ud[1], ud[2], ud[3])
    
    defaults = [ min( [ v[0] for k, v in units_data.iteritems() ] ),
                 min( [ v[1] for k, v in units_data.iteritems() ] ),
                 min( [ v[2] for k, v in units_data.iteritems() ] ),
                 0.0 ]
    
    conv_range = lambda r: str(r[0]+1) if r[0] == r[1] else "%s-%s" % (r[0]+1, r[1]+1)
    fd.write("MAT_I_FLOW\n")
    for i in range(4):
        fd.write("%(material_type_index)d %(default_value)13e \"%(comment)s\"\n" %
                 { "material_type_index": mtidxs[i], "default_value": defaults[i], "comment": comments[i] } )
        # join same material values in a single dictionary.
        materials = {}
        for u in unit_indices:
            if units_data[u][i] != defaults[i]:
                v = units_data[u][i]
                if v in materials:
                    materials[v] = materials[v] + ranges[u]
                else:
                    materials[v] = ranges[u]
        
        for v in materials:
            rngs = sorted(materials[v])
            j = num_columns
            node_list = " ".join( map( conv_range, rngs[:num_columns] ) )
            fd.write("%(material_value)13e\t%(node_list)s\n" % { "material_value": v, "node_list": node_list } )
            while j < len(rngs):
                rng = rngs[j:j+num_columns]
                node_list = " ".join( map( conv_range, rng ) )
                fd.write("\t\t%(node_list)s\n" % { "node_list": node_list } )
                j += num_columns
    
    # TINI
    fd.write("TINI\n0\n")
    
    # STEPS2
    fd.write("STEPS2\n")
    fd.write("0,0.001,365\n")

    # GOBS
    fd.write("GOBS\n")
    fd.write("\n")

    # GRAVITY
    fd.write("GRAVITY\n")
    fd.write("%(gx)3d %(gy)2d %(gz)2d\n" % { "gx": 0, "gy": 0, "gz": -1})
    
    # PROJGRAVITY
    fd.write("PROJGRAVITY\n")
    fd.write("%(refid)2d %(used)2d\n" % {"refid": -1,"used": 0})
    
    # ERROR_NORM
    fd.write("ERROR_NORM\n")
    fd.write("%(norm)1d\n" % {"norm": 2})
    
    # PHYSICS
    fd.write("PHYSICS\n")
    first = "%(cond_well_bore)1d,%(mass_matrices)1d,%(process_boussinesq)1d,%(process_viscosity)1d, " % {"cond_well_bore": 1, 
                                                                                                         "mass_matrices": 0, 
                                                                                                         "process_boussinesq": 0, 
                                                                                                         "process_viscosity": 1} 
    second = "%(fs_constraint.bot)1d,%(flags.quadr)1d,%(fs_constraint.top)1d,%(max_itera)d, " % {"fs_constraint.bot": 1, 
                                                                                                  "flags.quadr": 0, 
                                                                                                  "fs_constraint.top": 1, 
                                                                                                  "max_itera": 10} 
    third = "%(dry_delta_head)e,%(flags.revff)1d,%(delta_hysteresis)e,%(fs_constraint.bot_spec)1d, " % {"dry_delta_head": 1.00000e-03, 
                                                                                                               "flags.revff": 0, 
                                                                                                               "delta_hysteresis": 1.000000e-03, 
                                                                                                               "fs_constraint.bot_spec": 0 }
    fourth = "%(fs_constraint.top_spec)1d,%(process_densratio)1d,%(scale_ss_by_pseudosaturation)1d" % {"fs_constraint.top_spec": 0, 
                                                                                                           "process_densratio": 0, 
                                                                                                           "scale_ss_by_pseudosaturation": 1}
    fd.write( first + second + third + fourth + "\n")
    
    # OPTIONS
    fd.write("OPTIONS\n")
    
    options = [("SkipMeshFill", "false"),
               ("SkipMeshDraw", "false"),
               ("BackingStore", "false"),
               ("OpaqueFringeMode", "false"),
               ("LegendViewMode", "false"),
               ("VelocityApproximationType", "0"),
               ("EquationSolverType", "11"),
               ("IterSymmType", "0"),
               ("IterNonsymmType", "4"),
               ("DirectSolverType", "1"),
               ("PreconditionSymmType", "0"),
               ("MaxNumbOrthogonal", "5"),
               ("ReorderingMethod", "2"),
               ("RecordCPU", "false"),
               ("HideVelocity", "false"),
               ("VelocityType", "0"),
               ("GeoCS", ""),
               ("SubdivisionCurvedEdges", "32"),
               ("ShowOverview", "true"),
               ("ShowMeshWindow", "true"),
               ("TimeUnit", "0"),
               ("SpatialIndexing", "false"),
               ("MinimalSliceDistance", "0.1"),
               ("UseUnsmoothVelocityField", "false"),
               ("UnsatFractureMode", "1"),
               ("ComputeIterativeSolverResiduals", "0"),
               ("SkipOutputSteps", "1")]
    
    for opt in options:
        fd.write("  %(name)s=%(value)s\n" % {"name": opt[0], "value": opt[1]} )  
    
    # SYMSOLVE ?
    # NONSYMSOLVE ?
    
    # NODALSETS
    fd.write("NODALSETS\n")
    for k in node_lists:
        rngs = node_lists[k]
        j = num_columns
        node_list = " ".join( map( conv_range, rngs[:j] ) )
        fd.write("  %(setname)s %(setlist)s\n" % { "setname": k, "setlist": node_list })
        while j < len(rngs):
            rng = rngs[j:j+num_columns]
            node_list = " ".join( map( conv_range, rng ) )
            fd.write("\t\t%(setlist)s\n" % { "setname": k, "setlist": node_list })
            j += num_columns
    
    # ELEMENTALSETS
    fd.write("ELEMENTALSETS\n")
    for u in ranges:
        rngs = ranges[u]
        j = num_columns
        u = unicodedata.normalize('NFKD', u).encode('ascii','ignore')
        u = u.replace(" ", "_")
        node_list = " ".join( map( conv_range, rngs[:j] ) )
        fd.write("  %(setname)s %(setlist)s\n" % { "setname": u, "setlist": node_list })
        while j < len(rngs):
            rng = rngs[j:j+num_columns]
            node_list = " ".join( map( conv_range, rng ) )
            fd.write("\t\t%(setlist)s\n" % { "setname": u, "setlist": node_list })
            j += num_columns
    
    # END
    fd.write("END\n")
    
    output = {'num_layers': layers}
    return(output)
