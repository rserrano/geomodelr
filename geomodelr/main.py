#! /usr/bin/env python

import argparse
import json
import fileinput
import sys
from model import GeologicalModel

class ParametersException(Exception):
    pass

def prep_model(geojson):
    model = GeologicalModel(json.loads(geojson.read()))
    model.calc_cache()
    if not model.has_interpolation():
        print >> sys.stderr, ">> no match in file, calculating match"
        model.calc_interpolation()
        geojson.close()
        fp = open(geojson.name, "w")
        fp.write(json.dumps(model.geojson))
        print >> sys.stderr, ">> matching created"
    
    return model

def query_coordinates(geojson, grid=False):
    model = prep_model(geojson)

    line = sys.stdin.readline()
    while line:
        point = map(float, line.split())
        model.query_point(point)
        line = sys.stdin.readline()

def query_grid(geojson, grid=False):
    model = prep_model(geojson)

    line = sys.stdin.readline()
    while line:
        args = line.split()
        try:
            if len(args) < 9:
                raise ParametersException("wrong number of parameters for grid")
            
            mnx, mny, mnz, mxx, mxy, mxz = map(float, args[:6])
            nx, ny, nz = map(int, args[6:])
            
            if nx < 1:
                raise ParametersException("nx is not positive")
            if ny < 1:
                raise ParametersException("ny is not positive")
            if nz < 1:
                raise ParametersException("nz is not positive")
            
            dx = (mxx - mnx)/(nx-1) if nx > 1 else 0
            dy = (mxy - mny)/(ny-1) if ny > 1 else 0
            dz = (mxz - mnz)/(nz-1) if nz > 1 else 0
            
            for i in xrange(nz):
                for j in xrange(ny):
                    for k in xrange(nx):
                        point = (mnx + dx*k, mny + dy*j, mnz + dz*i)
                        sys.stdout.write(model.query_point(point).replace(" ", "_") + " ")
                    sys.stdout.flush()
        except ParametersException as e:
            print >> sys.stderr, ">>", e
        except Exception as e:
            raise
        line = sys.stdin.readline()
    

def get_information(geojson, verbose):
    # Show map, cross sections, polygons, etc.
    model = GeologicalModel(json.loads(geojson.read()))
    model.print_information(verbose)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""GeoModelR helps you to manage the versions 
                                                    of your geological models and use them to 
                                                    generate meshes or make any calculation
                                                    you want.""")
    
    group = parser.add_mutually_exclusive_group()
    
    group.add_argument("-q", "--query", const=query_coordinates, action='store_const',
                       help = """ reads lines with 3 coordinates: x y z 
                                  and returns a line per coordinate with 
                                  the formation at that spatial point. """)
    
    group.add_argument("-g", "--grid", const=query_grid, action='store_const',  
                       help = """ grid reads lines with 9 values: mnx mny mnz mxz mxy mxz nz ny nz, 
                                  the minimum and maximum coordinates of the box to query plus
                                  the number of points to query in each direction. 
                                  Returns a line with nx*ny*nz values in the following order: 
                                  f(mnx,mny,mnz) f(mnx+dx,mny,mnz) f(mnx+2*dx,mny,mnz) ... f(mxx,mny,mnz) f(mnx,mny+dy,mnz) ... f(mxx,mxy,mnz) f(mnx,mny,mnz+dz) ... f(mxx,mxy,mxz)
                                  where and dx=(mxx-mnx)/(nx-1), ... and f is the formation at that spatial point """)
    
    group.add_argument("-i", "--info", const=get_information, action='store_const',  
                       help = """ gets information from the geological 
                                  model and presents it to the user. """)
    
    parser.add_argument("-v", "--verbose", action="store_true",
                        help = """shows more information to the user.""")
    
    parser.add_argument("-p", "--profile", action="store_true",
                        help = """profiles geomodelr.""")
    
    parser.add_argument("model", nargs="?", type=argparse.FileType('r'))
    args = parser.parse_args()
    
    # Profile GeoModelR
    if args.profile:
        import cProfile, pstats, StringIO
        pr = cProfile.Profile()
        pr.enable()
    
    # Process.
    if args.query:
        args.query(args.model, args.verbose)
    elif args.grid:
        args.grid(args.model, args.verbose)
    elif args.info:
        args.info(args.model, args.verbose)
    else:
        parser.print_help()
    
    # Print profiling of GeoModelR.
    if args.profile:
        pr.disable()
        s = StringIO.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print >> sys.stderr, s.getvalue()

