#! /usr/bin/env python

"""
	Geomodelr query tool. Tools for using a geomodelr.com geological model.
	Copyright (C) 2016 Geomodelr, Inc.
	
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU Affero General Public License as
	published by the Free Software Foundation, either version 3 of the
	License, or (at your option) any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU Affero General Public License for more details.
	
	You should have received a copy of the GNU Affero General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import json
import fileinput
import sys
import shared
from model import GeologicalModel

class ParametersException(Exception):
    pass

def prep_model(geojson, save=False):
    model = GeologicalModel(json.loads(geojson.read()))
    changed = False
    
    if changed and save:
        geojson.close()
        fp = open(geojson.name, "w")
        fp.write(json.dumps(model.geojson))
    
    return model

def query_coordinates(geojson, verbose=False):
    model = prep_model(geojson)
    line  = sys.stdin.readline()
    while line:
        try:
            point = map(float, line.split())
        except ValueError:
            raise ParametersException("three numerical values are required per line")
        sys.stdout.write(shared.force_encode(model.closest(point)[0]).replace(" ", "_") + "\n")
        line = sys.stdin.readline()

def query_grid(geojson, verbose=False):
    model = prep_model(geojson)
    line = sys.stdin.readline()
    while line:
        args = line.split()
        try:
            if len(args) < 9:
                raise ParametersException("wrong number of parameters for grid")
            try:
                mnx, mny, mnz, mxx, mxy, mxz = map(float, args[:6])
                nx, ny, nz = map(int, args[6:])
            except ValueError:
                raise ParametersException("nine numerical values are required per line")
            
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
                        sys.stdout.write(shared.force_encode(model.closest(point)[0]).replace(" ", "_") + " ")
                    sys.stdout.flush()
        except ParametersException as e:
            print >> sys.stderr, ">>", e
        except Exception as e:
            raise
        line = sys.stdin.readline()

def intersect_plane( geojson, verbose=False ):
    model = prep_model(geojson)
    line = sys.stdin.readline()
    while line:
        args = line.split()
        try:
            if len(args) < 12:
                raise ParametersException("wrong number of parameters for grid")
            try:
                plane = []
                for i in range(4):
                    point = []
                    for j in range(3):
                        point.append(float(args[i*3+j]))
                    plane.append(point)
            except ValueError:
                raise ParametersException("nine numerical values are required per line")
            print model.intersect_plane(plane)
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
    
    group.add_argument("-ip", "--intersect_plane", const=intersect_plane, action="store_const",
                       help = """ intersects a plane with the faults of the model and 
                                  returns the lines using the coordinate system dictated by the plane. """)
    
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
    elif args.intersect_plane:
        args.intersect_plane(args.model, args.verbose)
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

