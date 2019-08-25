# coding=utf

# Geomodelr query tool. Tools for using a geomodelr.com geological model.
# Copyright (C) 2016 Geomodelr, Inc.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function
from __future__ import divison

from .shared import ModelException

def get_feature_fault_names( cs ):
    faults = set()
    for feature in cs['features']:
        if feature['geometry']['type'] == 'LineString':
            if 'name' in feature['properties']:
                faults.add(feature['properties']['name'])
    return faults

def get_all_fault_names( model ):
    faults = set()
    for idx, feature in enumerate(model['features']):
        faults |= set(get_feature_fault_names( feature ))
    return list(faults)

def to_version_0_1_6( model ):
    """
    This is the first version, that has a versioning number.
    """
    fnames = get_all_fault_names( model )
    model['properties']['lines'] = { n: "FAULT" for n in fnames }

def to_version_0_1_13( model ):
    """
    This version changes a units dict from name: color to name: { color: color, pattern: pattern }
    """
    units = model['properties']['units']
    for unit in units:
        units[unit] = { 'color': units[unit], 'pattern': None }

UPGRADES = [ ( (0, 1, 6), to_version_0_1_6 ), ( ( 0, 1, 13), to_version_0_1_13 ) ]

def upgrade_model( model ):
    if not 'geomodelr' in model['properties']:
        version = "0.1.5"
    else:
        version = model['properties']['geomodelr']
    try:
        version = tuple(map(int, version.split(".")))
    except:
        raise ModelException("Error with geomodelr version %s" % version )
    for up, fun in UPGRADES:
        if version < up:
            fun(model)

