from shared import ModelException

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
        if feature['geology_type'] == 'section' and 'interpolation' in feature['properties'] and feature['properties']['interpolation']:
            faults |= set(get_feature_fault_names( feature ))
    return list(faults)

def to_version_0_1_6( model ):
    """
    This is the first version, that has a versioning number.
    """
    fnames = get_all_fault_names( model )
    model['properties']['lines'] = { n: "FAULT" for n in fnames }

UPGRADES = [ ( (0, 1, 6), to_version_0_1_6 ) ]

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

