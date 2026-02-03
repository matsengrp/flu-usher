import pandas as pd

def get_simplified_host_group(host):
    if pd.isna(host):
        return 'other'

    host = host.lower()

    # Unknown/Environmental
    if host in [
        'host', 'mammals', 'other mammals', 'feces', 'environment', 'water sample', 'other environment',
        'animal', 'surface swab', 'air sample', 'unknown', 'insect'
    ]:
        return 'other'

    # Laboratory derived
    if 'laboratory' in host:
        return 'other'

    # Human
    if host in ['human']:
        return 'human'
    
    # Swine
    if host in ['swine', 'pig']:
        return 'swine'
    if 'sus scrofa' in host or 'sus ' in host:
        return 'swine'
    
    # Equine
    if host in ['equine', 'equus caballus', 'horse', 'donkey']:
        return 'other'

    # Bovine
    if host in ['dairy cow', 'bovine', 'cow', 'cattle']:
        return 'bovine'

    # Canine
    if host in ['canine', 'dog']:
        return 'other'
    if 'canis' in host:
        return 'other'

    # Feline
    if host in ['feline', 'cat']:
        return 'other'
    if 'felis' in host:
        return 'other'

    # Marine mammals
    if host in ['seal', 'dolphin', 'whale', 'sea lion', 'walrus']:
        return 'other'

    # Other mammals
    if host in [
        'ferret', 'mink', 'mouse', 'mus musculus', 'rodent', 'bat', 'meerkat', 'panda', 'camel',
        'primate'
    ]:
        return 'other'

    # ===== AVIAN CLASSIFICATION =====
    # All bird-related hosts should be classified as 'avian'
    
    # Check for any bird-related terms first
    bird_indicators = [
        'avian', 'bird', 'fowl', 'poultry', 'waterfowl',
        'duck', 'goose', 'swan', 'teal', 'mallard', 'pintail', 'wigeon', 'shoveler',
        'chicken', 'turkey', 'quail', 'pheasant', 'partridge', 'grouse', 'peafowl', 'peacock',
        'gull', 'tern', 'eagle', 'hawk', 'falcon', 'kestrel', 'harrier', 'buzzard',
        'crow', 'raven', 'magpie', 'pigeon', 'dove', 'sparrow', 'finch', 'swallow',
        'penguin', 'cormorant', 'pelican', 'heron', 'ibis', 'stork', 'crane', 'bustard',
        'sandpiper', 'turnstone', 'knot', 'dunlin', 'sanderling', 'curlew', 'plover',
        'shearwater', 'albatross', 'petrel', 'grebe', 'coot', 'rail', 'ostrich', 'great tit',
        'passerine'
    ]
    
    for indicator in bird_indicators:
        if indicator in host:
            return 'avian'
    
    # Check scientific names (genus names) for birds
    bird_genera = [
        'gallus', 'meleagris', 'coturnix', 'coturnic', 'numida', 'phasianus', 'phasanius',
        'alectoris', 'pavo', 'perdix', 'cyrtonyx', 'oreortyx', 'francolinus', 'polyplectron',
        'lophura', 'anas', 'anser', 'cygnus', 'branta', 'chen', 'cairina', 'aythya',
        'bucephala', 'tadorna', 'aix', 'mareca', 'melanitta', 'somateria', 'clangula',
        'oxyura', 'mergus', 'lophodytes', 'dendrocygna', 'netta', 'larus', 'sterna',
        'chroicocephalus', 'leucophaeus', 'gelochelidon', 'chlidonias', 'hydroprogne',
        'rissa', 'xema', 'larosterna', 'rynchops', 'calidris', 'arenaria', 'numenius',
        'tringa', 'scolopax', 'catoptrophorus', 'himantopus', 'gallinago', 'falco',
        'buteo', 'accipiter', 'haliaeetus', 'halietus', 'circus', 'parabuteo', 'nisaetus',
        'morphnus', 'necrosyrtes', 'gyps', 'corvus', 'pica', 'columba', 'streptopelia',
        'pygoscelis', 'spheniscus', 'podiceps', 'tachybaptus', 'gallinula', 'anthropoides',
        'passer', 'hirundo', 'zosterops', 'copsychus', 'garrulax', 'gracula', 'stercorarius',
        'falsco', 'anseriformes'
    ]
    
    # Check if the genus matches the first word in the host name
    for genus in bird_genera:
        if genus == host.split()[0]:
            return 'avian'

    # Raise error for unrecognized hosts
    raise ValueError(f"Unrecognized host: {host}")