from . import utils

def overwrite_config(config_pkl, args):
    config = utils.load_pkl(config_pkl)
    if args.norm is not None:
        config['dataset_args']['norm'] = args.norm
    v = vars(args)
    if 'D' in v and args.D is not None:
        config['lattice_args']['D'] = args.D + 1
    if 'l_extent' in v and args.l_extent is not None:
        config['lattice_args']['extent'] = args.l_extent
    for arg in ('qlayers','qdim','zdim','encode_mode','players','pdim','enc_mask','pe_type','pe_dim','domain'):
        if arg == 'pe_dim' and arg not in config['model_args']:
            assert v[arg] is None
            config['model_args']['pe_dim'] = None
            continue # maintain backwards compatibility 
        if v[arg] is not None:
            config['model_args'][arg] = v[arg]
    return config


