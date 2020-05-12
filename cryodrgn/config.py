from . import utils

def load_config(config_pkl, args):
    config = utils.load_pkl(config_pkl)
    if args.norm is None:
        args.norm = config['dataset_args']['norm']
    v = vars(args)
    if 'D' in v and args.D is None:
        args.D = config['lattice_args']['D'] - 1
    if 'l_extent' in v and args.l_extent is None:
        args.l_extent = config['lattice_args']['extent']
    for arg in ('qlayers','qdim','zdim','encode_mode','players','pdim','enc_mask','pe_type','pe_dim','domain'):
        if arg == 'pe_dim' and arg not in config['model_args']:
            assert v[arg] is None
            continue # maintain backwards compatibility 
        if v[arg] is None:
            v[arg] = config['model_args'][arg]
    return args


