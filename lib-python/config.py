import utils

def load_config(config_pkl, args):
    config = utils.load_pkl(config_pkl)
    if args.norm is None:
        args.norm = config['dataset_args']['norm']
    if args.D is None:
        args.D = config['lattice_args']['D'] - 1
    if args.l_extent is None:
        args.l_extent = config['lattice_args']['extent']
    v = vars(args)
    for arg in ('qlayers','qdim','zdim','encode_mode','players','pdim','enc_mask','pe_type','domain'):
        if v[arg] is None:
            v[arg] = config['model_args'][arg]
    return args


