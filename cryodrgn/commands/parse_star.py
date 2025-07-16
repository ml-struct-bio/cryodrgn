"""Parse CTF and pose parameters from a RELION .star file"""

import argparse
import os
import logging
from cryodrgn import utils

logger = logging.getLogger(__name__)

def add_args(parser):
    parser.add_argument(
        'input', help='Input RELION .star file containing metadata for parsing')
    parser.add_argument(
        '-o', '--outdir', type=os.path.abspath, default='.', help='Output directory for ctf.pkl and pose.pkl files (default: current directory)')
    parser.add_argument(
        '--overwrite', action='store_true', help='Overwrite existing output files if they exist')

def main(args):
    input_star = args.input
    os.makedirs(args.outdir, exist_ok=True)
    ctf_out = os.path.join(args.outdir, 'ctf.pkl')
    pose_out = os.path.join(args.outdir, 'pose.pkl')

    if not args.overwrite:
        if os.path.exists(ctf_out):
            raise RuntimeError(f"{ctf_out} already exists. Use --overwrite to overwrite.")
        if os.path.exists(pose_out):
            raise RuntimeError(f"{pose_out} already exists. Use --overwrite to overwrite.")

    # Run parse_ctf_star
    logger.info(f"Parsing CTF parameters from {input_star}...")
    cmd = f'cryodrgn parse_ctf_star {input_star} -o {ctf_out}'
    logger.info(f"Running command: {cmd}")
    utils.run_command(cmd)

    # Run parse_pose_star
    logger.info(f"Parsing pose parameters from {input_star}...")
    cmd = f'cryodrgn parse_pose_star {input_star} -o {pose_out}'
    logger.info(f"Running command: {cmd}")
    utils.run_command(cmd)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    args = parser.parse_args()
    main(args)

