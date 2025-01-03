import argparse
import subprocess
import os

def add_args(parser):
    parser.add_argument(
        'input_star', help='Input RELION STAR file containing metadata for parsing')
    parser.add_argument(
        '-o', '--output_dir', help='Output directory (default: current directory)')
    parser.add_argument(
        '--overwrite', action='store_true', help='Overwrite existing output files if they exist')

def run_command(command, description):
    try:
        print(f"Running {description}...")
        subprocess.run(command, check=True)
        print(f"{description} completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error while running {description}: {e}")
        raise

def main(args):
    input_star = args.input_star
    output_dir = os.path.abspath(args.output_dir) if args.output_dir else os.getcwd()
    os.makedirs(output_dir, exist_ok=True)
    ctf_output = os.path.join(output_dir, 'ctf.pkl')
    pose_output = os.path.join(output_dir, 'pose.pkl')

    if not args.overwrite:
        if os.path.exists(ctf_output):
            print(f"Error: {ctf_output} already exists. Use --overwrite to overwrite.")
            return
        if os.path.exists(pose_output):
            print(f"Error: {pose_output} already exists. Use --overwrite to overwrite.")
            return

    # Run parse_ctf_star
    run_command([
        'cryodrgn', 'parse_ctf_star', input_star, '-o', ctf_output
    ], "parse_ctf_star")

    # Run parse_pose_star
    run_command([
        'cryodrgn', 'parse_pose_star', input_star, '-o', pose_output
    ], "parse_pose_star")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse RELION STAR file into CTF and pose pkl files')
    add_args(parser)
    args = parser.parse_args()
    main(args)

