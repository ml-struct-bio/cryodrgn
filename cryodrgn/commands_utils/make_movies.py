"""Make MP4 movies of .mrc volumes produced by cryodrgn analyze* commands.

You must install ChimeraX under the alias `chimerax` before running this command, see:
https://www.cgl.ucsf.edu/chimerax/download.html

Example usage
-------------
# Latent k-means and PCA movies
$ cryodrgn_utils make_movies spr_runs/07/out 19 latent --iso=31 --camera="-0.03377,-0.97371,0.22528,545.75,0.89245,-0.13085,-0.43175,87.21,0.44988,0.18647,0.87341,1039"

# Volume PCA movies
$ cryodrgn_utils make_movies spr_runs/07/out 19 volume --name=front --iso=210 --camera="0.12868,-0.9576,0.25778,95.4,-0.85972,-0.23728,-0.45231,15.356,0.4943,-0.16341,-0.8538,-33.755"

"""
import argparse
import os
from datetime import datetime as dt
import logging
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)


def generate_movie_prologue(width: int, height: int) -> list[str]:
    return [
        "set bgColor white",
        "graphics silhouettes true",
        f"windowsize {width} {height}",
    ]


def generate_movie_epilogue(
    cam_matrix, iso_threshold, num_vols, directory, name, all_cornfl, framerate
):

    epilogue = list()
    if all_cornfl:
        epilogue = ["vol all color cornfl"]

    epilogue += [f"view matrix camera {cam_matrix}"]

    if iso_threshold:
        epilogue += [f"vol all level {iso_threshold}"]

    epilogue += [
        "surface dust all size 10",
        "lighting soft",
        "mov record",
        "mseries all pause 1 step 1",
        f"wait {num_vols}",
        f"mov encode {directory}/{name}.mp4 framerate {framerate}",
        "exit",
    ]

    return epilogue


def add_args(parser: argparse.ArgumentParser) -> None:
    """The command-line arguments for use with `cryodrgn_utils make_movies`."""

    parser.add_argument(
        "workdir", type=os.path.abspath, help="Directory with cryoDRGN results"
    )
    parser.add_argument(
        "epoch",
        type=str,
        help="Epoch number N to analyze "
        "(1-based indexing, corresponding to z.N.pkl, weights.N.pkl)",
    )
    parser.add_argument(
        "type",
        type=str,
        help="Analysis type to generate movies for ('latent' or 'volume')",
    )
    parser.add_argument(
        "--camera",
        required=True,
        type=str,
        help="Camera matrix string for the movies",
    )
    parser.add_argument(
        "--iso",
        type=str,
        help="Isosurface threshold for the movies (default: ChimeraX default level)",
    )
    parser.add_argument(
        "--frame",
        type=int,
        help="Frame rate (fps) for the movies (default: 3)",
    )
    parser.add_argument(
        "--name",
        type=str,
        help="Video name (default: 'movie')",
    )
    parser.add_argument(
        "--width",
        type=int,
        help="Video width in pixels (default: 600)",
    )
    parser.add_argument(
        "--height",
        type=int,
        help="Video height in pixels (default: 800)",
    )
    parser.add_argument(
        "--analysis-dir",
        type=os.path.abspath,
        help="Latent space analysis directory (default: [workdir]/analyze.[epoch])",
    )
    parser.add_argument(
        "--landscape-dir",
        type=os.path.abspath,
        help="Landscape analysis directory (default: [workdir]/landscape.[epoch])",
    )


def check_chimerax_installation() -> bool:
    """Checks if ChimeraX is installed."""

    command = "chimerax"
    try:
        subprocess.run(
            [command, "--version"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        logger.info(f"{command} is installed.")
        return True
    except FileNotFoundError:
        logger.info(f"{command} is not installed, aborting.")
        return False
    except subprocess.CalledProcessError:
        logger.info(
            f"{command} is installed, but there was an issue executing it. Aborting."
        )
        return False


def find_subdirs(directory: str, keyword: str) -> list[Path]:
    directory_path = Path(directory)
    subdirs = [
        p for p in directory_path.iterdir() if p.is_dir() and p.name.startswith(keyword)
    ]

    values = [p.name.split(keyword)[-1] for p in subdirs]
    logger.info(f"{len(subdirs)} {keyword} directories were found: {', '.join(values)}")

    return subdirs


def get_vols(directory: Path, postfix_regex: str = "") -> list[str]:
    """Postfix regex string enables us to find the state_m_mean.mrc volumes"""
    vol_files = sorted(directory.glob(f"*{postfix_regex}.mrc"))
    logger.info(f"{len(vol_files)} {postfix_regex} volumes were found in {directory}")
    return [f"open {vol_file}" for vol_file in vol_files]


def record_movie(
    dir_list: list[Path],
    iso: str,
    cam_matrix: str,
    name: str,
    prologue: list[str],
    frame_rate: int,
    all_cornfl: bool = False,
    vol_postfix_regex: str = "",
):
    """Movie recording subroutine"""

    for directory in dir_list:

        vols = get_vols(directory, vol_postfix_regex)

        epilogue = generate_movie_epilogue(
            cam_matrix, iso, len(vols), directory, name, all_cornfl, frame_rate
        )

        movie_commands = prologue + vols + epilogue
        movie_script = "\n".join(movie_commands)

        script_path = f"{directory}/temp.cxc"

        with open(script_path, "w") as file:
            file.write(movie_script)

        logger.info("Running chimerax subprocess for movie making.")
        subprocess.run(
            ["chimerax", "--offscreen", script_path],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        logger.info(f"Movie saved in {directory}/{name}.mp4")
        os.remove(script_path)


def latent_movies(
    analysis_dir: str,
    iso: str,
    cam_matrix: str,
    name: str,
    width: int,
    height: int,
    frame_rate: int,
):
    """Record the movies for latent space analysis"""

    prologue = generate_movie_prologue(width, height)

    # 1) kmeans movies
    kmeans_dirs = find_subdirs(analysis_dir, "kmeans")
    record_movie(kmeans_dirs, iso, cam_matrix, name, prologue, frame_rate, False)

    # 2) pc movies
    pc_dirs = find_subdirs(analysis_dir, "pc")
    record_movie(pc_dirs, iso, cam_matrix, name, prologue, frame_rate, True)


def landscape_movies(
    landscape_dir: str,
    iso: str,
    cam_matrix: str,
    name: str,
    width: int,
    height: int,
    frame_rate: int,
):
    """Record the movies for volume space landscape analysis"""

    prologue = generate_movie_prologue(width, height)

    # 1) clustering movies
    clustering_dirs = find_subdirs(f"{landscape_dir}", "clustering")
    record_movie(
        clustering_dirs,
        iso,
        cam_matrix,
        name,
        prologue,
        frame_rate,
        False,
        vol_postfix_regex="mean",
    )

    # 2) pc movies
    vol_pc_dirs = find_subdirs(f"{landscape_dir}/vol_pcs", "pc")
    record_movie(vol_pc_dirs, iso, cam_matrix, name, prologue, frame_rate, True)


def main(args: argparse.Namespace) -> None:
    t1 = dt.now()
    logger.info(args)

    # parsing args
    E = args.epoch
    workdir = args.workdir
    analysis_type = args.type
    cam_matrix = args.camera
    iso = args.iso
    width = 600 if args.width is None else args.width
    name = "movie" if args.name is None else args.name
    height = 800 if args.height is None else args.height
    frame_rate = 3 if args.frame is None else args.frame

    # checking chimerax
    if not check_chimerax_installation():
        return

    if analysis_type != "latent" and analysis_type != "volume":
        logger.info("Analysis type unrecognized, aborting.")
        return

    analysis_dir = (
        f"{workdir}/analyze.{E}" if args.analysis_dir is None else args.analysis_dir
    )
    landscape_dir = (
        f"{workdir}/landscape.{E}" if args.landscape_dir is None else args.landscape_dir
    )

    if analysis_type == "latent":
        logger.info(f"Working in {analysis_dir}")
        latent_movies(analysis_dir, iso, cam_matrix, name, width, height, frame_rate)
    else:
        logger.info(f"Working in {landscape_dir}")
        landscape_movies(
            landscape_dir, iso, cam_matrix, name, width, height, frame_rate
        )

    td = dt.now() - t1
    logger.info(f"Finished in {td}")
