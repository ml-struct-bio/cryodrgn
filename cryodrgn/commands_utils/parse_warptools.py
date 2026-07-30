"""Parse WarpTools subtomogram STAR files into cryoDRGN-compatible 2D rows.

WarpTools/RELION5 2D stacks are often CTF-premultiplied
(``rlnCtfDataAreCtfPremultiplied``); disable CTF correction in cryoDRGN when
that flag is set.

Example usage
-------------
cryodrgn_utils parse_warptools -t tomograms.star -p particles.star \\
    --tilt-dim 5760 4092 -o particles_2d.star

"""
import argparse
from ast import literal_eval

import numpy as np
import pandas as pd
import starfile
from scipy.spatial.transform import Rotation


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-t",
        "--tomograms",
        required=True,
        help="Path to tomograms.star file.",
    )
    parser.add_argument(
        "-p",
        "--particles",
        required=True,
        help="Path to main particles.star file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="particles_2d.star",
        help="Output name for expanded 2D star file (default: %(default)s)",
    )
    parser.add_argument(
        "--tilt-dim",
        nargs=2,
        type=int,
        required=True,
        help="Tilt image dimensions in pixels, e.g. --tilt-dim 5760 4092",
    )
    parser.add_argument(
        "--image-pixel-size",
        type=float,
        help=(
            "Pixel size of the cropped particle images in A/px. Defaults to "
            "rlnImagePixelSize from the particles optics table."
        ),
    )
    parser.add_argument(
        "--image-size",
        type=int,
        help=(
            "Box size of the cropped particle images in pixels. Defaults to "
            "rlnImageSize from the particles optics table."
        ),
    )


def _get_optics_value(optics_row, field, override=None, option_name=None):
    if override is not None:
        return override

    value = optics_row.get(field, None)
    if value is None or pd.isna(value):
        if option_name is None:
            raise ValueError(f"Could not determine {field} from the optics table.")
        raise ValueError(
            f"Could not determine {field}; provide {option_name} explicitly."
        )

    return value


def _get_tomogram_pixel_size(global_row, optics_row):
    pixel_size_ang = global_row.get("rlnTomoTiltSeriesPixelSize", None)
    if pixel_size_ang is None or pd.isna(pixel_size_ang):
        pixel_size_ang = _get_optics_value(optics_row, "rlnTomoTiltSeriesPixelSize")

    return float(pixel_size_ang)


def _parse_star_list(value):
    """Parse a STAR list cell that may already be a list/array or a string."""
    if isinstance(value, str):
        return literal_eval(value)
    return list(value)


def _warp_euler_matrix_from_row(row, fields):
    angles = [float(row.get(field, 0.0)) for field in fields]
    return Rotation.from_euler("ZYZ", angles, degrees=True).as_matrix()


def _particle_orientation_matrix(row):
    return _warp_euler_matrix_from_row(
        row,
        ("rlnAngleRot", "rlnAngleTilt", "rlnAnglePsi"),
    )


class Tomogram:
    """
    Represents one tilt-series' geometry & defocus data. Precomputes a 4x4 transform
    (rotation + translation) that maps tomogram coordinates to each 2D tilt image.
    """

    def __init__(
        self,
        tilt_image_dims,
        pixel_size,
        defocus_u_array,
        defocus_v_array,
        defocus_angle_array,
        hand,
        projection_matrices,
    ):
        self.tilt_image_dims = tilt_image_dims
        self.pixel_size = pixel_size
        self.defocus_u_array = defocus_u_array
        self.defocus_v_array = defocus_v_array
        self.defocus_angle_array = defocus_angle_array
        self.hand = hand

        # Pose and coordinate projection use the WarpTools projection matrices
        # directly. Handedness only changes the defocus Z correction below.
        self.projection_matrices = {
            i: m.astype(float) for i, m in enumerate(projection_matrices)
        }
        self.rotation_matrices = {
            i: m.astype(float)[:3, :3] for i, m in enumerate(projection_matrices)
        }
        self.n_tilts = len(self.projection_matrices)

    def project_point(self, point_3d_px, i_tilt):
        """
        Project centered tomogram-pixel coordinates into 2D tilt-image pixel coordinates.
        """
        pt_homog = np.append(point_3d_px, 1.0)
        M = self.projection_matrices[i_tilt]

        xy = (M @ pt_homog)[:2]

        # WarpTools projection matrices do not include the 2D image center.
        xy += np.array(
            [
                self.tilt_image_dims[0] / 2.0,
                self.tilt_image_dims[1] / 2.0,
            ]
        )

        return xy

    def calculate_local_defocus_uv(self, i_tilt, point_3d_px):
        """
        Compute local defocus by projecting centered tomogram-pixel coordinates
        along the beam/depth direction, then converting depth from pixels to Angstrom.
        """
        # The depth along the beam is the third row of the same projection
        # rotation used for the 2D coordinates; `hand` flips it for tomograms
        # reconstructed with the opposite handedness (matches parse_relion).
        depth_offset_px = (self.rotation_matrices[i_tilt] @ point_3d_px.astype(float))[
            2
        ]
        depth_offset_ang = depth_offset_px * self.pixel_size * self.hand

        loc_u = self.defocus_u_array[i_tilt] + depth_offset_ang
        loc_v = self.defocus_v_array[i_tilt] + depth_offset_ang
        loc_angle = self.defocus_angle_array[i_tilt]

        return loc_u, loc_v, loc_angle

    def expand_particle_to_2drows(
        self,
        point_3d_px,
        original_image_name,
        tilt_series_df,
        tomo_name,
        group_name,
        base_orientation_matrix=None,
    ):

        coords_2d = []
        defocusU_list = []
        defocusV_list = []
        defocusAngle_list = []
        final_zyz_list = []

        for i in range(self.n_tilts):
            coords_2d.append(self.project_point(point_3d_px, i))
            lu, lv, la = self.calculate_local_defocus_uv(i, point_3d_px)
            defocusV_list.append(lv)
            defocusU_list.append(lu)
            defocusAngle_list.append(la)

            if base_orientation_matrix is not None:
                # Match parse_relion / WarpTools expanded GT: R_particle @ R_tilt,
                # exported as RELION ZYZ (not cryoDRGN-internal R_to_relion_scipy).
                final_matrix = base_orientation_matrix @ self.rotation_matrices[i]
                a, b, c = Rotation.from_matrix(final_matrix).as_euler(
                    "ZYZ", degrees=True
                )
                final_zyz_list.append([a, b, c])

        coords_2d = np.array(coords_2d)
        final_zyz_list = np.array(final_zyz_list)

        # CTF scale factor from tilt_series_df or estimate from tilt angles.
        if "rlnCtfScalefactor" in tilt_series_df.columns:
            ctf_scale = tilt_series_df["rlnCtfScalefactor"].to_numpy()
        elif "rlnTomoYTilt" in tilt_series_df.columns:
            ctf_scale = np.cos(np.deg2rad(tilt_series_df["rlnTomoYTilt"].to_numpy()))
        else:
            ctf_scale = np.ones(len(tilt_series_df))

        # Make 2D dataframe
        n_tilts = self.n_tilts
        image_names_2d = [f"{i+1:06d}@{original_image_name}" for i in range(n_tilts)]

        if "rlnMicrographName" in tilt_series_df.columns:
            mic_names = tilt_series_df["rlnMicrographName"].values
        else:
            mic_names = np.array(
                [f"{tomo_name}_tilt_{i+1:06d}" for i in range(self.n_tilts)]
            )

        df_2d = pd.DataFrame(
            {
                "rlnMagnification": 10000.0,  # placeholder
                "rlnDefocusU": defocusU_list,
                "rlnDefocusV": defocusV_list,
                "rlnDefocusAngle": defocusAngle_list,
                "rlnImageName": image_names_2d,
                "rlnMicrographName": mic_names,
                "rlnCoordinateX": coords_2d[:, 0],
                "rlnCoordinateY": coords_2d[:, 1],
                "rlnCtfBfactor": 0.0,  # placeholder - possibly remove
                "rlnCtfScalefactor": ctf_scale,
                "rlnGroupName": group_name,
                "rlnTiltName": mic_names,
            }
        )
        if "rlnMicrographPreExposure" in tilt_series_df.columns:
            df_2d["rlnMicrographPreExposure"] = tilt_series_df[
                "rlnMicrographPreExposure"
            ].values
        if "rlnTomoYTilt" in tilt_series_df.columns:
            df_2d["rlnTomoYTilt"] = tilt_series_df["rlnTomoYTilt"].values

        if base_orientation_matrix is not None and len(final_zyz_list) > 0:
            df_2d["rlnAngleRot"] = final_zyz_list[:, 0]
            df_2d["rlnAngleTilt"] = final_zyz_list[:, 1]
            df_2d["rlnAnglePsi"] = final_zyz_list[:, 2]
        else:
            df_2d["rlnAngleRot"] = 0.0
            df_2d["rlnAngleTilt"] = 0.0
            df_2d["rlnAnglePsi"] = 0.0

        # Sort only this particle's expanded tilt-image rows. The caller appends
        # these per-particle frames in particle order, matching parse_relion.
        if "rlnMicrographPreExposure" in df_2d.columns:
            df_2d = df_2d.sort_values(
                "rlnMicrographPreExposure", ascending=True, kind="stable"
            ).reset_index(drop=True)

        return df_2d


def _projection_matrices_from_df(ts_df: pd.DataFrame) -> list[np.ndarray]:
    """Parse WarpTools ProjX/Y/Z/W list columns into 4x4 matrices (once per table)."""
    pX = ts_df["rlnTomoProjX"].apply(_parse_star_list)
    pY = ts_df["rlnTomoProjY"].apply(_parse_star_list)
    pZ = ts_df["rlnTomoProjZ"].apply(_parse_star_list)
    pW = ts_df["rlnTomoProjW"].apply(_parse_star_list)
    return [np.vstack([x, y, z, w]) for x, y, z, w in zip(pX, pY, pZ, pW)]


def main(args: argparse.Namespace) -> None:

    # Load star files
    tomo_star = starfile.read(args.tomograms, always_dict=True)
    particles_star = starfile.read(args.particles, always_dict=True)

    global_df = tomo_star["global"]
    particles_df = particles_star["particles"]
    optics_df = particles_star["optics"]

    # Index global / optics rows once; Warp Proj parsing dominates runtime and is
    # identical for all particles that share a tomogram (+ visible-frame mask).
    global_by_tomo = {row["rlnTomoName"]: row for _, row in global_df.iterrows()}
    optics_by_group = {row["rlnOpticsGroup"]: row for _, row in optics_df.iterrows()}
    optics_resolved = {}
    for group, optics_row in optics_by_group.items():
        optics_resolved[group] = {
            "row": optics_row,
            "image_pixel_size": float(
                _get_optics_value(
                    optics_row,
                    "rlnImagePixelSize",
                    args.image_pixel_size,
                    "--image-pixel-size",
                )
            ),
            "image_size": int(
                _get_optics_value(
                    optics_row,
                    "rlnImageSize",
                    args.image_size,
                    "--image-size",
                )
            ),
        }

    tomo_meta_cache = {}
    tomogram_cache = {}
    has_visible_frames = "rlnTomoVisibleFrames" in particles_df.columns

    all_2d_rows = []
    ps = 0

    for idx, row in particles_df.iterrows():
        tomo_name = row["rlnTomoName"]
        optics_group = row["rlnOpticsGroup"]
        optics_info = optics_resolved[optics_group]
        optics_row = optics_info["row"]

        if tomo_name not in tomo_meta_cache:
            global_row = global_by_tomo[tomo_name]
            ts_df = tomo_star[tomo_name]
            pixel_size_ang = _get_tomogram_pixel_size(global_row, optics_row)
            handedness = global_row.get("rlnTomoHand", 1)
            tomo_meta_cache[tomo_name] = {
                "global_row": global_row,
                "ts_df": ts_df,
                "proj_mats": _projection_matrices_from_df(ts_df),
                "pixel_size_ang": pixel_size_ang,
                "hand": -1 if float(handedness) == -1 else 1,
                "tomo_center_px": np.array(
                    [
                        global_row["rlnTomoSizeX"] / 2.0,
                        global_row["rlnTomoSizeY"] / 2.0,
                        global_row["rlnTomoSizeZ"] / 2.0,
                    ],
                    dtype=float,
                ),
            }

        meta = tomo_meta_cache[tomo_name]

        if has_visible_frames:
            vis = _parse_star_list(row["rlnTomoVisibleFrames"])
            visible_indices = tuple(i for i, val in enumerate(vis) if int(val) == 1)
        else:
            visible_indices = None

        cache_key = (tomo_name, visible_indices)
        if cache_key not in tomogram_cache:
            if visible_indices is None:
                sub_ts_df = meta["ts_df"]
                proj_mats = meta["proj_mats"]
            else:
                sub_ts_df = meta["ts_df"].iloc[list(visible_indices)]
                proj_mats = [meta["proj_mats"][i] for i in visible_indices]

            tomogram_cache[cache_key] = (
                Tomogram(
                    tilt_image_dims=args.tilt_dim,
                    pixel_size=meta["pixel_size_ang"],
                    defocus_u_array=sub_ts_df["rlnDefocusU"].to_numpy(),
                    defocus_v_array=sub_ts_df["rlnDefocusV"].to_numpy(),
                    defocus_angle_array=sub_ts_df["rlnDefocusAngle"].to_numpy(),
                    hand=meta["hand"],
                    projection_matrices=proj_mats,
                ),
                sub_ts_df,
            )

        tomogram, sub_ts_df = tomogram_cache[cache_key]

        particle_orientation = _particle_orientation_matrix(row)

        point_px = np.array(
            [
                row["rlnCoordinateX"],
                row["rlnCoordinateY"],
                row["rlnCoordinateZ"],
            ],
            dtype=float,
        )

        origin_px = (
            np.array(
                [
                    row.get("rlnOriginXAngst", 0.0),
                    row.get("rlnOriginYAngst", 0.0),
                    row.get("rlnOriginZAngst", 0.0),
                ],
                dtype=float,
            )
            / meta["pixel_size_ang"]
        )

        point_3d = point_px - meta["tomo_center_px"] - origin_px

        df_2d = tomogram.expand_particle_to_2drows(
            point_3d_px=point_3d,
            original_image_name=row["rlnImageName"],
            tilt_series_df=sub_ts_df,
            tomo_name=tomo_name,
            group_name=row["rlnTomoParticleName"],
            base_orientation_matrix=particle_orientation,
        )

        image_pixel_size = optics_info["image_pixel_size"]
        image_size = optics_info["image_size"]

        df_2d["rlnOriginalParticle"] = idx + 1
        df_2d["rlnRandomSubset"] = row.get("rlnRandomSubset", 1)
        df_2d["rlnImagePixelSize"] = image_pixel_size
        df_2d["rlnImageSize"] = image_size
        df_2d["rlnDetectorPixelSize"] = image_pixel_size
        df_2d["rlnVoltage"] = optics_row["rlnVoltage"]
        df_2d["rlnSphericalAberration"] = optics_row["rlnSphericalAberration"]
        df_2d["rlnAmplitudeContrast"] = optics_row["rlnAmplitudeContrast"]
        df_2d["rlnPhaseShift"] = ps
        df_2d["rlnOriginX"] = 0
        df_2d["rlnOriginY"] = 0

        all_2d_rows.append(df_2d)

    final_2d_df = pd.concat(all_2d_rows, ignore_index=True)
    print(f"Total 2D rows in output: {len(final_2d_df)}")

    out_star_dict = {"": final_2d_df}
    starfile.write(out_star_dict, args.output, overwrite=True)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
