import argparse
import numpy as np
import pandas as pd
import starfile
from ast import literal_eval
from scipy.spatial.transform import Rotation as R


"""
Parse 3D particle, tomogram, and tilt series metadata from RELION v5 .star files to calculate 2D particle (subtilt image) geometry.

Example:
cryodrgn_utils parse_relion -t tomograms.star -p particles.3D.star -D 5760 4092 -o particles.2D.star
"""


# minimum required columns for tomograms .star file
TOMOGRAMS_REQUIRED_COLUMNS = {
    "rlnTomoName",
    "rlnTomoHand",
    "rlnTomoTiltSeriesStarFile",
    "rlnTomoTiltSeriesPixelSize",
}

# minimum required columns for 3D particles .star file
PARTICLES_3D_REQUIRED_COLUMNS = {
    "rlnImageName",
    "rlnTomoName",
    "rlnTomoParticleName",
    "rlnTomoVisibleFrames",
    "rlnAngleRot",
    "rlnAngleTilt",
    "rlnAnglePsi",
    "rlnCenteredCoordinateXAngst",
    "rlnCenteredCoordinateYAngst",
    "rlnCenteredCoordinateZAngst",
    # "rlnTomoSubtomogramRot", # optional
    # "rlnTomoSubtomogramTilt", # optional
    # "rlnTomoSubtomogramPsi", # optional
    "rlnVoltage",
    "rlnSphericalAberration",
    "rlnAmplitudeContrast",
    "rlnImageSize",
    "rlnImagePixelSize",
}

# minimum required columns for tilt series .star file
TILT_SERIES_REQUIRED_COLUMNS = {
    "rlnMicrographName",
    "rlnTomoXTilt",
    "rlnTomoYTilt",
    "rlnTomoZRot",
    "rlnTomoXShiftAngst",
    "rlnTomoYShiftAngst",
    "rlnMicrographPreExposure",
    # "rlnCtfScalefactor" # optional
    "rlnDefocusU",
    "rlnDefocusV",
    "rlnDefocusAngle",
    "rlnPhaseShift",
}

# output columns for 2D particles .star file (ordered)
PARTICLES_2D_OUTPUT_COLUMNS = [
    "_image_index",
    "_particle_index",
    "_tilt_index",
    "rlnImageName",
    "rlnGroupName", # equivalent to rlnTomoParticleName
    "rlnCoordinateX",
    "rlnCoordinateY",
    "rlnAngleRot",
    "rlnAngleTilt",
    "rlnAnglePsi",
    "rlnMicrographPreExposure",
    "rlnCtfScalefactor",
    "rlnDefocusU",
    "rlnDefocusV",
    "rlnDefocusAngle",
    "rlnPhaseShift",
    "rlnVoltage",
    "rlnSphericalAberration",
    "rlnAmplitudeContrast",
    "rlnImageSize",
    "rlnImagePixelSize",
]


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-t", "--tomograms",
        required=True,
        help="Path to the input tomograms .star file.",
    )
    parser.add_argument(
        "-p", "--particles",
        required=True,
        help="Path to the input 3D particles .star file.",
    )
    parser.add_argument(
        "-D", "--tilt-dims",
        nargs=2,
        type=int,
        metavar=("WIDTH", "HEIGHT"),
        required=True,
        help="Tilt image dimensions in pixels (example: 5760 4092).",
    )
    parser.add_argument(
        "-F", "--full-data",
        action="store_true",
        help="Include all intermediate metadata columns in output 2D particles .star file. By default, only necessary data columns are included.",
    )
    parser.add_argument(
        "-o", "--output",
        default="particles.2D.star",
        help="Path to the output 2D particles .star file (default: %(default)s).",
    )


class Tomogram:
    """
    Single tilt series with geometry and defocus data.
    Precomputes a 4x4 transform (rotation + translation) that maps tomogram coordinates to each 2D tilt image.
    """

    def __init__(
        self,
        tilt_image_dims,
        pixel_size,
        defocus_u_arr,
        defocus_v_arr,
        defocus_angle_arr,
        tomo_x_tilt_arr,
        tomo_y_tilt_arr,
        tomo_z_rot_arr,
        tomo_x_shift_angst_arr,
        tomo_y_shift_angst_arr,
        tomo_hand,
    ):
        self.tilt_image_dims = tilt_image_dims
        self.pixel_size = pixel_size
        self.defocus_u_arr = defocus_u_arr
        self.defocus_v_arr = defocus_v_arr
        self.defocus_angle_arr = defocus_angle_arr
        self.tomo_x_tilt_arr = tomo_x_tilt_arr
        self.tomo_y_tilt_arr = tomo_y_tilt_arr
        self.tomo_z_rot_arr = tomo_z_rot_arr
        self.tomo_x_shift_angst_arr = tomo_x_shift_angst_arr
        self.tomo_y_shift_angst_arr = tomo_y_shift_angst_arr
        self.tomo_hand = tomo_hand

        self.n_tilts = len(self.tomo_x_tilt_arr)

        self.projection_matrices = {}
        self._build_projection_matrices()

    def _translation_matrix(self, shift_3D_pix):
        """
        Build 4x4 translation matrix for "shift_3D_pix" in pixels.
        """
        mat = np.eye(4)
        mat[:3, 3] = shift_3D_pix
        return mat

    def _rotation_matrix(self, axis, angle_deg):
        """
        Build a 4x4 rotation about "axis" by "angle_deg" degrees.
        """
        rot = R.from_rotvec(np.deg2rad(angle_deg) * np.array(axis))
        mat = np.eye(4)
        mat[:3, :3] = rot.as_matrix()
        return mat

    def _build_projection_matrices(self):
        global_center_pix = np.zeros(3)

        for i in range(self.n_tilts):
            s_0 = self._translation_matrix(-1*global_center_pix)
            # print('s_0:', s_0)

            r_0 = self._rotation_matrix([1, 0, 0], self.tomo_x_tilt_arr[i])
            # print('r_0:', r_0)

            r_1 = self._rotation_matrix([0, 1, 0], self.tomo_y_tilt_arr[i])
            # print('r_1:', r_1)

            r_2 = self._rotation_matrix([0, 0, 1], self.tomo_z_rot_arr[i])
            # print('r_2:', r_2)

            shift_3D_pix = np.array([
                self.tomo_x_shift_angst_arr[i] / self.pixel_size,
                self.tomo_y_shift_angst_arr[i] / self.pixel_size,
                0.0
            ])
            s_1 = self._translation_matrix(shift_3D_pix)
            # print('s_1:', s_1)

            tilt_center_pix = np.array([
                self.tilt_image_dims[0] / 2.0,
                self.tilt_image_dims[1] / 2.0,
                0.0
            ])
            s_2 = self._translation_matrix(tilt_center_pix)
            # print('s_2:', s_2)

            R_zyx = (
                R.from_matrix(r_2[:3, :3])
                * R.from_matrix(r_1[:3, :3])
                * R.from_matrix(r_0[:3, :3])
            )
            # print('R_zyx:', R_zyx.as_matrix())

            R_zyx_inv = np.eye(4)
            R_zyx_inv[:3, :3] = R_zyx.inv().as_matrix()
            # print('R_zyx_inv:', R_zyx_inv)

            M = s_2 @ s_1 @ R_zyx_inv @ s_0
            # print('M:', M)

            # if self.tomo_hand == -1:
            #     flip_z_4x4 = np.eye(4, dtype=np.float32)
            #     flip_z_4x4[2, 2] = -1
            #     M = flip_z_4x4 @ M

            self.projection_matrices[i] = M

    def project_point(self, point_3d, i_tilt):
        """
        Apply the 4x4 transformation for tilt i_tilt to project a 3D coordinate
        (in tomogram voxels) to 2D tilt coords.
        """
        pt_pix = point_3d / self.pixel_size
        pt_homog = np.append(pt_pix, 1.0)
        M = self.projection_matrices[i_tilt]
        return (M @ pt_homog)[:2]

    def calculate_local_defocus(self, i_tilt, point_3d):
        R_x = R.from_euler("x", self.tomo_x_tilt_arr[i_tilt], degrees=True)
        R_y = R.from_euler("y", self.tomo_y_tilt_arr[i_tilt], degrees=True)
        R_z = R.from_euler("z", self.tomo_z_rot_arr[i_tilt], degrees=True)
        R_zyx = (R_z * R_y * R_x).as_matrix()

        proj_mat = np.eye(4)
        proj_mat[:3, :3] = R_zyx

        proj_mat[0, 3] = self.tomo_x_shift_angst_arr[i_tilt]
        proj_mat[1, 3] = self.tomo_y_shift_angst_arr[i_tilt]

        proj_pos = proj_mat @ np.append(point_3d, 1.0)
        proj_center = proj_mat @ np.array([0.0, 0.0, 0.0, 1.0])

        depth_offset = (proj_pos[2] - proj_center[2]) * self.tomo_hand

        loc_u = self.defocus_u_arr[i_tilt] + depth_offset
        loc_v = self.defocus_v_arr[i_tilt] + depth_offset
        loc_angle = self.defocus_angle_arr[i_tilt]

        return loc_u, loc_v, loc_angle

    def expand_3D_particle_to_2D_particle(self, point_3d, image_name, base_orientation_zyz):
        coords_2d = []
        defocusU_list = []
        defocusV_list = []
        defocusAngle_list = []
        final_zyz_list = []

        for i in range(self.n_tilts):
            coords_2d.append(self.project_point(point_3d, i))
            lu, lv, la = self.calculate_local_defocus(i, point_3d)
            defocusU_list.append(lu)
            defocusV_list.append(lv)
            defocusAngle_list.append(la)

            if base_orientation_zyz is not None:
                tilt_rot_3x3 = self.projection_matrices[i][:3, :3]
                R_final = R.from_matrix(base_orientation_zyz.as_matrix() @ tilt_rot_3x3)
                a, b, c = R_final.as_euler("ZYZ", degrees=True)
                final_zyz_list.append([a, b, c])

        coords_2d = np.array(coords_2d)
        final_zyz_list = np.array(final_zyz_list)

        # make 2D dataframe
        n_tilts = self.n_tilts
        image_names_2d = [f"{i+1:06d}@{image_name}" for i in range(n_tilts)]

        df_2d = pd.DataFrame(
            {
                "rlnImageName": image_names_2d,
                "rlnDefocusU": defocusU_list,
                "rlnDefocusV": defocusV_list,
                "rlnDefocusAngle": defocusAngle_list,
                "rlnCoordinateX": coords_2d[:, 0],
                "rlnCoordinateY": coords_2d[:, 1],
            }
        )

        if base_orientation_zyz is not None and len(final_zyz_list) > 0:
            df_2d["rlnAngleRot"] = final_zyz_list[:, 0]
            df_2d["rlnAngleTilt"] = final_zyz_list[:, 1]
            df_2d["rlnAnglePsi"] = final_zyz_list[:, 2]
        else:
            df_2d["rlnAngleRot"] = 0.0
            df_2d["rlnAngleTilt"] = 0.0
            df_2d["rlnAnglePsi"] = 0.0

        return df_2d


def main(args):
    # load tomograms .star file
    tomograms_star = starfile.read(args.tomograms, always_dict=True)
    tomograms_df = tomograms_star["global"]

    print(f"loaded {len(tomograms_df)} tomograms from [{args.tomograms}]")

    missing_columns = TOMOGRAMS_REQUIRED_COLUMNS - set(tomograms_df.columns)
    if missing_columns:
        raise ValueError(
            f"tomograms .star file [{args.tomograms}] is missing required columns: {missing_columns}"
        )
    
    # load particles .star file
    particles_3D_star = starfile.read(args.particles, always_dict=True)
    try:
        keys = set(particles_3D_star.keys())

        if {"optics", "particles"}.issubset(keys):
            particles_3D_df = particles_3D_star["particles"]
            optics_3D_df = particles_3D_star["optics"]

            particles_3D_df = particles_3D_df.merge(
                optics_3D_df,
                on="rlnOpticsGroup",
                how="left"
            )
            
        elif len(keys) == 1:
            k = next(iter(keys))
            particles_3D_df = particles_3D_star[k]

        else:
            raise ValueError(
                f"unexpected  .star file format, found: {keys}, expected [optics, particles] or single block"
            )

    except Exception as e:
        raise RuntimeError(f"failed to parse .star file: {e}")

    print(f"loaded {len(particles_3D_df)} particles from [{args.particles}]")

    missing_columns = PARTICLES_3D_REQUIRED_COLUMNS - set(particles_3D_df.columns)
    if missing_columns:
        raise ValueError(
            f"3D particles .star file [{args.particles}] is missing required columns: {missing_columns}"
        )

    # begin processing
    tilt_series_dict = {}
    particles_2D_list = []
    particles_2D_index = 0

    for particles_3D_i_index, particles_3D_i in particles_3D_df.iterrows():
        # parse particles_3D_i
        particles_3D_i_tomo_name = particles_3D_i["rlnTomoName"]
        particles_3D_i_angle_rot = particles_3D_i["rlnAngleRot"]
        particles_3D_i_angle_tilt = particles_3D_i["rlnAngleTilt"]
        particles_3D_i_angle_psi = particles_3D_i["rlnAnglePsi"]
        particles_3D_i_centered_coordinate_x_angst = particles_3D_i["rlnCenteredCoordinateXAngst"]
        particles_3D_i_centered_coordinate_y_angst = particles_3D_i["rlnCenteredCoordinateYAngst"]
        particles_3D_i_centered_coordinate_z_angst = particles_3D_i["rlnCenteredCoordinateZAngst"]
        particles_3D_i_tomo_particle_name = particles_3D_i["rlnTomoParticleName"]
        particles_3D_i_tomo_visible_frames = particles_3D_i["rlnTomoVisibleFrames"]
        particles_3D_i_image_name = particles_3D_i["rlnImageName"]
        particles_3D_i_tomo_subtomogram_rot = particles_3D_i.get("rlnTomoSubtomogramRot", 0.0)
        particles_3D_i_tomo_subtomogram_tilt = particles_3D_i.get("rlnTomoSubtomogramTilt", 0.0)
        particles_3D_i_tomo_subtomogram_psi = particles_3D_i.get("rlnTomoSubtomogramPsi", 0.0) 
        particles_3D_i_voltage = particles_3D_i["rlnVoltage"]
        particles_3D_i_spherical_aberration = particles_3D_i["rlnSphericalAberration"]
        particles_3D_i_amplitude_contrast = particles_3D_i["rlnAmplitudeContrast"]
        particles_3D_i_image_size = particles_3D_i["rlnImageSize"]
        particles_3D_i_image_pixel_size = particles_3D_i["rlnImagePixelSize"]

        # find tomograms_i (matched to particles_3D_i by rlnTomoName)
        tomograms_filter = tomograms_df["rlnTomoName"] == particles_3D_i_tomo_name
        assert tomograms_filter.sum() == 1, f"expected 1 row for tomogram [{particles_3D_i_tomo_name}] in tomograms .star file, found {tomograms_filter.sum()}"
        
        tomograms_i = tomograms_df[tomograms_filter].iloc[0]
        tomograms_i_index = tomograms_df.index[tomograms_filter][0]

        tomograms_i_tomo_name = tomograms_i["rlnTomoName"]
        tomograms_i_tomo_hand = -1 if tomograms_i.get("rlnTomoHand", 1) == -1 else 1
        tomograms_i_tilt_series_pixel_size = tomograms_i["rlnTomoTiltSeriesPixelSize"]
        tomograms_i_tilt_series_star_file = tomograms_i["rlnTomoTiltSeriesStarFile"]

        # load tilt series .star file for tomograms_i
        if tomograms_i_tomo_name not in tilt_series_dict:
            tilt_series_i = starfile.read(tomograms_i_tilt_series_star_file, always_dict=True)
            tilt_series_i_df = next(iter(tilt_series_i.values()))

            print(f"loaded {len(tilt_series_i_df)} tilts for tomogram [{tomograms_i_tomo_name}] from tilt series .star file [{tomograms_i_tilt_series_star_file}]")

            missing_columns = TILT_SERIES_REQUIRED_COLUMNS - set(tilt_series_i_df.columns)
            if missing_columns:
                raise ValueError(
                    f"tilt series .star file [{tomograms_i_tilt_series_star_file}] is missing required columns: {missing_columns}"
                )
            
            tilt_series_dict[tomograms_i_tomo_name] = tilt_series_i_df
        
        else:
            tilt_series_i_df = tilt_series_dict[tomograms_i_tomo_name]

        # filter tilt_series_i to the subset of consecutive visible frames/tilts
        visible_frames = literal_eval(particles_3D_i_tomo_visible_frames)
        visible_frames_filter = []
        if visible_frames and visible_frames[0] == 1:
            for i, val in enumerate(visible_frames):
                if val == 1:
                    visible_frames_filter.append(i)
                else:
                    break
            print(f" ● selected {len(visible_frames_filter)} visible frames for particle [{particles_3D_i_tomo_particle_name}]")
        else:
            visible_frames_filter = []
            print(f" ● dropped particle [{particles_3D_i_tomo_particle_name}] with no visible frames")
            continue
        
        tilt_series_i_df_visible = tilt_series_i_df.iloc[visible_frames_filter].copy()
        tilt_series_i_index_arr = tilt_series_i_df_visible.index.to_numpy()

        # build Tomogram object for this particle's tomogram and tilt series
        _tomogram = Tomogram(
            tilt_image_dims=args.tilt_dims,
            pixel_size=tomograms_i_tilt_series_pixel_size,
            defocus_u_arr=tilt_series_i_df_visible["rlnDefocusU"].to_numpy(),
            defocus_v_arr=tilt_series_i_df_visible["rlnDefocusV"].to_numpy(),
            defocus_angle_arr=tilt_series_i_df_visible["rlnDefocusAngle"].to_numpy(),
            tomo_x_tilt_arr=tilt_series_i_df_visible["rlnTomoXTilt"].to_numpy(),
            tomo_y_tilt_arr=tilt_series_i_df_visible["rlnTomoYTilt"].to_numpy(),
            tomo_z_rot_arr=tilt_series_i_df_visible["rlnTomoZRot"].to_numpy(),
            tomo_x_shift_angst_arr=tilt_series_i_df_visible["rlnTomoXShiftAngst"].to_numpy(),
            tomo_y_shift_angst_arr=tilt_series_i_df_visible["rlnTomoYShiftAngst"].to_numpy(),
            tomo_hand=tomograms_i_tomo_hand,
        )

        # particle orientation
        R_particle = R.from_euler(
            "ZYZ", [
                particles_3D_i_angle_rot,
                particles_3D_i_angle_tilt,
                particles_3D_i_angle_psi
            ], degrees=True
        )

        # box orientation
        R_subtomo = R.from_euler(
            "ZYZ", [
                particles_3D_i_tomo_subtomogram_rot,
                particles_3D_i_tomo_subtomogram_tilt,
                particles_3D_i_tomo_subtomogram_psi
            ], degrees=True
        )

        # combined orientation from local to tomogram
        R_base = R_particle * R_subtomo
        point_3d_rotated = np.array([
            particles_3D_i_centered_coordinate_x_angst,
            particles_3D_i_centered_coordinate_y_angst,
            particles_3D_i_centered_coordinate_z_angst
        ])

        # expand particles rows into tilt rows
        _particle = _tomogram.expand_3D_particle_to_2D_particle(
            point_3d=point_3d_rotated,
            image_name=particles_3D_i_image_name,
            base_orientation_zyz=R_base,
        )

        # add additional metadata columns
        _particle["rlnGroupName"] = particles_3D_i_tomo_particle_name
        _particle["rlnImageSize"] = particles_3D_i_image_size
        _particle["rlnImagePixelSize"] = particles_3D_i_image_pixel_size
        _particle["rlnMicrographPreExposure"] = tilt_series_i_df_visible["rlnMicrographPreExposure"].to_numpy()
        _particle["rlnPhaseShift"] = tilt_series_i_df_visible["rlnPhaseShift"].to_numpy()
        _particle["rlnVoltage"] = particles_3D_i_voltage
        _particle["rlnSphericalAberration"] = particles_3D_i_spherical_aberration
        _particle["rlnAmplitudeContrast"] = particles_3D_i_amplitude_contrast

        _particle["rlnCtfScalefactor"] = (
            tilt_series_i_df_visible["rlnCtfScalefactor"].to_numpy()
            if "rlnCtfScalefactor" in tilt_series_i_df_visible.columns
            else np.cos(np.deg2rad(tilt_series_i_df_visible["rlnTomoYTilt"].to_numpy()))
        )

        if args.full_data:
            _particle["_particles_3D_index"] = particles_3D_i_index # index of the 3D particle in the input 3D particles .star file
            for col, val in particles_3D_i.items():
                _particle[f"_particles_3D_{col}"] = val

            _particle["_tomograms_index"] = tomograms_i_index # index of the tomogram in the input tomograms .star file
            for col, val in tomograms_i.items():
                _particle[f"_tomograms_{col}"] = val
            
            _particle["_tilt_series_index"] = tilt_series_i_index_arr # index of the tilt within the tilt series .star file for this tomogram
            for col in tilt_series_i_df_visible.columns:
                _particle[f"_tilt_series_{col}"] = tilt_series_i_df_visible[col].values

        _particle = _particle.sort_values("rlnMicrographPreExposure", ascending=True).reset_index(drop=True)
        
        _particle["_particle_index"] = particles_2D_index # particle index
        particles_2D_index += 1
        
        _particle["_tilt_index"] = tilt_series_i_index_arr # tilt index within each particle

        particles_2D_list.append(_particle)

    # combine all particles_2D_i into single particles_2D_df dataframe
    particles_2D_df = pd.concat(particles_2D_list, ignore_index=True)
    particles_2D_df["_image_index"] = particles_2D_df.index # unique index for each 2D particle row

    # reorder columns
    if args.full_data:
        particles_2D_df = particles_2D_df[PARTICLES_2D_OUTPUT_COLUMNS + [c for c in particles_2D_df.columns if c not in PARTICLES_2D_OUTPUT_COLUMNS]]
    else:
        particles_2D_df = particles_2D_df[PARTICLES_2D_OUTPUT_COLUMNS]

    print(f"generated metadata for {len(particles_2D_df)} images from {particles_2D_index}/{len(particles_3D_df)} 3D particles in {len(tomograms_df)} tomograms")

    # save output
    output_df = {"": particles_2D_df}
    starfile.write(output_df, args.output, overwrite=True)
    print(f"saved [{args.output}]")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse 3D particle, tomogram, and tilt series metadata from RELION v5 .star files to calculate 2D particle (subtilt image) geometry."
    )

    add_args(parser)
    main(parser.parse_args())
