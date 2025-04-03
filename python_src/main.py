#!/bin/python3

import numpy as np
import subprocess
from typing import List
import sys
import os
import glob


class MagneticField:
    def __init__(self, field_strength: float, lateral_a: float, field_max_loc: float):
        # in T
        self.field_strength = field_strength
        # in mm
        self.lateral_a = lateral_a
        # in mm
        self.field_max_loc = field_max_loc


class LensParameters:
    def __init__(
        self,
        sampling_points: int,
        magnetic_fields: List[MagneticField],
        optical_axis_length: float,
        acc_voltage: int,
    ):
        assert len(magnetic_fields) in [1, 2]
        self.sampling_points = sampling_points
        self.magnetic_fields = magnetic_fields
        # in mm
        self.optical_axis_length = optical_axis_length
        # in V
        self.acc_voltage = acc_voltage

    def to_ndarray(self) -> np.ndarray:
        if len(self.magnetic_fields) == 1:
            return np.array([
                self.sampling_points,
                self.magnetic_fields[0].field_strength,
                self.magnetic_fields[0].lateral_a,
                self.magnetic_fields[0].field_max_loc,
                self.optical_axis_length,
                self.acc_voltage,
            ])
        elif len(self.magnetic_fields) == 2:
            return np.array([
                self.sampling_points,
                self.magnetic_fields[0].field_strength,
                self.magnetic_fields[0].lateral_a,
                self.magnetic_fields[0].field_max_loc,
                self.magnetic_fields[1].field_strength,
                self.magnetic_fields[1].lateral_a,
                self.magnetic_fields[1].field_max_loc,
                self.optical_axis_length,
                self.acc_voltage,
            ])
        else:
            assert False

    @classmethod
    def from_ndarray(cls, vec: np.ndarray):
        if len(vec) == 6:
            return cls(
                sampling_points=vec[0],
                magnetic_fields=[
                    MagneticField(
                        field_strength=vec[1], lateral_a=vec[2], field_max_loc=vec[3]
                    )
                ],
                optical_axis_length=vec[4],
                acc_voltage=vec[5],
            )
        elif len(vec) == 9:
            return cls(
                sampling_points=vec[0],
                magnetic_fields=[
                    MagneticField(
                        field_strength=vec[1], lateral_a=vec[2], field_max_loc=vec[3]
                    ),
                    MagneticField(
                        field_strength=vec[4], lateral_a=vec[5], field_max_loc=vec[6]
                    ),
                ],
                optical_axis_length=vec[7],
                acc_voltage=vec[8],
            )
        else:
            assert False


class LensInput:
    def __init__(
        self, output_ps_file: str, output_png_file: str, lens_parameters: LensParameters
    ):
        self.output_ps_file = output_ps_file
        self.output_png_file = output_png_file
        self.lens_parameters = lens_parameters

    def to_input_str(self) -> str:
        output = f"{self.output_ps_file}\n"
        output += f"{int(self.lens_parameters.sampling_points)}\n"
        output += f"{len(self.lens_parameters.magnetic_fields)}\n"
        for magnetic_field in self.lens_parameters.magnetic_fields:
            output += f"{magnetic_field.field_strength:.10f}\n"
            output += f"{magnetic_field.lateral_a:.10f}\n"
            output += f"{magnetic_field.field_max_loc:.10f}\n"
        output += f"{self.lens_parameters.optical_axis_length:.10f}\n"
        output += f"{int(self.lens_parameters.acc_voltage)}\n"
        return output


def execute_lenses(inputs: List[LensInput]):
    input_strs = list(map(lambda i: i.to_input_str(), inputs))
    lens_procs: List[subprocess.Popen[bytes]] = []
    for input, input_str in zip(inputs, input_strs):
        print(f"running lens for {input.output_ps_file}")
        # lens should be in PATH
        ps_proc = subprocess.Popen(
            "lens",
            stdin=subprocess.PIPE,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        ps_proc.communicate(input=str.encode(input_str))
        lens_procs.append(ps_proc)

    gs_procs: List[subprocess.Popen[bytes]] = []
    for ps_proc, input, input_str in zip(lens_procs, inputs, input_strs):
        if ps_proc.wait() != 0:
            print("lens failed")
            print("input:")
            print(input_str)
            sys.exit(1)
        if not os.path.isfile(input.output_ps_file):
            print("lens failed to produce an output")
            print("input:")
            print(input_str)
            sys.exit(1)
        print(f"running ghostscript for {input.output_png_file}")
        gs_proc = subprocess.Popen(
            [
                "gs",
                "-dSAFER",
                "-dEPSCrop",
                "-r100",
                "-sDEVICE=pngalpha",
                "-o",
                input.output_png_file,
                input.output_ps_file,
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        gs_procs.append(gs_proc)

    for gs_proc, input in zip(gs_procs, inputs):
        if gs_proc.wait() != 0:
            print(f"ghostscript failed for {input.output_png_file}")
            sys.exit(1)
        if not os.path.isfile(input.output_ps_file):
            print(f"ghostscript failed to produce output for {input.output_png_file}")
            sys.exit(1)


def delete_old_frames(out_dir: str):
    files = glob.glob(f"{out_dir}/*.png") + glob.glob(f"{out_dir}/*.mp4")
    for file in files:
        os.remove(file)


def create_lens_frames(
    out_dir: str, first_frame: LensParameters, last_frame: LensParameters, frames: int
) -> List[LensInput]:
    t_vals = np.linspace(0, 1, frames)
    lens_parameters_vec_list = (
        1 - t_vals[:, None]
    ) * first_frame.to_ndarray() + t_vals[:, None] * last_frame.to_ndarray()
    lens_parameters_list = [
        LensParameters.from_ndarray(lens_parameters_vec)
        for lens_parameters_vec in lens_parameters_vec_list
    ]
    lens_inputs = [
        LensInput(
            output_ps_file=f"{out_dir}/{n:03d}.ps",
            output_png_file=f"{out_dir}/{n:03d}.png",
            lens_parameters=lens_parameters,
        )
        for n, lens_parameters in enumerate(lens_parameters_list)
    ]
    # for lens_input in lens_inputs:
    #     print(lens_input.to_input_str())
    execute_lenses(lens_inputs)
    return lens_inputs


# def create_file_list(lens_inputs: List[LensInput], file_list_path: str):
#     with open(file_list_path, "w") as file:
#         file.writelines([
#             f"file '{lens_input.output_png_file}'\n" for lens_input in lens_inputs
#         ])


def convert_to_video(lens_inputs: List[LensInput], out_dir: str):
    # this doesn't work for some reason :<
    # this needs to be in the current working directory because ffmpeg uses paths relative to the location of this file
    # file_list_path = "filelist.txt"
    # create_file_list(lens_inputs, file_list_path)

    output_video = f"{out_dir}/lens.mp4"
    proc = subprocess.Popen(
        [
            "ffmpeg",
            "-framerate",
            "30",
            "-i",
            f"{out_dir}/%03d.png",
            "-vcodec",
            "libx264",
            "-f",
            "mp4",
            "-vb",
            "1024k",
            "-preset",
            "slow",
            output_video,
        ],
        # stdout=subprocess.DEVNULL,
        # stderr=subprocess.DEVNULL,
    )
    if proc.wait() != 0:
        print("ffmpeg failed")
        sys.exit(1)


def main():
    print("CTEMsoft lens")
    OUT_DIR = "out"
    os.makedirs(OUT_DIR, exist_ok=True)
    delete_old_frames(OUT_DIR)
    lens_inputs = create_lens_frames(
        OUT_DIR,
        first_frame=LensParameters(
            sampling_points=1000,
            magnetic_fields=[
                MagneticField(field_strength=1, lateral_a=1, field_max_loc=0)
            ],
            optical_axis_length=10,
            acc_voltage=100000,
        ),
        last_frame=LensParameters(
            sampling_points=1000,
            magnetic_fields=[
                MagneticField(field_strength=4, lateral_a=1, field_max_loc=0)
            ],
            optical_axis_length=10,
            acc_voltage=100000,
        ),
        frames=300,
    )
    convert_to_video(lens_inputs, OUT_DIR)
    print("All done")
    print("Have a very safe and productive day.")


if __name__ == "__main__":
    main()
