#!/bin/python3

import numpy as np
import subprocess
from typing import List
import sys
import os


class MagneticField:
    def __init__(self, field_strength: float, lateral_a: float, field_max_loc: float):
        # in T
        self.field_strength = field_strength
        # in mm
        self.lateral_a = lateral_a
        # in mm
        self.field_max_loc = field_max_loc


class LensInput:
    def __init__(
        self,
        output_ps_file: str,
        output_png_file: str,
        sampling_points: int,
        magnetic_fields: List[MagneticField],
        optical_axis_length: float,
        acc_voltage: int,
    ):
        assert len(magnetic_fields) in [1, 2]
        self.output_ps_file = output_ps_file
        self.output_png_file = output_png_file
        self.sampling_points = sampling_points
        self.magnetic_fields = magnetic_fields
        # in mm
        self.optical_axis_length = optical_axis_length
        # in V
        self.acc_voltage = acc_voltage


def create_lens_input_str(lens_input: LensInput) -> str:
    output = f"{lens_input.output_ps_file}\n"
    output += f"{lens_input.sampling_points}\n"
    output += f"{len(lens_input.magnetic_fields)}\n"
    for magnetic_field in lens_input.magnetic_fields:
        output += f"{magnetic_field.field_strength:.1f}\n"
        output += f"{magnetic_field.lateral_a:.1f}\n"
        output += f"{magnetic_field.field_max_loc:.1f}\n"
    output += f"{lens_input.optical_axis_length:.1f}\n"
    output += f"{lens_input.acc_voltage}\n"
    return output


def execute_lenses(inputs: List[LensInput]):
    input_strs = list(map(create_lens_input_str, inputs))
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
                "-r600",
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


def main():
    print("CTEMsoft lens")
    OUT_DIR = "out"
    os.makedirs(OUT_DIR, exist_ok=True)
    execute_lenses([
        LensInput(
            output_ps_file=f"{OUT_DIR}/001.ps",
            output_png_file=f"{OUT_DIR}/001.png",
            sampling_points=1000,
            magnetic_fields=[
                MagneticField(field_strength=2, lateral_a=1, field_max_loc=0)
            ],
            optical_axis_length=10,
            acc_voltage=100000,
        )
    ])
    print("All done")
    print("Have a very safe and productive day.")


if __name__ == "__main__":
    main()
