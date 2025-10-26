"""Pipeline for generating spline-based track models from spline data."""

import json
import subprocess
import tempfile
import os
import sys
from typing import Dict, List

from gmsh_rail_generator import generate_rail_sweep
from sleeper_generator import generate_sleeper_sweep


class SplineTrackPipeline:
    """Pipeline class to manage the entire spline-based track generation process."""
    
    def __init__(
        self,
        spline_data_file: str,
        output_dir: str = "Wheelset.output",
        sleeper_length: float = 0.29,
        rail_burn_in: float = 0.205,
        sleeper_spacing: float = 0.7,
    ):
        """Initialize the pipeline.
        
        Args:
            spline_data_file (str): Path to the spline data .txt file
            output_dir (str): Directory for output files
            sleeper_length (float): Length of the sleeper
            rail_burn_in (float): Burn-in length of the rail
            sleeper_spacing (float): Spacing between sleepers
        """
        self.spline_data_file = spline_data_file
        self.output_dir = output_dir
        self.step_files = {}
        self.json_config = {}

        self.sleeper_length = sleeper_length
        self.rail_burn_in = rail_burn_in
        self.sleeper_spacing = sleeper_spacing
        
        self.left_half_points = [
            (1.244588, 283.628E-03),
            (1.224584, 283.628E-03),
            (1.204789, 285.628E-03),
            (1.204789, 291.404E-03),
            (1.208504, 295.394E-03),
            (1.224584, 296.545E-03),
            (1.244588, 302.192E-03),
            (1.263218, 308.967E-03),
            (1.267464, 313.325E-03),
            (1.269523, 322.022E-03),
            (1.271539, 343.928E-03),
            (1.271539, 359.538E-03),
            (1.271131, 385.814E-03),
            (1.268339, 403.458E-03),
            (1.265682, 410.873E-03),
            (1.26189, 414.324E-03),
            (1.25177, 418.175E-03),
            (1.244818, 420.82E-03),
            (1.243191, 429.775E-03),
            (1.243773, 441.367E-03),
            (1.246324, 448.472E-03),
            (1.253665, 453.342E-03),
            (1.269539, 455.453E-03),
        ]
        
        self.center_top = [(1.279789, 455.628E-03)]
        self.center_bottom = [(1.279789, 283.628E-03)]
        
        self.sleeper_profile_points = [
            (-0.756813, 0.172471),
            (-0.756813, 0.061628),
            (-0.312048, 0.061628),
            (-0.312048, 0.031628),
            (-0.222048, 0.031628),
            (-0.132048, 0.031628),
            (-0.132048, 0.061628),
            (1.189789, 0.061628),
            (1.189789, 0.031628),
            (1.269789, 0.031628),
            (1.369789, 0.031628),
            (1.369789, 0.061628),
            (1.843187, 0.061628),
            (1.843187, 0.161628),
            (1.543187, 0.261628),
            (1.369789, 0.261628),
            (1.369789, 0.281628),
            (1.269789, 0.281628),
            (1.189789, 0.281628),
            (1.189789, 0.261628),
            (1.010661, 0.261628),
            (0.710661, 0.161628),
            (0.310661, 0.161628),
            (0.010661, 0.261628),
            (-0.132048, 0.261628),
            (-0.132048, 0.281628),
            (-0.222048, 0.281628),
            (-0.312048, 0.281628),
            (-0.312048, 0.261628),
            (-0.489339, 0.261628),
        ]
    
    def reflect_rail_point(self, point: tuple) -> tuple:
        """Reflect a point across the vertical centerline.
        
        Args:
            point (tuple): Point (x, y)
            
        Returns:
            tuple: Reflected point
        """
        x, y = point
        center_x = self.center_top[0][0]
        new_x = 2 * center_x - x
        return (new_x, y)
    
    def generate_rail_profile_points(self) -> List[tuple]:
        """Generate the complete symmetric rail profile.
        
        Returns:
            List[tuple]: Complete rail profile points
        """
        right_half_points = [self.reflect_rail_point(p) for p in reversed(self.left_half_points)]
        rail_profile_points = self.left_half_points + self.center_top + right_half_points + self.center_bottom
        
        y_coords = [y for y, _ in rail_profile_points]
        z_coords = [z for _, z in rail_profile_points]
        y_center = (min(y_coords) + max(y_coords)) / 2
        z_center = (min(z_coords) + max(z_coords)) / 2
        
        return [(y - y_center, z - z_center) for y, z in rail_profile_points]
    
    
    def generate_step_files(self) -> bool:
        txt_file = os.path.abspath(os.path.join(self.output_dir, "Wheelset-Trk_Track.txt"))
        print(txt_file)
        rail_step = generate_rail_sweep(txt_file, self.generate_rail_profile_points())
        sleeper_step = generate_sleeper_sweep(self.sleeper_length, self.sleeper_profile_points)
        
        self.step_files = {
            'rail_step': rail_step,
            'sleeper_step': sleeper_step
        }
        
        return self.step_files
    
    def save_json_config(self) -> str:
        """Save configuration to JSON file.
        
        Returns:
            str: Path to the JSON file
        """
        config = {
            'rail_step_file': self.step_files['rail_step'],
            'sleeper_step_file': self.step_files['sleeper_step'],
            'base_cae_file': "SleeperMeshRefinement.cae",
            'SleeperLength': self.sleeper_length,
            'RailBurnIn': self.rail_burn_in,
            'SleeperSpacing': self.sleeper_spacing,
        }
        
        json_path = os.path.join(self.output_dir, "spline_track_config.json")
        with open(json_path, 'w') as f:
            json.dump(config, f, indent=4)
        
        self.json_config = config
    
    def modify_abaqus_script(self, script_path: str):
        """Modify the JSON file reference in the Abaqus script.
        
        Args:
            script_path (str): Path to the Abaqus script
        """
        with open(script_path, 'r') as f:
            lines = f.readlines()
        
        json_config_path = os.path.join(self.output_dir, "spline_track_config.json")
        for i, line in enumerate(lines):
            if 'json_file="spline_track_config.json"' in line:
                lines[i] = f'        json_file="{json_config_path}"\n'
                break
        
        with open(script_path, 'w') as f:
            f.writelines(lines)
        
    def run_abaqus_script(self, script_path: str):
        """Run the Abaqus script.
        
        Args:
            script_path (str): Path to the Abaqus script
            
        Returns:
            bool: True if successful
        """
        self.modify_abaqus_script(script_path)
        cmd = f'abaqus cae noGUI={script_path}'
        subprocess.run(cmd, shell=True, check=True, capture_output=False, text=True)
            
    def run(self, run_abaqus: bool = False, abaqus_script: str = "spline_based_track.py"):
        """Run the complete pipeline.
        
        Args:
            run_abaqus (bool): Whether to run the Abaqus script after generation
            abaqus_script (str): Path to the Abaqus script
            
        Returns:
            tuple: (step_files_dict, config_dict)
        """
        self.generate_step_files()
        self.save_json_config()
        
        if run_abaqus:
            self.run_abaqus_script(abaqus_script)
            
        return self.step_files, self.json_config


def main():
    """Main function to run the pipeline."""
    pipeline = SplineTrackPipeline(
        spline_data_file="Wheelset.output/Wheelset-Trk_Track.txt",
        output_dir="."
    )
    
    step_files, config = pipeline.run()
    print(f"\nGenerated STEP files:")
    print(f"  Rail: {step_files['rail_step']}")
    print(f"  Sleeper: {step_files['sleeper_step']}")
    print(f"\nNow you can run the Abaqus script with: abaqus cae noGUI=spline_based_track.py")


if __name__ == "__main__":
    main()

