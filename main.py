"""Script to run the spline-based track pipeline."""

from spline_track_pipeline import SplineTrackPipeline

def main():
    """Main function to run the pipeline."""
    pipeline = SplineTrackPipeline(
        spline_data_file="Wheelset-Trk_Track.txt",
        output_dir="Wheelset.output"
    )
    
    pipeline.run(run_abaqus=True)


if __name__ == "__main__":
    main()