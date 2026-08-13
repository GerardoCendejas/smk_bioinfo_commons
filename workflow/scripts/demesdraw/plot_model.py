#!/usr/bin/env module
import argparse
import sys
import matplotlib.pyplot as plt
import demes
import demesdraw
import demesdraw.utils

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot demes models for tubes and size history representations."
    )
    parser.add_argument(
        "-i", "--input", 
        required=True, 
        help="Path to the input model.yaml file"
    )
    parser.add_argument(
        "-t", "--tubes", 
        required=True, 
        help="Output path for the tubes plot (e.g., tubes.png)"
    )
    parser.add_argument(
        "-s", "--size-history", 
        required=True, 
        help="Output path for the size history plot (e.g., size_history.png)"
    )
    return parser.parse_args()

def main():
    args = parse_args()

    try:
        # Load the demes graph
        graph = demes.load(args.input)
    except Exception as e:
        print(f"Error loading YAML file: {e}", file=sys.stderr)
        sys.exit(1)

    # Use demesdraw heuristics to determine optimal scales
    log_time = demesdraw.utils.log_time_heuristic(graph)
    log_size = demesdraw.utils.log_size_heuristic(graph)

    print(f"Generating plots for: {args.input}")

    # 1. Generate and save the size history plot
    # demesdraw creates its own figure and returns the Axes object
    ax_size = demesdraw.size_history(
        graph,
        invert_x=True,
        log_time=log_time,
        log_size=log_size,
        title="Population Size History"
    )
    # Target the figure object and save it to the specified path
    ax_size.figure.savefig(args.size_history, dpi=300, bbox_inches="tight")
    plt.close(ax_size.figure)
    print(f"Saved size history plot to: {args.size_history}")

    # 2. Generate and save the tubes plot
    ax_tubes = demesdraw.tubes(
        graph,
        log_time=log_time,
        title="Demographic Tubes Graph"
    )
    ax_tubes.figure.savefig(args.tubes, dpi=300, bbox_inches="tight")
    plt.close(ax_tubes.figure)
    print(f"Saved tubes plot to: {args.tubes}")

if __name__ == "__main__":
    main()
