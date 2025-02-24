#! /usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
This script is used to quickly visualise the data of a .bm (benchmarking data text file)
output from a Commec Screen run. It generates an HTML file of provided name. It is run
automatically at the end of a Commec Screen, but can also be run as main on a valid .bm file.
"""
import os
import sys
import pandas as pd
import plotly.express as px

def time_str_to_seconds(time_str : str) -> float:
    """ 
    Convert string format 
    from "HH:MM:SS.SSS" 
    to seconds.milliseconds (float).
    """
    h, m, s = time_str.split(':')
    return int(h) * 3600 + int(m) * 60 + float(s)

def capitalize_and_concatenate(input_str):
    """ Helper function to make function names a little nicer to read."""
    tokens = input_str.replace('_', ' ').split()
    capitalized_tokens = [token.capitalize() for token in tokens]
    return ' '.join(capitalized_tokens)

def visualize_data(filename : os.PathLike, outfile : os.PathLike = None):
    """ 
    Takes the output benchmarking log from commec, 
    and produces a pretty html visual using plotly. 
    """

    if not os.path.exists(filename):
        sys.stderr.write("\t... input benchmarking file does not exist!\n" + filename)
        return

    # Parse the benchmark file:
    entries = []
    with open(filename, 'r', encoding="utf-8") as file:
        for line in file:
            parts = line.strip().split()

            # Skip any bad lines.
            if len(parts) < 5:
                continue

            y_pos = int(parts[0])                             # Stack depth
            label = capitalize_and_concatenate(parts[1]) # Function name
            duration_str = parts[2]                      # Duration string (for hover info)
            start_time_str = parts[3]                    # Start time as string
            end_time_str = parts[4]                      # End time as string

            start_time = time_str_to_seconds(start_time_str)
            end_time = time_str_to_seconds(end_time_str)
            duration = time_str_to_seconds(duration_str)

            entries.append({
                "y_pos": y_pos,
                "label": label,
                "duration": duration,
                "start_time": start_time,
                "end_time": end_time
            })

    df = pd.DataFrame(entries)

    # Plotly timelines REQUIRE a date string as the x-axis, so we start counting from 0... i.e. 1970.
    base_date = pd.Timestamp("1970-01-01")
    df['start_datetime'] = base_date + pd.to_timedelta(df['start_time'], unit='s')
    df['end_datetime'] = base_date + pd.to_timedelta(df['end_time'], unit='s')
    
    fig = px.timeline(
        df,
        x_start="start_datetime",
        x_end="end_datetime",
        y="y_pos",
        color="label",              # Unique colours for unique functions.
        hover_name="label",
        hover_data={"Duration": df["duration"]},
        title="Benchmark Visualization"
    )

    fig.update_layout(
        title="Benchmark Visualization",
        xaxis=dict(
            title="Time (seconds)",
            showgrid=True,
            zeroline=False,
            tickformat="%s"
        ),
        yaxis=dict(
            title="Call Stack",
            tickmode="linear",
            showgrid=True,
            zeroline=False
        )
    )

    # Customize hover to show only label and duration
    fig.update_traces(
        hovertemplate="<b>%{hovertext}</b><br>Duration: %{customdata[0]}<extra></extra>"
    )

    root, _ext = os.path.splitext(filename)
    output_filename = outfile or f"{root}_benchmarking.html"
    fig.write_html(output_filename, full_html = False, include_plotlyjs='cdn')

# Usage:
if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_filename = sys.argv[1]  # Take the first argument as the filename
    else:
        print("Error: No filename provided.")
        print("Usage: python script.py <filename>")
        sys.exit(1)  # Exit if no filename is provided

    visualize_data(input_filename)
