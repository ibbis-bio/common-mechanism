"""
Used to convert a json object created from Commec Screen, or a ScreenData object,
into a visual HTML representation of the Commec output. Which can be embedded into
any other HTML document as appropriate.
"""

from enum import StrEnum
import argparse
import plotly.graph_objects as go
import pandas as pd
import plotly.express as px
import datetime
from mako.template import Template

from commec.config.json_io import (
    get_screen_data_from_json,
    ScreenData,
    QueryData,
    HitDescription,
    CommecScreenStep,
)

class CommecPalette():
    """ 
    Enum of the colour palette used with Commec 
    """
    WHITE = [255,255,255]
    DK_BLUE = [35,42,88]
    LT_BLUE = [66,155,185]
    ORANGE = [241,80,36]
    # Yellow and Red are not official Commec colours,
    # however I like these to compliment the above.
    YELLOW = [241,80,36]
    RED = [207,27,81]

    @staticmethod
    def mod(modulate, alpha: int = 255, multiplier : float = 1.0) -> str:
        m : float = multiplier
        return str(f'rgba({round(modulate[0]*m)},{round(modulate[1]*m)},{round(modulate[2]*m)},{alpha})')

    rgba_WHITE = 'rgba(255,255,255,255)'
    rgba_DK_BLUE = 'rgba(35,42,88,255)'
    rgba_LT_BLUE = 'rgba(66,155,185,255)'
    rgba_ORANGE = 'rgba(241,80,36,255)'
    # Yellow and Red are not official Commec colours,
    # however I like these to compliment the above.
    rgba_YELLOW = 'rgba(241,80,36,255)'
    rgba_RED = 'rgba(207,27,81,255)'

def color_from_hit(hit : HitDescription):
    """ Convert a Screen step into an associated Colour."""
    if hit.recommendation.from_step == CommecScreenStep.BIORISK:
        return CommecPalette.RED
    if hit.recommendation.from_step == CommecScreenStep.BENIGN:
        return CommecPalette.LT_BLUE
    if hit.recommendation.from_step == CommecScreenStep.TAXONOMY_AA:
        return CommecPalette.ORANGE
    if hit.recommendation.from_step == CommecScreenStep.TAXONOMY_NT:
        return CommecPalette.YELLOW
    return CommecPalette.DK_BLUE

def generate_html_from_screen_data(input_data : ScreenData, output_file : str):
    """
    Interpret the ScreenData from Commec Screen output as a visualisation, in 
    the form on an HTML output.
    """

    # Create the plot
    fig = go.Figure()
    query_to_draw = input_data.queries[0]

    #fig = make_subplots(
    #    rows=len(input_data.queries),
    #    cols=1,
    #    shared_xaxes=True,  # Share the x-axis across subplots
    #    vertical_spacing=0.05  # Adjust vertical space between subplots
    #)

    stacks = draw_query_to_plot(fig, query_to_draw)

    # Add an invisible trace to the secondary x-axis, else it wont show.
    #fig.add_trace(go.Scatter(
    #    x=[0, query_to_draw.length],  # The x-values; these can be any value
    #    y=[0, 0],  # The y-values; these can be zeros or any dummy values
    #    name='Invisible Trace',  # Optional name
    #    xaxis='x2',  # Assign to the secondary x-axis
    #    mode='lines',  # Use lines to create the trace
    #    line=dict(color='rgba(255, 255, 255, 0)'),  # Fully transparent line
    #    hoverinfo='none',  # No hover information
    #    showlegend=False,
    #))

    # Update layout to display X-axis on top and hide Y-axis labels
    fig.update_layout(
        height = 100 * (stacks + 1),
        #barmode='stack',
        title="Query: " + query_to_draw.query + "  (" + str(query_to_draw.length)+" b.p.)",
        xaxis_title="",
        yaxis_title="",
        yaxis=dict(
            showticklabels=False,
            autorange='reversed',
            fixedrange=True,
            tickmode="linear",
            zeroline=False,
            range=[-0.5,stacks+1+0.5]
        ),
        bargap=0.0,

        template='plotly_white',
        
        # First x-axis (default)
        xaxis=dict(
            showline=True,
            #range=[convert_basepairs_to_datetime(0), convert_basepairs_to_datetime(query_to_draw.length)],
            constrain="domain",
            ticks='outside',
            showgrid=True,
            #fixedrange=True,
            zeroline=False,
            tickformat="%s",  # Display seconds as a numeric value
            type="date"
            #tickmode="linear",  # Ensure ticks are linearly spaced
        ),


    #fig.update_layout(
    #    title="Benchmark Visualization",
    #    xaxis=dict(
    #        title="Time (seconds)",
    #        showgrid=True,
    #        zeroline=False,
    #        tickformat="%s"
    #    ),
    #    yaxis=dict(
    #        title="Call Stack",
    #        tickmode="linear",
    #        showgrid=True,
    #        zeroline=False
    #    )
    #)

        # Second x-axis (xaxis2) at the top
        # the issue with having two axis is they don't necessarily scroll together.
        # We could try draw the same thing to both axis, and hope for the best.
        #xaxis2=dict(
        #    overlaying='x',  # Overlay it on the same domain as xaxis
        #    side='top',  # Place it on top
        #    showline=True,  # Show the X-axis line
        #    constrain="domain",  # Constrain the x-axis to its domain
        #    #range=[convert_basepairs_to_datetime(0), convert_basepairs_to_datetime(query_to_draw.length)],  # Same range as xaxis (or modify if needed)
        #    ticks='outside',  # Show ticks outside
        #    showgrid=False,  # Hide gridlines
        #    #fixedrange=True,
        #    #tickformat="%s",  # Display seconds as a numeric value
        #    #tickmode="linear",  # Ensure ticks are linearly spaced
        #)
    )

    # Show the plot
    fig.write_html(output_file.strip()+".html")

def convert_basepairs_to_datetime(bp : int) -> datetime:
    """ 
    We can't draw timelines with any dtype, so we convert bps to seconds, 
    and use a datetime object which parses as a string. Hurrah for hacky work-arounds.
    """
    return datetime.datetime(1970, 1, 1) + datetime.timedelta(seconds=bp)

def draw_query_to_plot(fig : go.Figure, query_to_draw : QueryData):
    """ 
    Write the data from a single query into the figure for plotly. 
    """
    print(query_to_draw.query, " with ", str(len(query_to_draw.hits)), " hits.")

    graph_data = [
        {"label": query_to_draw.query, "outcome" : "", "start": 0, "stop": query_to_draw.length, "color" : CommecPalette.DK_BLUE, "stack" : 0},
    ]

    write_stack_start_points = [query_to_draw.length]

    for hit in query_to_draw.hits:
        for match in hit.ranges:

            # Find the best spot.
            stack_write = None # Query is 0, so 0 is also invalid default.
            for i, start_point in enumerate(write_stack_start_points):
                # Check if we fit after any existing stacks.
                if start_point < match.query_start:
                    print(i)
                    write_stack_start_points[i-1] = match.query_end
                    stack_write = i
                    break
            
            # We don't fit, create a new stack...
            if not stack_write:
                stack_write = len(write_stack_start_points)
                write_stack_start_points.append(match.query_end)
            
            graph_data.append(
                {
                    "label" : hit.recommendation.from_step + " " + hit.description[:25] + "...",
                    "outcome" : hit.recommendation.outcome + ": " + hit.description[:25] + "...",
                    "start" : match.query_start,
                    "stop" : match.query_end,
                    "color" : color_from_hit(hit),
                    "stack" : stack_write
                }
            )

    for r in graph_data:
        print(r)

    # Add the actual data ranges with individual colors and text inside the bars
    write_stack : int = len(write_stack_start_points)

    df = pd.DataFrame(graph_data)

    # Convert start and stop times from seconds to datetime strings
    epoch = datetime.datetime(1970, 1, 1)
    df['start'] = df['start'].apply(lambda x: epoch + datetime.timedelta(seconds=x))
    df['stop'] = df['stop'].apply(lambda x: epoch + datetime.timedelta(seconds=x))

    print(df['start'])
    print(df['stop'])


    # Convert RGB colors to hex format
    def rgb_to_hex(rgb):
        return '#{:02x}{:02x}{:02x}'.format(rgb[0], rgb[1], rgb[2])

    df['color'] = df['color'].apply(rgb_to_hex)

    # Generate hover text
    df['hovertext'] = df.apply(lambda row: f"{row['label']}<br>({row['start']}-{row['stop']})", axis=1)

    print(df)

    # From benchmarking.
    #fig = px.timeline(
    #    df,
    #    x_start="start_datetime",
    #    x_end="end_datetime",
    #    y="y_pos",
    #    color="label",              # Unique colours for unique functions.
    #    hover_name="label",
    #    hover_data={"Duration": df["duration"]},
    #    title="Benchmark Visualization"
    #)


    # Plot the timeline with Plotly
    temp_fig = px.timeline(
        df, 
        x_start="start",
        x_end="stop",
        y="stack",
        color="color",
        hover_name="hovertext",
        labels={"stack": "Stack"}
    )

    # Transfer the traces from temp_fig to the existing fig
    for trace in temp_fig.data:
        fig.add_trace(trace)

    # Update the layout for better visualization
    #fig.update_layout(
    #    showlegend=False,
    #    title="Timeline Plot",
    #    xaxis_title="Time",
    #    yaxis_title="Stack",
    #)

    return write_stack

def generate_rounded_rect(x0,y0,x1,y1,rx,ry):
    """ 
    Creates an SVG path for a rounded rectangle given the following dimensions:
    xy0: Top Left Corner
    xy1: Bottom Right Corner.
    rxy: Radius of curvature for each axis.
    """
    # Create an SVG path for the rounded rectangle:
    path = (
        f'M {x0 + rx},{y0} '
        f'L {x1 - rx},{y0} '
        f'Q {x1},{y0} {x1},{y0 + ry} '
        f'L {x1},{y1 - ry} '
        f'Q {x1},{y1} {x1 - rx},{y1} '
        f'L {x0 + rx},{y1} '
        f'Q {x0},{y1} {x0},{y1 - ry} '
        f'L {x0},{y0 + ry} '
        f'Q {x0},{y0} {x0 + rx},{y0}'
        f'Z'
    )
    return path

def generate_html_from_screen_json(input_file : str, output_file : str):
    """ 
    Wrapper for input filepath, rather than screen data object.
    """
    input_data : ScreenData = get_screen_data_from_json(input_file)
    generate_html_from_screen_data(input_data, output_file)

def main():
    '''
    Convert a JSON output from Commec Screen, into a HTML data visualisation.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", dest="in_file",
        required=True, help="Input json file path")
    parser.add_argument("-o","--output", dest="out_file",
        required=True, help="Output html filepath, not including .html extension.")
    args = parser.parse_args()
    generate_html_from_screen_json(args.in_file, args.out_file)

if __name__ == "__main__":
    main()
