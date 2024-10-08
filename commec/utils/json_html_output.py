"""
Used to convert a json object created from Commec Screen, or a ScreenData object,
into a visual HTML representation of the Commec output. Which can be embedded into
any other HTML document as appropriate.
"""

from enum import StrEnum
import argparse
import plotly.graph_objects as go
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
    fig.add_trace(go.Scatter(
        x=[0, query_to_draw.length],  # The x-values; these can be any value
        y=[0, 0],  # The y-values; these can be zeros or any dummy values
        name='Invisible Trace',  # Optional name
        xaxis='x2',  # Assign to the secondary x-axis
        mode='lines',  # Use lines to create the trace
        line=dict(color='rgba(255, 255, 255, 0)'),  # Fully transparent line
        hoverinfo='none',  # No hover information
        showlegend=False,
    ))

    # Update layout to display X-axis on top and hide Y-axis labels
    fig.update_layout(
        height= 100 * (stacks + 1),
        barmode='stack',
        title="Query: " + query_to_draw.query + "  (" + str(query_to_draw.length)+" b.p.)",
        xaxis_title="",
        yaxis_title="",
        yaxis=dict(
            showticklabels=False,
            autorange='reversed',
            fixedrange=True,
        ),
        template='plotly_white',
        
        # First x-axis (default)
        xaxis=dict(
            showline=True,
            range=[0, query_to_draw.length],
            constrain="domain",
            ticks='outside',
            showgrid=True,
            fixedrange=True,
        ),
        
        # Second x-axis (xaxis2) at the top
        # the issue with having two axis is they don't necessarily scroll together.
        # We could try draw the same thing to both axis, and hope for the best.
        xaxis2=dict(
            overlaying='x',  # Overlay it on the same domain as xaxis
            side='top',  # Place it on top
            showline=True,  # Show the X-axis line
            constrain="domain",  # Constrain the x-axis to its domain
            range=[0, query_to_draw.length],  # Same range as xaxis (or modify if needed)
            ticks='outside',  # Show ticks outside
            showgrid=False,  # Hide gridlines
            fixedrange=True,
        )
    )

    # Show the plot
    fig.write_html(output_file.strip()+".html")


def draw_query_to_plot(fig : go.Figure, query_to_draw : QueryData):
    print(query_to_draw.query, " with ", str(len(query_to_draw.hits)), " hits.")

    graph_data = [
        {"label": query_to_draw.query, "outcome" : "", "start": 0, "stop": query_to_draw.length, "color" : CommecPalette.DK_BLUE},
    ]

    for hit in query_to_draw.hits:
        for match in hit.ranges:
            graph_data.append(
                {
                    "label" : hit.recommendation.from_step + " " + hit.description[:25]+"...",
                    "outcome" : hit.recommendation.outcome + ": " + hit.description[:25] + "...",
                    "start" : match.query_start,
                    "stop" : match.query_end,
                    "color" : color_from_hit(hit),
                }
            )

    # Prepare the data for the stacked bar plot
    labels = [f'{r["label"]}_{i}' for i, r in enumerate(graph_data)]
    #outcomes = [r["outcome"] for r in graph_data]
    data_ranges = [r["stop"] - r["start"] for r in graph_data]  # The data ranges (stop - start)
    colors = [r["color"] for r in graph_data]  # The data ranges (stop - start)

    for r in graph_data:
        print(r)

    # Add the actual data ranges with individual colors and text inside the bars
    write_stack : int = 0
    current_pos : int = 0
    for i, r in enumerate(graph_data):

        # Naive approach adds to the same column, if range comes after.
        if r["start"] < current_pos:
            write_stack+=1
            current_pos=0
        
        # Calculate the gap size, and add it if required.
        # (This is secretly a stacked horizontal bar graph)
        gap_size = r['start'] - current_pos
        if gap_size > 0:
            add_gap(fig, gap_size, write_stack)

        current_pos = r['stop']

        fill_colour =CommecPalette.mod(colors[i], multiplier=1.0)
        line_colour = CommecPalette.mod(colors[i], multiplier=0.5)

        # Add the rounded rectangle as a path shape
        #fig.add_shape(
        #    layer='below',
        #    type="path",
        #    path=generate_rounded_rect(r["start"], i-0.4, r["stop"], i+0.4, 25, 0.3),
        #    line=dict(color=line_colour),
        #    fillcolor=fill_colour,
        #)

        fig.add_trace(go.Bar(
            x=[data_ranges[i]],  # Length of the data ranges (on x-axis)
            y=[write_stack],  # Labels for each range
            orientation='h',  # Horizontal bars
            name=r["label"],
            #marker={"color":fill_colour},  # Assign different color for each range
            marker={"color":fill_colour,
                            "line": {
                                "color": line_colour,  # Make border line invisible
                                "width": 1  # Set line width (optional, set to 0 if no border is desired)
                                    }
                    },  # Make gaps invisible
            hoverinfo='text',
            hovertext=str(labels[i]) + "<br>("+str(r["start"])+"-"+str(r["stop"])+")",
            #hoverinfo='x+y+name',  # Show hover info
            #hovertemplate="%{name} <br>("+str(r["start"])+"-"+str(r["stop"])+")",
            text=r["outcome"],  # Text to display inside the bar
            textposition='inside',  # Position text inside the bar
            insidetextanchor='start', # Align text to the start of the bar (left-aligned)
            #hovertext="",
        ),
        )
    return write_stack

def add_gap(fig : go.Figure, gap_size : int, y_axis):
    """ Add an invisible stack to the bar graph, emulating a gap. """
    fig.add_trace(go.Bar(
        x=[gap_size],  # Gap length (on x-axis)
        y=[y_axis],  # Labels for each range (on y-axis, which will be horizontal due to orientation)
        orientation='h',  # Horizontal bars
        name="Gap",
        marker={"color":'rgba(255, 255, 255, 0)',
                "line": {"color": 'rgba(0, 0, 0, 0)',
                        "width": 0}},
        hoverinfo='skip',  # Don't show hover info for the gaps
        showlegend=False
    ))

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
