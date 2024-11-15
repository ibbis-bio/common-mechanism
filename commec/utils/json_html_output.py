"""
Used to convert a json object created from Commec Screen, or a ScreenData object,
into a visual HTML representation of the Commec output. Which can be embedded into
any other HTML document as appropriate.
"""

import os
import argparse
import plotly.graph_objects as go
import pandas as pd
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

def color_from_hit(hit : HitDescription) -> CommecPalette:
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
    If Commec Screen handled multiple Queries, then they are combined in the HTML output.
    """

    plot_filenames = []
    figures_html = []

    # Render each query as its own Plotly HTML visualisation:
    for i, query in enumerate(input_data.queries):
        fig = go.Figure()
        vertical_stack_count : int = draw_query_to_plot(fig, query)
        update_layout(fig, query, vertical_stack_count)
        file_name = output_file.strip()+"_"+str(i)+".html"
        html = fig.to_html(file_name, full_html = False, include_plotlyjs='cdn')
        figures_html.append(html)
        #plot_filenames.append(file_name)
        
    # Read each Plotly HTML file content
    #for filename in plot_filenames:
    #    with open(filename, "r", encoding = "utf-8") as file:
    #        figures_html.append(file.read())

    # Construct the composite HTML
    template = Template(filename="commec/utils/template.html")
    rendered_html = template.render(figures_html=figures_html)

    # Save the combined HTML output
    output_filename = output_file.strip()+".html"
    with open(output_filename, "w", encoding = "utf-8") as output_file:
        output_file.write(rendered_html)

    # Remove intermediaries
    #for file in plot_filenames:
        #os.remove(file)

def update_layout(fig, query_to_draw, stacks):
    """ 
    Applies some default settings to the plotly figure,
    also adjusts the figure height based on the number of vertically stacking bars.
    """

    figure_base_height = 180
    figure_stack_height = 30

    # Update layout to display X-axis on top and hide Y-axis labels for specified subplot
    fig.update_layout({
        # General layout properties
        'height': figure_base_height + (figure_stack_height * stacks),
        'title': f"Query: {query_to_draw.query}  ({query_to_draw.length} b.p.)",
        'barmode': 'overlay',
        'template': 'plotly_white',
        'plot_bgcolor': 'rgba(0,0,0,0)',  # Transparent plot area
        'paper_bgcolor': 'rgba(0,0,0,0)',  # Transparent outer area
        # Update specific xaxis and yaxis only for the specified subplot
        "xaxis": dict(
            title="Query Basepairs (bp)",
            showline=True,
            constrain="domain",
            ticks='outside',
            showgrid=True,
            fixedrange=True,
            zeroline=False,
            # Uncomment if time format and linear mode are needed
            #tickformat="%s",
            #type="date",
            #tickmode="linear",
        ),
        "yaxis": dict(
            title="",
            showticklabels=False,
            autorange='reversed',
            range=[-0.5, stacks + 0.5],
            fixedrange=True,
            tickmode="linear",
            zeroline=False,
        ),
        'bargap': 0.0,
    })

def draw_query_to_plot(fig : go.Figure, query_to_draw : QueryData):
    """ 
    Write the data from a single query into the figure for plotly. 
    """
    # Interpret the QueryData into bars for the plot.
    graph_data = [
        {"label": query_to_draw.query, "outcome" : "", "start": 0, "stop": query_to_draw.length, "color" : CommecPalette.DK_BLUE, "stack" : 0},
    ]

    # Keep track of how many vertical stacks this image has.
    n_stacks = 1

    for hit in query_to_draw.hits:
        for match in hit.ranges:
            # Find the best vertical position to reduce collisions, and fill all space.
            collision_free = False
            stack_write = 0
            while not collision_free:
                stack_write += 1
                collision_free = True
                for entry in graph_data:
                    if entry["stack"] == stack_write:
                        collision_free = (collision_free and
                                            (match.query_start > entry["stop"] 
                                            or match.query_end < entry["start"])
                                            )

            n_stacks = max(n_stacks, stack_write + 1)

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

    df = pd.DataFrame(graph_data)

    # Convert RGB colors to hex format
    def rgb_to_hex(rgb):
        return '#{:02x}{:02x}{:02x}'.format(rgb[0], rgb[1], rgb[2])
    df['color'] = df['color'].apply(rgb_to_hex)

    # Generate hover text
    df['hovertext'] = df.apply(lambda bar_data: f"{bar_data['label']}<br>({bar_data['start']}-{bar_data['stop']})", axis=1)
    df['clicktext'] = df.apply(lambda bar_data: f"{bar_data['label']}<br>({bar_data['start']}-{bar_data['stop']})<br>"
                               "This is some special custom text to be shown when you click on something.", 
                               axis=1)

    # Add each bar to the existing figure
    for _i, bar_data in df.iterrows():
        fig.add_trace(
            go.Bar(
                x=[bar_data['stop'] - bar_data['start']],  # Width of bar as timedelta
                y=[bar_data['stack']],  # Stack position
                base=bar_data['start'],  # Starting point on x-axis
                orientation='h',
                marker=dict(color=bar_data['color']),
                hovertext=bar_data['hovertext'],
                hoverinfo='text',
                customdata=df['clicktext'],
                name=bar_data['label'],
                width = 1.0,
            ),
        )

    return n_stacks

def generate_rounded_rect(x0,y0,x1,y1,rx,ry):
    """ 
    Creates an SVG path for a rounded rectangle given the following dimensions:
    xy0: Top Left Corner
    xy1: Bottom Right Corner.
    rxy: Radius of curvature for each axis.

    Currently unused.
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
