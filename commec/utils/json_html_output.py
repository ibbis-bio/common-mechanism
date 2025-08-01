"""
Used to convert a json object created from Commec Screen, or a ScreenResult object,
into a visual HTML representation of the Commec output. Which can be embedded into
any other HTML document as appropriate.
"""

import textwrap
import argparse
import plotly.graph_objects as go
import pandas as pd
import importlib
from mako.template import Template

from commec.config.json_io import get_screen_data_from_json
from commec.config.result import (
    ScreenResult,
    QueryResult,
    HitResult,
    ScreenStep,
    ScreenStatus,
)

class CommecPalette():
    """ 
    Static Container for colours used in Commec.
    """
    WHITE = [255,255,255]
    DK_BLUE = [35,42,88]
    LT_BLUE = [66,155,185]
    ORANGE = [241,80,36]
    # Yellow and Red are not official Commec colours.
    YELLOW = [241,168,29]
    RED = [207,27,81]
    BLACK = [0,0,0]

    @staticmethod
    def mod(modulate, alpha: int = 255, multiplier : float = 1.0) -> str:
        m : float = multiplier
        return str(f'rgba({round(modulate[0]*m)},{round(modulate[1]*m)},{round(modulate[2]*m)},{alpha})')

    rgba_WHITE = 'rgba(255,255,255,255)'
    rgba_DK_BLUE = 'rgba(35,42,88,255)'
    rgba_LT_BLUE = 'rgba(66,155,185,255)'
    rgba_ORANGE = 'rgba(241,80,36,255)'
    # Yellow and Red are not official Commec colours,
    rgba_YELLOW = 'rgba(241,80,36,255)'
    rgba_RED = 'rgba(207,27,81,255)'

def color_from_status(status : ScreenStatus) -> CommecPalette:
    """ 
    Convert a Screen step into an associated Colour.
    Previously, colours were based on step, now they are based on outcome.
    """
    if status >= ScreenStatus.FLAG:
        return CommecPalette.RED
    if status == ScreenStatus.WARN:
        return CommecPalette.ORANGE
    return CommecPalette.LT_BLUE

def color_from_hit(hit : HitResult) -> CommecPalette:
    """ 
    Convert a Screen step into an associated Colour.
    Previously, colours were based on step, now they are based on outcome.
    """
    return color_from_status(hit.recommendation.status)

    if hit.recommendation.from_step == ScreenStep.BIORISK:
        return CommecPalette.RED
    if (hit.recommendation.from_step == ScreenStep.LOW_CONCERN_PROTEIN or
        hit.recommendation.from_step == ScreenStep.LOW_CONCERN_RNA or
        hit.recommendation.from_step == ScreenStep.LOW_CONCERN_DNA):
        return CommecPalette.LT_BLUE
    if hit.recommendation.from_step == ScreenStep.TAXONOMY_AA:
        return CommecPalette.ORANGE
    if hit.recommendation.from_step == ScreenStep.TAXONOMY_NT:
        return CommecPalette.YELLOW
    return CommecPalette.DK_BLUE

def constrast_color_from_hit(hit : HitResult):
    """ 
    Convert a Screen step into an associated Colour.
    This is now always white, which works well with all outcome colours.
    As this colour is meant for text, and not box objects, its associated
    allowed values are restricted to hex, or built-in strings.
    """
    return "#FFFFFF"
    if (hit.recommendation.from_step == ScreenStep.TAXONOMY_NT):
        return "#000000"
    return "#FFFFFF"

def overlay_text_from_hit(hit : HitResult):
    outcome = str(hit.recommendation.status)
    outcome = outcome.replace("Warning", "Warn")
    outcome = outcome.replace("(Cleared)","")
    step = hit.recommendation.from_step
    match step:
        case ScreenStep.BIORISK:
            return "Biorisk "+outcome
        case ScreenStep.TAXONOMY_AA:
            return "Protein "+outcome
        case ScreenStep.TAXONOMY_NT:
            return "Nucl. "+outcome
        case ScreenStep.LOW_CONCERN_PROTEIN:
            return "Benign Protein"
        case ScreenStep.LOW_CONCERN_RNA:
            return "Benign RNA"
        case ScreenStep.LOW_CONCERN_DNA:
            return "Benign DNA"

def generate_html_from_screen_data(input_data : ScreenResult, output_file : str):
    """
    Interpret the ScreenResult from Commec Screen output as a visualisation, in 
    the form on an HTML output.
    If Commec Screen handled multiple Queries, then they are combined in the HTML output.
    """

    figures_html = []

    # Render each query as its own Plotly HTML visualisation:
    for i, query in enumerate(input_data.queries.values()):
        fig = go.Figure()
        vertical_stack_count : int = draw_query_to_plot(fig, query)
        update_layout(fig, query, vertical_stack_count)
        file_name = output_file.strip()+"_"+str(i)+".html"
        html = fig.to_html(file_name, full_html = False, include_plotlyjs='cdn')
        figures_html.append(html)

    # Additional template information
    n_query = input_data.query_info.number_of_queries
    plural = "y" if n_query == 1 else "ies"
    html_title = f"Commec Screen Summary: {n_query} Quer{plural}."

    # Construct the composite HTML
    template_path = str(importlib.resources.files("commec").joinpath("utils").joinpath("template.html"))
    template = Template(filename = template_path)
    rendered_html = template.render(figures_html=figures_html, page_title=html_title)

    # Save the combined HTML output
    output_filename = output_file.strip()+".html"
    with open(output_filename, "w", encoding = "utf-8") as output_file:
        output_file.write(rendered_html)


def update_layout(fig, query_to_draw : QueryResult, stacks):
    """ 
    Applies some default settings to the plotly figure,
    also adjusts the figure height based on the number of vertically stacking bars.
    """

    figure_base_height = 180
    figure_stack_height = 30

    r,g,b = color_from_status(query_to_draw.status.screen_status)
    css_color = f"rgb({r},{g},{b})"

    fig.update_layout(showlegend=False)
    # Update layout to display X-axis on top and hide Y-axis labels for specified subplot
    fig.update_layout({
        # General layout properties
        'height': figure_base_height + (figure_stack_height * stacks),
        'title' : (
            f"<span style='background-color:{css_color};padding:6px 10px;"
            f"border-radius:6px;color:{css_color};font-weight:bold;'>"
            f"{query_to_draw.status.screen_status}</span> : {query_to_draw.query} ({query_to_draw.length} b.p.) <br>{query_to_draw.status.rationale}"
        ),
        'barmode': 'overlay',
        'template': 'plotly_white',
        'plot_bgcolor': 'rgba(0,0,0,0)',  # Transparent plot area
        'paper_bgcolor': 'rgba(0,0,0,0)',  # Transparent outer area
        "xaxis": dict(
            #title="Query Basepairs (bp)",
            showline=True,
            constrain="domain",
            ticks='outside',
            showgrid=True,
            fixedrange=True,
            zeroline=False,
            side="top",
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
        'bargap': 0.01,
    })


def generate_outcome_string(query : QueryResult, hit : HitResult) -> str:
    """
    Takes a Query, and associated hit, and formats a human readable output string,
    handling associated construction logic.
    """
    if "Biorisk" in hit.recommendation.from_step:
        if "regulated" in hit.annotations:
            return hit.annotations["regulated"][0]
        return ""
    
    if "Benign" in hit.recommendation.from_step:
        return ""

    # If Blast NR or NT result:
    if "Taxonomy" in hit.recommendation.from_step:
        output_string = ""

        n_regulated_eukaryotes = 0
        n_regulated_bacteria = 0
        n_regulated_viruses = 0
        regulated_taxids = []
        non_regulated_taxids = []
        regulated_species = []

        if "regulated_taxonomy" in hit.annotations:
            reg = hit.annotations["regulated_taxonomy"]
            for r in reg:
                n_regulated_eukaryotes += int(r["regulated_eukaryotes"])
                n_regulated_bacteria += int(r["regulated_bacteria"])
                n_regulated_viruses += int(r["regulated_viruses"])
                regulated_taxids.extend(r["regulated_taxids"])
                non_regulated_taxids.extend(r["non_regulated_taxids"])
                regulated_species.extend(r["regulated_species"])
            
            domain_string = " across multiple domains."
            if n_regulated_bacteria > 0 and n_regulated_viruses == 0 and n_regulated_eukaryotes == 0:
                domain_string = " bacteria"
            elif n_regulated_eukaryotes > 0 and n_regulated_viruses == 0:
                domain_string = " eukaryotes"
            elif n_regulated_viruses > 0:
                domain_string = " viruses"

            total_hits : int = len(regulated_taxids) + len(non_regulated_taxids)
            peptide_type = "protein" if "Protein" in hit.recommendation.from_step else "nucleotides"

            if len(non_regulated_taxids) > 0:
                output_string += "Mix of regulated and non-regulated" + domain_string + ".<br>"
            else:
                output_string += "Best match to regulated" + domain_string + ".<br>"

            output_string += str(total_hits) + " Best match hits to " if total_hits > 1 else "Best hit to "
            output_string += peptide_type + " found in " + str(len(regulated_species)) + " species, with regulated pathogen taxId in "
            output_string += "lineage " if len(regulated_species) == 1 else "lineages "

            species_list = textwrap.fill(
                ", ".join(regulated_species), 100
            ).replace("\n", "<br>")
            
            output_string += "<br>" + species_list

            return output_string
            #Best match to regulated viruses. 2 best match hits to nucleotides found in 1 species, 2 with regulated pathogen taxId in lineage (Influenza A)
        return "No Annotations."

#@pytest.mark.filterwarnings("ignore")
def draw_query_to_plot(fig : go.Figure, query_to_draw : QueryResult):
    """ 
    Write the data from a single query into the figure for plotly. 
    """
    # Interpret the QueryResult into bars for the plot.
    graph_data = [
        {"label": query_to_draw.query, 
         "label_verbose": query_to_draw.query, 
         "outcome" : f"Commec Recommendation for this query: {query_to_draw.status.screen_status}",
         "outcome_verbose":f"{query_to_draw.status.rationale}",
         "start": 1, "stop": query_to_draw.length,
         "color" : CommecPalette.DK_BLUE, 
         "stack" : 0,
         "text" : query_to_draw.query[:25],
         "text_color" : "#FFFFFF",}
    ]

    # Keep track of how many vertical stacks this image has.
    n_stacks = 1

    for hit in query_to_draw.hits.values():
        for match in hit.ranges:
            # Find the best vertical position to reduce collisions, and fill all space.
            collision_free = False
            stack_write = 1 # 1 gives a bit of space between the query and the data.
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
                    "label" : hit.description[:25] + "...",
                    "label_verbose" : hit.description[:],
                    "outcome" : f"{hit.recommendation.status} from {hit.recommendation.from_step}",
                    "outcome_verbose" : generate_outcome_string(query_to_draw, hit),
                    "start" : match.query_start,
                    "stop" : match.query_end,
                    "color" : color_from_hit(hit),
                    "stack" : stack_write,
                    "text" : overlay_text_from_hit(hit),
                    "text_color" : constrast_color_from_hit(hit),
                }
            )

    df = pd.DataFrame(graph_data)

    # Convert RGB colors to hex format
    def rgb_to_hex(rgb):
        return '#{:02x}{:02x}{:02x}'.format(rgb[0], rgb[1], rgb[2])
    df['color'] = df['color'].apply(rgb_to_hex)


    df['hovertext'] = df.apply(lambda bar_data: 
                               f"{bar_data["outcome"]}<br>"
                               f"bases {bar_data['start']}-{bar_data['stop']}<br>"
                               f"{bar_data['label_verbose']}<br>"
                               f"{bar_data["outcome_verbose"]}", axis=1)

    # Add each bar to the existing figure
    for _i, bar_data in df.iterrows():

        marker_dictionary = {"color" : bar_data["color"]}
        if "Cleared" in bar_data["outcome"]:
            marker_dictionary["pattern"] = dict(shape="/")  # Options: "/", "\\", "x", "-", "|", ".", "+"

        fig.add_trace(
            go.Bar(
                x=[bar_data['stop'] - (bar_data['start'] - 1)],  # Width of bar as timedelta
                y=[bar_data['stack']],  # Stack position
                base=bar_data['start']-1,  # Starting point on x-axis
                orientation='h',
                marker=marker_dictionary,
                hovertext=bar_data['hovertext'],
                hoverinfo='text',
                #customdata=df['clicktext'],
                name=bar_data['label'],
                width = 1.0,
            ),
        )

        # Compute x center for horizontal bar
        bar_center_x = (bar_data['start']-1) + (bar_data['stop'] - (bar_data['start']-1)) / 2
        bar_y = bar_data['stack']

        fig.add_annotation(
            x=bar_center_x,
            y=bar_y,
            text=f"<b>{bar_data['text']}</b>",
            showarrow=False,
            font=dict(
                color=bar_data['text_color'],
                size=13,
            ),
            xanchor='center',
            yanchor='middle',
            align='center',
            xref='x',
            yref='y',
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
    input_data : ScreenResult = get_screen_data_from_json(input_file)
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
