#! /usr/bin/env python
import math
import pandas as pd
import networkx as nx

from bokeh.plotting import figure
from bokeh.models import Circle, Text, Row, Column, ColumnDataSource, TapTool, CustomJS, MultiLine
from bokeh.models.graphs import from_networkx
from bokeh.palettes import viridis
from bokeh.transform import factor_cmap
from bokeh.io import output_file, save
from bokeh.models.widgets import Div

import click


def get_palette(size):
    return viridis(max(2, size))


def make_cluster_plot(data):
    data = data.copy()
    data["phased"] = data.apply(lambda x: 0 if x["phase"] == -1 else 1, axis=1)
    factors = [str(c) for c in sorted(list(set(data["cluster"])))]
    data["cluster"] = data["cluster"].astype(str)

    TOOLTIPS = [
        ("Passes", "@np"),
        ("length", "@len"),
        ("qual", "@rq"),
    ]

    p=figure(tooltips=TOOLTIPS)

    p.circle(
        x="x",
        y="y",
        source=ColumnDataSource(data),
        size="np",
        color=factor_cmap("cluster", palette=get_palette(len(factors)), factors=factors),
        fill_alpha=0.5,
        line_alpha="phased",
        legend="cluster",
        name="sequences"
    )


    p.toolbar.logo = None
    return p


def make_graph(data):
    counts = data.groupby("cluster").count()
    color_map = get_palette(counts.shape[0])

    scale = 25 / counts.shape[0]

    G = nx.DiGraph()
    G.add_node(
        "root",
        n=data.shape[0],
        size=scale * math.sqrt(data.shape[0] / math.pi),
        name="root\n({})".format(data.shape[0]),
        color="#C0C0C0",
        cluster=-1,
        phase=-1
    )

    for cluster in counts.index:
        node = "cluster{}".format(cluster)
        n = counts.loc[cluster].values[0]
        G.add_node(
            node,
            n=n,
            size=scale * math.sqrt(n / math.pi),
            name="{}\n({})".format(node, n),
            color=color_map[cluster],
            cluster=cluster,
            phase=-1
        )
        G.add_edge("root", node)

        phasing = data[data["cluster"]==cluster].groupby("phase").count()
        phases = [x for x in phasing.index if x >= 0]
        for phase in phases:
            child = "{}.phase{}".format(node, phase)
            n=phasing.loc[phase].values[0]
            G.add_node(
                child,
                n=n,
                size=scale * math.sqrt(n / math.pi),
                name="phase{}\n({})".format(phase, n),
                color=color_map[cluster],
                cluster=cluster,
                phase=phase
            )
            G.add_edge(node, child)

    return G


def get_limits(layout):
    xs = [l[0] for _, l in layout.items()]
    ys = [l[1] for _, l in layout.items()]

    x_lim = (min(xs), max(xs))
    y_lim = (min(ys), max(ys))

    x_range = max(1, x_lim[1] - x_lim[0])
    y_range = max(1, y_lim[1] - y_lim[0])

    x_pad = 0.2 * x_range
    y_pad = 0.2 * y_range

    x_lim = (x_lim[0] - x_pad, x_lim[1] + x_pad)
    y_lim = (y_lim[0] - y_pad, y_lim[1] + y_pad)

    return (x_lim, y_lim)


def make_phase_plot(data):
    G = make_graph(data)
    layout = nx.nx_agraph.graphviz_layout(G, prog="twopi", root="root")
    (x_lim, y_lim) = get_limits(layout)

    p = figure(x_range=x_lim, y_range=y_lim)

    graph = from_networkx(G, layout, center=(0,0))
    graph.node_renderer.glyph = Circle(size="size", fill_color="white", line_color="white")
    graph.edge_renderer.glyph = MultiLine(line_color="black", line_alpha=0.5, line_width=1)
    graph.name="base_graph"
    p.renderers.append(graph)

    graph2 = from_networkx(G, layout, center=(0,0))
    graph2.node_renderer.glyph = Circle(size="size", fill_color="color", fill_alpha=0.5, line_color="color")
    graph2.edge_renderer.visible=False
    graph2.name="nodes"
    p.renderers.append(graph2)

    graph3 = from_networkx(G, layout, center=(0,0))
    graph3.node_renderer.glyph = Text(text="name", text_align="center", text_baseline="middle", text_font_size="9pt")
    graph3.edge_renderer.visible=False
    graph3.name="labels"
    p.renderers.append(graph3)

    p.tools = []
    p.toolbar.logo = None
    p.xaxis.visible = False
    p.yaxis.visible = False
    p.xgrid.visible = False
    p.ygrid.visible = False

    return p


def add_phaseplot_callback(clusterplot, phaseplot):
    clusterplot_datasource = clusterplot.select_one("sequences").data_source
    phaseplot_datasource = phaseplot.select_one("nodes").node_renderer.data_source

    callback = CustomJS(args=dict(cp_source=clusterplot_datasource, pp_source=phaseplot_datasource), code="""
        var inds = pp_source.selected.indices;
        var p_data = pp_source.data;
        var c_data = cp_source.data;
        cp_source.selected.indices = [];

        if (inds.length > 0)
        {
            for(var i=0; i<c_data['cluster'].length; i++)
            {
                var clusterA = c_data['cluster'][i];
                var phaseA = c_data['phase'][i];
                for (var j=0; j<inds.length; j++)
                {
                    clusterB = p_data['cluster'][inds[j]];
                    if (clusterA == clusterB)
                    {
                        phaseB = p_data['phase'][inds[j]];

                        if ((phaseB != -1 && phaseA == phaseB) || (phaseB == -1))
                        {
                            cp_source.selected.indices.push(i);
                            break;
                        }
                    }
                }
            }
        }

        cp_source.change.emit();
        """)
    phaseplot.add_tools(TapTool(renderers=[phaseplot.select_one("nodes")]))
    phaseplot_datasource.selected.js_on_change('indices', callback)


def make_plot(data, title):
    clusterplot = make_cluster_plot(data)
    phaseplot = make_phase_plot(data)
    add_phaseplot_callback(clusterplot, phaseplot)
    title_div = Div(text=title)

    return Column(
        title_div,
        Row(clusterplot, phaseplot)
    )


@click.command(context_settings=dict(
    ignore_unknown_options=True,
    ),
    short_help="Generate an html plot of the clustering"    
)
@click.option("--title", "-t", type=str, default="Clusters",
              help="Title for the plot")
@click.argument("embedding_file", type=click.Path(exists=True))
@click.argument("info_file", type=click.Path(exists=True))
@click.argument("plot_file", type=click.Path())
def cli_handler(title, embedding_file, info_file, plot_file):
    output_file(plot_file, title="Cluster report", mode="inline")
    try:
        embeddings = pd.read_csv(embedding_file, sep="\t", header=None)
        bam_info = pd.read_csv(info_file, sep="\t")
    except pd.errors.EmptyDataError:
        save(Div(text="No Data"))
        return

    embeddings.columns = ["x", "y"]
    data = bam_info.join(embeddings).groupby("cluster").filter(lambda x: len(x) > 1).copy()
    if data.shape[0] > 0:
        save(make_plot(data, title))
    else:
        save(Div(text="No Data"))


if __name__ == '__main__':
    cli_handler()
