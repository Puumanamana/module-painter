from pathlib import Path
from math import pi
import numpy as np
import pandas as pd

import numpy as np
from bokeh.models import HoverTool
from bokeh.palettes import Turbo256, linear_palette

import holoviews as hv
from holoviews import dim
import colorcet as cc
hv.extension('bokeh')


def get_palette(n):
    palette = linear_palette(Turbo256, n)
    palette = np.random.choice(palette, n, replace=False).tolist()

    return palette

def custom_cmap():
    return dict(
        A="#F1C40F", # yellow/orange
        B="#40E0D0", # light blue
        C="#008080", # teal
        D="#FF0000", # red
        E="#808000", # olive
        F="#FF7F50", # orange
        G="#E98DFF", # purple
        H="#008000", # green
        I="#0000FF", #  
        J="#3498DB", # blue
        K="#641E16", # dark red
        L="#800000", # brown
        M="#C911F6" # pink
    )

def normalize(x, start, end):
    return 100 * (x-start)/(end-start)

def display_genomes(graphs, clusters=None, norm=True, outdir=None):
    graphs = {graph["sacc"]: graph for graph in graphs}
    graphs = [graphs[sacc] for cluster in clusters for sacc in cluster]
    
    cols = ["sstart", "send", "parent", "sacc", "slen"]

    data = pd.DataFrame([
        vals for g in graphs
        for vals in zip(*[g.vs[col] for col in cols])
    ], columns=cols)

    if norm:
        data["extra"] = data.send - data.slen

        if any(data.extra > 0):
            data.loc[data.extra > 0, "send"] = data.slen
            to_add = data[data.extra > 0].copy()
            to_add["sstart"] = 0
            to_add["send"] = data.extra

            data = pd.concat([data, to_add]).sort_values(by=['sacc', 'sstart', 'send']).reset_index(drop=True)

    hover = HoverTool(tooltips=[(x, f"@{x.lower()}") for x in ["parent", "sacc", "sstart", "send"]])
    
    legend_opts = dict(
        glyph_height=15,
        label_height=15,
     )

    # max_cluster_size = max(len(c) for c in clusters)
    plot_opts = dict(
        width=900, #height=max(50*max_cluster_size, 200),
        xlabel='Position', ylabel='Phage',
        gridstyle=dict(ygrid_line_color='gray', xgrid_line_alpha=0, ygrid_line_dash=[4, 4]),
        show_grid=True,
        show_legend=False,
        legend_position="right",
        legend_limit=50,
        legend_offset=(30, 0),
        legend_opts=legend_opts,
    )

    cmap = custom_cmap()
    
    for parent in data.parent.unique():
        if 'NA' in parent:
            cmap[parent] = 'gray'
    remaining_parents = set(data.parent.unique()).difference(cmap.keys())
    colors = get_palette(len(remaining_parents))
    cmap.update(dict(zip(remaining_parents, colors)))

    subplots = []

    for cluster in clusters:
        # if len(cluster) == 1:
        #     continue
        data_c = data[data.sacc.isin(cluster)].copy()
        data_c['sacc_loc'] = pd.Categorical(data_c.sacc).codes

        parents_all = {p for p in data_c.parent if "NA" not in p}
       
        subplot_layer = []

        for i, p in enumerate(parents_all):
            data_c_p = data_c[data_c.parent==p].copy()

            # Add an offset to better see breakpoint overlap
            data_c_p.sacc_loc = data_c_p.sacc_loc + 0.05*i
            
            subplot_layer.append(
                hv.Segments(
                    data_c_p, ['sstart', 'sacc_loc', 'send', 'sacc_loc'], label=p
                )
                .opts(line_width=15, color=cmap[p],
                      tools=[hover],
                      height=len(cluster)*100 + 50*len(clusters),
                      default_tools=["box_zoom", "reset"])
            )

        yaxis_pos = range(data_c.sacc.nunique())
        yaxis_ticks = data_c.sacc.unique()
        # yaxis_ticks = [x if len(x) < 20 else x[:20] + "..." for x in data_c.sacc.unique()]
        
        overlay = (
            hv.Overlay(subplot_layer)
            .opts(**plot_opts)
            .relabel("Phage").opts(yticks=list(zip(yaxis_pos, yaxis_ticks)))
        )

        subplots.append(overlay)
        
    seg = hv.Layout(subplots).opts(shared_axes=True).cols(1)

    hv.save(seg, Path(outdir, 'paintings.html'))
