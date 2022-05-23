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
    graphs = {graph["ref"]: graph for graph in graphs}
    graphs = [graphs[ref] for cluster in clusters for ref in cluster]
    
    cols = ["start", "end", "parent", "ref"]

    data = pd.DataFrame([
        vals for g in graphs
        for vals in zip(*[g.vs[col] for col in cols])
    ], columns=cols)

    genome_sizes = {g["ref"]: g["size"] for g in graphs}
    
    if norm:
        to_add = []

        for ref in data.ref.unique():
            data_t = data[data.ref == ref].sort_values(by=["start", "end"])
            last_idx = data_t.index[-1]
            extra = data_t.end.iloc[-1] - genome_sizes[ref]

            if extra > 0:
                data.loc[last_idx, 'end'] = genome_sizes[ref]
                
                to_add.append((0, extra,
                               data_t.loc[last_idx, 'parent'],
                               data_t.loc[last_idx, 'ref']))

        offset = data.index.max() + 1
        for i, entry in enumerate(to_add):
            data.loc[i+offset] = entry
                
        data = data.sort_values(by=['ref', 'start', 'end']).reset_index(drop=True)

    hover = HoverTool(tooltips=[('Parent', '@parent')])
    
    legend_opts = dict(
        glyph_height=15,
        label_height=15,
     )

    plot_opts = dict(
        width=700, height=max(50*data.ref.nunique(), 200),
        xlabel='Position', ylabel='Phage',
        gridstyle=dict(ygrid_line_color='gray', xgrid_line_alpha=0, ygrid_line_dash=[4, 4]),
        show_grid=True,
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
        data_c = data[data.ref.isin(cluster)].copy()
        data_c['ref_loc'] = pd.Categorical(data_c.ref).codes

        parents_all = {p for p in data_c.parent if "NA" not in p}
       
        subplot_layer = []

        for i, p in enumerate(parents_all):
            data_c_p = data_c[data_c.parent==p].copy()

            # Add an offset to better see breakpoint overlap
            data_c_p.ref_loc = data_c_p.ref_loc + 0.05*i
            
            subplot_layer.append(
                hv.Segments(
                    data_c_p, ['start', 'ref_loc', 'end', 'ref_loc'], label=p
                )
                .opts(line_width=15, color=cmap[p],
                      tools=[hover],
                      default_tools=["box_zoom", "reset"])
            )

        yaxis_pos = range(data_c.ref.nunique())
        yaxis_ticks = data_c.ref.unique()
        
        overlay = (
            hv.Overlay(subplot_layer)
            .opts(**plot_opts)
            .relabel("Phage").opts(yticks=list(zip(yaxis_pos, yaxis_ticks)))
        )

        subplots.append(overlay)
        
    seg = hv.Layout(subplots).opts(shared_axes=True).cols(1)

    hv.save(seg, Path(outdir, 'paintings.html'))
