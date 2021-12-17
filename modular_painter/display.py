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


def access(L, i):
    if isinstance(L, list):
        if i < len(L):
            return L[i]
        return None
    return L

def get_palette(n):
    palette = linear_palette(Turbo256, n)
    palette = np.random.choice(palette, n, replace=False).tolist()

    return palette

def normalize(x, start, end):
    return 100 * (x-start)/(end-start)

def display_genomes(genomes, clusters=None, norm=True):
    data = pd.concat({
        genome.target: genome.to_pandas()
        for genome in genomes
    }).reset_index(level=1, drop=True).rename_axis(index='target').reset_index()
    
    if norm:
        to_add = []

        for target in data.target.unique():
            data_t = data[data.target == target]
            first_idx = data_t.index[0]
            extra = abs(data.loc[first_idx, 'start'] - 1)
            # extra = data.loc[last_idx, 'end'] - genomes[target].mod

            if extra > 0:
                data.loc[first_idx, 'start'] = 1
                to_add += [(data_t.loc[first_idx, 'target'],
                            data_t.loc[first_idx, 'parent'],
                            genomes[target].mod-extra+1,
                            genomes[target].mod)]

        first_idx = data.index.max() + 1
        for i, entry in enumerate(to_add):
            data.loc[i+first_idx] = entry
                
        data = data.sort_values(by=['target', 'start', 'end']).reset_index(drop=True)
            
    hover = HoverTool(tooltips=[('Parent', '@parent')])
    
    legend_opts = dict(
        # label_text_font_size="12px",
        # label_text_line_height=20, # doesnt seem to be doing anything
        # spacing=0,
        glyph_height=15,
        label_height=15,
     )

    plot_opts = dict(
        width=700, height=40*data.target.nunique(),
        xlabel='Position', ylabel='Phage',
        legend_position="right",
        legend_limit=50,
        legend_offset=(30, 0),
        legend_opts=legend_opts,
    )

    cmap = {}
    i = 0
    
    for parent in data.parent.unique():
        if 'NoCov' in parent:
            cmap[parent] = 'gray'
        else:
            cmap[parent] = cc.glasbey_light[i]
            i += 1

    subplots = []

    for cluster in clusters:
        data_c = data[data.target.isin(cluster)]
        data_c['target_loc'] = pd.Categorical(data_c.target).codes

        parents_all = data_c.parent.unique()
        
        subplot_layer = []

        for i ,s in enumerate(sorted(data_c.parent.unique())):
            data_c_s = data_c[data_c.parent==s]

            # Add an offset to better see breakpoint overlap
            data_c_s.target_loc = data_c_s.target_loc + 0.03*i
            
            subplot_layer.append(
                hv.Segments(
                    data_c_s, ['start', 'target_loc', 'end', 'target_loc'], label=s
                )
                .opts(line_width=15, color=cmap[s],
                      tools=[hover],
                      default_tools=["box_zoom", "reset"])
            )

        yaxis_pos = range(data_c.target.nunique())
        yaxis_ticks = data_c.target.unique()
        
        overlay = (
            hv.Overlay(subplot_layer)
            .opts(**plot_opts)
            .relabel("Phage").opts(yticks=list(zip(yaxis_pos, yaxis_ticks)))
        )

        subplots.append(overlay)
        
    seg = hv.Layout(subplots).opts(shared_axes=True).cols(2)

    hv.save(seg, '/tmp/cedric/modular_painting_tests/paintings.html')
