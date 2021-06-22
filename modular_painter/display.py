from math import pi
import numpy as np
import pandas as pd

import numpy as np
from bokeh.models import HoverTool

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

def display_genomes(genomes):
    data = pd.concat({
        name: genome.to_pandas()
        for (name, genome) in genomes.items()
    }).reset_index(level=1, drop=True).rename_axis(index='target').reset_index()

    starts = data.groupby('target').start.min()
    ends = data.groupby('target').end.max()
    
    data.start = data.groupby('target').start.transform(lambda x: normalize(x, starts[x.name], ends[x.name]))
    data.end = data.groupby('target').end.transform(lambda x: normalize(x, starts[x.name], ends[x.name]))

    hover = HoverTool(tooltips=[('Parent', '@source')])
    
    # seg = hv.Segments(data, ['start', 'target', 'end', 'target'])
    # seg.opts(color='source', line_width=30, width=800, height=800, cmap='tab10',
    #          xlabel='Position', ylabel='Phage', legend_limit=50,
    #          tools=[hover])

    legend_opts = dict(
        # label_text_font_size="12px",
        # label_text_line_height=20, # doesnt seem to be doing anything
        # spacing=0,
        glyph_height=15,
        label_height=15,
     )

    plot_opts = dict(
        width=1000, height=60*data.target.nunique(),
        xlabel='Position', ylabel='Phage',
        legend_position="right",
        legend_limit=50,
        legend_offset=(30, 0),
        legend_opts=legend_opts
    )

    parents = sorted([x for x in data.source.unique() if 'NoCov' not in x])
    parents_nocov = sorted([x for x in data.source.unique() if 'NoCov' in x], key=lambda x: int(x.split('-')[-1]))

    cmap_parents = cc.glasbey_light[:len(parents)]
    cmap_parents_nocov = cc.gray[0:len(cc.gray):(len(cc.gray)//len(parents_nocov))]
    
    seg = hv.Overlay([
        hv.Segments(data[data.source==s], ['start', 'target', 'end', 'target'], label=s)
        .opts(line_width=30, color=c, tools=[hover])
        for s, c in zip(parents + parents_nocov, cmap_parents + cmap_parents_nocov)
    ]).opts(**plot_opts)

    hv.save(seg, '/tmp/cedric/modular_painting_tests/paintings.html')
