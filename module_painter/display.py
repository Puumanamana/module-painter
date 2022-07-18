from pathlib import Path
from math import pi
import numpy as np
import pandas as pd

import numpy as np
from bokeh.models import HoverTool
from bokeh.palettes import Colorblind, Turbo256, linear_palette
import seaborn as sns
import matplotlib.pyplot as plt

import holoviews as hv
from holoviews import dim
import colorcet as cc
hv.extension("bokeh")


def get_cmap(factor):
    categories = set(factor)
    n = len(categories)
    if n <= 8:
        palette = Colorblind[max(3, n)]
    else:
        palette = linear_palette(Turbo256, n)
        palette = np.random.choice(palette, n, replace=False).tolist()

    cmap = dict(zip(categories, palette))

    return cmap

def normalize(x, start, end):
    return 100 * (x-start)/(end-start)

def remove_overlap(df):
    cur_end = df.send
    next_start = np.roll(df.sstart, -1)
    
    overlap = (cur_end - next_start).values
    # Wrapping condition. Rolling doesnt quite work
    overlap[-1] -= df.slen.iloc[0]
    # Cannot use boolean indexing for shifting
    overlap_cur_idx = df.index[overlap > 0]
    overlap_next_idx = np.roll(df.index, -1)[overlap>0]
    overlap_vals = overlap[overlap > 0]
    
    df.loc[overlap_cur_idx, "send"] -= overlap_vals // 2
    df.loc[overlap_next_idx, "sstart"] += overlap_vals // 2

    return df

def prepare_data(graphs, norm=True, clusters=None, disjoint=False, remove_na=True):
    graphs = {graph["sacc"]: graph for graph in graphs}
    graphs = [graphs[sacc] for cluster in clusters for sacc in cluster]
    
    cols = ["sstart", "send", "parent", "sacc", "slen"]

    data = pd.DataFrame([
        vals for g in graphs
        for vals in zip(*[g.vs[col] for col in cols])
    ], columns=cols)

    if remove_na:
        data = data[~data.parent.str.contains("@")]
    if disjoint:
        data = data.groupby("sacc").apply(remove_overlap)
    if norm:
        data["extra"] = data.send - data.slen

        if any(data.extra > 0):
            data.loc[data.extra > 0, "send"] = data.slen
            to_add = data[data.extra > 0].copy()
            to_add["sstart"] = 0
            to_add["send"] = data.extra

            data = pd.concat([data, to_add]).sort_values(by=["sacc", "sstart", "send"]).reset_index(drop=True)
    return data
    
def display_phages_hv(graphs, clusters=None, norm=True, outdir=None, fmt="html"):
    data = prepare_data(graphs, norm=norm, clusters=clusters)

    hover = HoverTool(tooltips=[(x, f"@{x.lower()}") for x in ["parent", "sacc", "sstart", "send"]])
    
    legend_opts = dict(
        glyph_height=15,
        label_height=15,
     )

    # max_cluster_size = max(len(c) for c in clusters)
    plot_opts = dict(
        width=900, #height=max(50*max_cluster_size, 200),
        xlabel="Position", ylabel="Phage",
        gridstyle=dict(ygrid_line_color="gray", xgrid_line_alpha=0, ygrid_line_dash=[4, 4]),
        show_grid=True,
        show_legend=True,
        legend_position="right",
        legend_limit=50,
        legend_offset=(30, 0),
        legend_opts=legend_opts,
    )

    cmap = get_cmap(data.parent)

    subplots = []

    for cluster in clusters:
        # if len(cluster) == 1:
        #     continue
        data_c = data[data.sacc.isin(cluster)].copy()
        data_c["sacc_loc"] = pd.Categorical(data_c.sacc).codes

        parents_all = {p for p in data_c.parent}
       
        subplot_layer = []

        for i, p in enumerate(parents_all):
            data_c_p = data_c[data_c.parent==p].copy()

            # Add an offset to better see breakpoint overlap
            data_c_p.sacc_loc = data_c_p.sacc_loc + 0.05*i
            
            subplot_layer.append(
                hv.Segments(
                    data_c_p, ["sstart", "sacc_loc", "send", "sacc_loc"], label=p
                )
                .opts(line_width=15, color=cmap[p],
                      tools=[hover],
                      height=max(2, len(cluster))*100,
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
        
    seg = hv.Layout(subplots).opts(shared_axes=False).cols(1)

    hv.save(seg, Path(outdir, "paintings.html"))
        
def display_phages_plt(graphs, clusters=None, norm=True, outdir=None, fmt="pdf"):
    data = prepare_data(graphs, norm=norm, clusters=clusters, disjoint=True)
    cluster_map = {p: i for (i, phages) in enumerate(clusters) for p in phages}
    data["cluster"] = [cluster_map[sacc] for sacc in data.sacc]

    colormap = get_cmap(data.parent)

    sns.set(style='ticks', font_scale=1)
    g = sns.FacetGrid(data=data, row="cluster", aspect=3, sharey=False, sharex=False)
    g.map_dataframe(plt_draw_intervals, x1="sstart", x2="send", y="sacc", hue="parent",
                    cmap=colormap, data=data)
    g.add_legend()
    g.set(ylabel="")
    g.savefig(Path(outdir, f"paintings.{fmt}"))
        
def plt_draw_intervals(data=None, x1="sstart", x2="send", y="sacc", hue="parent",
                       cmap=None, lw=15, **kw):
    
    y_pos = {sacc:i for (i, sacc) in enumerate(data.sacc.unique())}

    grouped = data.groupby(hue)[[x1,x2,y]].agg(list).T.to_dict()

    vbar_h = 0.09 * (data[y].nunique() - 1/data[y].nunique())
    for (parent, itv) in grouped.items():
        color = cmap[parent]
        y_order = [y_pos[sacc] for sacc in itv[y]]
        (vbar_min, vbar_max) = zip(*[(pos-vbar_h, pos+vbar_h) for pos in y_order])
        
        plt.hlines(y_order, itv[x1], itv[x2], color=color, lw=lw, label=parent, alpha=0.8)
        plt.vlines(itv[x1], vbar_min, vbar_max, color, lw=2)
        plt.vlines(itv[x2], vbar_min, vbar_max, color, lw=2)
    
    plt.yticks(list(y_pos.values()), list(y_pos.keys()))
    # plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=14)
