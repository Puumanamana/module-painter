from math import pi
import numpy as np
import pandas as pd

from bokeh.io import output_file, save
from bokeh.plotting import figure
from bokeh.transform import cumsum
from bokeh.palettes import Turbo256, linear_palette
from bokeh.models import Legend

TOOLS = "box_zoom,reset,hover"
TOOLTIPS = [("source", "@source"),
            ("identity", "@identity{0.00%}"),
            ("interval", "@start-@end"),
            ("position", "$x{0,0}")]

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

def display_genomes(df, genomes):
    return

def display_linear_genome(cov, output, title="", **plot_kw):
    data = cov.to_pandas()

    end_after_mod = data.end > cov.mod

    if sum(end_after_mod) > 0:
        data_wrapped = data.loc[end_after_mod].copy()
        data_wrapped.start = 0
        data_wrapped.end = data_wrapped.end % cov.mod

        data.loc[end_after_mod, 'end'] = cov.mod

        data = pd.concat([data_wrapped, data])

    # Set up different colors per source bacteriophage
    sources = data.source.unique()
    data.identity /= data.end - data.start
    data = data.groupby(['start', 'end'], as_index=False).agg(list)

    palette = pd.Series(dict(zip(sources, get_palette(len(sources)))))

    # Number of levels for painting
    # (when multiple bacteriophages align on the same fragment)
    max_depth = data.source.apply(len).max()

    # Additional columns for plotting
    data['mid'] = (data.start + data.end) / 2
    data['fs'] = data.end - data.start
    data['color'] = data.source.apply(lambda x: palette.loc[x].tolist())

    # Adjust figure limits when plotting linear fragments
    layer_height = 2

    plot_kw = {
        'min_border': 200,
        'plot_width': 1200, 'plot_height': 800,
        'y_range': (-layer_height*10, 5*max_depth*layer_height)
    }

    p = figure(title=title, tools=TOOLS, tooltips=TOOLTIPS, **plot_kw)

    # Plot one layer at a time
    for depth in range(max_depth):
        # data_i is the i^th layer
        data_i = data.applymap(lambda x: access(x, depth))
        data_i.fillna('None', inplace=True)
        data_i.loc[data_i.source == 'NoCov', 'alpha'] = 0

        # Hide the uncovered sections
        data_i['alpha'] = (data_i.source != 'None').astype(int) * 0.6

        p.rect(x='mid', y=(max_depth - depth)*3, source=data_i,
               width='fs', height=layer_height,
               fill_color='color', fill_alpha='alpha',
               line_alpha='alpha', line_width=1, line_color='black',
               legend_group='source')

    leg_items = {}

    for item in p.legend.items:
        label = item.label['value']
        if label not in leg_items and label not in {'NoCov', 'None'}:
            leg_items[label] = item
            
    p.legend.items.clear()

    legend = Legend(items=list(leg_items.values()), glyph_height=20, location=(20, 0))
    p.add_layout(legend, 'right')

    p.legend.label_text_font_size = '10pt'
    p.ygrid.visible = False
    p.yaxis.visible = False
    output_file(output)
    save(p)
    

def display_circular_genome(cov, output, title="", **plot_kw):
    data = cov.to_pandas()

    end_after_mod = data.end > cov.mod

    if sum(end_after_mod) > 0:
        data_wrapped = data.loc[end_after_mod].copy()
        data_wrapped.start = 0
        data_wrapped.end = data_wrapped.end % cov.mod

        data.loc[end_after_mod, 'end'] = cov.mod

        data = pd.concat([data_wrapped, data])

    # Set up different colors per source bacteriophage
    sources = data.source.unique()
    data.identity /= data.end - data.start
    data = data.groupby(['start', 'end'], as_index=False).agg(list)

    palette = pd.Series(dict(zip(sources, get_palette(len(sources)))))

    # Number of levels for painting
    # (when multiple bacteriophages align on the same fragment)
    max_depth = data.source.apply(len).max()

    # Additional columns for plotting
    data['mid'] = (data.start + data.end) / 2
    data['fs'] = data.end - data.start
    data['angle'] = data.fs / data.fs.sum() * 2*pi
    data['color'] = data.source.apply(lambda x: palette.loc[x].tolist())

    p = figure(title=title, tools=TOOLS, tooltips=TOOLTIPS, **plot_kw)

    # Plot one layer at a time
    for depth in range(max_depth):
        # data_i is the i^th layer
        data_i = data.applymap(lambda x: access(x, depth))
        data_i.fillna('None', inplace=True)
        data_i.loc[data_i.source == 'NoCov', 'alpha'] = 0

        # Hide the uncovered sections
        data_i['alpha'] = (data_i.source != 'None').astype(int) * 0.6

        outer_r = (max_depth - depth)*0.2
        inner_r = (max_depth - depth - 1)*0.2

        p.annular_wedge(source=data_i, x=0, y=1, inner_radius=inner_r, outer_radius=outer_r,
                        start_angle=cumsum('angle', include_zero=True), end_angle=cumsum('angle'),
                        line_color="white", fill_color='color', legend_group='source')

    leg_items = {}

    for item in p.legend.items:
        label = item.label['value']
        if label not in leg_items and label not in {'NoCov', 'None'}:
            leg_items[label] = item
            
    p.legend.items.clear()

    legend = Legend(items=list(leg_items.values()), glyph_height=20, location=(20, 0))
    p.add_layout(legend, 'right')

    p.legend.label_text_font_size = '10pt'
    p.ygrid.visible = False
    p.yaxis.visible = False
    output_file(output)
    save(p)
