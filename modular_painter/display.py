from math import pi
import pandas as pd

from bokeh.io import output_file, save
from bokeh.plotting import figure
from bokeh.transform import cumsum
from bokeh.palettes import Set1, Set3

TOOLS = "pan,wheel_zoom,box_zoom,reset,save,hover"
TOOLTIPS = "@source - identity=@identity{0.00%}"
PALETTE = Set1[9] + Set3[12]

def access(L, i):
    if isinstance(L, list):
        if i < len(L):
            return L[i]
        return None
    return L

def donut_display(cov, output, title="", circular=True, **plot_kw):

    # Set up different colors per source bacteriophage
    sources = cov.data.source.unique()
    cov.data.identity /= cov.data.tend - cov.data.tstart
    data = cov.data.groupby(['tstart', 'tend'], as_index=False).agg(list)

    color_map = pd.Series(dict(
        [(s, PALETTE[i % len(PALETTE)]) for (i, s) in enumerate(sources)]
        + [('NC', '#FFFFFF')]
        + [('None', '#FFFFFF')]        
    ))

    # Number of levels for painting
    # (when multiple bacteriophages align on the same fragment)
    max_depth = data.source.apply(len).max()

    # Additional columns for plotting
    data['mid'] = (data.tstart + data.tend) / 2
    data['fs'] = data.tend - data.tstart
    data['angle'] = data.fs / data.fs.sum() * 2*pi
    data['color'] = data.source.apply(lambda x: color_map.loc[x].tolist())

    # Adjust figure limits when plotting linear fragments
    if not circular:
        layer_height = 2

        plot_kw = {
            'plot_width': 1000, 'plot_height': 400,
            'y_range': (-layer_height*10, 5*max_depth*layer_height)
        }

    p = figure(title=title, tools=TOOLS, tooltips=TOOLTIPS, **plot_kw)

    # Plot one layer at a time
    for depth in range(max_depth):
        # data_i is the i^th layer
        data_i = data.applymap(lambda x: access(x, depth))
        data_i.fillna('None', inplace=True)

        # Hide the uncovered sections
        data_i['alpha'] = (data_i.source != 'None').astype(int)

        if circular:
            outer_r = (max_depth - depth)*0.2
            inner_r = (max_depth - depth - 1)*0.2

            p.annular_wedge(source=data_i, x=0, y=1, inner_radius=inner_r, outer_radius=outer_r,
                            start_angle=cumsum('angle', include_zero=True), end_angle=cumsum('angle'),
                            line_color="white", fill_color='color', legend_group='source')
        else:
            p.rect(x='mid', y=(max_depth - depth)*3, source=data_i,
                   width='fs', height=layer_height,
                   fill_color='color', line_alpha=0, line_width=0, fill_alpha='alpha',
                   legend_group='source')

    p.axis.axis_label=None
    p.axis.visible=False
    p.grid.grid_line_color = None

    lkup = {x.label['value']: x for x in p.legend.items
            if x.label['value'] not in ['None', 'NC']}
    p.legend.items.clear()

    for itm in lkup.values():
        p.legend.items.append(itm)
    
    output_file(output)
    save(p)
