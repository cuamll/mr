#!/usr/bin/env python3
import numpy as np
import warnings

# colours
purp = "#832591"
blu = "#278BE8"
rd = "#99001A"
grn = "#145E00"
yl = "#F2D33C"
orng = "#E88427"
levels = 128

def s_p(x, y, G_x, G_y):
    '''
        This is basically the function f from Steve's correlation notes.
        We can think of it as a kind of "projection" from the tensor 
        (\chi / S)^{\alpha \beta}, the theoretical object, onto the
        perpendicular component we get from neutron scattering.
        Without superposing this function we don't see the pinch point;
        it's the combination of the \chi/S tensor and f which gives us
        the pinch point.

        NB: the reason for catch warnings is that the zone centre throws
        a warning at one index of each array, where x = G_x and y = G_y.
        But I can't find an easy way of using meshgrid below and
        only picking out the zone centre to skip; it complains about
        truthiness being ambiguous for arrays. Hence, ignore the warning.
        Bit of a hack, but whatever.

    '''
    # the zone centre throws a warning because x = G_x and y = G_y,
    # but I can't find an easy way of using meshgrid below and
    # only picking out the zone centre to skip; it complains about
    # truthiness being ambiguous for arrays. hence, ignore the warning.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = (((G_x - x)*x + (G_y - y)*y)**2)/((x**2 + y**2)*((x - G_x)**2 + (y - G_y)**2))

    return res

def do_plots(no_plots, plot_titles, output_file, dots, X, Y, Z):
    ''' 
        I have a few different places where I want to make side-by-side
        contour plots of some data. This is just a wrapper to do that for me.
        NB: need to set plot rc somewhere in the file first?
    '''
    import matplotlib.pyplot as plt
    import matplotlib.ticker as tck
    import colormaps as cm
    fig, axes = plt.subplots(ncols=int(no_plots), nrows=1, figsize=(5 * int(no_plots),4))
    # chonk
    for chonk in range(int(no_plots)):
        axes[chonk].xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
        axes[chonk].xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
        axes[chonk].yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
        axes[chonk].yaxis.set_major_locator(tck.MultipleLocator(base=1.0))
        ax = axes[chonk]

        ax.set_title(plot_titles[chonk])
        cs = ax.contourf(X[chonk], Y[chonk], Z[chonk], levels, cmap=cm.inferno)
        fig.colorbar(cs, ax=ax)

    fig.tight_layout()
    fig.savefig(output_file, format='eps', dpi=dots)
    plt.close(fig)

