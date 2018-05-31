#!/usr/bin/env python3
import numpy as np
import warnings

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
