import io
from contextlib import redirect_stdout
import numpy as _np
import matplotlib.pyplot as _plt
from mdciao.utils.residue_and_atom import shorten_AA as _shortenAA
from mdciao.plots import compare_groups_of_contacts as _compare_groups_of_contacts
import ipywidgets
def indv_plot(argmap, iCG, im_wdg=None, capture=True, site=False):
    r"""
    Generate or update an individual plot from an argmap and a dictionary of ContactGroups

    The plot is transformed to bitmap internally

    If :obj:`img_wdg` is passed, then its :obj:`value` is updated with the bitmap.

    Else, a new :obj:`ipywidgets.Image` is instantiated, plotted-onto, and returned.

    Wraps around :obj:`lambda_compare_neighborhoods`.

    Parameters
    ----------
    argmap : dict
        Will be passed to :obj:`argmap2kwargs_compare_groups_of_contacts`
    iCG : dict
        Dictionary of :obj:`~mdciao.contacts.ContactGroup` objects
    im_wdg : :obj:`ipywidgets.Image` or None
        Pass an existing image widget here to be updated,
        else one will be instantiated and returned
    capture : bool, default is True
        Capture output
    site : bool, default is False
        Plot as site and not neighborhood

    Returns
    -------
    img_wdg : :obj:`ipywidgets.Image`

    """

    if len(iCG) > 0:
        b = io.StringIO()
        if capture:
            with redirect_stdout(b):
                ifig = lambda_compare_neighborhoods(iCG, argmap, site=site)
        else:
            ifig = lambda_compare_neighborhoods(iCG, argmap, site=site)
        b.close()
        ifig: _plt.Figure
        width_px, height_px = [_np.round(ii / 1.25) for ii in ifig.get_size_inches() * ifig.get_dpi()]
        b = io.BytesIO()
        ifig.savefig(b, format="png", bbox_inches="tight")
        _plt.close()
        if im_wdg is None:
            return ipywidgets.Image(value=b.getvalue(),
                                      width=width_px,
                                      height=height_px)
        else:
            im_wdg.value = b.getvalue()
            im_wdg.width = width_px,
            im_wdg.height = height_px

def lambda_compare_neighborhoods(CG_dict, argmap, site=False):
    r"""
    Wrap around :obj:`~mdciao.plots.compare_groups_of_contacts` using :obj:`argmap`:

    Parameters
    ----------
    CG_dict : dict
        The ContactGroups
    argmap : dict
        The control dictionary containing needed arguments.
        This object is a WIP. It's currently a dictionary
        so it can be updated by whatever method downstream
        using new values of the widgets. Will be passed to
        :obj:`argmap2kwargs_compare_groups_of_contacts`
    site : bool, default is None

    Returns
    -------

    """
    """
    anchor_str = _shortenAA(str(list(CG_dict.values())[0]._contacts[0].residues.anchor_residue),
                            substitute_fail="long",
                            keep_index=True)
    """

    kwargs = argmap2kwargs_compare_groups_of_contacts(argmap)
    anchor = None
    if site:
        title_str = str(list(CG_dict.values())[0].name)
    else:
        title_str = list(CG_dict.values())[0].anchor_res_and_fragment_str_short
        if argmap["kwargs"]["anchor"].value:
            anchor=title_str

        #if anchor is not None:
        #    anchor_str = list(CG_dict.values())[0]._contacts[0]._attribute_neighborhood_names.anchor_residue_short
    return _compare_groups_of_contacts(CG_dict,
                                       anchor = anchor,
                                       title='%s ' % str(title_str.replace("@", "^")),
                                       **kwargs,
                                       )[0]

def argmap2kwargs_compare_groups_of_contacts(argmap=None):
    r"""
    Map an :obj:`argmap`-dict to a kwargs dict

    The argmap is a dictionary containing :obj:`ipywidgets`, with
    keys that make sense in :obj:`mdcviewer` but not in :obj:`mdciao`.
    These :obj:`ipywidgets` have a :obj:`value` property,
    which is turned/mapped into a keyword argument
    for calling :obj:`~mdciao.plots.compare_groups_of_contacts`.

    Parameters
    ----------
    argmap : dict

    Returns
    -------
    kwargs : dict
    """


    kwargs = {
    "ctc_cutoff_Ang" : float(argmap["frequency"].value.split(" ")[0]),
    "remove_identities" : argmap["kwargs"]["remove_identities"].value,
    "identity_cutoff" : argmap["kwargs"]["identity_cutoff"].value,
    "assign_w_color" : argmap["kwargs"]["assign_w_color"].value,
    "lower_cutoff_val" : argmap["kwargs"]["lower_cutoff_val"].value * int(argmap["tgl_freqs_below"].value),
    "fontsize" : argmap["kwargs"]["fontsize"].value,
    "colors" : _center_colors,
    }

    if argmap["tgl_consensus"].value:
        kwargs.update({"defrag":None})
    else:
        kwargs.update({"defrag": "@"})
    return kwargs


_center_colors = {'A': '#7d4dfbff',
                  'B': '#ff58ffff',
                  'C': '#8cffffff'}