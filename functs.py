import io
from contextlib import redirect_stdout
import numpy as _np
import matplotlib.pyplot as _plt
from mdciao.utils.residue_and_atom import shorten_AA as _shortenAA
from mdciao.plots import compare_groups_of_contacts as _compare_groups_of_contacts
import ipywidgets
def indv_plot(argmap, iCG, im_wdg=None):
    if len(iCG) > 0:
        b = io.StringIO()
        with redirect_stdout(b):
            ifig = lambda_compare_neighborhoods(iCG, argmap)
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

def lambda_compare_neighborhoods(CGs,argmap):
    anchor_str = _shortenAA(str(list(CGs.values())[0]._contacts[0].residues.anchor_residue),
                            substitute_fail="long",
                            keep_index=True)
    if argmap["kwargs"]["anchor"].value:
        anchor=anchor_str
    else:
        anchor=None
    if argmap["tgl_consensus"].value:
        defrag=" "
    else:
        defrag="@"
    return _compare_groups_of_contacts(CGs,
                                       ctc_cutoff_Ang=argmap["frequency"].value,
                                       anchor=anchor,
                                       remove_identities=argmap["kwargs"]["remove_identities"].value,
                                       identity_cutoff=argmap["kwargs"]["identity_cutoff"].value,
                                       assign_w_color=argmap["kwargs"]["assign_w_color"].value,
                                       lower_cutoff_val=argmap["kwargs"]["lower_cutoff_val"].value * int(
                                           argmap["tgl_freqs_below"].value),
                                       fontsize=argmap["kwargs"]["fontsize"].value,
                                       defrag=defrag,
                                       colors=_center_colors,
                                       title='%s' % anchor_str
                                       )[0]


_center_colors = {'A': '#7d4dfbff',
                  'B': '#ff58ffff',
                  'C': '#8cffffff'}