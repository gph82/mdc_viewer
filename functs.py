import io
from contextlib import redirect_stdout
import numpy as _np
import matplotlib.pyplot as _plt
from mdciao.utils.residue_and_atom import shorten_AA as _shortenAA
from mdciao.plots import compare_groups_of_contacts as _compare_groups_of_contacts

def indv_plot(argmap, iCG, im_wdg):
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
        # if imgs is not None:
        #    imgs.append(widgets.Image(value=b.getvalue(),
        #                          width=width_px,
        #                          height=height_px))
        # else:
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
                                                   title='%s'% anchor_str
                                            )[0]

    def run(res_idxs,argmap, individual_controls):

        _plt.close("all")
        plot_widgets = []
        shown_CGs=[]
        for rr in res_idxs["res"]:
            rr = int(rr)
            Errors.clear_output()
            for key, iarch in indict.items():
                if rr in iarch.keys():
                    try:
                        CGs[rr][key]
                    except KeyError:
                        CGs[rr][key] = mdcv_io.CGdict2CG(iarch[rr], top=top)
                else:
                    with Errors:
                        Errors.clear_output(wait="True")
                        print("%s was found in the topology but not in the indict"%top.residue(rr))
            iCG: mdciao.contacts.ContactGroup
            if len(CGs[rr])>0:
                shown_CGs.append(CGs[rr])
                b = io.StringIO()
                with redirect_stdout(b):
                    ifig = functs.lambda_compare_neighborhoods(CGs[rr],argmap)
                b.close()
                ifig: _plt.Figure
                width_px, height_px = [_np.round(ii / 1.25) for ii in ifig.get_size_inches() * ifig.get_dpi()]
                b = io.BytesIO()
                ifig.savefig(b, format="png", bbox_inches="tight")
                _plt.close()
                img_wdg = widgets.Image(value=b.getvalue(),
                                          width=width_px,
                                          height=height_px)


                first_run["res"]=True

            if individual_controls.value:
                fs = freq_slider()
                fs.description = "cutoff"

                indv_options_box, indv_argmap = general_options(horizontal=False)
                indv_argmap["frequency"]=fs

                indv_options_box = _VBox([fs]+list(indv_options_box.children))
                [setattr(child,"layout",{"width":"100%"}) for child in indv_options_box.children]

                indv_options_box.layout = {"width":"31%",
                              #'border': '1px solid black'
                                           }

                [setattr(val, "_res_idx", rr) for key, val in indv_argmap.items() if key!="kwargs"]
                [setattr(val, "_res_idx", rr) for key, val in indv_argmap["kwargs"].items()]
                per_fig_controls[rr] = indv_argmap
                image_widgets[rr] = img_wdg
                plot_widgets.append(_HBox([_HBox([img_wdg], layout={"width": "68%"}),
                                                    indv_options_box],
                                                   layout={"width": "99%"}))
            else:
                plot_widgets.append(img_wdg)


            img_box.children = plot_widgets

_center_colors = {'A': '#7d4dfbff',
                  'B': '#ff58ffff',
                  'C': '#8cffffff'}