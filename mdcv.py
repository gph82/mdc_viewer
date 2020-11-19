import ipywidgets as widgets
import numpy as _np
import mdciao
import ipywidgets
from mdciao.fragments.fragments import _allowed_fragment_methods
from IPython.display import display
from contextlib import redirect_stdout
import pandas
from matplotlib import pyplot as _plt
import io
import mdcv_io
from collections import defaultdict as _defdict
import mdtraj as _md

def show():
    top = _md.load("data/top.pdb.gz").top
    data = mdcv_io.load_data(verbose=0, decompress_here=False)
    screen3(data,
             top,
             initial_value="R131,GDP"
            );

def fragmentation_Hbox(top):
    r"""
    Returns interactive viewer for the fragmentation heuristics

    Parameters
    ----------
    top

    Returns
    -------

    """
    Heuristics_Dropdown = widgets.Dropdown(description="Choose fragmentation heuristic",
                                           options=_allowed_fragment_methods, value="lig_resSeq+",
                                           layout={"width": "max-content"},
                                           style={"description_width": "initial"})


    Evaluated_Output= widgets.Output(layout={'border': '1px solid black'})


    fragments={}
    def fragment_run():
        with Evaluated_Output:
            Evaluated_Output.clear_output(wait="True")
            fragments["res"] = mdciao.fragments.get_fragments(top,
                                                              method=Heuristics_Dropdown.value)
    Heuristics_Dropdown.observe(lambda _ : fragment_run(), names="value")

    fragment_run()

    return ipywidgets.VBox([
        Heuristics_Dropdown,
        Evaluated_Output
    ]), fragments


def residue_selection(top, fragments=None, initial_value="R131,GDP*,L393-L394"):
    r"""
    Returns interactive viewer for the residue selection

    Parameters
    ----------
    top
    fragments
    initial_value

    Returns
    -------

    """

    residue_input = ipywidgets.Text(description="Choose residue(s):",
                                    style={"description_width":"auto",
                                           },
                                    value=initial_value,
                                    continuous_update=False,
                                    layout={"width":"500px"}
                                    )
    residue_output = ipywidgets.Textarea(placeholder="'unpacked' list of target residues will appear here",
                                         layout={"width":"100%"})
    evaluate_residues = widgets.Button(description="preview",
                                       layout={"width":"auto"})

    residue_input.observe(lambda _ : eval_residue_selection(), names="value")
    if fragments is None:
        fragments = [_np.arange(top.n_residues)]

    residue_idxs = {}
    def eval_residue_selection():
        if len(residue_input.value)>0:
            try:
                b = io.StringIO()
                with redirect_stdout(b):
                    residue_idxs["res"] = mdciao.utils.residue_and_atom.rangeexpand_residues2residxs(residue_input.value.replace(" ","").strip(","),
                                                                                      fragments,
                                                                                      top )
                b.close()
                istr = ", ".join([str(top.residue(rr)) for rr in residue_idxs["res"]])
            except ValueError as e:
                residue_idxs["res"] = []
                istr = str(e)

            residue_output.value = istr


    evaluate_residues.on_click(lambda _ : eval_residue_selection())
    return widgets.VBox([widgets.HBox([residue_input, evaluate_residues], layout={"width": "100%"}),
                         residue_output]),residue_idxs

def args_box():
    desc = ipywidgets.Button(description="AA-Label-Options:",
                             layout={"width": "100%"},
                             )

    """
    tgl_short_AAs = ipywidgets.ToggleButton(value=True,
                                              description="short AA names",
                                              icon="check",
                                              layout={"width": "50%"})
    """
    tgl_hide_anchor = ipywidgets.ToggleButton(value=True,
                                     description="use anchor AA",
                                     icon="check",
                                     layout={"width": "50%"})

    tgl_consensus_labs = ipywidgets.ToggleButton(value=True,
                                                    description="consensus labels",
                                                    icon="check",
                                                    layout={"width": "50%"})
    tgl_color_hint = ipywidgets.ToggleButton(value=False,
                                         description="color hint",
                                         layout={"width": "50%"})
    fontsize = ipywidgets.IntText(value=16, description="fontsize",
                                  layout={"width": "50%"})

    AA_labels = ipywidgets.VBox([desc,
                                 ipywidgets.HBox([
                                     #tgl_short_AAs,
                                                  tgl_hide_anchor,
                                                  tgl_consensus_labs]),
                                 ipywidgets.HBox([tgl_color_hint, fontsize])],
                                layout={"width": "35%"})

    tgl_freqs_above = ipywidgets.ToggleButton(value=True,
                                              description="Don't show freqs >=",
                                              style={"description_width": "max-content"},
                                              layout={"width": "70%"},
                                              icon="check"
                                              )
    thrs_above = ipywidgets.FloatText(1.0,
                         layout={"width": "30%"})
    tgl_freqs_below = ipywidgets.ToggleButton(value=True,
                                              description="Don't show freqs <=",
                                              layout={"width": "70%"},
                                              icon="check"
                                              )
    thrs_below = ipywidgets.FloatText(0.2,
                         style={"description_width": "max-content"},
                         layout={"width": "30%"})

    colors = color_pickers()


    identity = ipywidgets.VBox([
        ipywidgets.Button(description="Bar-options:",
                          layout={"width": "100%"}),
        ipywidgets.HBox([tgl_freqs_above, thrs_above]),
        ipywidgets.HBox([tgl_freqs_below, thrs_below]),
        # ipywidgets.HBox(colors)
        ],
        layout={"width": "30%"})


    for tgl in [tgl_freqs_above,tgl_freqs_below,
                #tgl_short_AAs,
                tgl_hide_anchor, tgl_consensus_labs, tgl_color_hint]:
        tgl.observe(lambda traits : change_icon(traits),names="value")

    argmap = {"kwargs":{"remove_identities": tgl_freqs_above,
                        "identity_cutoff" :  thrs_above,
                        "assign_w_color" : tgl_color_hint,
                        "lower_cutoff_val" : thrs_below,
                        "fontsize" : fontsize,
                        "anchor":tgl_hide_anchor,
                        },
              "tgl_freqs_below" : tgl_freqs_below,
              "tgl_consensus" : tgl_consensus_labs,
              #"tgl_short_AAs" : tgl_short_AAs,
              }

    return ipywidgets.HBox([
        identity,
        AA_labels
    ],
        layout={"width": "100%"}
    ), argmap

def cutoff_selection_HBox():
    clear_button = widgets.Button(description="clear neighborhoods", button_style="danger",
                                  layout={"width": "auto"})

    FreqSlider = widgets.FloatSlider(
        value=3.5,
        min=0,
        max=5,
        step=0.1,
        description='Neighbor Distance cutoff',
        style={"description_width": "auto"},
        readout_format='.1f',
        continuous_update=False,
    )

    run_button = widgets.Button(description="show neighborhoods", button_style="success",
                                layout={"width": "auto"}
                                )

    return widgets.HBox([FreqSlider, run_button, clear_button])

def screen2(archive_dict,top, **kwargs):


    residue_selection_box, res_idxs = residue_selection(top,**kwargs)





    img_box = widgets.VBox([])

    CGs = {}

    options_wdg, argmap = args_box()
    accordion = ipywidgets.Accordion(
        [options_wdg],
        value="options",
        layout={"width": "100%"})
    accordion.set_title(0, "Advanced Options")

    out_VBox = widgets.VBox([
        residue_selection_box,
        accordion,
        img_box,
        Errors,
    ])

    display(out_VBox)

    _plt.ioff()
    def run(lambda_of_CG, res_idxs, argmap):
        _plt.close("all")
        imgs = []
        for rr in res_idxs["res"]:
            rr = int(rr)
            if rr in archive_dict.keys():
                Errors.clear_output()
                if rr not in CGs.keys():
                    CGs[rr] = mdcv_io.CGdict2CG(archive_dict[rr], top=top)
                iCG: mdciao.contacts.ContactGroup
                ifig = lambda_of_CG(CGs[rr],argmap)
                ifig: _plt.Figure
                width_px, height_px = [_np.round(ii / 1.5) for ii in ifig.get_size_inches() * ifig.get_dpi()]
                b = io.BytesIO()
                ifig.savefig(b, format="png", bbox_inches="tight")
                _plt.close()
                imgs.append(widgets.Image(value=b.getvalue(),
                                          width=width_px,
                                          height=height_px))
            else:
                with Errors:
                    Errors.clear_output(wait="True")
                    print("%s was found in the topology but not in the archive"%top.residue(rr))

        img_box.children = imgs

    lambda_plot_neighborhood_freqs = lambda CG : CG.plot_neighborhood_freqs(FreqSlider.value,4,display_sort=True).figure
    lambda_compare_neighborhoods = lambda CG, argmap : mdciao.plots.compare_groups_of_contacts([CG],
                                                                                               ctc_cutoff_Ang=FreqSlider.value,
                                                                                               anchor=mdciao.utils.residue_and_atom.shorten_AA(str(CG._contacts[0].residues.anchor_residue), keep_index=True),
                                                                                               remove_identities=argmap["kwargs"]["remove_identities"].value,
                                                                                               identity_cutoff=argmap["kwargs"]["identity_cutoff"].value,
                                                                                               assign_w_color=argmap["kwargs"]["assign_w_color"].value,
                                                                                               lower_cutoff_val=argmap["kwargs"]["lower_cutoff_val"].value*int(argmap["tgl_freqs_below"].value),
                                                                                               fontsize=argmap["kwargs"]["fontsize"].value)[0]
    #run_lambda = lambda out_dict : run(lambda_plot_neighborhood_freqs, res_idxs)
    run_lambda = lambda out_dict : run(lambda_compare_neighborhoods, res_idxs,argmap)

    FreqSlider.observe(lambda _: run_lambda(None), names="value")
    run_button.on_click(lambda _ : residue_selection_box.children[0].children[-1].click())
    run_button.on_click(run_lambda)
    clear_button.on_click(lambda _ :   [img.close() for img in img_box.children])
    return out_VBox


_center_colors = {'A': '#7d4dfbff',
                  'B': '#ff58ffff',
                  'C': '#8cffffff'}

def color_pickers(centers_colors=None):
    if centers_colors is None:
        centers_colors = _center_colors
    cps = []
    width = "%s%%"%_np.floor(100/len(centers_colors)*.90)
    for key, value in centers_colors.items():
        iwdg = ipywidgets.ColorPicker(value = value[:-2],
                                      description=key,
                                      layout={"width": width},
                                      style={"description_width": "initial"},
                                      #ayout={"width":width}
        )
        cps.append(iwdg)

    return cps

def screen3(indict,top, **kwargs):

    residue_selection_box, res_idxs = residue_selection(top,**kwargs)

    residue_selection_acc = ipywidgets.Accordion([residue_selection_box],
                                                 layout={"width": "100%"})
    residue_selection_acc.set_title(0, "Residue Selection")

    cutoff_selection = cutoff_selection_HBox()
    FreqSlider, run_button, clear_button = cutoff_selection.children

    Errors =  widgets.Output(layout={'border': '1px solid black'})

    img_box = widgets.VBox([])

    CGs = _defdict(dict)

    output_acc = widgets.Accordion([widgets.VBox([img_box,Errors])],
                                   layout={"width":"100%"})
    output_acc.set_title(0, "Output")

    options_wdg, argmap = args_box()
    accordion = ipywidgets.Accordion(
        [options_wdg],
        layout={"width": "100%"},
        selected_index=None)

    accordion.set_title(0, "Plot Options")

    out_VBox = widgets.VBox([
        residue_selection_acc,
        accordion,
        cutoff_selection,
       output_acc,
    ])

    display(out_VBox)

    _plt.ioff()
    first_run={"res":False}
    def run(res_idxs,argmap):
        _plt.close("all")
        imgs = []
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
                b = io.StringIO()
                with redirect_stdout(b):
                    ifig = lambda_compare_neighborhoods(CGs[rr],argmap)
                b.close()
                ifig: _plt.Figure
                width_px, height_px = [_np.round(ii / 1.25) for ii in ifig.get_size_inches() * ifig.get_dpi()]
                b = io.BytesIO()
                ifig.savefig(b, format="png", bbox_inches="tight")
                _plt.close()
                imgs.append(widgets.Image(value=b.getvalue(),
                                          width=width_px,
                                          height=height_px))
                first_run["res"]=True

        img_box.children = imgs


    def lambda_compare_neighborhoods(CGs,argmap):
        if argmap["kwargs"]["anchor"].value:
            anchor=mdciao.utils.residue_and_atom.shorten_AA(str(list(CGs.values())[0]._contacts[0].residues.anchor_residue),
                                                            substitute_fail="long",
                                                            keep_index=True)
        else:
            anchor=None
        if argmap["tgl_consensus"].value:
            defrag=" "
        else:
            defrag="@"
        return mdciao.plots.compare_groups_of_contacts(CGs,
                                                ctc_cutoff_Ang=FreqSlider.value,
                                                anchor=anchor,
                                                remove_identities=argmap["kwargs"]["remove_identities"].value,
                                                identity_cutoff=argmap["kwargs"]["identity_cutoff"].value,
                                                assign_w_color=argmap["kwargs"]["assign_w_color"].value,
                                                lower_cutoff_val=argmap["kwargs"]["lower_cutoff_val"].value * int(
                                                    argmap["tgl_freqs_below"].value),
                                                fontsize=argmap["kwargs"]["fontsize"].value,
                                                defrag=defrag,
                                                       colors=_center_colors,
                                                )[0]

    run_lambda = lambda out_dict : run(res_idxs, argmap)

    FreqSlider.observe(lambda _: run_lambda(None), names="value")
    for wdg in list(argmap["kwargs"].values())+[val for key, val in argmap.items() if key!="kwargs"]:
        wdg.observe( lambda first_run :  run_lambda(None), names="value")

    run_button.on_click(lambda _ : residue_selection_box.children[0].children[-1].click())
    run_button.on_click(run_lambda)
    clear_button.on_click(lambda _ :   [img.close() for img in img_box.children])
    return out_VBox

def change_icon(traitlets):
    #assert traitlets["owner"].value and traitlets["owner"].icon=="check",traitlets
    if not traitlets["owner"].value:
        traitlets["owner"].icon=""
    else:
        traitlets["owner"].icon="check"