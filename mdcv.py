import ipywidgets as widgets
import numpy as _np
import mdciao
import ipywidgets
from mdciao.fragments.fragments import _allowed_fragment_methods
from mdciao.utils.residue_and_atom import rangeexpand_residues2residxs as _rangeexpand_residues2residxs
from ipywidgets import HBox as _HBox, VBox as _VBox
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
    data = mdcv_io.load_data(verbose=0,
                             #decompress_here=False
                             )
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

    return _VBox([
        Heuristics_Dropdown,
        Evaluated_Output
    ]), fragments

def figure_options():

    title = ipywidgets.Button(description="Figure Options",
                              layout={"width":"100%"})
    zoom = ipywidgets.FloatSlider(
        value=100,
        min=50,
        max=150,
        step=5,
        description='Zoom',
        disabled=False,
        continuous_update=False,
        #orientation='vertical',
        readout=True,
        readout_format='.1f',
    )

    return _VBox([title,zoom])

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

    residue_input = ipywidgets.Text(description="Input manually",
                                    style={"description_width":"auto",
                                           },
                                    value=initial_value,
                                    continuous_update=False,
                                    layout={"width":"500px"}
                                    )
    residue_output = ipywidgets.Textarea(placeholder="'unpacked' list of target residues will appear here",
                                         layout={"width":"100%"})

    residue_list = ipywidgets.Select(
        options=[str(rr) for rr in top.residues],
        rows=1,
        description='or pick from list:',
        # disabled=False,
        layout={"width": "max-content", },
        style={"description_width": "auto",
               },
    )

    pick_list = ipywidgets.Button(description="add",
                                  layout={"width":"auto"})

    def add_list_residue():
        res = [ii for ii in residue_input.value.split(",") if len(ii)>0]
        print(res)
        if residue_list.value not in res:
            residue_input.value = ",".join(res+[residue_list.value])


    pick_list.on_click(lambda _ : add_list_residue())


    evaluate_residues = ipywidgets.Button(description="preview",
                                       layout={"width":"auto"},
                                       button_style="info")

    clear = ipywidgets.Button(description="clear",
                              layout={"width":"auto"},
                              button_style="danger")
    clear.on_click(lambda _ : setattr(residue_input,"value",""))
    clear.on_click(lambda _ : setattr(residue_output,"value",""))


    residue_input.observe(lambda _ : eval_residue_selection(), names="value")
    if fragments is None:
        fragments = [_np.arange(top.n_residues)]

    residue_idxs = {}
    def eval_residue_selection():
        if len(residue_input.value)>0:
            try:
                b = io.StringIO()
                with redirect_stdout(b):
                    residue_idxs["res"] = _rangeexpand_residues2residxs(residue_input.value.replace(" ","").strip(","),
                                                                                      fragments,
                                                                                      top )
                b.close()
                istr = ", ".join([str(top.residue(rr)) for rr in residue_idxs["res"]])
            except ValueError as e:
                residue_idxs["res"] = []
                istr = str(e)

            residue_output.value = istr


    evaluate_residues.on_click(lambda _ : eval_residue_selection())
    box = ipywidgets.VBox([ipywidgets.HBox([residue_input, residue_list, pick_list, evaluate_residues, clear],
                                      layout={"width": "100%"}),
                         residue_output])
    widgets={"residue_input":residue_input,
             "preview":evaluate_residues,
             "clear":clear}
    return box, residue_idxs, widgets

def AA_Label_Options():
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

    AA_labels = _VBox([desc,
                                 _HBox([
                                     # tgl_short_AAs,
                                     tgl_hide_anchor,
                                     tgl_consensus_labs]),
                                 _HBox([tgl_color_hint, fontsize])],
                                layout={"width": "35%"})
    argmap = {"kwargs": {
        "assign_w_color": tgl_color_hint,
        "fontsize": fontsize,
        "anchor": tgl_hide_anchor,
    },
        "tgl_consensus": tgl_consensus_labs
    }
    return AA_labels, argmap

def Bar_Options():
    tgl_freqs_above = ipywidgets.ToggleButton(value=True,
                                              description="Don't show freqs >=",
                                              style={"description_width": "max-content"},
                                              layout={"width": "70%"},
                                              icon="check"
                                              )
    thrs_above = ipywidgets.FloatText(1.0,
                                      layout={"width": "30%"},
                                      step=.05)
    tgl_freqs_below = ipywidgets.ToggleButton(value=True,
                                              description="Don't show freqs <=",
                                              layout={"width": "70%"},
                                              icon="check"
                                              )
    thrs_below = ipywidgets.FloatText(0.2,
                                      style={"description_width": "max-content"},
                                      layout={"width": "30%"})


    identity = _VBox([
        ipywidgets.Button(description="Bar-options:",
                          layout={"width": "100%"}),
        _HBox([tgl_freqs_above, thrs_above]),
        _HBox([tgl_freqs_below, thrs_below]),
        # _HBox(colors)
    ],
        layout={"width": "30%"})

    argmap = {"kwargs": {
        "remove_identities": tgl_freqs_above,
        "identity_cutoff": thrs_above,
        "lower_cutoff_val": thrs_below,
    },
        "tgl_freqs_below": tgl_freqs_below,
    }

    return identity, argmap

def args_box():

    AA_labels, argmap1 = AA_Label_Options()

    Bars, argmap2 = Bar_Options()

    colors = color_pickers()

    argmap1["kwargs"].update(argmap2["kwargs"])
    argmap1.update({key:val for key, val in argmap2.items() if key!="kwargs"})


    for tgl in [argmap1["kwargs"]["remove_identities"],
                argmap1["tgl_freqs_below"],
                argmap1["kwargs"]["anchor"],
                argmap1["tgl_consensus"],
                argmap1["kwargs"]["assign_w_color"]]:
        tgl.observe(lambda traits : change_icon(traits),names="value")



    return _HBox([Bars,AA_labels],
                 layout={"width": "100%"}
                 ), argmap1

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

    green_button = widgets.Button(description="show neighborhoods", button_style="success",
                                layout={"width": "auto"}
                                )

    box = widgets.HBox([FreqSlider, green_button, clear_button])
    wdgs = {
        "FreqSlider":FreqSlider,
        "green_button":green_button,
        "clear_button":clear_button
    }
    return box, wdgs

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

def screen3(indict,top, individual_controls=False, **kwargs):

    residue_selection_box, res_idxs, res_wdg = residue_selection(top,**kwargs)

    residue_selection_acc = ipywidgets.Accordion([residue_selection_box],
                                                 layout={"width": "100%"})
    residue_selection_acc.set_title(0, "Residue Selection")

    cutoff_selection_box, co_wdgs = cutoff_selection_HBox()
    FreqSlider, run_button, clear_button = [co_wdgs[key] for key in ["FreqSlider", "green_button", "clear_button"]]

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

    accordion.set_title(0, "General Options")
    print("AAA",res_idxs)

    out_VBox = widgets.VBox([
        residue_selection_acc,
        accordion,
        cutoff_selection_box,
       output_acc,
    ])

    display(out_VBox)


    _plt.ioff()
    print(res_idxs)
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

        children=[]
        for ii in imgs:
            if individual_controls:
                FOs = figure_options()
                BOs, __ = Bar_Options()
                BOs.layout={"width":"90%"}
                AAs, __ = AA_Label_Options()
                AAs.layout={"width":"90%"}
                children.append(_VBox([_HBox([ii,
                                              _VBox([FOs, BOs,AAs]),
                                              ]),
                                       ipywidgets.Label(value="Separator"),
                                       ]
                                      ))
            else:
                children.append(ii)

        img_box.children=children

                            


    def lambda_compare_neighborhoods(CGs,argmap):
        anchor_str = mdciao.utils.residue_and_atom.shorten_AA(
            str(list(CGs.values())[0]._contacts[0].residues.anchor_residue),
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
                                                       title='%s neighborhood (4 bonded excluded)'% anchor_str
                                                )[0]
    # First consequence of clicking green : evaluate expression and update res_idx dictionary
    run_button.on_click(lambda _: res_wdg["preview"].click())
    # Now we can generate a lambda
    print("lambda",res_idxs)
    run_lambda = lambda out_dict : run(res_idxs, argmap)
    # No the lambda gets updated
    run_button.on_click(run_lambda, res_idxs)


    FreqSlider.observe(lambda _: run_lambda(None), names="value")
    for wdg in list(argmap["kwargs"].values())+[val for key, val in argmap.items() if key!="kwargs"]:
        wdg.observe( lambda first_run :  run_lambda(None), names="value")


    clear_button.on_click(lambda _ :   [img.close() for img in img_box.children])
    return out_VBox

def change_icon(traitlets):
    #assert traitlets["owner"].value and traitlets["owner"].icon=="check",traitlets
    if not traitlets["owner"].value:
        traitlets["owner"].icon=""
    else:
        traitlets["owner"].icon="check"