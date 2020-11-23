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
import functs
# TODO use interact instead of the argmaps!
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
                              layout={"width":"99%"})
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

    return _VBox([title], layout={"width":"99%"})

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
                                         layout={"width":"99%"})

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
        if residue_list.value not in res:
            residue_input.value = ",".join(res+[residue_list.value])


    pick_list.on_click(lambda _ : add_list_residue())


    evaluate_residues = ipywidgets.Button(description="preview",
                                          tooltip="Evaluate your residue selection before computing the neighborhoods",
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
                                      layout={"width": "99%"}),
                         residue_output])
    widgets={"residue_input":residue_input,
             "preview":evaluate_residues,
             "clear":clear}
    return box, residue_idxs, widgets

def AA_Label_Options():
    desc = ipywidgets.Button(description="AA-Label-Options:",
                             layout={"width": "99%"},
                             )

    """
    tgl_short_AAs = ipywidgets.ToggleButton(value=True,
                                              description="short AA names",
                                              icon="check",
                                              layout={"width": "50%"})
    """
    tgl_hide_anchor = ipywidgets.ToggleButton(value=True,
                                              tooltip="The AA shared by all contact pairs "
                                                      "is called 'anchor'. You can remove "
                                                      "it from the contact label.",
                                              description="toggle 'anchor'",
                                              icon="check",
                                              layout={"width": "50%"})

    tgl_consensus_labs = ipywidgets.ToggleButton(value=True,
                                                 tooltip="toggle consensus labels",
                                                 description="consensus labels",
                                                 icon="check",
                                                 layout={"width": "50%"})
    tgl_color_hint = ipywidgets.ToggleButton(value=False,
                                             tooltip="Use colors to hint if a single dataset is either "
                                                     "entirely missing or the only one present."
                                                     "E.g., if D30-K40 is the only"
                                                     " missing in dataset A (colored in purple), then you "
                                                     "will see the label '-D30-K40' in purple. Conversely, "
                                                     "if it is the only one present, you will see '+D30-K40"
                                                     " in puruple. Notice the signs '-' and '+' denoting "
                                                     "absence or presence.",
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
                                              description="hide freqs >=",
                                              tooltip="The pairs where all datasets show freqs >= this number "
                                                      "will not be plotted.\n"
                                                      "Their total sum will still appear "
                                                      "in the legend as 'a' (a='above threshold').\n"
                                                      "For example, if the first 'n' contacts are always"
                                                      " formed you can hide them to"
                                                      " unclutter the plot and focus on the differences.",
                                              style={"description_width": "max-content"},
                                              layout={"width": "70%"},
                                              icon="check"
                                              )
    thrs_above = ipywidgets.FloatText(1.0,
                                      layout={"width": "30%"},
                                      step=.05)
    tgl_freqs_below = ipywidgets.ToggleButton(value=True,
                                              tooltip="The pairs where all datasets show freqs <= this number "
                                                      "will not be plotted.\n"
                                                      "Their total sum will still appear "
                                                      "in the legend as 'b' (b='below threshold').\n"
                                                      "This unclutters the plot by truncating low frequencies.",

                                              description="hide freqs <=",
                                              layout={"width": "70%"},
                                              icon="check"
                                              )
    thrs_below = ipywidgets.FloatText(0.2,
                                      step=.05,
                                      style={"description_width": "max-content"},
                                      layout={"width": "30%"})


    identity = _VBox([
        ipywidgets.Button(description="Contact Frequency Options:",
                          layout={"width": "99%"}),
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

def general_options(horizontal=True):

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


    if horizontal:
        funct=_HBox
    else:
        funct=_VBox
    return funct([Bars,AA_labels],
                 layout={"width": "99%"}
                 ), argmap1

def freq_slider(value=3.5):
    return ipywidgets.FloatSlider(
        value=value,
        min=0,
        max=5,
        step=0.1,
        description='Neighbor Distance cutoff',
        style={"description_width": "auto"},
        readout_format='.1f',
        continuous_update=False,
    )


def cutoff_selection_HBox():
    clear_button = widgets.Button(description="clear neighborhoods", button_style="danger",
                                  layout={"width": "auto"})

    FreqSlider = freq_slider()

    green_button = widgets.Button(description="show neighborhoods", button_style="success",
                                layout={"width": "auto"}
                                )


    collap_button = widgets.ToggleButton(description="toggle neighborhoods", button_style="info",
                                layout={"width": "auto"}
                                )

    box = widgets.HBox([FreqSlider, green_button, clear_button,
                        collap_button
                        ])
    wdgs = {
        "FreqSlider":FreqSlider,
        "green_button":green_button,
        "clear_button":clear_button,
        "collapse_button":collap_button
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

def screen3(indict,top, individual_controls=False,
            start=False,
            **kwargs):

    residue_selection_box, res_idxs, res_wdg = residue_selection(top,**kwargs)

    residue_selection_acc = ipywidgets.Accordion([residue_selection_box],
                                                 layout={"width": "99%"})
    residue_selection_acc.set_title(0, "Residue Selection")

    runbox, co_wdgs = cutoff_selection_HBox()
    FreqSlider, run_button, clear_button, collapse_button = [co_wdgs[key] for key in ["FreqSlider",
                                                                                      "green_button",
                                                                                      "clear_button",
                                                                                      "collapse_button"]]

    Errors =  widgets.Output(layout={'border': '1px solid black'})

    img_box = widgets.VBox([],
                           layout={"width":"99%"})

    CGs = _defdict(dict)

    output_acc = widgets.Accordion([widgets.VBox([img_box,Errors])],
                                   layout={"width":"99%"})
    output_acc.set_title(0, "Output")

    options_wdg, main_argmap = general_options()
    tgl_per_figure = ipywidgets.ToggleButton(description={True:"on",
                                                          False:"off"}[individual_controls],
                                             value=individual_controls,
                                             icon={True:"check",
                                                   False:""}[individual_controls])

    options_wdg = _HBox(list(options_wdg.children)+[_VBox([ipywidgets.Button(description="Per-Figure-Options"),
                                                           tgl_per_figure,
                                                           ])])
    main_argmap["frequency"]=FreqSlider
    gen_opts_acc = ipywidgets.Accordion(
        [options_wdg],
        layout={"width": "99%"},
        selected_index=None)

    prog = ipywidgets.IntProgress(
    value=0,
    min=0,
    max=10,
    step=1,
    description='Loading:',
    bar_style='', # 'success', 'info', 'warning', 'danger' or ''
    orientation='horizontal'
)
    gen_opts_acc.set_title(0, "General Options")

    out_VBox = widgets.VBox([
        residue_selection_acc,
        gen_opts_acc,
        runbox,
        #prog,
       output_acc,
    ],
    layout={"width":"95%"}
    )

    display(out_VBox)


    _plt.ioff()

    per_fig_controls={}
    image_widgets = {}
    def run(res_idxs,argmap, individual_controls):

        _plt.close("all")
        plot_widgets = []
        prog.max=len(res_idxs)
        prog.value=0
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
            if len(CGs[rr])>0:

                img_wdg = functs.indv_plot(argmap,CGs[rr])

                fs = freq_slider()
                fs.description = "cutoff"

                indv_options_box, indv_argmap = general_options(horizontal=False)
                indv_argmap["frequency"]=fs

                indv_options_box = _VBox([fs]+list(indv_options_box.children)
                                         #+[FileChooser()]
                                         )
                [setattr(child,"layout",{"width":"100%"}) for child in indv_options_box.children]

                indv_options_box.layout = {"width":"31%",
                              #'border': '1px solid black'
                                           }

                [setattr(val, "_res_idx", rr) for key, val in indv_argmap.items() if key!="kwargs"]
                [setattr(val, "_res_idx", rr) for key, val in indv_argmap["kwargs"].items()]
                for key, val in argmap.items():
                    if key!="kwargs":
                        indv_argmap[key].value = val.value
                    else:
                        for key2,val2 in val.items():
                            indv_argmap["kwargs"][key2].value=val2.value
                per_fig_controls[rr] = indv_argmap
                image_widgets[rr] = img_wdg
            if individual_controls.value:
                acc_children = [_HBox([_HBox([img_wdg], layout={"width": "68%"}),
                                                        indv_options_box],
                                                       layout={"width": "99%"}),
                                                 ]
            else:
                acc_children = [_HBox([img_wdg], layout={"width": "68%"})]

            acc_wdig = ipywidgets.Accordion(acc_children)

            acc_wdig.set_title(0, str(top.residue(rr)))
            plot_widgets.append(acc_wdig)



            for rr, indv_argmap in per_fig_controls.items():
                indv_argmap["frequency"].observe(lambda change: update_lambda(change), names="value")
                for wdg in list(indv_argmap["kwargs"].values()) + [val for key, val in indv_argmap.items() if
                                                                   key != "kwargs"]:
                    wdg.observe(lambda change: update_lambda(change), names="value")

            img_box.children = plot_widgets
            prog.value+=1
        prog.bar_style="danger"

    # First consequence of clicking green : evaluate expression and update res_idx dictionary
    run_button.on_click(lambda _: res_wdg["preview"].click())
    # Now we can generate a lambda
    run_lambda = lambda __ : run(res_idxs, main_argmap, tgl_per_figure)
    # No the lambda gets updated
    run_button.on_click(run_lambda, res_idxs)

    def tgl_per_figure_effect(traits):
        tgl = traits["owner"]
        print(tgl)
        run_lambda(None)

    tgl_per_figure.observe(lambda traits: change_icon(traits), names="value")
    tgl_per_figure.observe(lambda traits: change_on_off(traits), names="value")
    tgl_per_figure.observe(lambda traits : tgl_per_figure_effect(traits), names="value")


    FreqSlider.observe(lambda __ : run_lambda(None), names="value")
    for wdg in list(main_argmap["kwargs"].values())+[val for key, val in main_argmap.items() if key!="kwargs"]:
        wdg.observe( lambda __ :  run_lambda(None), names="value")
    clear_button.on_click(lambda _ :   [img.close() for img in img_box.children])
    collapse_button.observe(lambda _ : [setattr(img,"selected_index",{0:None,
                                                                      None:0}[img.selected_index]) for img in img_box.children],
                            names="value")


    if start:
        run_button.click()

    def update_lambda(change):
        rr = change["owner"]._res_idx
        functs.indv_plot(per_fig_controls[rr],CGs[rr], im_wdg=image_widgets[rr])
        print()


    return out_VBox

def change_icon(traitlets):
    #assert traitlets["owner"].value and traitlets["owner"].icon=="check",traitlets
    if not traitlets["owner"].value:
        traitlets["owner"].icon=""
    else:
        traitlets["owner"].icon="check"

def change_on_off(traitlets):
    #assert traitlets["owner"].value and traitlets["owner"].icon=="check",traitlets
    if traitlets["owner"].value:
        traitlets["owner"].description="on"
    else:
        traitlets["owner"].description="off"