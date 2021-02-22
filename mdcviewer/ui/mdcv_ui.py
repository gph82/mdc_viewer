import ipywidgets as widgets
import numpy as _np
import mdciao
import ipywidgets
from mdciao.fragments.fragments import _allowed_fragment_methods
from mdciao.utils.residue_and_atom import rangeexpand_residues2residxs as _rangeexpand_residues2residxs
from ipywidgets import HBox as _HBox, VBox as _VBox
from IPython.display import display
from contextlib import redirect_stdout
from collections import namedtuple
from matplotlib import pyplot as _plt
import io
from mdcviewer import io as _mdcvio
from collections import defaultdict as _defdict
import mdtraj as _md
from mdcviewer.plots import plots
from typing import Tuple

# TODO use interact instead of the argmaps!
def show():
    top = _md.load("data/top.pdb.gz").top
    data = _mdcvio.load_data(verbose=0,
                             #decompress_here=False
                             )
    screen_neighborhoods(data,
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
                                           style={"description_width": "auto"})


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

def residue_selection(top, fragments=None, initial_value="R131,GDP*,L393-L394") -> Tuple[namedtuple, dict]:
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

    widgets = namedtuple("residue_selection_box",["itself","residue_input","preview","clear"])(box, residue_input, evaluate_residues, clear)
    return widgets, residue_idxs

class Indexed_Combobox(ipywidgets.Combobox):

    def __init__(self,*args,**kwargs):
        opts = kwargs.get("options")
        kwargs["options"]=["%04u: %s"%(ii,istr) for ii, istr in enumerate(opts)]
        ipywidgets.Combobox.__init__(self,*args,**kwargs)
        #self.value=kwargs["options"][0]

    @property
    def index(self):
        try:
            return int(self.value.split(":")[0])
        except:
            return None


def connected_bond_lists(top, bond_dict, l1_desc='or pick from list:') -> namedtuple:
    r""" Create two widgets containing residue pairs

    If you select residue "i" on the first widget,
    the second widget shows the set of {j,k,l...}
    residues for which a ContactPair can be created

    Parameters
    ----------
    top
    bond_dict

    Returns
    -------
    bond_from_list : namedtuple

    """

    bond_from_list_NamedTuple = namedtuple("bond_from_list", ["list1", "list2", "add_button"])


    residue_list_1 =Indexed_Combobox(

        options=[str(rr) for rr in top.residues],
        #rows=1,
        description=l1_desc,
        # disabled=False,
        layout={"width": "max-initial"},
        style={"description_width": "auto",
               },
        placeholder="Select and/or find residues"
    )

    """
    residue_list_1 =ipywidgets.Select(
        options=[str(rr) for rr in top.residues],
        rows=1,
        description=l1_desc,
        # disabled=False,
        layout={"width": "max-content", },
        style={"description_width": "auto",
               },
        #placeholder="Select and/or find residues"
    )
    """

    residue_list_2 = ipywidgets.Select(
        options=[],
        rows=1,
        description='',
        # disabled=False,
        layout={"width": "max-content", },
        style={"description_width": "auto",
               },
    )

    # Populate the second residue with bond-partners
    def update_residue_list_2():
        if residue_list_1.index is None:
            idxs=[None]
            opts = idxs
        else:
            idxs = bond_dict[residue_list_1.index]
            opts = [str(top.residue(ii)) for ii in idxs]

        residue_list_2.options=tuple(opts)
    update_residue_list_2()

    residue_list_1.observe(lambda _ : update_residue_list_2(), names="value")
    add_bond_from_list = ipywidgets.Button(description="add",
                                           layout={"width": "auto"})
    return bond_from_list_NamedTuple( residue_list_1, residue_list_2, add_bond_from_list)

def bond_selection(top, bond_dict, fragments=None,
                   #initial_value="R131-E392"
                   ) -> Tuple[namedtuple, dict]:
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

    manual_bond_input = ipywidgets.Textarea(description="Input manually",
                                    style={"description_width":"auto",
                                           },
                                    #value=initial_value,
                                    continuous_update=False,
                                    layout={"width":"500px"}
                                    )

    evaluated_bonds = ipywidgets.Textarea(placeholder="Evaluated list of bonds will appear here",
                                         layout={"width":"99%"})

    bond_from_list_nt = connected_bond_lists(top, bond_dict)
    #residue_list_1, residue_list_2, add_bond_from_list =

    residue_idxs = {}
    residue_idxs["res"] = []

    def add_bond_from_lists():
        ii, jj = bond_from_list_nt.list1.index, bond_dict[bond_from_list_nt.list1.index][bond_from_list_nt.list2.index]
        if None not in [ii, jj]:
            evaluated_bonds.value += "%s[%u]-%s[%u]\n" % (bond_from_list_nt.list1.value, ii,
                                                          bond_from_list_nt.list2.value, jj)

            residue_idxs["res"].append([ii, jj])

    bond_from_list_nt.add_button.on_click(lambda _: add_bond_from_lists())



    evaluate_bonds = ipywidgets.Button(description="preview",
                                          tooltip="Evaluate your bond selection before computing the sites",
                                       layout={"width":"auto"},
                                       button_style="info")

    clear = ipywidgets.Button(description="clear",
                              layout={"width":"auto"},
                              button_style="danger")
    clear.on_click(lambda _ : setattr(manual_bond_input,"value",""))
    clear.on_click(lambda _ : setattr(evaluated_bonds,"value",""))

    #manual_bond_input.observe(lambda _ : eval_bond_selection(), names="value")
    if fragments is None:
        fragments = [_np.arange(top.n_residues)]

    def eval_manual_bonds():
        if len(manual_bond_input.value)>0:
            try:
                b = io.StringIO()
                with redirect_stdout(b):
                    for line in manual_bond_input.value.splitlines():
                        for pair in line.split(","):
                            ii, jj = [int(kk) for kk in pair.split("-")]
                            evaluated_bonds.value += "%s[%u]-%s[%u]\n" % (top.residue(ii), ii,
                                                                          top.residue(jj), jj)
                            residue_idxs["res"].append([ii,jj])

            except Exception as e:
                #residue_idxs["res"] = []
                istr = str(e)
                evaluated_bonds.value = istr

    evaluate_bonds.on_click(lambda _ : eval_manual_bonds())
    box = ipywidgets.VBox([ipywidgets.HBox([manual_bond_input, bond_from_list_nt.list1, bond_from_list_nt.list2, bond_from_list_nt.add_button,
                                            evaluate_bonds,
                                            clear],
                                      layout={"width": "99%"}),
                         evaluated_bonds])
    bond_selection_box = namedtuple("bond_selection_box",["itself", "manual_bond_input", "evaluate_bonds", "clear_input"])(box, manual_bond_input,evaluate_bonds,clear)

    return bond_selection_box, residue_idxs

AA_labels_NamedTuple = namedtuple("AA_labels",
                                  ["itself", "tgl_hide_anchor", "tgl_consensus_labs",
                                   "tgl_color_hint", "fontsize", "argmap"])
def AA_Label_Options() -> AA_labels_NamedTuple:
    r""""
    Widget for controlling AA-labeling options:

    The options are:
     * toggle anchor
     * consensus labels
     * color hint
     * fontsize

    Returns : widget, argmap
    widget : _HBox
    argmap : dict
    """
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
                                              layout={"width": "49%"})

    tgl_consensus_labs = ipywidgets.ToggleButton(value=True,
                                                 tooltip="toggle consensus labels",
                                                 description="consensus labels",
                                                 icon="check",
                                                 layout={"width": "49%"})
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
                                             layout={"width": "49%"})
    fontsize = ipywidgets.IntText(value=16, description="fontsize",
                                  layout={"width": "49%"})

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

    AA_labels_nt = AA_labels_NamedTuple(AA_labels, tgl_hide_anchor, tgl_consensus_labs, tgl_color_hint, fontsize, argmap)

    return AA_labels_nt


Bars_NamedTuple = namedtuple("Bars_Options", ["itself", "tgl_freqs_above", "thrs_above",
                                              "tgl_freqs_below", "thrs_below", "freqslider",
                                              "argmap"
                                              ])
def Bar_Options(freqslider=None) -> Bars_NamedTuple:
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

    if freqslider is None:
        identity = _VBox([
            ipywidgets.Button(description="Contact Frequency Options:",
                              layout={"width": "99%"}),
            _HBox([tgl_freqs_above, thrs_above]),
            _HBox([tgl_freqs_below, thrs_below]),
            # _HBox(colors)
        ],
            layout={"width": "30%"})
    else:
        identity = _VBox([
            ipywidgets.Button(description="Contact Frequency Options:",
                              layout={"width": "99%"}),
            _HBox([freqslider]),
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

    Bars_nt = Bars_NamedTuple(identity,tgl_freqs_above, thrs_above,tgl_freqs_below, thrs_below, freqslider, argmap)
    return Bars_nt

options_NamedTuple = namedtuple("option_box", ["itself", "AA_labels", "Bars", "Freqslider", "tgl_per_figure", "per_fig_positon", "argmap"])

def general_options(freqslider=None,horizontal=True, individual_controls=False,
                    is_general=True) -> options_NamedTuple:
    r"""

    Instantiate a widget offering options to control a frequency plot

    Parameters
    ----------
    freqslider : :obj:`~ipywidgets.Slider`
        Passed to :obj:`Bar_Options`, which
        will create one if None is passed.
    horizontal : bool, default is True
        Whether to place the widgets for
        AA_label-control and Bar-control
        next to each other (default) or
        on top of each other

    Returns
    -------
    nt : namedtuple
        A :obj:`~collections.namedtuple` with
        the fields "itself" and others (TODO)
    """

    AA_labels_nt = AA_Label_Options()

    Bars_nt = Bar_Options(freqslider=freqslider)

    colors = color_pickers()

    argmap1 = AA_labels_nt.argmap
    argmap1["kwargs"].update(Bars_nt.argmap["kwargs"])
    argmap1.update({key:val for key, val in Bars_nt.argmap.items() if key!="kwargs"})
    argmap1.update({"frequency":freqslider})

    for tgl in [argmap1["kwargs"]["remove_identities"],
                argmap1["tgl_freqs_below"],
                argmap1["kwargs"]["anchor"],
                argmap1["tgl_consensus"],
                argmap1["kwargs"]["assign_w_color"]]:
        tgl.observe(lambda traits : change_icon(traits),names="value")

    if horizontal:
        box =_HBox([Bars_nt.itself, AA_labels_nt.itself],
                 layout={"width": "99%"}
                 )
    else:
        box=_VBox([Bars_nt.itself,AA_labels_nt.itself],
                 layout={"width": "99%"}
                 )

    tgl_per_figure, per_fig_position, itself = None, None, box
    if is_general:
        tgl_per_figure = ipywidgets.ToggleButton(description={True: "on",
                                                              False: "off"}[individual_controls],
                                                 value=individual_controls,
                                                 icon={True: "check",
                                                       False: ""}[individual_controls])
        per_fig_position = ipywidgets.Dropdown(description="position:",
                                                    options=["left", "right", "bottom"],
                                                    style={"description_width": "initial"},
                                                    layout={"width": "max-content"},
                                                    value="right")
        itself = _HBox(
            list(box.children) + [_VBox([ipywidgets.Button(description="Per-Neighborhood-Options"),
                                                           tgl_per_figure,
                                                           per_fig_position,
                                                           ])])

    return options_NamedTuple(itself, AA_labels_nt, Bars_nt, freqslider, tgl_per_figure, per_fig_position, argmap1)

def freq_slider(value=3.5):

    return widgets.SelectionSlider(
        options=["%2.1f Ang"%ii for ii in _np.linspace(0,5,50)],
        value="3.5 Ang",
        description='Neighbor Distance cutoff',
        style={"description_width": "auto"},
        disabled=False,
        continuous_update=False,
        readout=True
    )

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

runbox = namedtuple("runbox", ["itself", "FreqSlider", "green_button", "clear_button", "collapse_button"])
def run_Hbox(CG_type)-> runbox:
    r"""
    Prepare the runbox, i.e. the box with the buttons that make it GUI run

    There's other GUI buttons, but only the green-button of this widget
    will make the program "run". By "run" we mean:

    - thing 1 (evaluate the Selection box options)
    - thing 2 (grab data, creating new CPs or CGs if necessary)
    - thing 3 (update the needed plots)


    Parameters
    ----------
    CG_type: str

    Returns
    -------
    runbox : namedtuple

    """
    clear_button = widgets.Button(description="clear %s"%CG_type, button_style="danger",
                                  layout={"width": "auto"})

    FreqSlider = freq_slider()

    green_button = widgets.Button(description="show %s"%CG_type, button_style="success",
                                layout={"width": "auto"}
                                )


    collap_button = widgets.ToggleButton(description="toggle %s"%CG_type, button_style="info",
                                layout={"width": "auto"}
                                )

    box = widgets.HBox([FreqSlider, green_button, clear_button,
                        collap_button
                        ])

    return runbox(
        box,
        FreqSlider,
        green_button,
        clear_button,
        collap_button)



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

def screen_neighborhoods(indict, top, individual_controls=False,
                         start=False,
                         **kwargs):

    residue_selection_box_nc, res_idxs = residue_selection(top,**kwargs)

    residue_selection_acc = ipywidgets.Accordion([residue_selection_box_nc.itself],
                                                 layout={"width": "99%"})
    residue_selection_acc.set_title(0, "Residue Selection")

    prepare_GUI_WIP(indict, top, individual_controls,
                    residue_selection_acc,
                    residue_selection_box_nc.preview, res_idxs, start)

def screen_sites(indict, top, individual_controls=False,
                 start=False,
                 **kwargs):

    bond_dict=get_bond_dict(indict,top)
    bond_selection_box_nc, res_idxs = bond_selection(top, bond_dict, **kwargs)
    bond_selection_box_nc.manual_bond_input.layout.width="25%"
    bond_selection_acc = ipywidgets.Accordion([bond_selection_box_nc.itself],
                                                 layout={"width": "99%"})
    bond_selection_acc.set_title(0, "Bond Selection")
    wid_dict = prepare_GUI_WIP(indict, top, individual_controls,
                               bond_selection_acc, bond_selection_box_nc.evaluate_bonds, res_idxs,
                               start=False,
                               sites=True)

def get_bond_dict(indict, top):
    bond_dict = _defdict(list)
    for key, val in indict.items():
        for r1, r1dict in val.items():
            bond_dict[r1].extend(list(r1dict.keys()))
    bond_dict = {key: _np.unique(val) for key, val in bond_dict.items()}
    bond_dict.update({key: [None] for key in range(top.n_residues) if key not in list(bond_dict.keys())})
    return bond_dict

def prepare_GUI_WIP(indict, top, individual_controls,
                    selection_accordion,
                    preview_button,
                    selection_dict,
                    start,
                    sites=False):
    r"""

    Parameters
    ----------
    indict
    top
    individual_controls
    selection_accordion
    preview_button
    selection_dict
    start
    sites

    Returns
    -------

    """
    runbox_nt = run_Hbox(CG_type={False: "neighborhoods",
                                  True:  "sites"}[sites])

    Errors = widgets.Output(layout={'border': '1px solid black'})

    img_box = widgets.VBox([],
                           layout={"width": "99%"})

    output_acc = widgets.Accordion([widgets.VBox([img_box, Errors])],
                                   layout={"width": "99%"})
    output_acc.set_title(0, "Output")

    output_acc_nt = namedtuple("output_accordion",["img_box","Errors"])(img_box,Errors)

    options_wdg_nt = general_options(freqslider=runbox_nt.FreqSlider, individual_controls=individual_controls)



    options_wdg_nt.argmap["frequency"] = runbox_nt.FreqSlider
    gen_opts_acc = ipywidgets.Accordion(
        [options_wdg_nt.itself],
        layout={"width": "99%"},
        selected_index=None)

    prog = ipywidgets.IntProgress(
        value=0,
        min=0,
        max=10,
        step=1,
        description='Loading:',
        bar_style='',  # 'success', 'info', 'warning', 'danger' or ''
        orientation='horizontal'
    )
    gen_opts_acc.set_title(0, "General Options")

    out_VBox = widgets.VBox([
        selection_accordion,
        gen_opts_acc,
        runbox_nt.itself,
        # prog,
        output_acc,
    ],
        layout={"width": "95%"}
    )

    display(out_VBox)

    _plt.ioff()

    GUI = namedtuple("GUI",["itself","runbox", "output_acc"])(out_VBox, runbox_nt, output_acc)

    CGs = _defdict(dict)
    per_fig_controls = {}
    image_widgets = {}
    if sites:
        run = prepare_run_function_sites(indict, prog, output_acc_nt, CGs, per_fig_controls,
                                         image_widgets,
                                         options_wdg_nt.per_fig_positon,
                                         top)
    else:
        run = prepare_run_function_neighborhoods(indict, prog, output_acc_nt, CGs, per_fig_controls,
                                                 image_widgets,
                                                 options_wdg_nt.per_fig_positon,
                                                 top,
                                                 )

    # First consequence of clicking green : evaluate expression and update res_idx dictionary
    runbox_nt.green_button.on_click(lambda _: preview_button.click())
    # Now we can generate a lambda
    run_lambda = lambda __: run(selection_dict, options_wdg_nt.argmap, options_wdg_nt.tgl_per_figure)
    # Now the lambda gets updated
    runbox_nt.green_button.on_click(run_lambda)

    options_wdg_nt.tgl_per_figure.observe(lambda traits: change_icon(traits), names="value")
    options_wdg_nt.tgl_per_figure.observe(lambda traits: change_on_off(traits), names="value")
    options_wdg_nt.tgl_per_figure.observe(lambda __: run_lambda(None), names="value")

    options_wdg_nt.per_fig_positon.observe(lambda __: run_lambda(None), names="value")

    runbox_nt.FreqSlider.observe(lambda __: run_lambda(None), names="value")
    for wdg in list(options_wdg_nt.argmap["kwargs"].values()) + [val for key, val in options_wdg_nt.argmap.items() if key != "kwargs"]:
        wdg.observe(lambda __: run_lambda(None), names="value")
    runbox_nt.clear_button.on_click(lambda _: [img.close() for img in img_box.children])
    runbox_nt.collapse_button.observe(lambda _: [setattr(img,
                                                         "selected_index", {0: None,
                                                                            None: 0}[img.selected_index]) for img
                                                 in
                                                 img_box.children],
                                      names="value")

    if start:
        runbox_nt.green_button.click()

    return GUI

def prepare_run_function_neighborhoods(indict, prog, output_acc_nt, CGs, per_fig_controls, image_widgets, per_figure_horizontal, top,
                                      ):

    def update_lambda(change):
        rr = change["owner"]._res_idx
        plots.indv_plot(per_fig_controls[rr], CGs[rr], im_wdg=image_widgets[rr])
        print()


    def funct(res_idxs, argmap, individual_controls):
        _plt.close("all")
        plot_widgets = []
        prog.max = len(res_idxs)
        prog.value = 0
        for rr in res_idxs["res"]:
            rr = int(rr)
            output_acc_nt.Errors.clear_output()
            for key, iarch in indict.items():
                if rr in iarch.keys():
                    try:
                        CGs[rr][key]
                    except KeyError:
                        CGs[rr][key] = _mdcvio.CGdict2CG(iarch[rr], top=top)
                else:
                    with output_acc_nt.Errors:
                        output_acc_nt.Errors.clear_output(wait="True")
                        print("%s was found in the topology but not in the indict" % top.residue(rr))
            if len(CGs[rr]) > 0:

                img_wdg = plots.indv_plot(argmap, CGs[rr])

                fs = freq_slider()
                fs.description = "cutoff"
                fs.layout = {"width": "99%"}

                indv_options_box_nt = general_options(freqslider=fs, horizontal=False)
                indv_argmap = indv_options_box_nt.argmap
                indv_argmap["frequency"] = fs

                [setattr(child, "layout", {"width": "100%"}) for child in indv_options_box_nt.itself.children]

                [setattr(val, "_res_idx", rr) for key, val in indv_argmap.items() if key != "kwargs"]
                [setattr(val, "_res_idx", rr) for key, val in indv_argmap["kwargs"].items()]
                for key, val in argmap.items():
                    if key != "kwargs":
                        indv_argmap[key].value = val.value
                    else:
                        for key2, val2 in val.items():
                            indv_argmap["kwargs"][key2].value = val2.value
                per_fig_controls[rr] = indv_argmap
                image_widgets[rr] = img_wdg
            if individual_controls.value:
                if per_figure_horizontal.value in ["right", "left"]:
                    indv_options_box = _VBox(indv_options_box_nt.itself.children)
                    # +[FileChooser()]
                    indv_options_box.layout = {"width": "31%",
                                               # 'border': '1px solid gray'
                                               }
                    for_HBox = [_HBox([img_wdg], layout={"width": "68%"}),
                                indv_options_box_nt.itself]
                    acc_children = [_HBox([for_HBox[ii] for ii in {"right": [0, 1],
                                                                   "left": [1, 0]}[per_figure_horizontal.value]],
                                          layout={"width": "99%"}),
                                    ]
                elif per_figure_horizontal.value == 'bottom':
                    acc_options = ipywidgets.Accordion([indv_options_box_nt.itself])
                    acc_options.set_title(0, "Options")
                    acc_children = [_VBox([_HBox([img_wdg], layout={"width": "99%"}),
                                           acc_options],
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

            output_acc_nt.img_box.children = plot_widgets
            prog.value += 1
        prog.bar_style = "danger"


    return funct

def ensure_CP_and_return(ii, jj, iarch, top):
    r"""
    Check if iarch[ii][jj] is ContactPair, if not overwrite in place, also for [jj][ii]
    Parameters
    ----------
    ii
    jj
    iarch
    top

    Returns
    -------

    """
    try:
        iarch[ii][jj] = _mdcvio.sCP2CP(iarch[ii][jj], top=top)
        iarch[jj][ii] = iarch[ii][jj]
        return iarch[ii][jj]
    except:
        return None


def prepare_run_function_sites(indict, prog, output_acc_nt, CGs,
                               per_fig_controls, image_widgets, per_figure_horizontal, top,
                                      ):
    def funct(res_idxs, argmap, individual_controls):
        _plt.close("all")
        prog.max = len(res_idxs)
        prog.value = 0
        output_acc_nt.Errors.clear_output()
        with output_acc_nt.Errors:
            for key, iarch in indict.items():
                try:
                    CPs = []
                    for ii, jj in res_idxs["res"]:
                        iCG = ensure_CP_and_return(ii, jj, iarch, top)
                        if iCG is not None:
                            CPs.append(iCG)

                    print(key, len(CPs))
                    #(print(id(iCG.top), iCG.top.__hash__()))
                    if len(CPs)>0:
                        CGs[0][key] = mdciao.contacts.ContactGroup(CPs,
                                                               neighbors_excluded=0,
                                                               top=top,
                                                               name="site")

                except Exception as e:
                    print(e)
                    print("There should be an exception message above this")
                    raise

            if len(CGs[0]) > 0:
                for key, val in CGs[0].items():
                    print(key)
                    print(val.frequency_table(5,None))
                if individual_controls.value:
                    _argmap=None
                else:
                    _argmap = argmap
                __, indv_acc_wdg_nt = CG2_accordion(CGs[0],
                                                    general_argmap=_argmap,
                                                    indv_argmap=None, top=top, indict=indict)

                output_acc_nt.img_box.children = [indv_acc_wdg_nt.itself]
                prog.value += 1
            prog.bar_style = "danger"


    return funct

indv_acc_wdg_nt1 = namedtuple("individual_accordion",["itself","indv_options_box_nt", "image_wdg", "CG_dict"])
def CG2_accordion(CG, indv_argmap=None, general_argmap=None, img_wdg=None, site=True, top=None, indict=None) -> Tuple[dict, indv_acc_wdg_nt1]:
    r"""
    ContactGroup to per-plot accordion, either with minimal or very informed input

    If no controls are given at all, assume it's an individual plot

    Parameters
    ----------
    CG  : dict
        Contains :obj:`~mdciao.contacts.ContactGroup`
    indv_argmap : dict, default is None
        Passed to :obj:`~mdcviewer.plots.indv_plot`, which
        calls :obj:`~mdciao.plots.compare_groups_of_contacts`.
        If None, one argmap will be created  by calling :obj:`general_options`
    general_argmap : dict, default is None
        If None, one will be created using
        :obj:`general_options`
    im_wdg : :obj:`ipywidgets.Image` or None
        Pass an existing image widget here to be updated,
        else one will be instantiated and returned
    site : bool, default is False
        Plot as site and not neighborhood

    Returns
    -------
    indv_argmap : dict
    indv_acc_wdg_nt : namedtuple

    """


    fs = freq_slider()
    fs.description = "cutoff"
    fs.layout = {"width": "99%"}

    if indv_argmap is None:
        if general_argmap in [None,False]:
            options_box_nt = general_options(freqslider=fs, is_general=False)
            individual_controls_value = True
            argmap = options_box_nt.argmap
            [setattr(child, "layout", {"width": "100%"}) for child in options_box_nt.itself.children]

        else:
            argmap = general_argmap
            options_box_nt = None
            individual_controls_value=False
        img_wdg = plots.indv_plot(argmap, CG, im_wdg=img_wdg, site=site)

    else:
        raise(NotImplementedError)

    if individual_controls_value:

        if True:  # per_figure_horizontal.value in ["right", "left"]:
            options_box = _VBox(options_box_nt.itself.children)
            options_box.layout = {"width": "98%",
                                          # 'border': '1px solid gray'
                                          }
            tab_children = [options_box]
            if top is not None and indict is not None:
                bond_dict = get_bond_dict(indict,top) #TODO this gets recalc'd every time
                bond_from_list_nt = connected_bond_lists(top, bond_dict, l1_desc="add pair")
                tab_children.append(_HBox([bond_from_list_nt.list1,bond_from_list_nt.list2, bond_from_list_nt.add_button],
                                          #layout={"width": "68%%"}
                                          ))
                bond_from_list_nt.add_button.on_click(lambda __ : update_CG())
                def update_CG():
                    for key in list(CG.keys()):
                        ii, jj = bond_from_list_nt.list1.index, bond_from_list_nt.list2.index
                        r1_idx, r2_idx = ii, bond_dict[ii][jj]
                        try:
                            ensure_CP_and_return(r1_idx, r2_idx, indict[key], top)
                            iCG = mdciao.contacts.ContactGroup(CG[key]._contacts + [indict[key][r1_idx][r2_idx]],
                                                               top=top,
                                                               neighbors_excluded=0)
                            CG[key]=iCG
                        except KeyError:
                            pass
                    plots.indv_plot(argmap, CG, im_wdg=img_wdg, site=True)

            # TODO make a named tuple out of this option box
            options_box = ipywidgets.Tab(tab_children, layout={"width":"33%"})
            options_box.set_title(0,"Options...")
            options_box.set_title(1,"Add...")
            for_HBox = [_HBox([img_wdg], layout={"width": "68%"}),
                        options_box]

            indv_acc_children = [_HBox([for_HBox[ii] for ii in {"right": [0, 1],
                                                                "left": [1, 0]}["right"]],
                                       # [per_figure_horizontal.value]],
                                       layout={"width": "99%"}),
                                 ]

    else:
        indv_acc_children = [_HBox([img_wdg], layout={"width": "68%"})]
        options_box = None

    indv_acc_wdg = ipywidgets.Accordion(indv_acc_children)
    indv_acc_wdg.set_title(0, "test")
    indv_acc_wdg_nt = indv_acc_wdg_nt1(indv_acc_wdg, options_box, img_wdg, CG)

    # This is where we link the widgets and the plot
    argmap["frequency"].observe(lambda __ : plots.indv_plot(argmap, CG, im_wdg=img_wdg, site=True), names="value")
    for wdg in list(argmap["kwargs"].values()) + [val for key, val in argmap.items() if
                                                       key != "kwargs"]:
        wdg.observe(lambda __: plots.indv_plot(argmap, CG, im_wdg=img_wdg, site=site), names="value")

    return  argmap, indv_acc_wdg_nt


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