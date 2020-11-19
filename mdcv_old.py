import ipywidgets as widgets
#from ipyfilechooser import FileChooser
from ipywidgets import Layout
from IPython.display import display
import numpy as _np


from mdciao.cli import residue_neighborhoods
def _def_widget(**kwargs):
    return widgets.Text(
        **kwargs
    )
def topology(description="Select topology",
             placeholder="Topology file",
             **kwargs):
    # Create and display a FileChooser widget
    upl = SelectFilesButton()
    return upl

    # TODO DECIDE WHETHER TO KEEP OR TOSS THIS CODE
    upl = FileChooser()
    #upl._select.description = description
    upl._select.description = "select"
    upl._label.value = ""

    upl = widgets.FileUpload(
        description=description,
        accept='',  # Accepted file extension e.g. '.txt', '.pdf', 'image/*', 'image/*,.pdf'
        multiple=False  # True to accept multiple files upload else False
    )
    return upl

def trajectories(description="trajectories",
                 placeholder="can be sppace or comma separated, e.g. traj1.xtc, traj2.xtc traj3.xtc\ntraj4.xtc",
                 **kwargs):
    upl = SelectFilesButton()
    return upl
    # TODO DECIDE WHETHER TO KEEP OR TOSS THIS CODE
    upl = widgets.FileUpload(
        description=description,
        accept='',  # Accepted file extension e.g. '.txt', '.pdf', 'image/*', 'image/*,.pdf'
        multiple=True  # True to accept multiple files upload else False
    )
    return upl

def resSeq_idxs(description="resSeq_idxs",
                 placeholder="resSeq_idxs (see tooltip)",
                 description_tooltip="Can be in a format 1, 2-6, 10, 20-25"):
    return _def_widget(description=description,
                       placeholder=placeholder,
                       description_tooltip=description_tooltip,
                       #placeholder_tooltip=description_tooltip
                       )

def ctc_cutoff_Ang(description="ctc_cutoff_Ang",
                   value=3,
                    description_tooltip="The cutoff distance between two residues for them to be considered in contact. Default is 3 Angstrom",
                   ):
    return widgets.FloatText(description=description, value=value,
                             step=.5,
                             description_tooltip=description_tooltip)

def stride():
    return widgets.IntText(description="stride",
                           value=1, step=1,
                           description_tooltip="Stride down the input trajectoy files by this factor. Default is 1.")

def n_ctcs():
    return widgets.IntText(description="n_ctcs",
                           description_tooltip="Only the first n_ctcs most-frequent contacts will be written to the ouput. Default is 5.",
                           value="5")
def n_nearest():
    return widgets.IntText(description="n_nearest",
                           description_tooltip="Ignore this many nearest neighbors when computing neighbor lists. 'Near' means 'connected by this many bonds'. Default is 4.",
                           value="4")

def chunksize_in_frames():
    return widgets.IntText(description="chunksize_in_frames",
                           description_tooltip="Trajectories are read in chunks of this size (helps with big files and memory problems). Default is 10000",
                           value="10000",
                           step=500)

def nlist_cutoff_Ang():
    return widgets.IntText(description="nlist_cutoff_Ang",
                           description_tooltip="Cutoff for the initial neighborlist. Only atoms that are within this distance in the original reference "
                                               "(the topology file) are considered potential neighbors of the residues in resSeq_idxs, s.t "
                                               "non-necessary distances (e.g. between N-terminus and G-protein) are not even computed. "
                                               "Default is 15 Angstrom.",
                           value=15)

def output_ext():
    return widgets.Dropdown(
        options= [".pdf",".svg", ".png", ".jpg"],
        description="output_ext",
        description_tooltip="Extension of the output graphics"
    )

def order():
    return widgets.Checkbox(
        indent=False,
        value=True,
        description='sort',
        disabled=False,
        tooltip="Sort the resSeq_idxs list. Defaut is True")

def pbc():
    return widgets.Checkbox(value=True,
                            indent=False,
                            description="use PBC",
                            description_tooltip="Consider periodic boundary conditions when computing distances.")

def fragments():
    return widgets.Checkbox(value=True,
                            indent=False,
                            description="auto-detect fragments in the peptide chain",
                            description_tooltip="Auto-detect fragments (i.e. breaks) in the peptide-chain. Default is true.")
def fragment_names(**kwargs):
    return widgets.Textarea(
        #description="fragment_names",
                            placeholder="Name of the fragments. Leave empty if you want them automatically named."
                                        " Otherwise, give a quoted list of strings separated by commas, e.g. "
                                        "'TM1, TM2, TM3,'",
                            description_tooltip="Name of the fragments. Leave empty if you want them automatically named."
                                        " Otherwise, give a quoted list of strings separated by commas, e.g. "
                                        "'TM1, TM2, TM3,'", **kwargs)

def contact_control():

    out_dict = {
        "resSeq_idxs": resSeq_idxs(),
        "order": order(),
        "ctc_cutoff_Ang" : ctc_cutoff_Ang(),
        "nlist_cutoff_Ang" : nlist_cutoff_Ang(),
        "n_nearest" : n_nearest(),
    }
    out_VBox = widgets.VBox([widgets.Label("Contact Control:"),
                             widgets.HBox([out_dict["resSeq_idxs"], out_dict["order"]]),
                             out_dict["ctc_cutoff_Ang"],
                             out_dict["nlist_cutoff_Ang"],
                             out_dict["n_nearest"]],
                            layout=Layout(width="450px")
                            )
    out_VBox._gui_GUI_dict = out_dict



    return out_VBox


def IO_control():
    out_VBox = widgets.VBox([widgets.Label("I/O Control:"),
                             stride(),
                             chunksize_in_frames(),
                             n_ctcs(),
                             output_ext()],
                            layout=Layout(
                                margin='0 0 0 30px',
                                width="450px"))

    out_VBox._gui_GUI_dict = {"stride":out_VBox.children[1],
                              "chunksize_in_frames":out_VBox.children[2],
                              "n_ctcs":out_VBox.children[3],
                              "output_ext":out_VBox.children[4]}
    return out_VBox

def fragment_control():
    out_VBox= widgets.VBox([widgets.Label("Fragment Control:"),
                         fragments(),
                         fragment_names(layout=Layout(height="100px"))],
                          #layout=Layout(
                              #margin='0 0 0 5px',
                          #    width="250px"
                                        )
                          # )

    out_VBox._gui_GUI_dict = {"fragments":out_VBox.children[1],
                              "fragment_names":out_VBox.children[2]}

    return out_VBox

def file_control():
    out_dict = {"topology":topology(),
                "pbc":pbc(),
                "trajectories":trajectories()}
    out_VBox = widgets.VBox([widgets.Label("File Control:"),
                             widgets.HBox([widgets.Label("topology:    ", layout=Layout(width="100px")), out_dict["topology"], out_dict["pbc"]]),
                             widgets.HBox([widgets.Label("trajectories:", layout=Layout(width="100px")), out_dict["trajectories"], out_dict["pbc"]])],
                            #layout=Layout(width="400px")
                            )

    out_VBox._gui_GUI_dict = out_dict
    return out_VBox

def screen():

    run_button = widgets.Button(description="Compute residue neighborhoods")

    ifile_control = file_control()
    out_dict = ifile_control._gui_GUI_dict

    icontact_control = contact_control()
    out_dict.update(icontact_control._gui_GUI_dict)

    iIO_control = IO_control()
    out_dict.update(iIO_control._gui_GUI_dict)

    ifragment_control = fragment_control()
    out_dict.update(ifragment_control._gui_GUI_dict)

    out_VBox = widgets.VBox([widgets.HBox([ifile_control, run_button]),
                             widgets.Label(),
                             widgets.HBox([icontact_control, iIO_control]),
                             fragment_control(),
                                       ])

    out_VBox._gui_GUI_dict = out_dict
    display(out_VBox)
    # We have to define the function here (I think...)

    def run():
        output = widgets.Output()
        display(output)
        a = get_current_widget_values(out_dict)
        #print(a.can_start)
        if a.can_start:
            with output:
                residue_neighborhoods(a)

    run_lambda = lambda out_dict : run()

        #a =
        #a.resSeqidxs
        #pass
    output = widgets.Output()
    run_button.on_click(run_lambda)

    return out_VBox

class get_current_widget_values(object):
    def __init__(self, indict):
        self.BW_file = str(None)
        self.CGN_PDB = str(None)
        self.ask = True
        self.can_start = True
        for key, val in indict.items():
            if key == "topology":
                try:
                    setattr(self,key, val.files[0])
                except IndexError:
                    print("No topology given!")
                    self.can_start = False
                    break
            elif key == "trajectories":
                if len(val.files)==0:
                    print("No trajectories given!")
                    self.can_start = False
                setattr(self,key, val.files)
            elif key=="fragments":
                setattr(self,"fragmentify", val.value)
            else:
                setattr(self,key, val.value)
        if len(self.resSeq_idxs)==0:
            print("You have to input at least on residue!")
            self.can_start = False

class SelectFilesButton(widgets.Button):
    """A file widget that leverages tkinter.filedialog."""

    def __init__(self, description="select"):
        import traitlets
        super(SelectFilesButton, self).__init__()
        # Add the selected_files trait
        self.add_traits(files=traitlets.traitlets.List())
        # Create the button.
        self.description = description
        self.icon = "square-o"
        #self.style.button_color = "orange"
        self.style.button_color = "Gainsboro"

        # Set on click behavior.
        self.on_click(self.select_files)

    @staticmethod
    def select_files(b):
        from tkinter import Tk, filedialog

        """Generate instance of tkinter.filedialog.

        Parameters
        ----------
        b : obj:
            An instance of ipywidgets.widgets.Button
        """
        # Create Tk root
        root = Tk()
        # Hide the main window
        root.withdraw()
        # Raise the root to the top of all windows.
        root.call('wm', 'attributes', '.', '-topmost', True)
        # List of selected fileswill be set to b.value
        b.files = filedialog.askopenfilename(multiple=True)

        b.description = "(%u) %s "%(len(b.files), "files selected")
        b.icon = "check-square-o"
        b.style.button_color = "lightgreen"
