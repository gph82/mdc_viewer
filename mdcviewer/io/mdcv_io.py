import mdciao
import os.path as _path
import numpy as _np
from joblib import Parallel, delayed
from tqdm import tqdm_notebook as _tqdm
import h5py

import numpy as _np


def load_data(
        files={
            "A": "data/neighborhoods.A.B2AR.zip.stride.10.hdf5",
            "B": "data/neighborhoods.B.B2AR.zip.stride.10.hdf5",
            "C": "data/neighborhoods.C.B2AR.zip.stride.10.hdf5"
        },
        verbose=0,
        database=False,
        **kwargs_hd5_2_archives,
       ):
    r"""
    Load pre-computed data in form of ``hdf5`` files.

    Parameters
    ----------
    files : dict
        Dictionary with the paths to the hd5f-files
        The keys will be used as state descriptors
    verbose : int, default is 1
        The verbosity of the joblib-Parallel read-in
    database : bool, default is False
        Whether to wrap around :obj:`hd5_2_database`
        or :obj:`hd5_2_dict_of_CGdicts` when reading
    kwargs_hd5_2_archives: dict
        Will be passed to the hd5_2_* functions

    Returns
    -------
    data : dict

    """
    if database:
        funct = hd5_2_database
    else:
        funct = hd5_2_dict_of_CGdicts

    Ns = Parallel(n_jobs=1, verbose=verbose)(delayed(funct)(ff,
                                                            **kwargs_hd5_2_archives
                                                            ) for ff in _tqdm(files.values(),
                                                                              desc="loading data"
                                                                              ))

    return {key:val for key, val in zip(files.keys(),Ns)}
def hd5_2_dict_of_CGdicts(obj, restrict_to_residxs=None, decompress_here=True):
    r"""

    Parameters
    ----------
    obj : string or :obj:`h5py.File` object
        If string, path to an .hdf5 file

    Returns
    -------
    neighborhoods : dict
        Dictionary keyed with residue indices
        and valued with :obj:`ContactGroups`
        representing the neighborhood of
        a residue with that index

    """
    if _path.exists(obj):
        h5py_fileobject = h5py.File(obj,"r")
    else:
        h5py_fileobject = obj

    if restrict_to_residxs is None:
        valid_res = lambda res : True
    else:
        valid_res = lambda res : res in restrict_to_residxs

    output_dict = {}
    if "compress" in h5py_fileobject.keys() and h5py_fileobject["compress"][()]:
        needs_decompression=True
        try:
            ref_t = h5py_fileobject["ref_t"][()]
            #TODO this is while we adapt the old hdf5 files
        except:
            ref_t = []
            for key in sorted([key for key in h5py_fileobject.keys() if key.startswith("ref_t")]):
                ref_t.append(h5py_fileobject[key][()])
    for key, CGdict in h5py_fileobject.items():
        if key.isdigit() and valid_res(int(key)):

            CG = {}

            serialized_CPs = list(
                [{key: val[()] for key, val in dict(val).items()} for key, val in CGdict.items() if key.isdigit()])

            CG["serialized_CPs"] = [decode_dict_values(idict) for idict in serialized_CPs]
            if needs_decompression:
                for sCP in serialized_CPs:
                    sCP["time_traces.time_trajs"] = ref_t
                    if decompress_here:
                        decompress_serialized_CP(sCP)
            iname = CGdict["name"][()]
            try:
                iname = iname.decode()
            except AttributeError:
                pass
            CG["name"] = [None if iname.lower()=="none" else iname][0]
            CG["interface_residxs"]=CGdict["interface_residxs"][()]
            CG["neighbors_excluded"]=CGdict["neighbors_excluded"][()]
            output_dict[int(key)] = CG

    return output_dict

def hd5_2_database(obj, restrict_to_residxs=None, decompress_here=True, database=False):
    r"""

    Return a per-contact datatabse. The database is a dictionary of dictionaries,
    keyed with residue indices and valued with residue-residue mindist values

    Parameters
    ----------
    obj : string or :obj:`h5py.File` object
        If string, path to an .hdf5 file

    Returns
    -------
    neighborhoods : dict
        Dictionary keyed with residue indices
        and valued with :obj:`ContactGroups`
        representing the neighborhood of
        a residue with that index

    """
    if _path.exists(obj):
        h5py_fileobject = h5py.File(obj,"r")
    else:
        h5py_fileobject = obj

    if restrict_to_residxs is None:
        valid_res = lambda res : True
    else:
        valid_res = lambda res : res in restrict_to_residxs

    output_dict = {}
    if "compress" in h5py_fileobject.keys() and h5py_fileobject["compress"][()]:
        needs_decompression=True
        ref_t = h5py_fileobject["ref_t"][()]
    outputCPs = {}
    for key, CGdict in h5py_fileobject.items():
        #print(key)
        if key.isdigit() and valid_res(int(key)):
            CG = {}

            serialized_CPs = list(
                [{key: val[()] for key, val in dict(val).items()} for key, val in CGdict.items() if key.isdigit()])

            for CP in serialized_CPs:
                ii, jj = CP["residues.idxs_pair"]
                a, b = _np.sort([ii, jj])
                if a not in outputCPs.keys():
                    outputCPs[a]={}
                if b not in outputCPs.keys():
                    outputCPs[b]={}
                if outputCPs[a].get(b) is None:
                    CP = decode_dict_values(CP)
                    if needs_decompression:
                        CP["time_traces.time_trajs"] = ref_t
                        if decompress_here:
                            decompress_serialized_CP(CP)
                    outputCPs[a][b] = CP#["time_traces.ctc_trajs"]
                    assert outputCPs[b].get(a) is None
                    outputCPs[b][a] = outputCPs[a][b]
    return outputCPs

from copy import deepcopy as _deepcopy
def decompress_serialized_CP(sCP,inplace=True):
    r"""
    Decompress a :obj:`~mdciao.contacts.ContactPair` object that's been compressed when serialized

      * ``time_traces.atom_pair_trajs`` get turned to pairs of 2d nd.arrays of ints (originally compressed to CSV-strings)
      * ``time_traces.ctc_trajs`` get turned to lists of floats and divided by 1000 (originally compressed to *1000- CSV-strings)

    Parameters
    ----------
    sCP : dict
        :obj:`~mdciao.contacts.ContactPair`
        serialized to a dict
    inplace : bool, default is True
        Decompress in place, else
        return a copy and leave
        :obj:`sCP` intact

    Returns
    -------
    dsCP : dict
        Only returns when :obj:`inplace`

    """

    if inplace:
        out = sCP
    else:
        out = _deepcopy(sCP)
    out["time_traces.atom_pair_trajs"]=[_np.vstack([[int(aa) for aa in pp.split()] for pp in tt.split(",")]) for tt in out["time_traces.atom_pair_trajs"]]
    out["time_traces.ctc_trajs"] = [[float(ff) / 1000 for ff in tt.split()] for tt in out["time_traces.ctc_trajs"]]
    if not inplace:
        return  out

def decode_dict_values(idict):
    r"""
    Decode strings stored as numpy arrays only
    for those values of the dict that need it

    The need for this is strings stored in hdf5 format
    Parameters
    ----------
    idict : dict
        Dictionary with some of its values being
        strings encoded as numpy arrays

    Returns
    -------
    idict : dict
        Same as dict with its its decodeable(?)
        strings decoded

    """
    for key, val in idict.items():
        if isinstance(val, (_np.ndarray, list)):
            idict[key] = [vv.decode() if hasattr(vv,"decode") else vv for vv in val]
        else:
            try:
                idict[key] = val.decode()
            except AttributeError:
                pass
    return idict


def CGdict2CG(filename, return_CPs=False, **cp_kwargs):
    r"""

    Parameters
    ----------
    filename
    cp_kwargs

    Returns
    -------

    """


    if isinstance(filename,str):
        a = _np.load(filename, allow_pickle=True)[()]
    elif isinstance(filename,dict):
        a = filename


    for sCP in a["serialized_CPs"]:
        try:
            decompress_serialized_CP(sCP,inplace=True)
        except AttributeError as e:
            pass
    contact_pairs = [sCP2CP(sCP, **cp_kwargs) for sCP in a["serialized_CPs"]]

    return mdciao.contacts.ContactGroup(contact_pairs, **{key:val for key,val in a.items() if key!="serialized_CPs"})

def sCP2CP(sCP,**cp_kwargs):
    r"""
    Turn a serialized ContactPair (dict) into an actual :obj:`~mdciao.contacts.ContactPair`-object

    Parameters
    ----------
    sCP : dict
    cp_kwargs : dict
        Will be passed to :obj:`~mdciao.contacts.ContactPair`

    Returns
    -------

    """
    cp_args = ["residues.idxs_pair", "time_traces.ctc_trajs", "time_traces.time_trajs"]
    _hd5f2CP_mapping = {
        'res_idxs_pair': 'residues.idxs_pair',
        'ctc_trajs': 'time_traces.ctc_trajs',
        'time_trajs': 'time_traces.time_trajs',
        'atom_pair_trajs': 'time_traces.atom_pair_trajs',
        'fragment_idxs': 'fragments.idxs',
        'fragment_names': 'fragments.names',
        'fragment_colors': 'fragments.colors',
        'anchor_residue_idx': 'residues.anchor_residue_index',
        'consensus_labels': 'residues.consensus_labels',
        'trajs' : "time_traces.trajs"
        }
    return mdciao.contacts.ContactPair(*[sCP[arg] for arg in cp_args],
                                       **{key: sCP[val] for key, val in _hd5f2CP_mapping.items() if
                                          _hd5f2CP_mapping[key] not in cp_args},
                                       **cp_kwargs)
def dict_of_CGs_2_hdf5(f, idict, compress=False,exclude=None, stride=1):
    if exclude is None:
        exclude=[]
    if compress:
        ref_t  = common_time_array_of_CG_dicts(idict)
        if ref_t is not None:
            ref_t = [t[::stride] for t in ref_t]
            exclude=["time_traces.time_trajs"]

    for key, cg in idict.items():
        try:
            g = f.create_group(str(key))
            iarch = cg.archive(exclude=exclude)
            #print(iarch.keys())
            # print(key, len(iarch["serialized_CPs"]))
            for ii, cp in enumerate(iarch["serialized_CPs"]):
                h = g.create_group(str(ii))
                for key2, val2 in cp.items():
                    if stride > 1 and key2 in ['time_traces.atom_pair_trajs',
                                               'time_traces.ctc_trajs',
                                               'time_traces.time_trajs']:

                        val2 = [tt[::stride] for tt in val2]

                    if compress and key2 =="time_traces.ctc_trajs":
                        val2 = [t for t in [" ".join([str(ff) for ff in _np.round(_np.array(tt)*1000).astype(_np.int32)]) for tt in val2]]
                    if compress and key2 =="time_traces.atom_pair_trajs":
                        val2 = [",".join(["%u %u"%tuple(pair) for pair in tt]) for tt in val2]
                    h.create_dataset(key2, data=val2)
            for key3 in ["interface_residxs","neighbors_excluded"]:
                g.create_dataset(key3, data=iarch[key3])
            g.create_dataset("name",data=str("name"))
        except ValueError as e:
            print("Neighborhood '%s' already exists in archive '%s'" % (str(key), f))
            print("Skipping")
            # raise(e)
    if "compress" in f.keys():
        assert f["compress"][()]==compress
    else:
        f.create_dataset(name="compress",data=compress)
    if "ref_t" in f.keys():
        _np.testing.assert_array_equal(_np.hstack(ref_t),_np.hstack(f["ref_t"]))
    else:
        if compress:
            for ii, itraj in enumerate(ref_t):
                key = "ref_t_%u"%ii
                if key not in f.keys():
                    f.create_dataset(name=key,data=itraj)
                else:
                    # we could assert with almost_equal but these should actually be EXACTLY the same
                    _np.testing.assert_array_equal(itraj, f[key],
                                                   err_msg="You cannot save ContactGroups with differing time-arrays"
                                                           "in the same hdf5, t")

def common_time_array_of_CG_dicts(CG_dict):
    ref_t = None
    for key, CG in CG_dict.items():
        for ii, cp in enumerate(CG._contacts):
            t = _np.hstack(cp.time_traces.time_trajs)
            if ref_t is None:
                ref_t = t
            try:
                _np.testing.assert_array_equal(t, ref_t)
            except:
                print("The CPs in the CGs of this CG_dict don't share their 'time_traces.time_trajs' attribute")
                return None
    return cp.time_traces.time_trajs
