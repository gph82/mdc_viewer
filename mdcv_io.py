import mdciao
import os.path as _path
import numpy as _np
from joblib import Parallel, delayed
from tqdm import tqdm_notebook as _tqdm
import numpy as _np
def load_data(
        files={
            "A": "data/neighborhoods.A.B2AR.zip.stride.10.hdf5",
            "B": "data/neighborhoods.B.B2AR.zip.stride.10.hdf5",
            "C": "data/neighborhoods.C.B2AR.zip.stride.10.hdf5"
        },
        verbose=0,
              **kwargs_hd5_2_archives):
   # _tqdm = lambda x : x


    Ns = Parallel(n_jobs=1, verbose=verbose)(delayed(hd5_2_dict_of_CGdicts)(ff,
                                                                        **kwargs_hd5_2_archives
                                                                        ) for ff in _tqdm(files.values(),
                                                                     desc="loading data"
                                                                     ))

    return {key:val for key, val in zip(files.keys(),Ns)}
def hd5_2_dict_of_CGdicts(obj, stride=1, restrict_to_residxs=None, decompress_here=True):
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
    import h5py
    if _path.exists(obj):
        data = h5py.File(obj,"r")
    else:
        data = obj

    if restrict_to_residxs is None:
        valid_res = lambda res : True
    else:
        valid_res = lambda res : res in restrict_to_residxs

    output_dict = {}
    if "compress" in data.keys() and data["compress"][()]:
        compress=True
        ref_t = data["ref_t"][()]
    for key, CG in data.items():
        # print(key)
        if key.isdigit() and valid_res(int(key)):
            pass

            archive = {}

            neighborhood_archs = list(
                [{key: val[()] for key, val in dict(val).items()} for key, val in CG.items() if key.isdigit()])

            neighborhood_archs = [decode_dict_values(idict) for idict in neighborhood_archs]
            archive["serialized_CPs"] = neighborhood_archs
            if compress:
                for iarch in neighborhood_archs:
                    iarch["time_traces.time_trajs"] = ref_t
                    if decompress_here:
                        decompress_serialized_CP(iarch)
            archive["name"] = [None if CG["name"][()].decode().lower()=="none" else CG["name"][()].decode()][0]
            archive["interface_residxs"]=CG["interface_residxs"][()]
            output_dict[int(key)] = archive

    return output_dict

from copy import deepcopy as _deepcopy
def decompress_serialized_CP(sCP,inplace=True):

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


def CGdict2CG(filename, **cp_kwargs):
    mapping = {
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

    if isinstance(filename,str):
        a = _np.load(filename, allow_pickle=True)[()]
    elif isinstance(filename,dict):
        a = filename

    cp_args = ["residues.idxs_pair", "time_traces.ctc_trajs", "time_traces.time_trajs"]

    for arch in a["serialized_CPs"]:
        try:
            decompress_arch(arch)
        except:
            pass
    contact_pairs = [mdciao.contacts.ContactPair(*[cp[arg] for arg in cp_args],
                                 **{key: cp[val] for key, val in mapping.items() if mapping[key] not in cp_args},
                                 **cp_kwargs,
                                 ) for cp in a["serialized_CPs"]]

    return mdciao.contacts.ContactGroup(contact_pairs, **{key:val for key,val in a.items() if key!="serialized_CPs"})

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
                    #print(key2)
                    h.create_dataset(key2, data=val2)
            g.create_dataset("interface_residxs", data=iarch["interface_residxs"])
            g.create_dataset("name", data=str(iarch["name"]))
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
            f.create_dataset(name="ref_t",data=ref_t)

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
