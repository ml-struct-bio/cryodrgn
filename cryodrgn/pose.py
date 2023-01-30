import pickle
from typing import Optional, Tuple, Union, List
import logging
import numpy as np
import torch
import torch.nn as nn
from torch import Tensor
from cryodrgn import lie_tools, utils

logger = logging.getLogger(__name__)


class PoseTracker(nn.Module):
    def __init__(
        self,
        rots_np: np.ndarray,
        trans_np: Optional[np.ndarray] = None,
        D: Optional[int] = None,
        emb_type: Optional[str] = None,
        device: Optional[torch.device] = None,
    ):
        super(PoseTracker, self).__init__()
        rots = torch.tensor(rots_np.astype(np.float32), device=device)
        trans = (
            torch.tensor(trans_np.astype(np.float32), device=device)
            if trans_np is not None
            else None
        )
        self.rots = rots
        self.trans = trans
        self.use_trans = trans_np is not None
        self.D = D
        self.emb_type = emb_type
        if emb_type is None:
            pass
        else:
            if trans is not None:
                trans_emb = nn.Embedding(trans.shape[0], 2, sparse=True)
                trans_emb.weight.data.copy_(trans)
                self.trans_emb = trans_emb.to(device)
            else:
                self.trans_emb = None
            if emb_type == "s2s2":
                rots_emb = nn.Embedding(rots.shape[0], 6, sparse=True)
                rots_emb.weight.data.copy_(lie_tools.SO3_to_s2s2(rots))
            elif emb_type == "quat":
                rots_emb = nn.Embedding(rots.shape[0], 4, sparse=True)
                rots_emb.weight.data.copy_(lie_tools.SO3_to_quaternions(rots))
            else:
                raise RuntimeError("Embedding type {} not recognized".format(emb_type))
            self.rots_emb = rots_emb.to(device)

    @classmethod
    def load(
        cls,
        infile: Union[str, List[str]],
        Nimg: int,
        D: int,
        emb_type: Optional[str] = None,
        ind: Optional[np.ndarray] = None,
        device: Optional[torch.device] = None,
    ):
        """
        Return an instance of PoseTracker

        Inputs:
            infile (str or list):   One or two files, with format options of:
                                    single file with pose pickle
                                    two files with rot and trans pickle
                                    single file with rot pickle
            Nimg:               Number of particles
            D:                  Box size (pixels)
            emb_type:           SO(3) embedding type if refining poses
            ind:                Index array if poses are being filtered
        """
        # load pickle
        if type(infile) is str:
            infile = [infile]
        assert len(infile) in (1, 2)
        if len(infile) == 2:  # rotation pickle, translation pickle
            poses = (utils.load_pkl(infile[0]), utils.load_pkl(infile[1]))
        else:  # rotation pickle or poses pickle
            poses = utils.load_pkl(infile[0])
            if type(poses) != tuple:
                poses = (poses,)

        # rotations
        rots = poses[0]
        if ind is not None:
            if len(rots) > Nimg:  # HACK
                rots = rots[ind]
        assert rots.shape == (
            Nimg,
            3,
            3,
        ), f"Input rotations have shape {rots.shape} but expected ({Nimg},3,3)"

        # translations if they exist
        if len(poses) == 2:
            trans = poses[1]
            if ind is not None:
                if len(trans) > Nimg:  # HACK
                    trans = trans[ind]
            assert trans.shape == (
                Nimg,
                2,
            ), f"Input translations have shape {trans.shape} but expected ({Nimg},2)"
            assert np.all(
                trans <= 1
            ), "ERROR: Old pose format detected. Translations must be in units of fraction of box."
            trans *= D  # convert from fraction to pixels
        else:
            logger.warning("WARNING: No translations provided")
            trans = None

        return cls(rots, trans, D, emb_type, device=device)

    def save(self, out_pkl: str) -> None:
        if self.emb_type == "quat":
            r = lie_tools.quaternions_to_SO3(self.rots_emb.weight.data).cpu().numpy()
        elif self.emb_type == "s2s2":
            r = lie_tools.s2s2_to_SO3(self.rots_emb.weight.data).cpu().numpy()
        else:
            r = self.rots.cpu().numpy()

        if self.use_trans:
            if self.emb_type is None:
                assert self.trans is not None
                t = self.trans.cpu().numpy()
            else:
                assert self.trans_emb is not None
                t = self.trans_emb.weight.data.cpu().numpy()
            t /= self.D  # convert from pixels to extent
            poses = (r, t)
        else:
            poses = (r,)

        pickle.dump(poses, open(out_pkl, "wb"))

    def get_pose(self, ind: Union[int, Tensor]) -> Tuple[Tensor, Optional[Tensor]]:
        if self.emb_type is None:
            rot = self.rots[ind]
            tran = self.trans[ind] if self.trans is not None else None
        else:
            if self.emb_type == "s2s2":
                rot = lie_tools.s2s2_to_SO3(self.rots_emb(ind))
            elif self.emb_type == "quat":
                rot = lie_tools.quaternions_to_SO3(self.rots_emb(ind))
            else:
                raise RuntimeError  # should not reach here
            tran = self.trans_emb(ind) if self.trans_emb is not None else None
        return rot, tran
