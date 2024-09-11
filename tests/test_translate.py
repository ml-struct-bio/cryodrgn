import pytest
import os
import shutil
import argparse
import numpy as np
import torch
from cryodrgn import fft
from cryodrgn.source import ImageSource
from cryodrgn.lattice import Lattice
from cryodrgn.commands_utils import translate_mrcs


def test_shifted_image():
    torch.manual_seed(15321)
    imgs = ImageSource.from_file(os.path.join(pytest.DATADIR, "hand.mrcs")).images()
    img = imgs[0]
    D = img.shape[0]
    ht = fft.ht2_center(img)
    ht = fft.symmetrize_ht(ht)
    D += 1

    lattice = Lattice(D)
    ht = ht.view(1, -1)  # type: ignore
    trans = torch.tensor([5.0, 10.0]).view(1, 1, 2)
    ht_shifted = lattice.translate_ht(ht, trans)
    ht_np = ht_shifted.view(D, D)[0:-1, 0:-1]

    img_shifted = fft.ihtn_center(ht_np)
    new_arr = torch.Tensor(np.load(os.path.join(pytest.DATADIR, "im_shifted.npy")))
    assert torch.allclose(new_arr, img_shifted, atol=1e-4)


@pytest.mark.parametrize("trans", ["toy"], indirect=True)
class TestTranslateStack:
    out_lbl = "trans-particles.mrcs"

    def get_outdir(self, tmpdir_factory, particles, trans):
        particles_lbl = particles if isinstance(particles, str) else particles.label
        dirname = os.path.join("TranslateStack", particles_lbl, trans.label)
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)
        return odir

    @pytest.mark.parametrize(
        "particles", ["toy.mrcs", "toy.txt", "toy.star"], indirect=True
    )
    def test_default_translate(self, tmpdir_factory, particles, trans):
        outdir = self.get_outdir(tmpdir_factory, particles, trans)
        out_fl = os.path.join(outdir, self.out_lbl)
        args = [particles.path, trans.path, "-o", out_fl]

        parser = argparse.ArgumentParser()
        translate_mrcs.add_args(parser)
        translate_mrcs.main(parser.parse_args(args))
        assert os.path.exists(out_fl)

    def test_filetype_consistency(self, tmpdir_factory, trans):
        outdir_mrcs = self.get_outdir(tmpdir_factory, "toy.mrcs", trans)
        outdir_txt = self.get_outdir(tmpdir_factory, "toy.txt", trans)
        outdir_star = self.get_outdir(tmpdir_factory, "toy.star", trans)

        imgs_mrcs = ImageSource.from_file(
            os.path.join(outdir_mrcs, self.out_lbl)
        ).images()
        imgs_txt = ImageSource.from_file(
            os.path.join(outdir_txt, self.out_lbl)
        ).images()
        assert np.allclose(imgs_mrcs, imgs_txt)

        imgs_star = ImageSource.from_file(
            os.path.join(outdir_star, self.out_lbl)
        ).images()
        assert np.allclose(imgs_mrcs, imgs_star)

    @pytest.mark.parametrize("particles", ["toy.mrcs"], indirect=True)
    @pytest.mark.parametrize("tscale", [-1, 1, 0.5])
    def test_tscales(self, tmpdir_factory, particles, trans, tscale):
        outdir = self.get_outdir(tmpdir_factory, particles, trans)
        out_fl = os.path.join(
            outdir, self.out_lbl.replace(".mrcs", f"_tscale.{tscale}.mrcs")
        )
        args = [particles.path, trans.path, "-o", out_fl, "--tscale", str(tscale)]

        parser = argparse.ArgumentParser()
        translate_mrcs.add_args(parser)
        translate_mrcs.main(parser.parse_args(args))
        assert os.path.exists(out_fl)

    @pytest.mark.parametrize("particles", ["toy.mrcs"], indirect=True)
    def test_tscale_consistency(self, tmpdir_factory, particles, trans):
        outdir_mrcs = self.get_outdir(tmpdir_factory, "toy.mrcs", trans)
        out_fl_orig = os.path.join(outdir_mrcs, self.out_lbl)
        out_fl_tscale = os.path.join(
            outdir_mrcs, self.out_lbl.replace(".mrcs", "_tscale.1.mrcs")
        )
        imgs_orig = ImageSource.from_file(out_fl_orig).images()
        imgs_tscale = ImageSource.from_file(out_fl_tscale).images()
        assert np.allclose(imgs_orig, imgs_tscale)

    @pytest.mark.parametrize(
        "particles", ["toy.mrcs", "toy.txt", "toy.star"], indirect=True
    )
    def test_png_output(self, tmpdir_factory, particles, trans):
        outdir = self.get_outdir(tmpdir_factory, particles, trans)
        out_fl_png = os.path.join(outdir, self.out_lbl.replace(".mrcs", "_png.mrcs"))
        out_png = os.path.join(outdir, "translated.png")
        args = [particles.path, trans.path, "-o", out_fl_png, "--out-png", out_png]

        parser = argparse.ArgumentParser()
        translate_mrcs.add_args(parser)
        translate_mrcs.main(parser.parse_args(args))
        assert os.path.exists(out_fl_png)
        assert os.path.exists(out_png)

        shutil.rmtree(outdir)
