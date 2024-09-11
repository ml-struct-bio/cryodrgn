import pytest
import os
import argparse
import torch
from cryodrgn.source import ImageSource
from cryodrgn.commands_utils import add_psize


@pytest.mark.parametrize("volume", ["toy", "hand", "spike", "empiar"], indirect=True)
@pytest.mark.parametrize("Apix", ["1", "1.7"])
def test_add_psize(tmpdir, volume, Apix):
    out_vol = os.path.join(tmpdir, "toy_projections_added_psize.mrc")
    parser = argparse.ArgumentParser()
    add_psize.add_args(parser)
    args = parser.parse_args([volume.path, "-o", out_vol])
    add_psize.main(args)

    # data should be unchanged
    old_data = ImageSource.from_file(volume.path).images()
    new_data = ImageSource.from_file(out_vol).images()
    assert torch.allclose(new_data, old_data)

    new_vol = str(out_vol).replace(".mrc", "_new.mrc")
    args = parser.parse_args([out_vol, "-o", new_vol, "--Apix", Apix])
    add_psize.main(args)

    # data should be unchanged
    old_data = ImageSource.from_file(out_vol).images()
    new_data = ImageSource.from_file(new_vol).images()
    assert torch.allclose(new_data, old_data)
