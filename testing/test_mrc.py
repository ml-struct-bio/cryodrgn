import numpy as np
from cryodrgn.source import ImageSource

data = ImageSource.from_file("data/toy_projections.mrcs", lazy=True)
data2 = ImageSource.from_file("data/toy_projections.mrcs", lazy=False).images()
data2 = np.array(data2)
data1 = np.array(data[:])
assert (data1 == data2).all()
print("ok")

data2 = np.array(ImageSource.from_file("data/toy_projections.star").images())
assert (data1 == data2).all()
print("ok")

data2 = np.array(ImageSource.from_file("data/toy_projections.txt").images())
assert (data1 == data2).all()
print("ok")

print("all ok")
