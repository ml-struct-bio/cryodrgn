import pickle
import sys

a = sys.argv[1]
b = sys.argv[2]
c = sys.argv[3]

with open(a, "rb") as f:
    a = pickle.load(f)
with open(b, "rb") as f:
    b = pickle.load(f)

if type(a) is tuple:
    diff_r = ((a[0] - b[0]) ** 2).sum()
    diff_t = ((a[1] - b[1]) ** 2).sum()
    assert diff_r < 1e-4
    assert diff_t < 1e-4
else:
    diff = ((a - b) ** 2).sum()
    print(f"{c} {1 - diff}")
    # assert diff < 1e-4
