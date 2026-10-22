# Removed from the published repository

The executor copied `model/mismatch_model.pth` here as input provenance for the reference run.
It is an **upstream OligoFormer artifact**, not something produced by this conversion, so it is not
redistributed here.

- Source: `model/mismatch_model.pth` in [lulab/OligoFormer](https://github.com/lulab/OligoFormer) at commit `e2f53ad63387bbe166bf123949151e2bc9bf6ec3`
- sha256: `fb32c9a31d436113ff5271edd138cc86fe4065b7b1082e7c1cadb3563d627043`
- size: 77734 bytes

`deps/setup_oligoformer_runtime.sh` checks out that commit, so the runtime it builds contains this
exact file. Nothing in the test suite reads this fixture path; the tests use the runtime copy.
