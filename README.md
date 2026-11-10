# aiamplitudes_common_public

Shared utilities for the AI for Scattering Amplitudes project: Phi2/Phi3
symbols, front and back bases, relations, and d-ray polynomials.

For usage examples, see `notebook/Amplitudes_101.ipynb` in the companion
`AIAmplitudes_common_dev` repository.

## Data

Data files are stored in `~/.local/AIAmplitudesData/` by default. To download
the public dataset into an empty data directory:

```python
from aiamplitudes_common_public.download_data import download_all_public

download_all_public()
```

Some basis files require access to `AIAmplitudes/data_public`. Use
`download_all_private(username, token)` from the same module with a GitHub
token that has read access to that repository.

## Bases and cached counts

Select compatible basis files by name in the data directory or by absolute
path. Use `"R"` for restrictive bases, `"P"` for permissive bases, and `"PD"`
for permissive bases dual to the restrictive basis. `"PD"` requires the
companion `aiamplitudes_common_dev` package.

```python
import aiamplitudes_common_public as amp
from aiamplitudes_common_public.fbspaces import (
    F_number, B_number, NF_rels, NB_rels,
)

amp.configure_basis(front_rest="multiinitial_E", back_rest="phi2multifinal_E")
front = amp.fspace(3, "R")
back = amp.bspace(3, "R")

n_front, n_back = F_number[3], B_number[3]
n_front_relations, n_back_relations = NF_rels[3], NB_rels[3]
```

Dimensions and relation counts come from the selected data and are cached on
first load. Repeated queries reuse them without rereading or recounting.
Available weights depend on the files. Relation counts include zero relations
and count stored entries, rather than matrix rank.

- `fr_rels(w)` and `br_rels(w)` load the cached relation dictionaries.
- `basis_size(...)` and `relation_count(...)` support explicit queries with
  `prefix=` for another file or `optrace="phi3"` for the Phi3 back space.
- `configure_basis` also accepts `front_perm` and `back_perm`;
  `F_perm_number[w]` and `B_perm_number[w]` expose their dimensions.
- Changing the selection clears basis-dependent caches. After editing a file
  in place, call `amp.clear_basis_cache()`. Treat cached dictionaries as read-only.

Switching files does not transform existing symbols or sewing matrices.
Reload previously loaded basis and sewing dictionaries after a switch.

## Main modules

All modules are under `src/aiamplitudes_common_public/`.

| Module | Purpose |
| --- | --- |
| `__init__.py` | Public loaders for symbols, bases, relations, and polynomials |
| `download_data.py`, `file_readers.py` | Download and parse data files |
| `fbspaces.py` | Front/back spaces, lookups, and cached counts |
| `rels_utils.py`, `commonclasses.py` | Homogeneous relations from *Bootstrapping a Stress Tensor...*, relation evaluation, and wrappers |
| `polynomial_utils.py` | d-ray polynomial utilities discussed in *Recurrent Features of Amplitudes...* |
| `lazy_symbol.py`, `uncompressor.py` | Lazy symbol access and decompression |
