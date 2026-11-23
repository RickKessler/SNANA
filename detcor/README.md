# detcor_map: Roman detcor map in C, built on SNANA tools (Task 2)

This code is developed **outside SNANA**, but it uses SNANA's own tools and links
against SNANA's compiled object files. Nothing in SNANA has to change.

| File | Role |
|---|---|
| `detcor_map.h`, `detcor_map.c` | The function. It reads the handshaker map into memory and interpolates R(λ; t). It uses SNANA's `SNTABLE_*` table reader, `errmsg`, `quickBinSearch` and `interp_1DFUN`. |
| `test_detcor_map.c` | Test program: random thicknesses, filter dumps, the handshaker system-call check, and timers. Filters are read with SNANA's `rd2columnFile`. |
| `Makefile`, `snana_vars.mk` | Compile with `-I$SNANA_DIR/src` and link `$SNANA_DIR/obj/sntools.o` and `sntools_output.o`. The compiler and libraries are read from SNANA's configured `src/Makefile`. |
| `roman_ice_map_0-300_10nm.dat` | The Task 1 grid (0–300 nm in 10 nm steps). |

## Run (Midway)

```bash
export SNANA_DIR=/path/to/SNANA           # already built (obj/sntools.o exists)
cd handshaker/snana/detcor_map
make show                                 # check CC / libs picked up from SNANA
make check                                # C vs handshaker, 5 random thicknesses
make timing                               # timers + dumps (FILTER_DIR=... for Roman filters)
```

- **`make check`:** for each random thickness it makes a system call to `python3 -m handshaker.icecor --thickness T` and compares the result with the C value. It passes if the difference is ≤ 2e-6. The handshaker `src/` must be importable; the default `HS_CMD` uses `../../src`.
- **`make timing`:** times reading the map, `get_detcor_ratio` per thickness and `apply_detcor_to_filter` per filter. It also writes R and every `--filter` at 5 random thicknesses to `detcor_map_test_dump.dat*`.
- **On a Mac:** add `EXTRA_CFLAGS=-DMACOS`, and SNANA has to be built there first.

## API

```c
#include "detcor_map.h"
DETCOR_MAP_DEF map;
read_detcor_map("roman_ice_map_0-300_10nm.dat", &map);          // aborts via errmsg on error
double R = interp_detcor_ratio(&map, thick_nm, wave_A);          // one wavelength
get_detcor_ratio(&map, thick_nm, R_out);                         // all map wavelengths
apply_detcor_to_filter(&map, thick_nm, NLAM, lam, T_nom, T_detcor); // T_detcor = T_nom * R
free_detcor_map(&map);
```

R is linear in thickness (clamped to the map range) and in wavelength, and R = 1
outside 4800–23000 Å. This is the same convention as handshaker's `ReleasedIceModel`.

## Later: SNANA calling it

When SNANA starts calling these functions, its build has to compile this
`detcor_map.c` (from here, or a copy) into a `detcor_map.o` and link it. The code
itself already runs on SNANA's tools.
