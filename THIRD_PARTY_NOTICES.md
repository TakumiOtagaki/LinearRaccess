# Third-party notices

The MIT license in `LICENSE` applies to original LinearRaccess code. It does
not replace the terms below.

## ViennaRNA parameter data

`src/energy_param.hpp` and `src/intloops.hpp` contain parameter data extracted
from ViennaRNA Package 2.5.1. The package's redistribution and attribution
terms are reproduced in `third_party/ViennaRNA-2.5.1-COPYING.txt`. Credit is
due to the ViennaRNA authors and the Institute for Theoretical Chemistry of
the University of Vienna. These files must not be redistributed for a fee
other than media costs under those terms.

Source: https://github.com/ViennaRNA/ViennaRNA/tree/v2.5.1

## CapR parameter data

`src/legacy_energy_param.hpp` and `src/legacy_intloops.hpp` were extracted
from CapR. CapR is distributed under the MIT license reproduced in
`third_party/CapR-LICENSE.txt` and acknowledges use of ViennaRNA Package 1.8.5
source. The conservative ViennaRNA notice above is retained for the bundled
thermodynamic data as well.

Source: https://github.com/fukunagatsu/CapR

## CONTRAfold numerical routine

The segmented polynomial in `include/linearraccess/miscs.hpp` used by the fast
log-sum-exp path was borrowed from CONTRAfold. Its BSD-style license is
reproduced in `third_party/CONTRAfold-LICENSE.txt`.

Source: https://github.com/csfoo/contrafold-se

## Optional Raccess adapter

`src/energy_raccess.cpp` and `src/energy_raccess.hpp` are original adapter code
and do not contain the Raccess implementation. Raccess is not bundled,
downloaded, or enabled by default because the version audited during
development prohibits redistribution. Enabling the adapter requires a
separately obtained include tree and does not grant any rights to Raccess.

This notice records the engineering release boundary and is not legal advice.
