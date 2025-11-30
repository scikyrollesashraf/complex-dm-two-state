# Complex-Field Two-State Dark Matter — Repository

**Title:** A Complex-Field Two-State Model for Dark Matter: Theory, Dynamics and Cosmological Constraints
**Author:** Kyrolles Ashraf Dawood
**Affiliation:** Department of Physics, Faculty of Science, Assiut University, Asyut, Egypt
**Contact:** [kyrls.dawd1231@science.aun.edu.eg](mailto:kyrls.dawd1231@science.aun.edu.eg)

## Overview

This repository contains the materials to reproduce the analytic derivations and numerical results for a model where a single complex scalar field $\Phi$ describes two effective sectors (real ↔ imaginary components) and their coherent conversions. The repository includes:

* `real_complex_DM_model.tex` — LaTeX manuscript (main text + appendices and derivations).
* `references.bib` — BibTeX reference file.
* `boltzmann_full_solver.py` — Scientific solver for the Boltzmann equation and plotting utilities.
* `LICENSE` — License for the code and data (choose and add appropriate text).
* `README.md` — This file.

## Requirements

* Python 3.8+
* TeX distribution (TeX Live, MiKTeX) with `pdflatex` and `bibtex` in PATH
* Python packages: `numpy`, `scipy`, `matplotlib`
  Install with:

```bash
pip install numpy scipy matplotlib
```

Optional (recommended): `pandas` for data handling, `seaborn` for enhanced plotting, and `joblib` or `multiprocessing` for parallel scans.

## Build the manuscript (locally)

From the repository root, compile the LaTeX source and bibliography:

```bash
pdflatex real_complex_DM_model.tex
bibtex real_complex_DM_model
pdflatex real_complex_DM_model.tex
pdflatex real_complex_DM_model.tex
```

Open `real_complex_DM_model.pdf` and verify figures, citations, and equations are correct.

## Run the Boltzmann solver (quick demo)

The script `boltzmann_full_solver.py` solves the Boltzmann equation for the comoving yield (Y(x)) and computes the relic abundance (\Omega_c h^2) using standard conversions. Run the demo:

```bash
python3 boltzmann_full_solver.py
```

This produces an `outputs/` folder containing:

* `Y_vs_x.csv` — computed yield vs x data,
* `Y_vs_x.pdf` — plot of yield vs x,
* `P_t.pdf` — example transition probability plot,
* `summary.txt` — run summary and numeric outputs.

## Using the solver programmatically

Import functions in a Jupyter notebook or another script:

```python
from boltzmann_full_solver import solve_Y, P_r_to_i_time, Omega_prefactor

# Example usage:
x, Y, params = solve_Y(m=10.0, sigma_v_cm3s=3e-26, Gamma_s_inv=0.0)
Y_inf = Y[-1]
omega_h2 = Omega_prefactor * 10.0 * Y_inf
```

## Input parameters and units

* `m` — dark matter mass in GeV.
* `sigma_v_cm3s` — thermally averaged annihilation cross section in cm³/s (the script converts to GeV⁻² internally). Typical thermal value: `3e-26 cm^3/s`.
* `Gamma_s_inv` — decay width (s⁻¹) if applicable (typically 0 for stable DM).
* `V` and `Delta` — mixing and splitting used for `P(t)` plots (in eV inside the plotting utility; conversion to natural units is done inside the code).
* `g_*` / `g_{*s}` — effective degrees of freedom: the script contains an approximate table for demonstration. For publication replace with high-resolution tabulations (e.g., from CLASS/PDG).


## Reproducibility & provenance

* The solver uses unit conversions and an approximate `g_*(T)` table. For publication-grade results, replace the approximate table with a validated `g_*(T)` and `g_{*s}(T)` dataset.
* For transparency, include a release of the code (GitHub release) and archive it on Zenodo to obtain a DOI. Add the DOI in the manuscript under "Data & Code Availability".

## License

Creative Commons Legal Code

Attribution 4.0 International (CC BY 4.0)

THE WORK (AS DEFINED BELOW) IS PROVIDED UNDER THE TERMS OF THIS CREATIVE
COMMONS PUBLIC LICENSE ("CCPL" OR "LICENSE"). THE WORK IS PROTECTED BY
COPYRIGHT AND/OR OTHER APPLICABLE LAW. ANY USE OF THE WORK OTHER THAN AS
AUTHORIZED UNDER THIS LICENSE OR COPYRIGHT LAW IS PROHIBITED.

BY EXERCISING ANY RIGHTS TO THE WORK PROVIDED HERE, YOU ACCEPT AND AGREE
TO BE BOUND BY THE TERMS OF THIS LICENSE. TO THE EXTENT THIS LICENSE MAY BE
A CONSUMER PROTECTION LAW THAT IS INCAPABLE OF BEING WAIVED BY AGREEMENT OF
THE PARTIES, NOTHING IN THIS LICENSE SHALL AFFECT YOUR RIGHTS UNDER THAT
LAW.

Section 1 — Definitions.

a. Adapted Material means material subject to Copyright and Similar Rights that
is derived from or based upon the Licensed Material and in which the Licensed
Material is translated, altered, arranged, transformed, or otherwise modified
in a manner requiring permission under the Copyright and Similar Rights held
by the Licensor. For purposes of this License, where the Licensed Material is
a musical work, performance, or sound recording, Adapted Material is always
inclusive of the sync rights necessary to use the Licensed Material in
combination with a moving image.

b. Adapter's License means the license You apply to Your Copyright and Similar
Rights in Your contributions to Adapted Material in accordance with the
requirements of this License.

c. BY-SA Compatible License means a license listed at
creativecommons.org that the Licensor has identified as essentially the same
as this License, including by allowing adaptations and requiring
derivative works to be licensed under the same terms.

d. Copyright and Similar Rights means copyright and/or similar rights closely
related to copyright including, without limitation, performance, publicity,
moral, neighboring and sui generis database rights, and rights to control
the collection, arrangement and presentation of data, but specifically
excluding Patents and Trademark rights.

e. Effective Technological Measures means those measures that, in the
absence of proper authority, may not be circumvented under laws fulfilling
obligations under Article 11 of the WIPO Copyright Treaty adopted on 20
December 1996, and similar international agreements.

f. Exceptions and Limitations means fair use, fair dealing or other
equivalents.

g. License Elements means the following high-level license attributes
as provided by the Licensor: attribution; license notice; disclaimer of
warranties and limitation of liability.

h. Licensor means the individual(s) or entity(ies) granting rights under this
License.

i. Licensed Material means the artistic or literary work, database, or other
material to which the Licensor applied this License.

j. Licensed Rights means the rights granted to You subject to the terms and
conditions of this License, which are limited to all Copyright and Similar
Rights that apply to Your use of the Licensed Material and that the Licensor
has authority to license.

k. Licensor's Designated Agent means the person or entity designated by the
Licensor to receive notices as required by Section 4(d).

l. Share means to provide material to the public by any means or process that
requires permission under the Licensed Rights, such as reproduction,
public display, public performance, distribution, dissemination, communication
or importation, and to make material available to the public including in
 ways that recipients may further make available or disseminate it.

m. Sui Generis Database Rights means rights other than copyright resulting
from Directive 96/9/EC of the European Parliament and of the Council of 11
March 1996 on the legal protection of databases, and similar national laws
that apply instead of, or in addition to, Directive 96/9/EC.

n. You means the individual or entity exercising the Licensed Rights under
this License. Your has a corresponding meaning.

Section 2 — Scope.

a. License grant.
Subject to the terms and conditions of this License, the Licensor hereby
grants You a worldwide, royalty-free, non-sublicensable, non-exclusive,
irrevocable license to exercise the Licensed Rights in the Licensed Material
to:

  i. reproduce and Share the Licensed Material, in whole or in part; and
  
  ii. produce, reproduce, and Share Adapted Material.

b. Exceptions and Limitations. For the avoidance of doubt, where Exceptions
and Limitations apply to Your use, this License does not apply and You do not
need to comply with its terms and conditions.

c. Term. The term of this License is specified in Section 6(a).

d. Media and formats; technical modifications allowed. The Licensor
authorizes You to exercise the Licensed Rights in all media and formats
whether now known or hereafter created, and to make technical modifications
necessary to do so. The Licensor waives and/or agrees not to assert any
right or authority to forbid You from making technical modifications necessary
to exercise the Licensed Rights, including technical modifications necessary
to circumvent Effective Technological Measures. For purposes of this
License, simply making modifications authorized by this Section 2(a) (such
as format-shifting) does not produce Adapted Material.

e. Downstream recipients. Those who receive Adapted Material from You
under the terms of this License will not have their Copyright and Similar
Rights in that Adapted Material revoked by virtue of your complying with
this License.

Section 3 — License Conditions.

Your exercise of the Licensed Rights is expressly made subject to the
following conditions.

a. Attribution. If You Share the Licensed Material (including in modified
form), You must:

  i. retain the following if it is supplied by the Licensor with the
  Licensed Material:

    A. identification of the creator(s) of the Licensed Material and
    any others designated to receive attribution, in any reasonable manner required
    by the Licensor; provided that, in all cases, attribution parties
    designated by the Licensor may be identified by the title or an
    identifier associated with the person(s) or entity(ies);

    B. a copyright notice;

    C. a notice that refers to this License;

    D. a notice that refers to the disclaimer of warranties;

    E. a URI or hyperlink to the Licensed Material to the extent
    reasonably practicable;

  ii. indicate if You modified the Licensed Material and retain an
  indication of any previous modifications; and

  iii. indicate the Licensed Material is licensed under this License, and
  include the text of, or the URI or hyperlink to, this License.

The above may be implemented in any reasonable manner; provided that
You do not imply the Licensor endorses You or Your use.

b. No additional restrictions. You may not apply legal terms or
technological measures that legally restrict others from doing anything the
License permits.

Section 4 — Sui Generis Database Rights.

Where the Licensed Rights include Sui Generis Database Rights that apply to
Your use of the Licensed Material, the Licensor waives and/or agrees not to
assert those rights to the limited extent necessary to allow You to exercise
the Licensed Rights as stated in Section 2.

Section 5 — Representations and Warranties.

Except for the limited warranty in Section 5(b), the Licensor offers the
Licensed Material as-is and makes no representations or warranties of any
kind concerning the Licensed Material.

a. To the extent possible, the Licensor asserts that to the best of their
knowledge they have the right to license the Licensed Material under this
License.

b. Limited warranty. Where the Licensor expressly states that the Licensed
Material is provided on terms that include a limited warranty, this warranty
will be stated in the Licensed Material or in the document provided with the
Licensed Material.

Section 6 — Term and Termination.

a. This License and the rights granted hereunder terminate automatically
upon any exercise of the Licensed Rights that is not expressly permitted
by this License and fails to meet any of its conditions.

b. Where You have received the Licensed Material under more than one
license, which includes this License and other terms, this License does not
terminate to the extent required by the other applicable license.

c. Subject to the above, the Licensor may expressly and specifically
place additional terms onto the Licensed Material notwithstanding the
terms of this License. In such cases the additional terms take precedence to
the extent they are compatible with this License.

Section 7 — Other Terms and Conditions.

a. The terms of this License are governed by the laws of the jurisdiction
where the Licensor resides, unless otherwise required by applicable law.

b. The Berne Convention for the Protection of Literary and Artistic Works,
the Universal Copyright Convention, and the World Intellectual Property
Organization Copyright Treaty explicitly form part of the terms of this
License.

c. No trademark license is granted under this License. This License does
not grant permission to use the trade names, trademarks, service marks,
or product names of the Licensor, except as required for reasonable and
customary use.

d. Where any part of this License is held to be invalid or unenforceable
under applicable law, that part shall be construed to the greatest extent
possible in a manner that results in the intentions of the Licensor being
fulfilled.

e. This License does not reduce, limit, or remove any of the rights
that others have under applicable law.

For further information, visit: https://creativecommons.org/licenses/by/4.0/legalcode


## Data & Code Availability (example text)

> The code and data used to generate results in this manuscript are available at `https://github.com/scikyrollesashraf/complex-dm-two-state` 

## Contact

Questions, issues, or requests for reproducibility assistance: Kyrolles Ashraf Dawood — [kyrls.dawd1231@science.aun.edu.eg](mailto:kyrls.dawd1231@science.aun.edu.eg)
