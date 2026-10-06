# Extended introduction — quantum chaos in brain data, for everyone

*No equations beyond one division. If you finished high school, you can read
this — and then actually use the code.*

## The big idea in one paragraph

Imagine a very long street where cars park one after another, and we only care
about the **gaps between neighbouring cars**. If drivers park at random, most
gaps are small and a few are huge — pure luck. But if every driver politely
keeps at least a car-length of space, tiny gaps become rare and the gaps look
evenly spaced. Physicists discovered that the *energy levels* of complicated
systems behave exactly like these parked cars: **random parking** (called
"Poisson statistics") when the system is simple and orderly, **polite parking**
(called "level repulsion", the signature of "quantum chaos") when the system is
strongly interconnected. This project asks a strange and fun question: when
your brain switches from daydreaming (rest) to solving a task, do the
"energy levels" of its activity-correlation matrix repark — from random-ish to
repelled-ish? And is your personal parking style as unique as a fingerprint?

## What is the brain data, exactly?

We use **fMRI** (functional magnetic resonance imaging). It does not measure
neurons directly; it measures blood oxygen, which rises where neurons were
recently active (this blood signal is called **BOLD** — blood-oxygen-level-
dependent contrast). A scanner takes a whole-brain picture roughly every
second while a person lies inside, sometimes resting, sometimes doing tasks
(remembering images, gambling games, moving fingers...).

The brain is divided into a few hundred regions (a **parcellation** — think of
it as a map's postal districts; we use the popular Glasser atlas, see the links
below). For each region we get one wiggly time series: activity up, activity
down, for several minutes.

## Step 1 — From wiggles to a matrix

For every pair of regions we compute the **correlation** of their time series:
+1 means "perfectly in sync", −1 "perfectly opposite", 0 "unrelated". With
200 regions this gives a 200×200 table of numbers, symmetric (region A vs B =
region B vs A). This table is called the **functional connectivity (FC)
matrix**. It is the single most studied object in human brain imaging.

## Step 2 — Eigenvalues: the matrix's "parking spots"

Every symmetric matrix hides a special set of numbers called **eigenvalues**
(German "eigen" = "its own"). You can think of them as the loudness of the
matrix's hidden vibration modes: push the matrix and it "rings" in a few
characteristic patterns, each with a strength — an eigenvalue. A 200×200 FC
matrix has 200 eigenvalues. Sorted from small to large and laid on a number
line, they are our "parked cars".

Why care? Because 70 years of physics (Wigner, Dyson, Mehta's textbook — see
links) proved that the **spacing pattern** of such eigenvalues tells you how
the underlying system is organised:

- **Random parking (Poisson)**: the system behaves like many *independent*
  parts. Gaps between neighbours are usually small, occasionally huge.
- **Polite parking (level repulsion, "GOE")**: the parts are strongly
  *coupled*; eigenvalues push each other apart, tiny gaps almost never happen.

GOE stands for **Gaussian Orthogonal Ensemble** — a fancy name for "fill a
symmetric matrix with random numbers from a bell curve". It is the standard
model of a strongly-coupled system, and its parking style is the anchor at one
end of our scale.

## Step 3 — One number to locate the parking style

Comparing whole histograms of gaps is fiddly, because every spectrum has a
different overall scale (like comparing parking in metres vs feet). Physicists
Oganesyan and Huse (2007) invented a trick that needs **no rescaling at all**:

1. Sort the eigenvalues, take each gap between neighbours: s1, s2, s3, ...
2. For each pair of *consecutive* gaps, divide the smaller by the larger:

   **r = min(s_i, s_{i+1}) / max(s_i, s_{i+1})**

   That is it — one division. r is always between 0 and 1.
3. Average all the r's. Call the average **⟨r⟩**.

The magic: the average lands on two known universal values,

- random parking (Poisson): **⟨r⟩ = 2·ln(2) − 1 ≈ 0.386**,
- polite parking (GOE): **⟨r⟩ ≈ 0.536**,

and anything in between means "partially coupled". Our code computes ⟨r⟩ in
one line. Try it:

```python
from spectral_chaos import poisson_spectrum, goe_matrix, sorted_eigenvalues, mean_spacing_ratio

print(mean_spacing_ratio(poisson_spectrum(2000)))                # ≈ 0.386
print(mean_spacing_ratio(sorted_eigenvalues(goe_matrix(2000))))  # ≈ 0.536
```

## The hypothesis (in plain words)

When you engage in a task, your brain regions coordinate more strongly —
the parts stop behaving independently. If that is true, the FC eigenvalues
should **repark toward the polite end**: ⟨r⟩ should rise above its resting
value. We will also check whether your ⟨r⟩ is stable across scan days (a
**fingerprint**) and whether people whose spectra shift more during tasks also
perform better. If none of that survives careful comparison with shuffled
"null" data (recordings with the same sound/spectrum but destroyed
coordination — see *phase randomisation* in the table), the idea is wrong, and
we will say so.

## Every keyword we use, explained

| Term | Plain meaning | Why we need it | Learn more |
|---|---|---|---|
| fMRI | Scanner measuring blood oxygen as a proxy for brain activity | Our raw data | https://en.wikipedia.org/wiki/Functional_magnetic_resonance_imaging |
| BOLD | The specific blood-oxygen signal fMRI records | What each region's time series contains | https://en.wikipedia.org/wiki/Blood-oxygen-level-dependent_imaging |
| Parcellation | Dividing the cortex into ~100–400 regions | Turns images into per-region time series | https://www.nature.com/articles/nature18933 (Glasser atlas) |
| Correlation | How in-sync two series are (−1 to +1) | Builds the FC matrix | https://en.wikipedia.org/wiki/Pearson_correlation_coefficient |
| Functional connectivity (FC) | The matrix of all pairwise correlations | The object we analyse | https://en.wikipedia.org/wiki/Functional_connectivity |
| Eigenvalue | Strength of a matrix's hidden vibration modes | Our "parked cars" | https://www.khanacademy.org/math/linear-algebra/alternate-bases/eigen-everything/v/linear-algebra-introduction-to-eigenvalues-and-eigenvectors |
| Random matrix theory (RMT) | Statistics of eigenvalues of random matrices | Gives the two anchor values | https://en.wikipedia.org/wiki/Random_matrix |
| GOE | "Bell-curve-filled symmetric matrix" ensemble; polite-parking anchor | One end of our scale | https://en.wikipedia.org/wiki/Random_matrix#Gaussian_ensembles |
| Poisson statistics | Random-parking spacing pattern (⟨r⟩ ≈ 0.386) | The other end | https://en.wikipedia.org/wiki/Poisson_distribution |
| Level repulsion | Eigenvalues avoid crowding; tiny gaps are rare | The chaos signature | https://en.wikipedia.org/wiki/Random_matrix#Level_repulsion |
| Spacing ratio ⟨r⟩ | Smaller-gap ÷ larger-gap, averaged | Our main coordinate | https://journals.aps.org/prb/abstract/10.1103/PhysRevB.75.155111 |
| Unfolding | Rescaling a spectrum to even density | Needed for some extra checks (not for ⟨r⟩) | https://en.wikipedia.org/wiki/Random_matrix#Spectral_rigidity |
| Spectral form factor | Fourier-space check of long-range spacing order | A second opinion | https://en.wikipedia.org/wiki/Spectral_form_factor |
| Wigner surmise | A one-line formula for the gap histogram | A third opinion (fits "repulsion strength" β) | https://en.wikipedia.org/wiki/Random_matrix#Wigner_surmise |
| Phase randomisation | Shuffle each series' phases: same wiggles, no sync | The null test protecting us from fooling ourselves | https://en.wikipedia.org/wiki/Surrogate_data_testing |
| HCP Young Adult | Public dataset: ~1 200 people, rest + 7 tasks | Where we test everything | https://www.humanconnectome.org/study/hcp-young-adult |
| Fingerprinting | Showing a person's brain data is re-identifiable | Motivates H2 | https://www.nature.com/articles/nn.4135 |
| ICC | Statistic for "same person, same number across days" | Quantifies the fingerprint | https://en.wikipedia.org/wiki/Intraclass_correlation |

## A tiny worked example

Five sorted eigenvalues: 1.0, 1.4, 2.0, 2.1, 3.0.

- Gaps: 0.4, 0.6, 0.1, 0.9.
- Consecutive gap pairs: (0.4, 0.6), (0.6, 0.1), (0.1, 0.9).
- Ratios: 0.4/0.6 = 0.667; 0.1/0.6 = 0.167; 0.1/0.9 = 0.111.
- ⟨r⟩ = (0.667 + 0.167 + 0.111)/3 ≈ 0.315 — closer to random parking
  (0.386) than to polite parking (0.536), as expected for such a lumpy,
  mostly-random little spectrum.

Real analyses use hundreds of eigenvalues and bootstrap error bars, but the
arithmetic never gets harder than this.

## What could make us wrong?

- The FC matrix's eigenvalues repel *even for pure noise* (we verified this in
  our tests — random correlation matrices are already GOE-like in the middle
  of their spectrum!). So the honest question is **"does ⟨r⟩ MOVE between rest
  and task?"**, not "is there any repulsion?".
- Task scans are shorter — fewer time points — and fewer points alone can
  shift spacing statistics. Our null models and effective-sample-size checks
  exist exactly for this trap.
- Blood flow, head motion, and heartbeat also create correlations. Standard
  cleaning pipelines (from the HCP) are applied first, and motion is included
  as a confound in every model.

## Where the pieces of this repository fit

- `src/spectral_chaos/spectral.py` — every statistic above, one function each.
- `tests/` — proves the code reproduces the two magic numbers (0.386 / 0.536)
  and the other anchors.
- `docs/METHODS.md` — the full pipeline on real HCP data.
- GitHub **Issues** — the project split into atomic tasks you can pick up.

## Further friendly resources

- 3Blue1Brown, *Eigenvectors and eigenvalues* (video):
  https://www.youtube.com/watch?v=PFDu9oVAE-g
- Khan Academy, correlation: https://www.khanacademy.org/math/statistics-probability/describing-relationships-quantitative-data
- HCP data access: https://db.humanconnectome.org/
- nilearn (Python fMRI library we build on): https://nilearn.github.io/
- Mehta's classic (advanced): *Random Matrices* — free previews at
  https://books.google.com/books?id=B7QgAQAAIAAJ
