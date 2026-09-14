# Discussion with the Diffsky team

## Project introduction

We use Diffsky LastJourney galaxies to build SNANA host libraries for DES
supernova simulations, comparing host distributions and SN-type-dependent
rate weighting against DES host data. We need consistent DES photometry,
stellar masses, star-formation rates, and galaxy profiles for host selection,
surface-brightness noise, and host association. Our current target extends to
z=1.8. We compute DECam grizY through your existing photometry engine and inject
the results into the SNANA HOSTLIB. Could this become native catalog photometry
alongside LSST and Roman?

## Priority questions

1. **Where should DECam integration live?** Is adding a transmission-curve
   collection to the existing photometry stage sufficient? What naming scheme
   and contribution format would you prefer? Can we reuse the same galaxy SED
   and random draws across LSST, Roman, and DES to preserve cross-band colors?
2. **Which DECam throughput definition should we use?** Are sedpy's `decam_*`
   curves appropriate for DES-SN5YR, including atmosphere and AB normalization?
   Do you have an authoritative filter set and validation reference?
3. **How do we reproduce native photometry exactly?** Which SSP/model release,
   dust prescription, emission-line treatment, redshift convention, lensing,
   and Milky Way extinction settings generated the existing LSST/Roman columns?
   Can you provide a small reference batch with expected magnitudes?
4. **How should we handle SSP scatter?** Our older environment could not read
   the compressed `delta_mag_ssp_scatter` column, so we set it to zero. What is
   the supported way to recover the correct values for resolved and synthetic
   cores? What biases in colors or selection could the zero-scatter approximation
   introduce? Is it valid to reuse scatter across filters and catalog releases?
5. **What is the completeness of this release?** Across 0.05 < z < 1.8, what
   stellar-mass limits and uncertainties apply to resolved and synthetic cores?
   Our current production cut is log10(M*/Msun) >= 8. Is that defensible for SN
   hosts, especially core-collapse hosts, or does it remove an important population?
   Earlier missing low-mass objects partly came from our reader excluding
   synthetic cores; how can we distinguish reader effects from true incompleteness?
6. **What is the stable row identity?** Is `gal_id` unique across patches,
   light-cone replicas, and releases? If not, which composite key should we use?
   Our join code can drop duplicate IDs; we need to establish whether those rows
   are distinct observable galaxies before accepting that behavior.

## SNANA-specific follow-ups

7. **Surface brightness and morphology:** Can you supply DECam disk, bulge,
   and knot fluxes as well as totals? Are component sizes half-light radii or
   profile scale lengths, in physical or comoving units? What are the angle and
   axis conventions? Is `bulge_to_total` a mass fraction or a band-specific light
   fraction? SNANA needs consistent flux fractions and Sersic profiles.
8. **Mass and SFR semantics:** What do `logsm_obs` and `logssfr_obs` represent,
   including IMF, averaging timescale, and scatter? How should quenched or
   effectively zero-SFR galaxies be represented? These choices affect our
   mass/SFR-based SN Ia and core-collapse rate weighting.
9. **Observed DES selection:** Are the native magnitudes total noiseless model
   magnitudes? How should we relate them to DES AUTO/aperture magnitudes and
   colors used by host-redshift efficiency curves? What external measurement
   and completeness model should remain in our simulation rather than the catalog?
10. **Redshift and velocity:** Which redshift should feed the photometry, and
    which should SNANA use as CMB-frame redshift? How are peculiar velocities
    and lensing already included, so we avoid applying them twice?
11. **Validation and supported software:** Can we agree on numerical tolerances,
    a redshift-grid convergence test, and a small maintained regression dataset
    covering low-mass, quenched, dusty, and synthetic-core galaxies? Which
    compatible package versions and catalog release should we target?

## Useful meeting outcome

Agree on the upstream insertion point, filter definitions, photometric
conventions, supported environment, stable IDs, and one reference batch.
Then validate DES photometry and quantify the host-selection impact before
regenerating the full HOSTLIB or interpreting SN-rate comparisons.
