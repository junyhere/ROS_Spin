import csv
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from spin_chemistry import (
    ACTIVATABLE_STATUSES,
    EVIDENCE_BACKED,
    REACTION_SQ_O2_ET,
    SENSITIVITY,
    USE_CONDITION_MATCHED,
    USE_LITERATURE_ARITHMETIC,
    EncounterParameters,
    I6,
    ModelPolicyError,
    O,
    P_DOUBLET,
    P_QUARTET,
    R,
    R_DOT_O,
    acid_base_fractions,
    associate_encounters,
    bulk_superoxide_formation_rate,
    clebsch_gordan_states,
    downstream_ros_species_resolved,
    evidence_backed_encounter_prediction,
    evaluate_bulk_superoxide_rate,
    execute_three_qubit_unitary_circuit,
    form_semiquinone,
    h2o2_loss_rate,
    hamiltonian,
    independent_circuit_capability,
    initial_density,
    load_execution_context,
    primary_superoxide_formation,
    propagate_encounter_independent,
    propagate_encounter_reference,
    sod_dismutation_loss_rate,
    spin_matrices,
    spontaneous_dismutation_rate,
    unitary_embedding_consistency_check,
    validate_authority_bundle,
)


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs" / "doxorubicin_parameters.json"
PROVENANCE = ROOT / "configs" / "parameter_provenance.csv"


class SpinAlgebraTests(unittest.TestCase):
    def test_spin_commutators_and_casimirs(self):
        for spin in (0.5, 1.0):
            sx, sy, sz = spin_matrices(spin)
            self.assertTrue(np.allclose(sx @ sy - sy @ sx, 1j * sz))
            self.assertTrue(
                np.allclose(sx @ sx + sy @ sy + sz @ sz,
                            spin * (spin + 1) * np.eye(sx.shape[0]))
            )

    def test_projector_algebra_and_dimensions(self):
        self.assertTrue(np.allclose(P_DOUBLET @ P_DOUBLET, P_DOUBLET))
        self.assertTrue(np.allclose(P_QUARTET @ P_QUARTET, P_QUARTET))
        self.assertTrue(np.allclose(P_DOUBLET @ P_QUARTET, 0))
        self.assertTrue(np.allclose(P_DOUBLET + P_QUARTET, I6))
        self.assertAlmostEqual(float(np.trace(P_DOUBLET).real), 2)
        self.assertAlmostEqual(float(np.trace(P_QUARTET).real), 4)

    def test_clebsch_gordan_states_resolve_projectors(self):
        states = clebsch_gordan_states()
        gram = np.array([[np.vdot(a, b) for b in states.values()] for a in states.values()])
        self.assertTrue(np.allclose(gram, np.eye(6)))
        doublet = sum(np.outer(state, state.conj()) for name, state in states.items()
                      if name.startswith("D"))
        quartet = sum(np.outer(state, state.conj()) for name, state in states.items()
                      if name.startswith("Q"))
        self.assertTrue(np.allclose(doublet, P_DOUBLET))
        self.assertTrue(np.allclose(quartet, P_QUARTET))

    def test_r_dot_o_analytic_eigenvalues(self):
        values = np.linalg.eigvalsh(R_DOT_O)
        self.assertTrue(np.allclose(values, [-1, -1, 0.5, 0.5, 0.5, 0.5]))


class InitialStateAndHamiltonianTests(unittest.TestCase):
    def test_initial_states_are_hermitian_positive_and_normalized(self):
        for scenario in ("unpolarized", "doublet", "quartet"):
            density = initial_density(scenario)
            self.assertTrue(np.allclose(density, density.conj().T), scenario)
            self.assertGreaterEqual(float(np.linalg.eigvalsh(density).min()), -1e-14)
            self.assertAlmostEqual(float(np.trace(density).real), 1.0)
        self.assertAlmostEqual(
            float(np.trace(P_DOUBLET @ initial_density("unpolarized")).real), 1 / 3
        )

    def test_variable_manifold_mixture_is_exposed_and_bounded(self):
        for p_doublet in (0.0, 0.2, 1 / 3, 0.8, 1.0):
            density = initial_density("mixture", p_doublet)
            expected = p_doublet * P_DOUBLET / 2 + (1 - p_doublet) * P_QUARTET / 4
            self.assertTrue(np.allclose(density, expected))
            self.assertAlmostEqual(float(np.trace(P_DOUBLET @ density).real), p_doublet)
        for invalid in (-0.1, 1.1, np.nan, np.inf):
            with self.assertRaises(ValueError):
                initial_density("mixture", invalid)
        with self.assertRaisesRegex(ValueError, "requires"):
            initial_density("mixture")

    def test_manifold_states_are_mixtures_not_pure_wavefunctions(self):
        doublet = initial_density("doublet")
        quartet = initial_density("quartet")
        self.assertAlmostEqual(float(np.trace(doublet @ doublet).real), 0.5)
        self.assertAlmostEqual(float(np.trace(quartet @ quartet).real), 0.25)

    def test_zero_and_exchange_hamiltonian_limits(self):
        self.assertTrue(np.allclose(hamiltonian(EncounterParameters()), 0))
        exchange = 3.25e6
        values = np.linalg.eigvalsh(
            hamiltonian(EncounterParameters(exchange_rad_s=exchange))
        )
        self.assertTrue(
            np.allclose(values, exchange * np.array([-1, -1, 0.5, 0.5, 0.5, 0.5]))
        )

    def test_common_zeeman_is_total_spin_generator(self):
        g_value, field = 2.0, (0.0, 0.0, 1e-4)
        parameters = EncounterParameters(field_t=field, g_radical=g_value, g_oxygen=g_value)
        expected = 8.79410005e10 * field[2] * g_value * (R[2] + O[2])
        self.assertTrue(np.allclose(hamiltonian(parameters), expected))
        self.assertTrue(np.allclose(hamiltonian(parameters) @ P_DOUBLET,
                                    P_DOUBLET @ hamiltonian(parameters)))


class EncounterDynamicsTests(unittest.TestCase):
    def test_unpolarized_state_is_invariant_under_arbitrary_unitary(self):
        parameters = EncounterParameters(
            field_t=(1e-4, -2e-4, 3e-4), g_radical=2.1, g_oxygen=1.8,
            exchange_rad_s=3e6, dipolar_rad_s=2e6, dipolar_axis=(1, 2, 3),
            oxygen_zfs_d_rad_s=1.1e6, oxygen_zfs_e_rad_s=-0.7e6,
            local_field_proxy_rad_s=(0.4e6, -0.8e6, 1.2e6),
        )
        output = propagate_encounter_reference(parameters, 2e-7, samples=31)
        self.assertLess(np.max(np.abs(output["density_matrices"] - I6 / 6)), 2e-12)

    def test_selective_loss_then_mixing_changes_unpolarized_reactive_yield(self):
        common = dict(k_doublet_s=3e6, k_quartet_s=1e5, k_escape_s=8e5)
        without_mixing = propagate_encounter_reference(
            EncounterParameters(**common), 3e-6, "unpolarized", 101
        )
        with_mixing = propagate_encounter_reference(
            EncounterParameters(
                **common, local_field_proxy_rad_s=(2e6, 0.0, 0.0)
            ),
            3e-6,
            "unpolarized",
            101,
        )
        self.assertGreater(
            abs(
                without_mixing["primary_superoxide_yield"]
                - with_mixing["primary_superoxide_yield"]
            ),
            1e-4,
        )

    def test_exchange_and_common_zeeman_do_not_mix_manifolds(self):
        cases = (
            EncounterParameters(exchange_rad_s=8e6),
            EncounterParameters(field_t=(2e-4, -1e-4, 3e-4),
                                g_radical=2.0, g_oxygen=2.0),
        )
        for parameters in cases:
            output = propagate_encounter_reference(
                parameters, 2e-7, "doublet", samples=101
            )
            self.assertLess(np.max(np.abs(output["p_quartet"])), 1e-10)

    def test_named_local_field_structure_can_mix_manifolds(self):
        parameters = EncounterParameters(local_field_proxy_rad_s=(4e6, 0.0, 0.0))
        output = propagate_encounter_reference(parameters, 1e-6, "doublet", 201)
        self.assertGreater(np.max(output["p_quartet"]), 0.1)

    def test_positive_trace_decreasing_and_yield_accounting(self):
        parameters = EncounterParameters(
            k_doublet_s=2e6, k_quartet_s=2e5, k_escape_s=1e6,
            local_field_proxy_rad_s=(3e5, 0.0, 0.0),
            radical_relaxation_s=4e5, oxygen_relaxation_s=2e5,
        )
        output = propagate_encounter_reference(parameters, 8e-6, "doublet", 801)
        minima = [
            np.linalg.eigvalsh((density + density.conj().T) / 2).min()
            for density in output["density_matrices"]
        ]
        self.assertGreaterEqual(float(np.min(minima)), -2e-11)
        self.assertTrue(np.all(np.diff(output["survival"]) <= 2e-12))
        accounted = (
            output["primary_superoxide_yield"] + output["escape_yield"]
            + output["unresolved_probability"]
        )
        self.assertAlmostEqual(accounted, 1.0, places=5)

    def test_relaxation_only_generator_preserves_trace_and_positivity(self):
        parameters = EncounterParameters(
            radical_relaxation_s=7e5,
            oxygen_relaxation_s=4e5,
        )
        for initial_state in ("unpolarized", "doublet", "quartet"):
            output = propagate_encounter_reference(
                parameters, 4e-6, initial_state, samples=81
            )
            traces = np.trace(
                output["density_matrices"], axis1=1, axis2=2
            ).real
            minima = [
                np.linalg.eigvalsh((density + density.conj().T) / 2).min()
                for density in output["density_matrices"]
            ]
            self.assertLess(np.max(np.abs(traces - 1.0)), 2e-12)
            self.assertGreaterEqual(float(np.min(minima)), -2e-12)
            self.assertAlmostEqual(output["probability_balance"], 1.0, places=12)

    def test_output_sampling_density_roundoff_only(self):
        # The reference step exponential is exact for the constant generator;
        # this checks output sampling/roundoff, not independent convergence.
        parameters = EncounterParameters(
            k_doublet_s=2e6, k_quartet_s=2e5, k_escape_s=1e6
        )
        coarse = propagate_encounter_reference(parameters, 5e-6, samples=501)
        fine = propagate_encounter_reference(parameters, 5e-6, samples=1001)
        for key in ("primary_superoxide_yield", "escape_yield"):
            self.assertLess(abs(coarse[key] - fine[key]), 2e-5, key)

    def test_analytic_reaction_and_escape_limits(self):
        kd, escape, duration = 2e6, 1e6, 2e-6
        output = propagate_encounter_reference(
            EncounterParameters(k_doublet_s=kd, k_escape_s=escape),
            duration, "doublet", 1001,
        )
        consumed = 1 - np.exp(-(kd + escape) * duration)
        self.assertAlmostEqual(output["doublet_reaction_yield"],
                               kd / (kd + escape) * consumed, places=6)
        self.assertAlmostEqual(output["escape_yield"],
                               escape / (kd + escape) * consumed, places=6)
        escape_only = propagate_encounter_reference(
            EncounterParameters(k_escape_s=escape), duration, samples=1001
        )
        self.assertAlmostEqual(escape_only["escape_yield"],
                               1 - np.exp(-escape * duration), places=6)
        self.assertEqual(escape_only["primary_superoxide_yield"], 0.0)

    def test_public_encounter_boundaries_reject_nonfinite_shapes_and_counts(self):
        with self.assertRaises(ValueError):
            propagate_encounter_reference(EncounterParameters(field_t=(0, 0)), 1e-6)
        with self.assertRaises(ValueError):
            propagate_encounter_reference(
                EncounterParameters(k_escape_s=np.nan), 1e-6
            )
        for samples in (1, 2.5, np.nan, True):
            with self.assertRaises(ValueError):
                propagate_encounter_reference(EncounterParameters(), 1e-6, samples=samples)


class IndependentEncounterValidationTests(unittest.TestCase):
    def test_full_trajectory_and_yields_match_reference(self):
        cases = (
            EncounterParameters(
                field_t=(1e-5, -2e-5, 3e-5), g_radical=2.0035, g_oxygen=2.0,
                exchange_rad_s=3e5, dipolar_rad_s=2e5,
                local_field_proxy_rad_s=(8e5, -2e5, 1e5),
                k_doublet_s=2e6, k_quartet_s=2e5, k_escape_s=7e5,
            ),
            EncounterParameters(
                oxygen_zfs_d_rad_s=5e5, oxygen_zfs_e_rad_s=2e5,
                local_field_proxy_rad_s=(9e5, 0, 0),
                radical_relaxation_s=4e5, oxygen_relaxation_s=3e5,
                k_doublet_s=1.5e6, k_quartet_s=4e5, k_escape_s=6e5,
            ),
        )
        for parameters in cases:
            reference = propagate_encounter_reference(
                parameters, 2e-6, "mixture", 41, 0.4
            )
            independent = propagate_encounter_independent(
                parameters, 2e-6, "mixture", 41, 0.4,
                atol=1e-12, rtol=1e-10,
            )
            self.assertLess(
                np.max(abs(reference["density_matrices"] - independent["density_matrices"])),
                2e-9,
            )
            for key in ("p_doublet", "p_quartet", "survival"):
                self.assertLess(np.max(abs(reference[key] - independent[key])), 3e-9)
            for key in (
                "doublet_reaction_yield", "quartet_reaction_yield",
                "primary_superoxide_yield", "escape_yield", "unresolved_probability",
            ):
                self.assertLess(abs(reference[key] - independent[key]), 3e-9, key)
            accounted = (
                independent["primary_superoxide_yield"]
                + independent["escape_yield"]
                + independent["unresolved_probability"]
            )
            self.assertAlmostEqual(accounted, 1.0, places=8)
            minimum = min(
                np.linalg.eigvalsh((density + density.conj().T) / 2).min()
                for density in independent["density_matrices"]
            )
            self.assertGreaterEqual(float(minimum), -2e-9)

    def test_independent_solver_analytic_loss_limit(self):
        rate, escape, duration = 2e6, 1e6, 2e-6
        output = propagate_encounter_independent(
            EncounterParameters(k_doublet_s=rate, k_escape_s=escape),
            duration, "doublet", 21,
        )
        consumed = 1 - np.exp(-(rate + escape) * duration)
        self.assertAlmostEqual(
            output["doublet_reaction_yield"], rate / (rate + escape) * consumed,
            places=8,
        )
        self.assertAlmostEqual(
            output["escape_yield"], escape / (rate + escape) * consumed, places=8
        )


class StagedChemistryTests(unittest.TestCase):
    def test_upstream_stages_are_classical_and_separate(self):
        semiquinone = form_semiquinone(10e-6, 0.4)
        encounters = associate_encounters(semiquinone, 2e-6, 0.5)
        primary = primary_superoxide_formation(encounters, 0.25)
        self.assertAlmostEqual(semiquinone, 4e-6)
        self.assertAlmostEqual(encounters, 1e-6)
        self.assertAlmostEqual(primary, 0.25e-6)

    def test_bulk_bimolecular_rate_remains_dimensionally_separate(self):
        self.assertAlmostEqual(
            bulk_superoxide_formation_rate(1e-6, 2e-4, 3.5e8), 0.07
        )
        parameters = EncounterParameters(k_doublet_s=3.5e8)
        self.assertEqual(parameters.k_doublet_s, 3.5e8)
        self.assertNotEqual("M^-1 s^-1", "s^-1")

    def test_acid_base_and_spontaneous_rate(self):
        fraction_ho2, fraction_o2minus = acid_base_fractions(4.88, 4.88)
        self.assertAlmostEqual(fraction_ho2, 0.5)
        self.assertAlmostEqual(fraction_o2minus, 0.5)
        rate = spontaneous_dismutation_rate(2e-6, 4.88, 4.88, 0.76e6, 8.5e7)
        self.assertAlmostEqual(rate, (0.76e6 + 8.5e7) * 1e-12)

    def test_sod_and_h2o2_loss_are_separate_rates(self):
        sod_loss = sod_dismutation_loss_rate(2e-6, 1e-7, 2.37e9, 7.4, 4.88)
        self.assertGreater(sod_loss, 0)
        self.assertEqual(h2o2_loss_rate(3e-6, None), 0.0)
        self.assertAlmostEqual(h2o2_loss_rate(3e-6, 2e3), 6e-3)

    def test_superoxide_h2o2_stoichiometry_without_loss(self):
        output = downstream_ros_species_resolved(
            1e-6, 2e-3, 7.4, 4.88, 0.76e6, 8.5e7,
            k_sod_m_inv_s=2.37e9, sod_m=1e-7, samples=1001,
        )
        balance = output["radical_pool_m"] + 2 * output["hydrogen_peroxide_m"]
        self.assertLess(np.max(np.abs(balance - 1e-6)), 2e-13)
        self.assertIn("H2O2_loss", output["omitted_interactions"])

    def test_species_are_reported_separately_and_balance_includes_h2o2_loss(self):
        output = downstream_ros_species_resolved(
            2e-6, 4e-3, 7.4, 4.88, 0.76e6, 8.5e7,
            k_sod_m_inv_s=2.37e9, sod_m=2e-7,
            k_o2minus_o2minus_m_inv_s=100.0,
            k_h2o2_loss_s=2e3, samples=101,
        )
        self.assertNotIn("superoxide_m", output)
        self.assertNotIn("superoxide_pool_m", output)
        self.assertTrue(np.allclose(
            output["ho2_m"] + output["o2minus_m"], output["radical_pool_m"]
        ))
        self.assertLess(output["max_radical_equivalent_balance_error_m"], 1e-16)
        self.assertGreater(output["accumulated_hydrogen_peroxide_loss_m"][-1], 0)
        self.assertTrue(any("fixed pH" in item for item in output["approximations"]))

    def test_optional_anion_anion_path_is_stable_even_when_extremely_fast(self):
        output = downstream_ros_species_resolved(
            1e-3, 1.0, 14.0, 4.88, 0.0, 0.0,
            k_o2minus_o2minus_m_inv_s=1e14, samples=2,
        )
        for key in ("radical_pool_m", "ho2_m", "o2minus_m", "hydrogen_peroxide_m"):
            self.assertTrue(np.all(np.isfinite(output[key])), key)
            self.assertGreaterEqual(float(np.min(output[key])), 0.0, key)
        self.assertLess(output["radical_pool_m"][-1], 1e-10)
        self.assertLess(output["max_radical_equivalent_balance_error_m"], 1e-15)

    def test_downstream_public_boundaries_reject_invalid_numerics(self):
        common = (1e-6, 1e-3, 7.4, 4.88, 0.76e6, 8.5e7)
        for bad_args in (
            (np.nan, *common[1:]),
            (common[0], np.inf, *common[2:]),
            (*common[:4], -1.0, common[5]),
        ):
            with self.assertRaises(ValueError):
                downstream_ros_species_resolved(*bad_args)
        for samples in (1, 3.5, True):
            with self.assertRaises(ValueError):
                downstream_ros_species_resolved(*common, samples=samples)

    def test_h2o2_loss_reduces_accumulation(self):
        common = dict(
            superoxide0_m=1e-6, duration_s=2e-3, pH=7.4, pKa_ho2=4.88,
            k_ho2_ho2_m_inv_s=0.76e6, k_ho2_o2minus_m_inv_s=8.5e7,
            k_sod_m_inv_s=2.37e9, sod_m=1e-7, samples=1001,
        )
        retained = downstream_ros_species_resolved(**common)
        lost = downstream_ros_species_resolved(**common, k_h2o2_loss_s=2e3)
        self.assertGreater(retained["hydrogen_peroxide_m"][-1],
                           lost["hydrogen_peroxide_m"][-1])


class AuthorityAndModeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw = validate_authority_bundle(CONFIG, PROVENANCE)

    def test_only_permitted_statuses_are_active(self):
        active = self.raw["active_model"]
        for block in (active["encounter"], active["downstream"]):
            for name, parameter in block.items():
                if parameter["enabled"]:
                    self.assertIn(parameter["status"], ACTIVATABLE_STATUSES, name)

    def test_unavailable_terms_are_disabled_not_physical_zero(self):
        encounter = self.raw["active_model"]["encounter"]
        unavailable = (
            "exchange_rad_s", "dipolar_rad_s", "oxygen_zfs_d_rad_s",
            "oxygen_zfs_e_rad_s", "local_field_proxy_rad_s",
            "radical_relaxation_s", "oxygen_relaxation_s",
        )
        for name in unavailable:
            self.assertFalse(encounter[name]["enabled"], name)
            self.assertEqual(encounter[name]["status"], "unavailable", name)
            self.assertIn("not", encounter[name]["limitations"].lower(), name)

    def test_bulk_constants_are_not_encounter_rates(self):
        encounter = self.raw["active_model"]["encounter"]
        for name in ("k_doublet_s", "k_quartet_s", "k_escape_s"):
            self.assertFalse(encounter[name]["enabled"])
            self.assertEqual(encounter[name]["unit"], "s^-1")
        for name, parameter in self.raw["parameter_records"].items():
            if parameter["reaction_id"] == REACTION_SQ_O2_ET:
                self.assertEqual(parameter["unit"], "M^-1 s^-1")

    def test_bulk_policy_rejects_same_unit_wrong_reaction(self):
        with self.assertRaisesRegex(ModelPolicyError, "not .*semiquinone"):
            evaluate_bulk_superoxide_rate(
                self.raw, "cuzn_sod_superoxide", 1e-6, 2e-4,
                use_scope=USE_LITERATURE_ARITHMETIC,
            )

    def test_condition_matched_bulk_requires_complete_exact_profile(self):
        parameter_id = "doxorubicin_semiquinone_plus_oxygen_pH6"
        with self.assertRaisesRegex(ModelPolicyError, "explicit condition"):
            evaluate_bulk_superoxide_rate(
                self.raw, parameter_id, 1e-6, 2e-4,
                use_scope=USE_CONDITION_MATCHED,
            )
        unmatched = dict(self.raw["condition_profiles"]["land_1985_pH6"])
        unmatched["pH"] = 7.4
        with self.assertRaisesRegex(ModelPolicyError, "pH"):
            evaluate_bulk_superoxide_rate(
                self.raw, parameter_id, 1e-6, 2e-4,
                use_scope=USE_CONDITION_MATCHED, condition_profile=unmatched,
            )
        matched = evaluate_bulk_superoxide_rate(
            self.raw, parameter_id, 1e-6, 2e-4,
            use_scope=USE_CONDITION_MATCHED,
            condition_profile=self.raw["condition_profiles"]["land_1985_pH6"],
        )
        self.assertAlmostEqual(matched["rate_m_s"], 0.07)
        self.assertEqual(matched["parameter"]["reaction_id"], REACTION_SQ_O2_ET)
        self.assertIn("condition-matched", matched["interpretation"])

    def test_unknown_literature_conditions_cannot_be_condition_matched(self):
        with self.assertRaisesRegex(ModelPolicyError, "does not permit"):
            evaluate_bulk_superoxide_rate(
                self.raw, "doxorubicin_semiquinone_plus_oxygen_1981",
                1e-6, 2e-4, use_scope=USE_CONDITION_MATCHED,
                condition_profile=self.raw["condition_profiles"]["land_1985_pH6"],
            )

    def test_complete_provenance_paths_are_resolved(self):
        with tempfile.TemporaryDirectory() as directory:
            temp = Path(directory)
            json_path = temp / CONFIG.name
            csv_path = temp / PROVENANCE.name
            json_path.write_text(CONFIG.read_text())
            with PROVENANCE.open(newline="") as source:
                rows = list(csv.DictReader(source))
            rows[0]["code_location"] = "active_model.encounter.no_such_parameter"
            with csv_path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaisesRegex(ValueError, "broken provenance path"):
                validate_authority_bundle(json_path, csv_path)

    def test_duplicated_parameter_value_consistency_is_checked(self):
        with tempfile.TemporaryDirectory() as directory:
            temp = Path(directory)
            raw = copy.deepcopy(self.raw)
            raw["active_model"]["encounter"]["g_radical"]["value"] = 2.1
            json_path = temp / CONFIG.name
            json_path.write_text(json.dumps(raw))
            with self.assertRaisesRegex(ValueError, "g_SQ|g_radical"):
                validate_authority_bundle(json_path, PROVENANCE)

    def test_validation_data_schema_and_inactive_status(self):
        for name, dataset in self.raw["validation_datasets"].items():
            self.assertTrue({"dataset_id", "unit", "dose", "method", "source", "enabled"}
                            .issubset(dataset), name)
            self.assertEqual(dataset["dataset_id"], name)
            self.assertFalse(dataset["enabled"], name)
            self.assertTrue("uncertainty" in dataset or "uncertainty_SE" in dataset)

    def test_evidence_mode_refuses_encounter_yield(self):
        context = load_execution_context(
            CONFIG, EVIDENCE_BACKED, provenance_path=PROVENANCE
        )
        self.assertFalse(context.non_predictive)
        self.assertIsNone(context.duration_s)
        with self.assertRaisesRegex(ModelPolicyError, "refuses encounter-level"):
            evidence_backed_encounter_prediction(context.parameters)

    def test_sensitivity_requires_opt_in_and_named_scenario(self):
        with self.assertRaisesRegex(ModelPolicyError, "explicit"):
            load_execution_context(
                CONFIG, SENSITIVITY,
                sensitivity_scenario="baseline_no_optional_interactions",
            )
        with self.assertRaisesRegex(ModelPolicyError, "named"):
            load_execution_context(CONFIG, SENSITIVITY, allow_sensitivity=True)

    def test_named_sensitivity_scenarios_are_non_predictive(self):
        for name in self.raw["sensitivity_scenarios"]:
            context = load_execution_context(
                CONFIG, SENSITIVITY, name, True, PROVENANCE
            )
            self.assertTrue(context.non_predictive, name)
            self.assertEqual(context.sensitivity_scenario, name)
            self.assertGreater(context.duration_s, 0)
            self.assertIsNotNone(context.resolved_inputs)
            overrides = self.raw["sensitivity_scenarios"][name]["encounter_overrides"]
            for parameter_name in overrides:
                self.assertEqual(
                    context.resolved_inputs[parameter_name]["origin"],
                    "sensitivity_override",
                )

    def test_cli_mode_enforcement_and_labels(self):
        evidence = subprocess.run(
            [sys.executable, "ROS.py"], cwd=ROOT, text=True, capture_output=True
        )
        self.assertEqual(evidence.returncode, 0, evidence.stderr)
        evidence_output = json.loads(evidence.stdout)
        self.assertEqual(evidence_output["execution_mode"], EVIDENCE_BACKED)
        self.assertEqual(evidence_output["encounter_level_yield"]["status"],
                         "refused_unavailable_parameters")
        self.assertIn("metadata", evidence_output)

        wrong_reaction = subprocess.run(
            [
                sys.executable, "ROS.py", "--bulk-sq-m", "1e-6",
                "--bulk-o2-m", "2e-4", "--bulk-rate-key", "cuzn_sod_superoxide",
                "--bulk-use-scope", USE_LITERATURE_ARITHMETIC,
            ],
            cwd=ROOT, text=True, capture_output=True,
        )
        self.assertNotEqual(wrong_reaction.returncode, 0)
        self.assertIn("not", wrong_reaction.stderr)

        matched = subprocess.run(
            [
                sys.executable, "ROS.py", "--bulk-sq-m", "1e-6",
                "--bulk-o2-m", "2e-4", "--bulk-use-scope", USE_CONDITION_MATCHED,
                "--bulk-condition-profile", "land_1985_pH6",
            ],
            cwd=ROOT, text=True, capture_output=True,
        )
        self.assertEqual(matched.returncode, 0, matched.stderr)
        self.assertAlmostEqual(json.loads(matched.stdout)["bulk_calculation"]["rate_m_s"], 0.07)

        denied = subprocess.run(
            [sys.executable, "ROS.py", "--mode", SENSITIVITY,
             "--sensitivity-scenario", "baseline_no_optional_interactions"],
            cwd=ROOT, text=True, capture_output=True,
        )
        self.assertNotEqual(denied.returncode, 0)
        self.assertIn("explicit", denied.stderr)

        allowed = subprocess.run(
            [sys.executable, "ROS.py", "--mode", SENSITIVITY,
             "--sensitivity-scenario", "baseline_no_optional_interactions",
             "--allow-sensitivity", "--samples", "101"],
            cwd=ROOT, text=True, capture_output=True,
        )
        self.assertEqual(allowed.returncode, 0, allowed.stderr)
        output = json.loads(allowed.stdout)
        self.assertTrue(output["non_predictive"])
        self.assertIn("NON-PREDICTIVE", output["output_label"])
        self.assertIn("resolved_inputs", output)
        self.assertIn("metadata", output)

        pipeline = subprocess.run(
            [sys.executable, "ROS.py", "--mode", SENSITIVITY,
             "--sensitivity-scenario", "full_pipeline_with_loss_probe",
             "--allow-sensitivity", "--samples", "101"],
            cwd=ROOT, text=True, capture_output=True,
        )
        self.assertEqual(pipeline.returncode, 0, pipeline.stderr)
        staged = json.loads(pipeline.stdout)["staged_pipeline"]
        self.assertIn("NON-PREDICTIVE", staged["output_label"])
        self.assertGreater(staged["semiquinone_formed_m"],
                           staged["associated_encounter_m"])
        self.assertGreaterEqual(staged["final_hydrogen_peroxide_m"], 0)
        self.assertGreaterEqual(staged["final_o2minus_m"], 0)


class CoherentValidationTests(unittest.TestCase):
    def test_unitary_embedding_consistency_only(self):
        parameters = EncounterParameters(
            field_t=(0, 0, 1e-4), g_radical=2.0035, g_oxygen=2.0023,
            exchange_rad_s=2e6, dipolar_rad_s=1e5,
            local_field_proxy_rad_s=(3e5, 0, 0),
        )
        self.assertLess(unitary_embedding_consistency_check(parameters, 1e-7), 1e-12)

    def test_real_three_qubit_circuit_executes_when_dependency_is_available(self):
        capability = independent_circuit_capability()
        self.assertIn("available", capability)
        self.assertIn("reason", capability)
        if not capability["available"]:
            self.skipTest(capability["reason"])
        parameters = EncounterParameters(
            field_t=(0, 0, 1e-4), g_radical=2.0035, g_oxygen=2.0023,
            exchange_rad_s=2e6, dipolar_rad_s=1e5,
            local_field_proxy_rad_s=(3e5, 0, 0),
        )
        result = execute_three_qubit_unitary_circuit(parameters, 1e-7)
        self.assertEqual(result["status"], "executed")
        self.assertIn("dense", result["gate_implementation"])
        self.assertLess(result["max_basis_ordering_error"], 1e-15)
        self.assertLess(result["max_statevector_error"], 1e-12)
        self.assertLess(result["doublet_observable_error"], 1e-12)
        self.assertLess(result["quartet_observable_error"], 1e-12)
        self.assertLess(result["physical_subspace_leakage_probability"], 1e-24)
        self.assertLess(result["max_doublet_quartet_observable_error"], 1e-12)
        self.assertIn("reaction", " ".join(result["does_not_validate"]))


class ReproducibleAnalysisTests(unittest.TestCase):
    def test_small_dimensionless_sweep_writes_csv_figure_and_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            completed = subprocess.run(
                [
                    sys.executable, "sensitivity_analysis.py",
                    "--reference-rate-s", "1e6", "--samples", "21",
                    "--output-dir", directory,
                ],
                cwd=ROOT, text=True, capture_output=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            paths = json.loads(completed.stdout)
            csv_path = Path(paths["csv"])
            svg_path = Path(paths["figure"])
            summary_path = Path(paths["summary"])
            self.assertTrue(csv_path.is_file())
            self.assertTrue(svg_path.is_file())
            summary = json.loads(summary_path.read_text())
            self.assertIn("source_fingerprint_sha256", summary["metadata"])
            self.assertEqual(summary["normalization"]["reference_rate_s^-1"], 1e6)
            self.assertTrue(summary["resolved_inputs"]["zero_mixing_control"])
            self.assertIn("kQ/kD=1", summary["resolved_inputs"]["spin_independent_null"])
            with csv_path.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), paths["rows"])
            self.assertTrue(any(row["kq_over_kd"] == "1.0" for row in rows))
            self.assertTrue(all(row["non_predictive"] == "True" for row in rows))


class DocumentationTraceabilityTests(unittest.TestCase):
    def test_original_topics_and_outputs_a_through_j_are_traced(self):
        review = (ROOT / "references" / "doxorubicin_parameter_review.md").read_text()
        for topic in (
            "g", "Hyperfine", "O2 parameters", "Spin selection", "ET rates",
            "Encounter lifetime / escape", "Relaxation", "Exchange, dipolar, geometry",
            "Dismutation", "Experimental ROS validation", "Spin-correlation evidence",
            "D/Q mixing interactions",
        ):
            self.assertIn(f"| {topic} |", review, topic)
        for label in "ABCDEFGHIJ":
            self.assertIn(f"| {label} —", review, label)

    def test_readme_contains_all_reproduction_commands_and_scope_warnings(self):
        readme = (ROOT / "README.md").read_text()
        for command_fragment in (
            "unittest discover", "condition_matched_prediction",
            "baseline_no_optional_interactions", "sensitivity_analysis.py",
            "--reference-rate-s", "--execute-circuit",
        ):
            self.assertIn(command_fragment, readme)
        self.assertIn("not concentrations or fluxes", readme)
        self.assertIn("not an independent Hamiltonian", readme)


if __name__ == "__main__":
    unittest.main()
