"""Evidence-gated CLI for doxorubicin semiquinone--triplet-oxygen chemistry."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from spin_chemistry import (
    EVIDENCE_BACKED,
    REACTION_SQ_O2_ET,
    SENSITIVITY,
    USE_CONDITION_MATCHED,
    USE_LITERATURE_ARITHMETIC,
    USE_SENSITIVITY_ONLY,
    ModelPolicyError,
    associate_encounters,
    downstream_ros_species_resolved,
    evaluate_bulk_superoxide_rate,
    evidence_backed_encounter_prediction,
    form_semiquinone,
    independent_circuit_capability,
    execute_three_qubit_unitary_circuit,
    load_execution_context,
    primary_superoxide_formation,
    propagate_encounter_reference,
    unitary_embedding_consistency_check,
)


ROOT = Path(__file__).resolve().parent


def _scalar_results(result):
    return {
        key: value
        for key, value in result.items()
        if not hasattr(value, "shape") and key != "omitted_interactions"
    }


def _run_metadata(config_path: Path, provenance_path: Path) -> dict:
    tracked_inputs = [
        ROOT / "spin_chemistry.py", ROOT / "ROS.py", config_path, provenance_path
    ]
    digest = hashlib.sha256()
    for path in tracked_inputs:
        try:
            identity = str(path.resolve().relative_to(ROOT))
        except ValueError:
            identity = str(path.resolve())
        digest.update(identity.encode())
        digest.update(path.read_bytes())
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
            text=True, capture_output=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain"], cwd=ROOT, check=True,
            text=True, capture_output=True,
        ).stdout.splitlines()
    except (OSError, subprocess.CalledProcessError):
        commit, status = None, ["git metadata unavailable"]
    return {
        "commit": commit,
        "dirty_tree": bool(status),
        "dirty_entries": status,
        "source_fingerprint_sha256": digest.hexdigest(),
        "config_identity": (
            str(config_path.resolve().relative_to(ROOT))
            if config_path.resolve().is_relative_to(ROOT)
            else str(config_path.resolve())
        ),
        "config_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
        "provenance_identity": (
            str(provenance_path.resolve().relative_to(ROOT))
            if provenance_path.resolve().is_relative_to(ROOT)
            else str(provenance_path.resolve())
        ),
        "provenance_sha256": hashlib.sha256(provenance_path.read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", default=str(ROOT / "configs" / "doxorubicin_parameters.json")
    )
    parser.add_argument(
        "--provenance",
        help="matching provenance CSV (default: parameter_provenance.csv beside config)",
    )
    parser.add_argument(
        "--mode", choices=(EVIDENCE_BACKED, SENSITIVITY), default=EVIDENCE_BACKED
    )
    parser.add_argument("--sensitivity-scenario")
    parser.add_argument(
        "--allow-sensitivity", action="store_true",
        help="explicitly opt in to non-predictive sensitivity output",
    )
    parser.add_argument(
        "--initial-state", choices=("unpolarized", "doublet", "quartet", "mixture"),
        default="unpolarized",
    )
    parser.add_argument("--p-doublet", type=float, help="required for mixture initial state")
    parser.add_argument("--duration", type=float, help="sensitivity-only duration override")
    parser.add_argument("--samples", type=int, default=1001)
    parser.add_argument("--bulk-sq-m", type=float)
    parser.add_argument("--bulk-o2-m", type=float)
    parser.add_argument(
        "--bulk-rate-key",
        default="doxorubicin_semiquinone_plus_oxygen_pH6",
        help="stable parameter_records ID for the bulk calculation",
    )
    parser.add_argument(
        "--bulk-use-scope",
        choices=(
            USE_LITERATURE_ARITHMETIC, USE_CONDITION_MATCHED, USE_SENSITIVITY_ONLY,
        ),
        help="required interpretation/policy scope when a bulk calculation is requested",
    )
    parser.add_argument(
        "--bulk-condition-profile",
        help="condition_profiles ID; required for condition_matched_prediction",
    )
    parser.add_argument("--check-unitary-embedding", action="store_true")
    parser.add_argument("--execute-circuit", action="store_true")
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    provenance = (
        Path(args.provenance).resolve()
        if args.provenance
        else config_path.with_name("parameter_provenance.csv")
    )
    try:
        context = load_execution_context(
            config_path,
            mode=args.mode,
            sensitivity_scenario=args.sensitivity_scenario,
            allow_sensitivity=args.allow_sensitivity,
            provenance_path=provenance,
        )
    except (ModelPolicyError, ValueError, OSError) as error:
        parser.error(str(error))

    if args.initial_state == "mixture" and args.p_doublet is None:
        parser.error("--initial-state mixture requires --p-doublet")
    if args.initial_state != "mixture" and args.p_doublet is not None:
        parser.error("--p-doublet is only valid with --initial-state mixture")

    if args.mode == EVIDENCE_BACKED:
        records = context.authority["parameter_records"]
        output = {
            "execution_mode": EVIDENCE_BACKED,
            "quantitative_prediction_supported": False,
            "encounter_level_yield": {
                "status": "refused_unavailable_parameters",
                "reason": str(_encounter_refusal()),
            },
            "available_bulk_bimolecular_parameters": {
                name: {
                    "parameter_id": value["parameter_id"],
                    "reaction_id": value["reaction_id"],
                    "value": value["value"],
                    "unit": value["unit"],
                    "allowed_uses": value["allowed_uses"],
                    "conditions": value["conditions"],
                    "uncertainty": value.get("uncertainty"),
                    "source": value["source"],
                    "limitations": value["limitations"],
                }
                for name, value in records.items()
                if value.get("reaction_id") == REACTION_SQ_O2_ET
            },
            "metadata": _run_metadata(config_path, provenance),
        }
        supplied = (args.bulk_sq_m is not None, args.bulk_o2_m is not None)
        if any(supplied) and not all(supplied):
            parser.error("--bulk-sq-m and --bulk-o2-m must be supplied together")
        if all(supplied):
            if args.bulk_use_scope is None:
                parser.error("bulk calculation requires explicit --bulk-use-scope")
            if args.bulk_use_scope == USE_SENSITIVITY_ONLY:
                parser.error(
                    "sensitivity-only bulk arithmetic is refused in evidence_backed mode"
                )
            profile = None
            if args.bulk_condition_profile is not None:
                profiles = context.authority.get("condition_profiles", {})
                if args.bulk_condition_profile not in profiles:
                    parser.error(
                        f"unknown --bulk-condition-profile {args.bulk_condition_profile!r}"
                    )
                profile = profiles[args.bulk_condition_profile]
            if args.bulk_use_scope == USE_CONDITION_MATCHED and profile is None:
                parser.error(
                    "condition_matched_prediction requires --bulk-condition-profile"
                )
            try:
                output["bulk_calculation"] = evaluate_bulk_superoxide_rate(
                    context.authority,
                    args.bulk_rate_key,
                    args.bulk_sq_m,
                    args.bulk_o2_m,
                    use_scope=args.bulk_use_scope,
                    condition_profile=profile,
                )
            except ModelPolicyError as error:
                parser.error(str(error))
        if args.check_unitary_embedding:
            output["unitary_embedding"] = {
                "status": "algebraic_consistency_only",
                "covers": "closed coherent Hamiltonian only",
                "does_not_cover": [
                    "Lindblad relaxation", "reaction", "escape",
                    "upstream kinetics", "downstream kinetics",
                ],
                "independent_circuit": independent_circuit_capability(),
            }
        if args.execute_circuit:
            capability = independent_circuit_capability()
            if capability["available"]:
                output["circuit_validation"] = execute_three_qubit_unitary_circuit(
                    context.parameters, 1e-7
                )
            else:
                output["circuit_validation"] = {
                    "status": "incomplete_dependency_blocker", **capability
                }
        print(json.dumps(output, indent=2))
        return

    duration = args.duration if args.duration is not None else context.duration_s
    if duration is None or duration <= 0:
        parser.error("sensitivity duration must be positive")
    result = propagate_encounter_reference(
        context.parameters, duration, args.initial_state, args.samples, args.p_doublet
    )
    output = {
        "execution_mode": SENSITIVITY,
        "sensitivity_scenario": context.sensitivity_scenario,
        "non_predictive": True,
        "output_label": "NON-PREDICTIVE SENSITIVITY OUTPUT",
        "quantitative_prediction_supported": False,
        "initial_state": args.initial_state,
        "p_doublet_requested": args.p_doublet,
        "resolved_inputs": context.resolved_inputs,
        "units": {
            "Hamiltonian_frequencies": "rad s^-1",
            "encounter_rates": "s^-1",
            "time": "s",
            "populations_and_per_encounter_yields": "dimensionless",
        },
        "scenario_limitations": context.authority["sensitivity_scenarios"][
            context.sensitivity_scenario
        ]["description"],
        "metadata": _run_metadata(config_path, provenance),
        "encounter": _scalar_results(result),
    }
    scenario = context.authority["sensitivity_scenarios"][context.sensitivity_scenario]
    if "pipeline_inputs" in scenario:
        pipeline = scenario["pipeline_inputs"]
        semiquinone_m = form_semiquinone(
            pipeline["quinone_m"], pipeline["semiquinone_formation_fraction"]
        )
        encounter_m = associate_encounters(
            semiquinone_m, pipeline["oxygen_m"],
            pipeline["encounter_association_fraction"],
        )
        primary_superoxide_m = primary_superoxide_formation(
            encounter_m, result["primary_superoxide_yield"]
        )
        evidence = context.authority["reference_evidence"]["downstream"]
        downstream = downstream_ros_species_resolved(
            primary_superoxide_m,
            pipeline["downstream_duration_s"],
            pipeline["pH"],
            evidence["HO2_pKa"]["value"],
            evidence["k_HO2_HO2"]["value"],
            evidence["k_HO2_O2minus"]["value"],
            k_sod_m_inv_s=evidence["k_CuZnSOD_overall"]["value"],
            sod_m=pipeline["sod_m"],
            k_h2o2_loss_s=pipeline["k_h2o2_loss_s"],
            samples=pipeline["samples"],
        )
        output["staged_pipeline"] = {
            "output_label": "NON-PREDICTIVE SENSITIVITY OUTPUT",
            "semiquinone_formed_m": semiquinone_m,
            "associated_encounter_m": encounter_m,
            "primary_superoxide_m": primary_superoxide_m,
            "final_radical_pool_m": float(downstream["radical_pool_m"][-1]),
            "final_ho2_m": float(downstream["ho2_m"][-1]),
            "final_o2minus_m": float(downstream["o2minus_m"][-1]),
            "final_hydrogen_peroxide_m": float(
                downstream["hydrogen_peroxide_m"][-1]
            ),
            "limitations": pipeline["limitations"],
            "kinetic_approximations": downstream["approximations"],
            "radical_equivalent_balance_error_m": downstream[
                "max_radical_equivalent_balance_error_m"
            ],
        }
    if args.check_unitary_embedding:
        output["unitary_embedding"] = {
            "status": "algebraic_consistency_only",
            "max_abs_error": unitary_embedding_consistency_check(
                context.parameters, duration
            ),
            "covers": "closed coherent Hamiltonian only",
            "does_not_cover": [
                "Lindblad relaxation", "reaction", "escape",
                "upstream kinetics", "downstream kinetics",
            ],
            "independent_circuit": independent_circuit_capability(),
        }
    if args.execute_circuit:
        capability = independent_circuit_capability()
        if capability["available"]:
            output["circuit_validation"] = execute_three_qubit_unitary_circuit(
                context.parameters, duration
            )
        else:
            output["circuit_validation"] = {
                "status": "incomplete_dependency_blocker", **capability
            }
    print(json.dumps(output, indent=2))


def _encounter_refusal():
    try:
        evidence_backed_encounter_prediction()
    except ModelPolicyError as error:
        return error
    raise AssertionError("evidence-backed encounter call unexpectedly succeeded")


if __name__ == "__main__":
    main()
