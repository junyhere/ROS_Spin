"""Evidence-gated CLI for doxorubicin semiquinone--triplet-oxygen chemistry."""
import argparse
import json
from pathlib import Path

from spin_chemistry import (
    EVIDENCE_BACKED,
    SENSITIVITY,
    ModelPolicyError,
    associate_encounters,
    bulk_superoxide_formation_rate,
    downstream_ros_species_resolved,
    evidence_backed_encounter_prediction,
    form_semiquinone,
    independent_circuit_capability,
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", default=str(ROOT / "configs" / "doxorubicin_parameters.json")
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
        "--initial-state", choices=("unpolarized", "doublet", "quartet"),
        default="unpolarized",
    )
    parser.add_argument("--duration", type=float, help="sensitivity-only duration override")
    parser.add_argument("--samples", type=int, default=1001)
    parser.add_argument("--bulk-sq-m", type=float)
    parser.add_argument("--bulk-o2-m", type=float)
    parser.add_argument(
        "--bulk-rate-key",
        default="doxorubicin_semiquinone_plus_oxygen_pH6",
        help="measured_parameters key for a bulk bimolecular calculation",
    )
    parser.add_argument("--check-unitary-embedding", action="store_true")
    args = parser.parse_args()

    provenance = ROOT / "configs" / "parameter_provenance.csv"
    try:
        context = load_execution_context(
            args.config,
            mode=args.mode,
            sensitivity_scenario=args.sensitivity_scenario,
            allow_sensitivity=args.allow_sensitivity,
            provenance_path=provenance,
        )
    except (ModelPolicyError, ValueError) as error:
        parser.error(str(error))

    if args.mode == EVIDENCE_BACKED:
        measured = context.authority["measured_parameters"]
        output = {
            "execution_mode": EVIDENCE_BACKED,
            "quantitative_prediction_supported": False,
            "encounter_level_yield": {
                "status": "refused_unavailable_parameters",
                "reason": str(_encounter_refusal()),
            },
            "available_bulk_bimolecular_parameters": {
                name: {"value": value["value"], "unit": value["unit"]}
                for name, value in measured.items()
                if "semiquinone_plus_oxygen" in name
            },
        }
        supplied = (args.bulk_sq_m is not None, args.bulk_o2_m is not None)
        if any(supplied) and not all(supplied):
            parser.error("--bulk-sq-m and --bulk-o2-m must be supplied together")
        if all(supplied):
            if args.bulk_rate_key not in measured:
                parser.error(f"unknown --bulk-rate-key {args.bulk_rate_key!r}")
            parameter = measured[args.bulk_rate_key]
            if parameter["status"] not in {"measured", "calculated", "fitted"}:
                parser.error("selected bulk rate is not evidence-backed")
            if parameter["unit"] != "M^-1 s^-1":
                parser.error("selected bulk rate is not bimolecular")
            output["bulk_superoxide_formation_rate_m_s"] = (
                bulk_superoxide_formation_rate(
                    args.bulk_sq_m, args.bulk_o2_m, parameter["value"]
                )
            )
            output["bulk_rate_parameter"] = args.bulk_rate_key
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
        print(json.dumps(output, indent=2))
        return

    duration = args.duration if args.duration is not None else context.duration_s
    if duration is None or duration <= 0:
        parser.error("sensitivity duration must be positive")
    result = propagate_encounter_reference(
        context.parameters, duration, args.initial_state, args.samples
    )
    output = {
        "execution_mode": SENSITIVITY,
        "sensitivity_scenario": context.sensitivity_scenario,
        "non_predictive": True,
        "output_label": "NON-PREDICTIVE SENSITIVITY OUTPUT",
        "quantitative_prediction_supported": False,
        "initial_state": args.initial_state,
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
            "final_superoxide_pool_m": float(downstream["superoxide_pool_m"][-1]),
            "final_hydrogen_peroxide_m": float(
                downstream["hydrogen_peroxide_m"][-1]
            ),
            "limitations": pipeline["limitations"],
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
    print(json.dumps(output, indent=2))


def _encounter_refusal():
    try:
        evidence_backed_encounter_prediction()
    except ModelPolicyError as error:
        return error
    raise AssertionError("evidence-backed encounter call unexpectedly succeeded")


if __name__ == "__main__":
    main()
