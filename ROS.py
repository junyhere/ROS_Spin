"""CLI for the doxorubicin semiquinone–triplet-oxygen spin model."""
import argparse
import json
from spin_chemistry import coherent_circuit_validation,downstream_ros,parameters_from_json,propagate_encounter

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config",default="configs/doxorubicin_parameters.json")
    parser.add_argument("--scenario",choices=("unpolarized","doublet","quartet"),default="unpolarized")
    parser.add_argument("--duration",type=float,help="Encounter time in seconds; overrides config")
    parser.add_argument("--validate-circuit",action="store_true")
    args=parser.parse_args(); params,active=parameters_from_json(args.config); duration=args.duration or active["duration_s"]
    result=propagate_encounter(params,duration,args.scenario); output={k:v for k,v in result.items() if not hasattr(v,"shape")}
    if args.validate_circuit: output["coherent_circuit_max_abs_error"]=coherent_circuit_validation(params,duration)
    downstream=active.get("downstream")
    if downstream and downstream.get("superoxide0_m") is not None:
        kinetics=downstream_ros(**downstream)
        output["downstream_final_superoxide_m"]=float(kinetics["superoxide_m"][-1]); output["downstream_final_hydrogen_peroxide_m"]=float(kinetics["hydrogen_peroxide_m"][-1])
    print(json.dumps(output,indent=2))

if __name__=="__main__": main()
