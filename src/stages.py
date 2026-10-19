"""Run one learning checkpoint and export its figures/data."""
from pathlib import Path
import argparse
import json
from experiments import basic,verify,geometry_sweeps,dielectric_sweep,screening

if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("stage",choices=["baseline","verify","design","dielectric","screening"])
    parser.add_argument("--output",type=Path,default=Path("reproduction/checkpoint"))
    parser.add_argument("--spacing",type=float,default=.5)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    action={"baseline":basic,"verify":verify,"design":geometry_sweeps,
            "dielectric":dielectric_sweep,"screening":screening}[args.stage]
    result=action(args.output,args.spacing)
    if result is not None:
        (args.output/(args.stage+".json")).write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(f"Checkpoint {args.stage} completed: {args.output.resolve()}")
