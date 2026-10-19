"""Run from any directory; defaults save beside this script."""
from pathlib import Path
import argparse
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent/"src"))
from experiments import run

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="2D depleted MOSFET electrostatic experiments")
    parser.add_argument("--output",type=Path,default=Path(__file__).resolve().parent/"figures")
    parser.add_argument("--spacing",type=float,default=0.5,help="Base mesh spacing (nm); verification refines to h/4")
    args=parser.parse_args()
    run(args.output,args.spacing)
